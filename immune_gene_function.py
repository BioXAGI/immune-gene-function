"""
AGENT NOTE: run `python immune_gene_function.py --describe` for a JSON schema of methods, output
columns and interpretation rules; every run also writes run_summary.json. See AGENT_GUIDE.md.

Infer the immune function of genes from bulk RNA-seq, in one cohort or pan-cancer.

Per cohort, for each gene X:
  1. Signature scores (mean z-score of marker genes): CTL, IFNG, exhaustion, stem-like T, Treg,
     NK, B cell, monocyte, M1-like, M2-like, macrophage, neutrophil, MDSC-like, cDC, CAF,
     endothelial, VEGF/angiogenesis, antigen presentation, leukocyte (PTPRC).
  2. Correlations (outcome-free; every tumor with RNA-seq is used). Signatures containing X are
     recomputed without X.
       rho_*        raw Spearman
       pcorr_*      partial Spearman, adjusted per --pcorr-adjust (default basic: leukocyte + CAF)
       pcorr_full_* partial Spearman adjusted for leukocyte, CAF, exhaustion, Treg, macrophage,
                    M2-like, monocyte, neutrophil, MDSC-like, VEGF
       (+ MSI-H in both when --msi and MSI data exist)
       --decompose  decomposition.csv: raw / each covariate alone / all / all-minus-one
  3. Survival (skipped with --no-survival). TIDE-style dysfunction test, Cox:
        hazard ~ S + X + S*X + stage + age + CAF [+ MSI-H] [+ --adjust-scores]
     S = CTL by default (--interaction-score to change), plus the hazard ratio of S within the
     X-low and X-high halves. Output columns keep the CTL/CTLxX names for whichever S is used.
  4. Benjamini-Hochberg FDR within each cohort; automatic interpretation label.
Optional:
  --score-prognosis  independent prognostic value of every signature score
                     (score + stage + age + CAF [+ MSI-H]), with and without CTL
  --km               Kaplan-Meier panels: S-high vs S-low inside X-low and X-high
Across cohorts:
  5. Fixed- and random-effects (DerSimonian-Laird) meta-analysis of log HR(S*X) and log HR(X).
  6. Heatmaps of interaction z-scores and partial correlations.

Data
  --source tcga   UCSC Xena TCGA hub, downloaded per cohort and cached (default)
  --source custom your own files
      --expr      genes x samples (tsv/csv, optionally .gz), first column = gene symbol,
                  log-scale values (e.g. log2(TPM+1))
      --clinical  optional; samples x columns with <time col>, <event col> and optional
                  age, stage, msi, cohort

Examples
  python immune_gene_function.py --cohorts COADREAD LUAD SKCM --genes CD274 IDO1 LRRC33 --msi
  python immune_gene_function.py --cohorts ALL --gene-file genes.txt --msi --out pancan_all
  python immune_gene_function.py --cohorts COADREAD --genes IDO1 TGFB1 --pcorr-adjust full \\
         --decompose CTL M2_LIKE CAF --score-prognosis --km
  python immune_gene_function.py --cohorts COADREAD --genes TGFB1 --interaction-score TREG \\
         --adjust-scores M2_LIKE NEUTROPHIL
  python immune_gene_function.py --source custom --expr expr.tsv.gz --genes IDO1 --no-survival
  python immune_gene_function.py --source custom --expr expr.tsv.gz --clinical clin.tsv \\
         --time-col OS.time --event-col OS --genes IDO1

Requires: numpy, pandas, scipy, lifelines, matplotlib
"""
import argparse
import gzip
import os
import sys
import textwrap
import urllib.parse
import urllib.request
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

XENA = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
TCGA_COHORTS = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COADREAD", "DLBC", "ESCA", "GBM", "HNSC",
                "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD",
                "PCPG", "PRAD", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]
# Tumor sample-type codes in order of preference (01 primary, 03 blood-derived primary, 06 metastatic)
TUMOR_CODES = ["01", "03", "06"]
MSI_COLUMNS = ["CDE_ID_3226963", "MSI_updated_Oct62011", "msi", "MSI", "msi_status"]
MIN_EVENTS = 20

SIGS = {
    "CTL": ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1"],
    "IFNG": ["IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA"],
    "EXHAUSTION": ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "CTLA4"],
    "TSTEM": ["SELL", "CCR7", "IL7R"],
    "TREG": ["FOXP3", "IL2RA", "IKZF2"],
    "NK": ["KLRD1", "NCR1", "KLRF1", "SH2D1B"],
    "BCELL": ["MS4A1", "CD79A", "CD19", "BANK1"],
    "MONOCYTE": ["CD14", "FCN1", "VCAN", "S100A12"],
    "M1_LIKE": ["CXCL9", "CXCL10", "IL1B", "CD80", "NOS2"],
    "M2_LIKE": ["CD163", "MRC1", "MSR1", "F13A1", "CD209"],
    "MACROPHAGE": ["CD68", "CSF1R", "C1QA", "C1QB", "C1QC"],
    "NEUTROPHIL": ["CSF3R", "FCGR3B", "CXCR2", "FPR1", "CMTM2"],
    "MDSC_LIKE": ["S100A8", "S100A9", "ARG1", "OLR1"],
    "CDC": ["CLEC9A", "XCR1", "CD1C", "LAMP3"],
    "CAF": ["FAP", "COL1A1", "COL1A2", "COL3A1", "PDGFRB", "THBS2"],
    "ENDOTHELIAL": ["PECAM1", "VWF", "CDH5", "KDR"],
    "VEGF_ANGIO": ["VEGFA", "VEGFC", "ESM1", "ANGPT2", "PGF"],
    "ANTIGEN_PRES": ["B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2"],
    "LEUKOCYTE": ["PTPRC"],
}
SCORES_REPORTED = [s for s in SIGS if s != "LEUKOCYTE"]
ADJUST_SETS = {
    "basic": ["LEUKOCYTE", "CAF"],
    "full": ["LEUKOCYTE", "CAF", "EXHAUSTION", "TREG", "MACROPHAGE", "M2_LIKE", "MONOCYTE",
             "NEUTROPHIL", "MDSC_LIKE", "VEGF_ANGIO"],
}
PROGNOSIS_SCORES = ["CTL", "IFNG", "EXHAUSTION", "TREG", "MONOCYTE", "MACROPHAGE", "M1_LIKE",
                    "M2_LIKE", "NEUTROPHIL", "MDSC_LIKE", "CAF", "VEGF_ANGIO"]
HEATMAP_SCORES = ["CTL", "EXHAUSTION", "TREG", "M2_LIKE", "NEUTROPHIL", "CAF", "VEGF_ANGIO"]
# Current symbol -> older symbols used in TCGA/Xena matrices
ALIASES = {
    "VSIR": ["C10orf54"], "NRROS": ["LRRC33"], "LRRC33": ["NRROS"], "SH2D1B": ["EAT2"],
    "CMTM2": ["CKLFSF2"], "PDCD1LG2": ["PDL2"], "HAVCR2": ["TIM3"], "TIGIT": ["VSIG9"],
    "KIRREL1": ["KIRREL"], "SEPTIN7": ["SEPT7"], "MARCHF1": ["MARCH1"],
}


# ============================================================================= data loading
def fetch(name, cache):
    os.makedirs(cache, exist_ok=True)
    path = os.path.join(cache, name.replace("/", "__"))
    if not os.path.exists(path):
        print(f"  downloading {name} ...", flush=True)
        tmp = path + ".part"
        with urllib.request.urlopen(XENA + urllib.parse.quote(name, safe=""), timeout=600) as r, open(tmp, "wb") as f:
            while chunk := r.read(1 << 20):
                f.write(chunk)
        os.replace(tmp, path)
    return path


def read_expr(path, genes_needed):
    """Stream a genes x samples matrix, keeping only the needed genes."""
    opener = gzip.open if path.endswith(".gz") else open
    sep = "," if ".csv" in path else "\t"
    rows = {}
    with opener(path, "rt") as f:
        cols = f.readline().rstrip("\n").split(sep)[1:]
        for line in f:
            g, rest = line.split(sep, 1)
            g = g.strip('"')
            if g in genes_needed:
                rows[g] = pd.to_numeric(pd.Series(rest.rstrip("\n").split(sep)), errors="coerce").values
    return pd.DataFrame(rows, index=[c.strip('"') for c in cols])


def genes_needed(genes):
    return set(genes) | {g for v in SIGS.values() for g in v} | {a for g in genes for a in ALIASES.get(g, [])}


def resolve_aliases(df, genes):
    rename = {}
    for g in genes:
        if g not in df:
            hit = next((a for a in ALIASES.get(g, []) if a in df), None)
            if hit:
                rename[hit] = g
    return df.rename(columns=rename)


def pick_tumor_samples(index):
    s = pd.Series(index, index=index)
    keep = []
    for _, grp in s.groupby(s.str[:12]):
        for code in TUMOR_CODES:
            hit = grp[grp.str[13:15] == code]
            if len(hit):
                keep.append(hit.iloc[0])
                break
    return keep


def parse_stage(series):
    st = series.astype(str).str.upper()
    return pd.Series(np.select([st.str.contains("IV"), st.str.contains("III"),
                                st.str.contains("II"), st.str.contains(r"\bI|STAGE I")],
                               [4, 3, 2, 1], np.nan), index=series.index)


def parse_msi(series):
    s = series.astype(str).str.upper()
    out = pd.Series(np.nan, index=series.index)
    out[s.str.contains("MSI-H|MSI_H|MSIH|^H$|HIGH", regex=True)] = 1.0
    out[s.str.contains("MSS|MSI-L|MSI_L|STABLE|^L$|LOW", regex=True)] = 0.0
    return out


def load_tcga(cohort, genes, cache):
    expr = read_expr(fetch(f"TCGA.{cohort}.sampleMap/HiSeqV2.gz", cache), genes_needed(genes))
    surv = pd.read_csv(fetch(f"survival/{cohort}_survival.txt", cache), sep="\t", index_col=0)
    clin = pd.read_csv(fetch(f"TCGA.{cohort}.sampleMap/{cohort}_clinicalMatrix", cache),
                       sep="\t", index_col=0, low_memory=False)
    expr = expr.loc[pick_tumor_samples(expr.index)]
    df = expr.join(surv.drop(columns=[c for c in ["_PATIENT", "Redaction"] if c in surv]), how="left")
    age_col = next((c for c in ["age_at_initial_pathologic_diagnosis", "age"] if c in clin), None)
    stage_col = next((c for c in ["pathologic_stage", "clinical_stage", "ajcc_pathologic_tumor_stage"] if c in clin), None)
    msi_col = next((c for c in MSI_COLUMNS if c in clin and clin[c].notna().sum() > 20), None)
    df["age"] = pd.to_numeric(clin[age_col], errors="coerce").reindex(df.index) if age_col else np.nan
    df["stage"] = parse_stage(clin[stage_col]).reindex(df.index) if stage_col else np.nan
    df["MSI_H"] = parse_msi(clin[msi_col]).reindex(df.index) if msi_col else np.nan
    return resolve_aliases(df, genes)


def load_custom(args, genes):
    expr = read_expr(args.expr, genes_needed(genes))
    if args.clinical:
        sep = "," if args.clinical.endswith(".csv") else "\t"
        df = expr.join(pd.read_csv(args.clinical, sep=sep, index_col=0), how="left")
        if args.endpoints:
            df = df.rename(columns={args.time_col: f"{args.endpoints[0]}.time", args.event_col: args.endpoints[0]})
    else:
        df = expr.copy()
    df["age"] = pd.to_numeric(df["age"], errors="coerce") if "age" in df else np.nan
    df["stage"] = parse_stage(df["stage"]) if "stage" in df else np.nan
    df["MSI_H"] = parse_msi(df["msi"]) if "msi" in df else np.nan
    df = resolve_aliases(df, genes)
    if "cohort" in df:
        return {c: df[df["cohort"] == c] for c in df["cohort"].dropna().unique()}
    return {"CUSTOM": df}


# ============================================================================= scores and statistics
def zscore(s):
    return (s - s.mean()) / s.std()


def add_scores(df):
    for name, gl in SIGS.items():
        present = [g for g in gl if g in df and df[g].std() > 0]
        df[name] = pd.concat([zscore(df[g]) for g in present], axis=1).mean(axis=1) if present else np.nan
    return df


def score_without(df, name, gene):
    """Signature score recomputed without `gene`, so a gene is never correlated with itself."""
    allg = [g for g in SIGS[name] if g in df and df[g].std() > 0]
    gl = [g for g in allg if g != gene]
    if not gl:
        return pd.Series(np.nan, index=df.index)
    if len(gl) == len(allg):
        return df[name]
    return pd.concat([zscore(df[g]) for g in gl], axis=1).mean(axis=1)


def gene_frame(df, gene, names):
    cols = {gene: df[gene]}
    for n in names:
        if n in SIGS:
            cols[n] = score_without(df, n, gene)
        elif n in df:
            cols[n] = df[n]
    return pd.DataFrame(cols)


def partial_spearman(df, x, y, covars):
    d = df[[x, y] + covars].dropna()
    if len(d) < 20:
        return np.nan, np.nan
    r = d.rank()
    C = np.column_stack([np.ones(len(r))] + [r[c] for c in covars])
    res = lambda v: v - C @ np.linalg.lstsq(C, v, rcond=None)[0]
    return stats.pearsonr(res(r[x].values), res(r[y].values))


def decompose(gf, gene, target, covars):
    """Raw, each covariate alone, all covariates, and all-minus-one (unique contribution)."""
    cov = [c for c in covars if c != target and c in gf]
    raw = partial_spearman(gf, gene, target, [])[0]
    full = partial_spearman(gf, gene, target, cov)[0]
    rows = [("raw", raw), ("all", full)]
    for c in cov:
        rows.append((f"only {c}", partial_spearman(gf, gene, target, [c])[0]))
        rows.append((f"all minus {c}", partial_spearman(gf, gene, target, [k for k in cov if k != c])[0]))
    return [{"adjustment": a, "pcorr": v, "change_from_raw": v - raw,
             "unique_contribution": (v - full) if a.startswith("all minus") else np.nan}
            for a, v in rows]


def usable_covariates(d, candidates, min_frac=0.7):
    return [c for c in candidates if c in d and d[c].notna().mean() >= min_frac and d[c].nunique() > 1]


def has_endpoint(df, ep):
    return ep in df and f"{ep}.time" in df and df[ep].notna().sum() > 0


def cox_interaction(df, gene, time, event, use_msi, extra=(), score="CTL"):
    from lifelines import CoxPHFitter
    base = ["stage", "age", "CAF"] + (["MSI_H"] if use_msi else []) + [e for e in extra if e != score]
    cov = usable_covariates(df, list(dict.fromkeys(base)))
    d = df[[time, event, score, gene] + cov].dropna().rename(columns={score: "CTL"})
    d = d[d[time] > 0].copy()
    out = {"n": len(d), "events": int(d[event].sum()) if len(d) else 0, "covariates": "+".join(cov)}
    if out["events"] < MIN_EVENTS or d[gene].std() == 0:
        return out
    d["CTL"], d["X"] = zscore(d["CTL"]), zscore(d[gene])
    d["CTLxX"] = d["CTL"] * d["X"]
    try:
        s = CoxPHFitter(penalizer=0.01).fit(d[[time, event, "CTL", "X", "CTLxX"] + cov], time, event).summary
        for term in ["X", "CTL", "CTLxX"]:
            out[f"logHR_{term}"] = s.loc[term, "coef"]
            out[f"se_{term}"] = s.loc[term, "se(coef)"]
            out[f"HR_{term}"] = s.loc[term, "exp(coef)"]
            out[f"p_{term}"] = s.loc[term, "p"]
        hi = d["X"] >= d["X"].median()
        for lab, sub in [("Xlow", d[~hi]), ("Xhigh", d[hi])]:
            ss = CoxPHFitter(penalizer=0.01).fit(sub[[time, event, "CTL"] + cov], time, event).summary
            out[f"HR_CTL_{lab}"] = ss.loc["CTL", "exp(coef)"]
            out[f"p_CTL_{lab}"] = ss.loc["CTL", "p"]
    except Exception as e:
        out["error"] = str(e)[:80]
    return out


def score_prognosis(name, df, endpoints, use_msi):
    """Independent prognostic value of each signature score, without and with CTL."""
    from lifelines import CoxPHFitter
    msi = ["MSI_H"] if use_msi else []
    rows = []
    for ep in endpoints:
        if not has_endpoint(df, ep):
            continue
        for s in PROGNOSIS_SCORES:
            for extra in ([], ["CTL"]):
                if s in extra:
                    continue
                cov = usable_covariates(df, ["stage", "age"] + ([] if s == "CAF" else ["CAF"]) + msi + extra)
                d = df[[f"{ep}.time", ep, s] + cov].dropna()
                d = d[d[f"{ep}.time"] > 0].copy()
                if d[ep].sum() < MIN_EVENTS:
                    continue
                d[s] = zscore(d[s])
                try:
                    r = CoxPHFitter(penalizer=0.01).fit(d, f"{ep}.time", ep).summary.loc[s]
                except Exception:
                    continue
                rows.append(dict(cohort=name, endpoint=ep, score=s, model="+CTL" if extra else "base",
                                 covariates="+".join(cov), n=len(d), events=int(d[ep].sum()),
                                 HR_per_SD=r["exp(coef)"], CI_low=r["exp(coef) lower 95%"],
                                 CI_high=r["exp(coef) upper 95%"], p=r["p"]))
    return rows


def bh(p):
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.sum():
        pv = p[ok]
        order = np.argsort(pv)
        ranked = np.minimum.accumulate((pv[order] * ok.sum() / (np.arange(ok.sum()) + 1))[::-1])[::-1]
        tmp = np.empty_like(pv)
        tmp[order] = np.minimum(ranked, 1)
        q[ok] = tmp
    return q


def meta(logs, ses):
    logs, ses = np.asarray(logs, float), np.asarray(ses, float)
    ok = ~(np.isnan(logs) | np.isnan(ses)) & (ses > 0)
    logs, ses = logs[ok], ses[ok]
    k = len(logs)
    if k < 2:
        return dict(k=k)
    w = 1 / ses**2
    fe = np.sum(w * logs) / w.sum()
    fe_se = np.sqrt(1 / w.sum())
    Q = np.sum(w * (logs - fe) ** 2)
    tau2 = max(0.0, (Q - (k - 1)) / (w.sum() - np.sum(w**2) / w.sum()))
    wr = 1 / (ses**2 + tau2)
    re = np.sum(wr * logs) / wr.sum()
    re_se = np.sqrt(1 / wr.sum())
    pz = lambda b, s: 2 * stats.norm.sf(abs(b / s))
    return dict(k=k, HR_fixed=np.exp(fe), p_fixed=pz(fe, fe_se), HR_random=np.exp(re),
                CI_low=np.exp(re - 1.96 * re_se), CI_high=np.exp(re + 1.96 * re_se),
                p_random=pz(re, re_se), I2=max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0, tau2=tau2)


# ============================================================================= per-cohort analysis
def analyze_cohort(name, df, genes, args, decomp_scores, decomp_rows):
    has_msi = args.msi and df["MSI_H"].notna().mean() >= 0.7 and df["MSI_H"].nunique() > 1
    msi = ["MSI_H"] if has_msi else []
    adj = list(ADJUST_SETS.get(args.pcorr_adjust) or args.pcorr_adjust.split(",")) + msi
    full = ADJUST_SETS["full"] + msi
    rows = []
    for g in genes:
        if g not in df or df[g].std() == 0:
            continue
        gf = gene_frame(df, g, list(dict.fromkeys(SCORES_REPORTED + ["LEUKOCYTE"] + msi)))
        r = {"cohort": name, "gene": g, "n_samples": len(df), "msi_adjusted": has_msi,
             "pcorr_adjusted_for": "+".join(adj), "interaction_score": args.interaction_score,
             "extra_covariates": "+".join(args.adjust_scores), "median_expr": df[g].median(),
             "in_signatures": "+".join(n for n, gl in SIGS.items() if g in gl)}
        for s in SCORES_REPORTED:
            r[f"rho_{s}"] = stats.spearmanr(gf[g], gf[s], nan_policy="omit")[0]
            r[f"pcorr_{s}"], r[f"pcorr_p_{s}"] = partial_spearman(gf, g, s, [c for c in adj if c != s])
            r[f"pcorr_full_{s}"], r[f"pcorr_full_p_{s}"] = partial_spearman(gf, g, s, [c for c in full if c != s])
            if decomp_rows is not None and s in decomp_scores:
                decomp_rows += [{"cohort": name, "gene": g, "score": s, **d} for d in decompose(gf, g, s, full)]
        for ep in args.endpoints:
            if has_endpoint(df, ep):
                for k, v in cox_interaction(df, g, f"{ep}.time", ep, has_msi,
                                            args.adjust_scores, args.interaction_score).items():
                    r[f"{ep}_{k}"] = v
        rows.append(r)
    res = pd.DataFrame(rows)
    for ep in args.endpoints:
        for term in ["CTLxX", "X"]:
            if f"{ep}_p_{term}" in res:
                res[f"{ep}_q_{term}"] = bh(res[f"{ep}_p_{term}"])
    if not res.empty:
        res["interpretation"] = res.apply(lambda r: interpret(r, args.endpoints, args.interaction_score), axis=1)
    return res, has_msi


def interpret(r, endpoints, score="CTL"):
    tags = []
    for ep in endpoints:
        p, hr = r.get(f"{ep}_p_CTLxX"), r.get(f"{ep}_HR_CTLxX")
        if pd.notna(p) and p < 0.05:
            tags.append(f"{'weakens' if hr > 1 else 'strengthens'} {score} benefit ({ep})"
                        + (" = dysfunction-like" if score == "CTL" and hr > 1 else ""))
        p, hr = r.get(f"{ep}_p_X"), r.get(f"{ep}_HR_X")
        if pd.notna(p) and p < 0.05:
            tags.append(f"{'adverse' if hr > 1 else 'favorable'} prognostic ({ep})")
    if r.get("pcorr_CTL", 0) <= -0.2 and r.get("pcorr_p_CTL", 1) < 0.05:
        tags.append("exclusion-like")
    top = sorted(((r.get(f"pcorr_{s}", np.nan), s) for s in SCORES_REPORTED if s != "CTL"),
                 key=lambda t: -np.nan_to_num(t[0], nan=-9))[:2]
    tags += [f"tracks {s} ({v:+.2f})" for v, s in top if pd.notna(v) and v >= 0.3]
    full_top = [s for s in SCORES_REPORTED if r.get(f"pcorr_full_{s}", 0) >= 0.3]
    if full_top:
        tags.append("specific after full adjustment: " + ", ".join(full_top))
    return "; ".join(tags) or "no clear immune signal"


# ============================================================================= plots
def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def heatmap(mat, title, path, vmax=None, fmt="{:.1f}"):
    if mat.empty or mat.isna().all().all():
        return
    plt = _plt()
    vmax = vmax or np.nanmax(np.abs(mat.values))
    fig, ax = plt.subplots(figsize=(max(4.5, 1.5 + 0.55 * mat.shape[1]), 1.8 + 0.35 * mat.shape[0]))
    im = ax.imshow(mat.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(mat.shape[1]), mat.columns, rotation=90)
    ax.set_yticks(range(mat.shape[0]), mat.index)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if pd.notna(mat.values[i, j]):
                ax.text(j, i, fmt.format(mat.values[i, j]), ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=ax, shrink=0.6)
    ax.set_title("\n".join(textwrap.wrap(title, 32)), fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def km_panel(df, cohort, gene, ep, score, path):
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test
    t = f"{ep}.time"
    d = df[[t, ep, score, gene]].dropna()
    d = d[d[t] > 0]
    if d[ep].sum() < MIN_EVENTS:
        return
    plt = _plt()
    hi = d[gene] >= d[gene].median()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, (lab, sub) in zip(axes, [(f"{gene} low", d[~hi]), (f"{gene} high", d[hi])]):
        c = sub[score] >= sub[score].median()
        for flag, nm, col in [(True, f"{score} high", "tab:red"), (False, f"{score} low", "tab:blue")]:
            g = sub[c == flag]
            KaplanMeierFitter().fit(g[t] / 365.25, g[ep], label=f"{nm} (n={len(g)})").plot_survival_function(
                ax=ax, ci_show=False, color=col)
        p = logrank_test(sub[c][t], sub[~c][t], sub[c][ep], sub[~c][ep]).p_value
        ax.set_title(f"{lab}: log-rank p={p:.3g}")
        ax.set_xlabel("Years")
    axes[0].set_ylabel(f"{ep} probability")
    fig.suptitle(f"{cohort}: does {gene} modify the effect of {score}? (unadjusted)")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# ============================================================================= self-description
SCHEMA = {
    "purpose": "Infer whether a gene is linked to T-cell dysfunction, T-cell exclusion, or a specific "
               "immune/stromal cell compartment, using bulk tumor RNA-seq (TCGA via UCSC Xena, or custom).",
    "method_summary": {
        "dysfunction_test": "TIDE-style (Jiang et al. Nat Med 2018): Cox model hazard ~ S + X + S*X + covariates. "
                            "S = CTL score (mean z of CD8A CD8B GZMA GZMB PRF1) unless --interaction-score. "
                            "Interaction HR > 1 means high X weakens the survival benefit of S.",
        "correlations": "Spearman raw (rho_*), partial Spearman (pcorr_*), and fully adjusted (pcorr_full_*). "
                        "Outcome-free. Signatures containing the gene are recomputed without it.",
        "meta_analysis": "DerSimonian-Laird random effects across cohorts on log HR.",
    },
    "interpretation_rules": [
        "interaction HR > 1 and p < 0.05: gene weakens CTL benefit (dysfunction-like).",
        "interaction HR < 1 and p < 0.05: gene strengthens CTL benefit.",
        "A single-cohort p near 0.05 is weak evidence: housekeeping genes reach |z| ~2 by chance. "
        "Prefer meta-analysis q_random and consistency across cohorts (I2).",
        "pcorr_full_* >= 0.3 means a link that survives adjustment for leukocyte, CAF, exhaustion, Treg, "
        "myeloid and VEGF scores; adjustment can over-correct (lower bound).",
        "Myeloid scores are NOT in the default Cox model (possible mediators); use --adjust-scores as sensitivity.",
        "Correlation is not causation; TCGA is untreated-cohort prognosis, not immunotherapy response.",
    ],
    "files": {
        "per_cohort_results.csv": "one row per cohort x gene",
        "meta_analysis.csv": "one row per gene x endpoint x term (CTLxX or X), random-effects HR",
        "decomposition.csv": "(--decompose) partial correlation raw / only-covariate / all / all-minus-one",
        "score_prognosis.csv": "(--score-prognosis) HR per SD of each signature score, base and +CTL",
        "heatmap_*.png, km_*.png": "figures",
        "run_summary.json": "machine-readable run metadata and top findings",
    },
    "columns": {
        "rho_<SCORE>": "raw Spearman correlation of gene with signature score",
        "pcorr_<SCORE>": "partial Spearman adjusted per pcorr_adjusted_for",
        "pcorr_full_<SCORE>": "partial Spearman adjusted for the full covariate set",
        "<EP>_HR_CTLxX / _p_CTLxX / _q_CTLxX": "interaction hazard ratio, p, BH-FDR within cohort (EP = OS or PFI)",
        "<EP>_HR_X / _p_X": "gene main effect on survival per SD",
        "<EP>_HR_CTL_Xlow / _Xhigh": "hazard ratio of S within X-low and X-high halves",
        "<EP>_n / _events / _covariates": "sample size, events, covariates actually used",
        "interpretation": "automatic text label",
        "in_signatures": "signatures that contain this gene",
    },
    "scores": None,  # filled at runtime from SIGS
    "endpoints": {"OS": "overall survival", "PFI": "progression-free interval", "DSS": "disease-specific survival",
                  "DFI": "disease-free interval"},
    "limits": [
        "TCGA events are few (tens to ~200 per cohort); power is low.",
        "Scores are marker-gene averages, not cell counts; bulk RNA cannot localize expression to a cell type.",
        "TIDE itself used only cohorts where CTL predicts better survival; this script does not filter on that.",
    ],
}


def describe():
    import json
    SCHEMA["scores"] = {k: v for k, v in SIGS.items()}
    SCHEMA["adjustment_sets"] = ADJUST_SETS
    print(json.dumps(SCHEMA, indent=2))


def write_summary(args, genes, res, meta_df, out):
    import json
    S = args.interaction_score
    summ = {"command": " ".join(sys.argv), "interaction_score": S, "genes": genes,
            "cohorts": sorted(res["cohort"].unique().tolist()), "endpoints": args.endpoints,
            "n_rows_per_cohort": res.groupby("cohort").size().to_dict(), "interaction_hits_p_lt_0.05": [],
            "meta_analysis_top": [], "specific_after_full_adjustment": []}
    for ep in args.endpoints:
        col = f"{ep}_p_CTLxX"
        if col in res:
            for _, r in res[res[col] < 0.05].sort_values(col).iterrows():
                summ["interaction_hits_p_lt_0.05"].append({
                    "endpoint": ep, "cohort": r["cohort"], "gene": r["gene"], "HR": round(float(r[f"{ep}_HR_CTLxX"]), 3),
                    "p": float(r[col]), "q_within_cohort": float(r[f"{ep}_q_CTLxX"])})
    if "p_random" in meta_df:
        for _, r in meta_df[meta_df["term"] == "CTLxX"].sort_values("p_random").head(20).iterrows():
            if pd.notna(r.get("p_random")):
                summ["meta_analysis_top"].append({"endpoint": r["endpoint"], "gene": r["gene"], "k": int(r["k"]),
                                                  "HR_random": float(r["HR_random"]), "p_random": float(r["p_random"]),
                                                  "q_random": float(r["q_random"]), "I2": float(r["I2"])})
    for _, r in res.iterrows():
        top = [s for s in SCORES_REPORTED if r.get(f"pcorr_full_{s}", 0) >= 0.3]
        if top:
            summ["specific_after_full_adjustment"].append({"cohort": r["cohort"], "gene": r["gene"], "scores": top})
    json.dump(summ, open(os.path.join(out, "run_summary.json"), "w"), indent=2)


# ============================================================================= main
def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_argument_group("data")
    g.add_argument("--source", choices=["tcga", "custom"], default="tcga")
    g.add_argument("--cohorts", nargs="+", default=["COADREAD"], help="TCGA codes, or ALL")
    g.add_argument("--expr", help="custom expression matrix (genes x samples)")
    g.add_argument("--clinical", help="custom clinical table (optional)")
    g.add_argument("--time-col", default="OS.time")
    g.add_argument("--event-col", default="OS")
    g.add_argument("--cache", default="xena_cache")
    g = ap.add_argument_group("genes and output")
    g.add_argument("--genes", nargs="*", default=[])
    g.add_argument("--gene-file", help="one gene symbol per line")
    g.add_argument("--out", default="immune_gene_results")
    g = ap.add_argument_group("correlations")
    g.add_argument("--msi", action="store_true", help="adjust for MSI-H where available")
    g.add_argument("--pcorr-adjust", default="basic",
                   help="'basic', 'full' or a comma list of scores for pcorr_* columns "
                        "(pcorr_full_* is always reported)")
    g.add_argument("--decompose", nargs="*", default=None, metavar="SCORE",
                   help="write decomposition.csv for these scores (flag alone = CTL IFNG MACROPHAGE M2_LIKE NEUTROPHIL)")
    g = ap.add_argument_group("survival")
    g.add_argument("--no-survival", action="store_true", help="correlations only")
    g.add_argument("--endpoints", nargs="+", default=["OS", "PFI"])
    g.add_argument("--interaction-score", default="CTL", choices=list(SIGS),
                   help="score whose survival effect is tested for modification by the gene")
    g.add_argument("--adjust-scores", nargs="*", default=[], choices=list(SIGS),
                   help="extra Cox covariates, e.g. M2_LIKE TREG (sensitivity; may be mediators)")
    g.add_argument("--score-prognosis", action="store_true",
                   help="also test each signature score as an independent prognostic factor")
    g.add_argument("--km", nargs="*", default=None, metavar="GENE",
                   help="Kaplan-Meier panels for these genes (flag alone = all genes)")
    g.add_argument("--describe", action="store_true",
                   help="print the JSON schema of methods, columns, files and interpretation rules, then exit")
    args = ap.parse_args()
    if args.describe:
        describe()
        sys.exit(0)
    if args.no_survival:
        args.endpoints = []
    return args


def main():
    args = parse_args()
    genes = list(args.genes)
    if args.gene_file:
        genes += [l.strip() for l in open(args.gene_file) if l.strip() and not l.startswith("#")]
    genes = list(dict.fromkeys(genes))
    if not genes:
        sys.exit("Provide --genes or --gene-file")
    os.makedirs(args.out, exist_ok=True)
    decomp_scores = (args.decompose or ["CTL", "IFNG", "MACROPHAGE", "M2_LIKE", "NEUTROPHIL"]) if args.decompose is not None else []
    decomp_rows = [] if decomp_scores else None
    km_genes = (args.km or genes) if args.km is not None else []

    if args.source == "custom":
        if not args.expr:
            sys.exit("--source custom needs --expr")
        if not args.clinical:
            args.endpoints = []
        args.endpoints = args.endpoints[:1]
        datasets = load_custom(args, genes)
    else:
        datasets = {}
        for c in (TCGA_COHORTS if args.cohorts == ["ALL"] else args.cohorts):
            print(f"[{c}] loading", flush=True)
            try:
                datasets[c] = load_tcga(c, genes, args.cache)
            except Exception as e:
                print(f"  skipped {c}: {e}")

    all_res, prog_rows = [], []
    for name, df in datasets.items():
        df = add_scores(df)
        missing = [g for g in genes if g not in df]
        print(f"[{name}] {len(df)} tumors" + (f"; genes not found: {missing}" if missing else ""), flush=True)
        res, has_msi = analyze_cohort(name, df, genes, args, decomp_scores, decomp_rows)
        if not res.empty:
            all_res.append(res)
        if args.score_prognosis:
            prog_rows += score_prognosis(name, df, args.endpoints, has_msi)
        for g in km_genes:
            if g in df and args.endpoints:
                km_panel(df, name, g, args.endpoints[0], args.interaction_score,
                         os.path.join(args.out, f"km_{name}_{g}_{args.endpoints[0]}.png"))
    if not all_res:
        sys.exit("No results")

    res = pd.concat(all_res, ignore_index=True)
    res.to_csv(os.path.join(args.out, "per_cohort_results.csv"), index=False, float_format="%.4g")
    if decomp_rows:
        pd.DataFrame(decomp_rows).to_csv(os.path.join(args.out, "decomposition.csv"), index=False, float_format="%.4g")
    if prog_rows:
        pd.DataFrame(prog_rows).to_csv(os.path.join(args.out, "score_prognosis.csv"), index=False, float_format="%.4g")

    meta_rows = []
    for g, grp in res.groupby("gene"):
        for ep in args.endpoints:
            for term in ["CTLxX", "X"]:
                if f"{ep}_logHR_{term}" in grp:
                    meta_rows.append({"gene": g, "endpoint": ep, "term": term,
                                      **meta(grp[f"{ep}_logHR_{term}"], grp[f"{ep}_se_{term}"])})
    meta_df = pd.DataFrame(meta_rows)
    if "p_random" in meta_df:
        meta_df["q_random"] = np.nan
        for _, idx in meta_df.groupby(["endpoint", "term"]).groups.items():
            meta_df.loc[idx, "q_random"] = bh(meta_df.loc[idx, "p_random"])
    if not meta_df.empty:
        meta_df.to_csv(os.path.join(args.out, "meta_analysis.csv"), index=False, float_format="%.4g")

    S = args.interaction_score
    for ep in args.endpoints[:1]:
        if f"{ep}_logHR_CTLxX" in res:
            z = res[f"{ep}_logHR_CTLxX"] / res[f"{ep}_se_CTLxX"]
            heatmap(res.assign(z=z).pivot(index="gene", columns="cohort", values="z"),
                    f"{ep}: {S} x gene interaction z (red = gene weakens {S} effect)",
                    os.path.join(args.out, f"heatmap_interaction_{ep}.png"), vmax=3)
    for s in HEATMAP_SCORES:
        for pref, lab in [("pcorr", args.pcorr_adjust), ("pcorr_full", "full")]:
            heatmap(res.pivot(index="gene", columns="cohort", values=f"{pref}_{s}"),
                    f"Partial Spearman with {s} (adjusted: {lab})",
                    os.path.join(args.out, f"heatmap_{pref}_{s}.png"), vmax=0.7, fmt="{:.2f}")

    # ---------------- console summary
    print("\nCorrelations (raw / basic-or-chosen / full) with CTL, exhaustion, M2-like, CAF:")
    for _, r in res.iterrows():
        print(f"  {r['cohort']:9s} {r['gene']:10s} " + "  ".join(
            f"{s}: {r[f'rho_{s}']:+.2f}/{r[f'pcorr_{s}']:+.2f}/{r[f'pcorr_full_{s}']:+.2f}"
            for s in ["CTL", "EXHAUSTION", "M2_LIKE", "CAF"]))
    if args.endpoints:
        print(f"\nGene modifies the survival effect of {S} (interaction p<0.05):")
        for ep in args.endpoints:
            col = f"{ep}_p_CTLxX"
            if col in res:
                for _, r in res[res[col] < 0.05].sort_values(col).head(20).iterrows():
                    print(f"  {ep:4s} {r['cohort']:9s} {r['gene']:10s} HR={r[f'{ep}_HR_CTLxX']:.2f} p={r[col]:.3g} "
                          f"q={r[f'{ep}_q_CTLxX']:.2g} | {S} HR in X-low {r.get(f'{ep}_HR_CTL_Xlow', np.nan):.2f}"
                          f" -> X-high {r.get(f'{ep}_HR_CTL_Xhigh', np.nan):.2f}")
    if "p_random" in meta_df:
        print(f"\nMeta-analysis (random effects) of {S} x gene interaction:")
        for _, r in meta_df[meta_df["term"] == "CTLxX"].sort_values("p_random").iterrows():
            if pd.notna(r.get("p_random")):
                print(f"  {r['endpoint']:4s} {r['gene']:10s} k={int(r['k'])} HR={r['HR_random']:.2f} "
                      f"({r['CI_low']:.2f}-{r['CI_high']:.2f}) p={r['p_random']:.3g} q={r['q_random']:.2g} I2={r['I2']:.0%}")
    if prog_rows:
        p = pd.DataFrame(prog_rows)
        sig = p[p["p"] < 0.05].sort_values("p")
        print("\nIndependent prognostic scores (p<0.05; HR per SD):")
        for _, r in sig.head(25).iterrows():
            print(f"  {r['endpoint']:4s} {r['cohort']:9s} {r['score']:11s} model {r['model']:5s} "
                  f"HR={r['HR_per_SD']:.2f} ({r['CI_low']:.2f}-{r['CI_high']:.2f}) p={r['p']:.3g}")
    write_summary(args, genes, res, meta_df, args.out)
    print(f"\nOutputs in {args.out}/ (see run_summary.json)")


if __name__ == "__main__":
    main()
