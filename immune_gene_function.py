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
import time
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
                vals = rest.rstrip("\n").split(sep)
                try:
                    rows[g] = np.array(vals, dtype=float)
                except ValueError:
                    rows[g] = pd.to_numeric(pd.Series(vals), errors="coerce").values
    return pd.DataFrame(rows, index=[c.strip('"') for c in cols])


class AnyGene:
    """Membership test that accepts every gene symbol (genome-wide scan)."""
    def __contains__(self, item):
        return True


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


def msi_labels(series):
    """Recognized labels only; unknown/ambiguous values never become MSS."""
    normalized = series.astype(str).str.strip().str.upper().str.replace("_", "-", regex=False)
    mapping = {"MSS": "MSS", "STABLE": "MSS", "MICROSATELLITE STABLE": "MSS",
               "MSI-L": "MSI-L", "MSIL": "MSI-L", "L": "MSI-L", "LOW": "MSI-L",
               "MSI-H": "MSI-H", "MSIH": "MSI-H", "H": "MSI-H", "HIGH": "MSI-H",
               "MICROSATELLITE INSTABILITY-HIGH": "MSI-H",
               "MICROSATELLITE INSTABILITY-LOW": "MSI-L"}
    return normalized.map(mapping).fillna("UNKNOWN")


def parse_msi(series):
    return msi_labels(series).map({"MSS": 0.0, "MSI-L": 0.0, "MSI-H": 1.0})


def select_msi_cohort(df, selection="pooled", msi_low="include", missing="keep"):
    """Filter before scores/models. Pooled defaults retain unknown annotations."""
    labels = df.get("MSI_status", pd.Series("UNKNOWN", index=df.index))
    if selection == "mss":
        keep = labels.isin(["MSS", "MSI-L"] if msi_low == "include" else ["MSS"])
    elif selection == "msi-h":
        keep = labels.eq("MSI-H")
    else:
        keep = pd.Series(True, index=df.index)
        if missing == "exclude":
            keep &= labels.ne("UNKNOWN")
    if missing == "error" and labels.eq("UNKNOWN").any():
        raise ValueError("Missing/unrecognized MSI annotations; use --msi-missing keep or exclude")
    out = df.loc[keep].copy()
    report = {"selection": selection, "msi_low_policy": msi_low,
              "missing_policy": missing, "n_input": len(df), "n_selected": len(out),
              "n_excluded": int((~keep).sum()), "n_unknown_input": int(labels.eq("UNKNOWN").sum()),
              "input_labels": labels.value_counts().to_dict(),
              "selected_labels": labels[keep].value_counts().to_dict(),
              "msi_source": df.attrs.get("msi_source", "custom:msi"),
              "note": "MSI-L grouped with MSS by default; MSI is not a direct pMMR assay."}
    return out, report


def endpoint_counts(df, endpoints):
    counts = {}
    for ep in endpoints:
        if has_endpoint(df, ep):
            valid = df[ep].isin([0, 1]) & np.isfinite(df[f"{ep}.time"]) & (df[f"{ep}.time"] > 0)
            counts[ep] = {"n": int(valid.sum()), "events": int(df.loc[valid, ep].sum())}
        else:
            counts[ep] = {"n": 0, "events": 0}
    return counts


# Clinical subsets usable as cohort names: name -> (parent TCGA cohort, selector on the clinical matrix)
def _tnbc(clin):
    neg = lambda col: clin[col].eq("Negative") if col in clin else pd.Series(False, index=clin.index)
    her2_neg = neg("lab_proc_her2_neu_immunohistochemistry_receptor_status") | (
        ~clin.get("lab_proc_her2_neu_immunohistochemistry_receptor_status", pd.Series(index=clin.index)).eq("Positive")
        & neg("lab_procedure_her2_neu_in_situ_hybrid_outcome_type"))
    return neg("breast_carcinoma_estrogen_receptor_status") & neg("breast_carcinoma_progesterone_receptor_status") & her2_neg


SUBSETS = {
    "BRCA_TNBC": ("BRCA", _tnbc),
    "BRCA_BASAL": ("BRCA", lambda clin: clin.get("PAM50Call_RNAseq", pd.Series(index=clin.index)).eq("Basal")),
}


def load_tcga(cohort, genes, cache, all_genes=False):
    parent, selector = SUBSETS.get(cohort, (cohort, None))
    expr = read_expr(fetch(f"TCGA.{parent}.sampleMap/HiSeqV2.gz", cache), AnyGene() if all_genes else genes_needed(genes))
    surv = pd.read_csv(fetch(f"survival/{parent}_survival.txt", cache), sep="\t", index_col=0)
    clin = pd.read_csv(fetch(f"TCGA.{parent}.sampleMap/{parent}_clinicalMatrix", cache),
                       sep="\t", index_col=0, low_memory=False)
    expr = expr.loc[pick_tumor_samples(expr.index)]
    if selector is not None:
        keep = set(clin.index[selector(clin).fillna(False).astype(bool)])
        expr = expr.loc[[i for i in expr.index if i in keep]]
    df = expr.join(surv.drop(columns=[c for c in ["_PATIENT", "Redaction"] if c in surv]), how="left")
    age_col = next((c for c in ["age_at_initial_pathologic_diagnosis", "age"] if c in clin), None)
    stage_col = next((c for c in ["pathologic_stage", "clinical_stage", "ajcc_pathologic_tumor_stage"] if c in clin), None)
    msi_col = max((c for c in MSI_COLUMNS if c in clin),
                  key=lambda c: parse_msi(clin[c]).notna().sum(), default=None)
    df["age"] = pd.to_numeric(clin[age_col], errors="coerce").reindex(df.index) if age_col else np.nan
    df["stage"] = parse_stage(clin[stage_col]).reindex(df.index) if stage_col else np.nan
    df["MSI_H"] = parse_msi(clin[msi_col]).reindex(df.index) if msi_col else np.nan
    df["MSI_status"] = msi_labels(clin[msi_col]).reindex(df.index).fillna("UNKNOWN") if msi_col else "UNKNOWN"
    out = resolve_aliases(df, genes)
    out.attrs["msi_source"] = msi_col
    out.attrs["genes"] = list(expr.columns)
    return out


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
    df["MSI_status"] = msi_labels(df["msi"]) if "msi" in df else "UNKNOWN"
    df = resolve_aliases(df, genes)
    if "cohort" in df:
        return {c: df[df["cohort"] == c] for c in df["cohort"].dropna().unique()}
    return {"CUSTOM": df}


# ============================================================================= external scores (e.g. ecotypes)
def _norm_patient(x):
    x = str(x).strip().strip('"').replace(".", "-")
    return x[:12] if x.upper().startswith("TCGA") else x


def load_score_file(path, transpose=False):
    """Table of per-sample scores (e.g. EcoTyper ecotype abundances). First column = sample id.
    TCGA barcodes in dotted or dashed form are reduced to the 12-character patient id; when a patient has
    several samples the primary tumor (sample code 01) is kept."""
    sep = "," if path.endswith(".csv") else "\t"
    t = pd.read_csv(path, sep=sep, index_col=0)
    if transpose:
        t = t.T
    t.index = t.index.map(lambda v: str(v).strip().strip('"'))
    code = pd.Series([i.replace(".", "-")[13:15] if i.upper().startswith("TCGA") else "01" for i in t.index], index=t.index)
    t = t.assign(_pat=t.index.map(_norm_patient), _pri=(code != "01").astype(int)).sort_values("_pri")
    t = t.drop_duplicates("_pat").set_index("_pat").drop(columns="_pri")
    return t.apply(pd.to_numeric, errors="coerce")


def attach_scores(df, score_df):
    if score_df is None:
        return df
    pat = pd.Series(df.index.map(_norm_patient), index=df.index)
    for c in score_df.columns:
        df[c] = pat.map(score_df[c]).values
    return df


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


def partial_spearman_details(df, x, y, covars):
    cols = list(dict.fromkeys([x, y] + covars))
    d = df.reindex(columns=cols).replace([np.inf, -np.inf], np.nan).dropna()
    out = dict(r=np.nan, p=np.nan, n=len(d), ci_low=np.nan, ci_high=np.nan,
               status="insufficient_samples", covariates="+".join(covars))
    if len(d) < 20:
        return out
    r = d.rank()
    C = np.column_stack([np.ones(len(r))] + [r[c] for c in covars])
    rank = np.linalg.matrix_rank(C)
    if rank < C.shape[1]:
        out["status"] = "singular_covariates"
        return out
    residual = lambda v: v - C @ np.linalg.lstsq(C, v, rcond=None)[0]
    a, b = residual(r[x].values), residual(r[y].values)
    if np.std(a) < 1e-10 or np.std(b) < 1e-10:
        out["status"] = "constant_residual"
        return out
    rho = float(np.clip(np.corrcoef(a, b)[0, 1], -1, 1))
    dof = len(d) - rank - 1
    if dof <= 1:
        return out
    t = rho * np.sqrt(dof / max(1e-15, 1 - rho*rho))
    z = np.arctanh(np.clip(rho, -1+1e-12, 1-1e-12))
    width = 1.96 / np.sqrt(len(d) - (rank - 1) - 3)
    out.update(r=rho, p=float(2*stats.t.sf(abs(t), dof)),
               ci_low=float(np.tanh(z-width)), ci_high=float(np.tanh(z+width)), status="ok")
    return out


def partial_spearman(df, x, y, covars):
    out = partial_spearman_details(df, x, y, covars)
    return out["r"], out["p"]


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


def cox_interaction(df, gene, time, event, use_msi, extra=(), score="CTL",
                    min_samples=40, min_events=20, events_per_parameter=10, fixed_covariates=None):
    """Penalized exploratory Cox; fail closed on sparse or unstable designs."""
    base = ["stage", "age", "CAF"] + (["MSI_H"] if use_msi else []) + [e for e in extra if e != score]
    cov = [c for c in usable_covariates(df, list(dict.fromkeys(base))) if c not in [gene, score]]
    if fixed_covariates is not None:
        cov = list(fixed_covariates) + (["MSI_H"] if use_msi else [])
    # Requested MSI adjustment must never silently disappear after complete-case filtering.
    if use_msi and "MSI_H" not in cov:
        return dict(n=0, events=0, status="msi_adjustment_unavailable")
    d = df[[time, event] + list(dict.fromkeys([score, gene] + cov))].replace([np.inf, -np.inf], np.nan).dropna()
    d = d[(d[time] > 0) & d[event].isin([0, 1])].copy()
    out = dict(n=len(d), events=int(d[event].sum()), covariates="+".join(cov),
               parameters=3+len(cov), status="not_fit", penalizer=0.01)
    if len(d) < min_samples:
        out["status"] = "insufficient_samples"
        return out
    if out["events"] < max(min_events, events_per_parameter*out["parameters"]):
        out["status"] = "insufficient_events"
        return out
    if any(d[c].nunique() < 2 for c in [score, gene] + cov):
        out["status"] = "constant_predictor"
        return out
    d["_S"], d["_X"] = zscore(d[score]), zscore(d[gene])
    d["_SxX"] = d["_S"] * d["_X"]
    design = d[["_S", "_X", "_SxX"] + cov]
    scaled = (design - design.mean()) / design.std()
    out["condition_number"] = float(np.linalg.cond(scaled))
    if np.linalg.matrix_rank(scaled) < scaled.shape[1] or out["condition_number"] > 1e8:
        out["status"] = "singular_design"
        return out
    try:
        from lifelines import CoxPHFitter
    except ImportError:
        out["status"] = "dependency_missing:lifelines"
        return out
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            fit = CoxPHFitter(penalizer=0.01).fit(d[[time, event] + list(design.columns)], time, event)
        messages = [str(w.message) for w in caught if "Convergence" in w.category.__name__]
        if messages:
            out.update(status="convergence_warning", error="; ".join(messages)[:300])
            return out
        summary = fit.summary
        if not np.isfinite(summary[["coef", "se(coef)"]].values).all() or (summary["se(coef)"] <= 0).any():
            out["status"] = "invalid_estimate"
            return out
        out["status"] = "ok"
        for term, lab in [("_X", "X"), ("_S", "CTL"), ("_SxX", "CTLxX")]:
            beta, se = summary.loc[term, ["coef", "se(coef)"]]
            out.update({f"logHR_{lab}": beta, f"se_{lab}": se,
                        f"HR_{lab}": np.exp(beta), f"p_{lab}": summary.loc[term, "p"],
                        f"CI_low_{lab}": np.exp(beta-1.96*se), f"CI_high_{lab}": np.exp(beta+1.96*se)})
        # Split-half models are descriptive only and get their own sample/event guard.
        hi = d["_X"] >= d["_X"].median()
        for lab, sub in [("Xlow", d[~hi]), ("Xhigh", d[hi])]:
            out[f"n_{lab}"], out[f"events_{lab}"] = len(sub), int(sub[event].sum())
            out[f"status_{lab}"] = "insufficient_samples_or_events"
            if len(sub) < min_samples or sub[event].sum() < max(min_events, events_per_parameter*(1+len(cov))):
                continue
            if any(sub[c].nunique() < 2 for c in ["_S"] + cov):
                out[f"status_{lab}"] = "constant_predictor"
                continue
            try:
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    ss = CoxPHFitter(penalizer=0.01).fit(sub[[time, event, "_S"] + cov], time, event).summary
                if any("Convergence" in w.category.__name__ for w in caught):
                    out[f"status_{lab}"] = "convergence_warning"
                    continue
                out[f"HR_CTL_{lab}"] = ss.loc["_S", "exp(coef)"]
                out[f"p_CTL_{lab}"] = ss.loc["_S", "p"]
                out[f"status_{lab}"] = "ok"
            except Exception as e:
                out[f"status_{lab}"] = "fit_failed"
                out[f"error_{lab}"] = str(e)[:150]
    except Exception as e:
        out.update(status="fit_failed", error=str(e)[:300])
    return out


def score_prognosis(name, df, endpoints, use_msi, extra_scores=()):
    """Independent prognostic value of each signature score, without and with CTL."""
    from lifelines import CoxPHFitter
    msi = ["MSI_H"] if use_msi else []
    rows = []
    for ep in endpoints:
        if not has_endpoint(df, ep):
            continue
        for s in list(PROGNOSIS_SCORES) + list(extra_scores):
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
        r = {"cohort": name, "gene": g, "n_samples": len(df), "msi_adjusted": has_msi, "msi_cohort": getattr(args, "msi_cohort", "pooled"),
             "msi_adjustment_status": "used" if has_msi else ("unavailable_or_constant" if args.msi else "not_requested"),
             "pcorr_adjusted_for": "+".join(adj), "interaction_score": args.interaction_score,
             "extra_covariates": "+".join(args.adjust_scores), "median_expr": df[g].median(),
             "in_signatures": "+".join(n for n, gl in SIGS.items() if g in gl)}
        for s in SCORES_REPORTED:
            r[f"rho_{s}"] = stats.spearmanr(gf[g], gf[s], nan_policy="omit")[0]
            for prefix, covars in [("pcorr", adj), ("pcorr_full", full)]:
                detail = partial_spearman_details(gf, g, s, [c for c in covars if c != s])
                r[f"{prefix}_{s}"] = detail.pop("r")
                for key, value in detail.items():
                    r[f"{prefix}_{key}_{s}"] = value
            if decomp_rows is not None and s in decomp_scores:
                decomp_rows += [{"cohort": name, "gene": g, "score": s, **d} for d in decompose(gf, g, s, full)]
        for ep in args.endpoints:
            if has_endpoint(df, ep):
                for k, v in cox_interaction(df, g, f"{ep}.time", ep, has_msi,
                                            args.adjust_scores, args.interaction_score,
                                            getattr(args, "min_model_samples", 40),
                                            getattr(args, "min_model_events", 20),
                                            getattr(args, "events_per_parameter", 10)).items():
                    r[f"{ep}_{k}"] = v
        rows.append(r)
    res = pd.DataFrame(rows)
    for ep in args.endpoints:
        for term in ["CTLxX", "X"]:
            if f"{ep}_p_{term}" in res:
                res[f"{ep}_q_{term}"] = bh(res[f"{ep}_p_{term}"])
    if not res.empty:
        # One BH family per cohort and adjustment set, across all genes x signatures.
        for pref in ["pcorr", "pcorr_full"]:
            cols = [f"{pref}_p_{score}" for score in SCORES_REPORTED]
            q = bh(res[cols].to_numpy().ravel()).reshape(len(res), len(cols))
            for j, score in enumerate(SCORES_REPORTED):
                res[f"{pref}_q_{score}"] = q[:, j]
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


# ============================================================================= genome-wide scan
CTL_GENES = set(SIGS["CTL"])
SCAN_CONTROLS_POS = ["IDO1", "TGFB1", "PDCD1LG2", "CD274", "SERPINB9", "IL10", "LGALS9", "CD47", "PTPN2", "VSIR"]
SCAN_CONTROLS_NEG = ["ACTB", "GAPDH", "TBP", "PGK1", "PPIA", "RPLP0", "HPRT1", "POLR2A"]


def _scan_chunk(task):
    """Worker: Cox interaction for a block of genes. Module-level so it can be pickled."""
    from lifelines import CoxPHFitter
    names, G, T, E, zc, C, cov_names, time_col, ev_col = task
    out = []
    base = pd.DataFrame(C, columns=cov_names)
    base[time_col], base[ev_col], base["CTL"] = T, E, zc
    for i, g in enumerate(names):
        x = G[:, i]
        sd = x.std()
        if not sd > 0:
            continue
        zx = (x - x.mean()) / sd
        d = base.copy()
        d["X"], d["CTLxX"] = zx, zc * zx
        try:
            sm = CoxPHFitter(penalizer=0.01).fit(d, time_col, ev_col).summary
        except Exception:
            continue
        out.append((g, float(np.median(x)), sm.loc["CTLxX", "coef"], sm.loc["CTLxX", "se(coef)"],
                    sm.loc["CTLxX", "p"], sm.loc["X", "coef"], sm.loc["X", "se(coef)"], sm.loc["X", "p"]))
    return out


def scan_cohort(name, df, ep, use_msi, min_expr, workers, score="CTL"):
    from concurrent.futures import ProcessPoolExecutor
    from lifelines import CoxPHFitter
    T = f"{ep}.time"
    cov = usable_covariates(df, ["stage", "age", "CAF"] + (["MSI_H"] if use_msi else []))
    base = df[[T, ep, score] + [c for c in cov if c != score]].dropna()
    cov = [c for c in cov if c != score]
    base = base[base[T] > 0]
    ev = int(base[ep].sum())
    if ev < MIN_EVENTS:
        return None, dict(cohort=name, n=len(base), events=ev, note="too few events")
    zc = zscore(base[score]).values
    ctl_hr = CoxPHFitter(penalizer=0.01).fit(base.assign(**{score: zc})[[T, ep, score] + cov], T, ep).summary.loc[score]
    skip = CTL_GENES if score == "CTL" else set()
    genes = [g for g in df.attrs["genes"] if g not in skip and g in df and df.loc[base.index, g].median() >= min_expr
             and df.loc[base.index, g].std() > 0]
    G = df.loc[base.index, genes].values.astype(float)
    C = base[cov].values.astype(float)
    nchunk = max(workers * 4, 1)
    tasks = [(genes[i::nchunk], G[:, i::nchunk], base[T].values, base[ep].values, zc, C, cov, T, ep)
             for i in range(nchunk) if len(genes[i::nchunk])]
    rows = []
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for r in ex.map(_scan_chunk, tasks):
                rows += r
    else:
        for t in tasks:
            rows += _scan_chunk(t)
    res = pd.DataFrame(rows, columns=["gene", "median_expr", "logHR_CTLxX", "se_CTLxX", "p_CTLxX",
                                      "logHR_X", "se_X", "p_X"])
    res["z_CTLxX"] = res["logHR_CTLxX"] / res["se_CTLxX"]
    res["HR_CTLxX"] = np.exp(res["logHR_CTLxX"])
    res["q_CTLxX"] = bh(res["p_CTLxX"])
    res.insert(0, "cohort", name)
    lam = float(np.median(res["z_CTLxX"] ** 2) / 0.4549)
    info = dict(cohort=name, score=score, n=len(base), events=ev, genes_tested=len(res), covariates="+".join(cov),
                CTL_HR=float(ctl_hr["exp(coef)"]), CTL_p=float(ctl_hr["p"]), lambda_inflation=lam,
                n_q_lt_0_1=int((res["q_CTLxX"] < 0.1).sum()), frac_q_lt_0_1=float((res["q_CTLxX"] < 0.1).mean()))
    return res, info


def run_scan(args):
    os.makedirs(args.out, exist_ok=True)
    workers = args.scan_workers or max(1, min(8, (os.cpu_count() or 2) - 2))
    cohorts = TCGA_COHORTS if args.cohorts == ["ALL"] else args.cohorts
    dfs, results, infos = {}, {}, []
    for ep in args.scan_endpoints:
        results[ep] = []
    for c in cohorts:
        print(f"[{c}] loading all genes", flush=True)
        try:
            df = attach_scores(add_scores(load_tcga(c, [], args.cache, all_genes=True)), SCORE_DF)
        except Exception as e:
            print(f"  skipped {c}: {e}")
            continue
        dfs[c] = df
        has_msi = args.msi and df["MSI_H"].notna().mean() >= 0.7 and df["MSI_H"].nunique() > 1
        for ep in args.scan_endpoints:
            t0 = time.time()
            res, info = scan_cohort(c, df, ep, has_msi, args.scan_min_expr, workers, args.interaction_score)
            info["endpoint"] = ep
            infos.append(info)
            if res is None:
                print(f"  {c} {ep}: skipped ({info['note']}, {info['events']} events)")
                continue
            res.to_csv(os.path.join(args.out, f"scan_{c}_{ep}.csv"), index=False, float_format="%.5g")
            results[ep].append(res)
            print(f"  {c} {ep}: {info['genes_tested']} genes, {info['events']} events, {args.interaction_score} HR={info['CTL_HR']:.2f} "
                  f"(p={info['CTL_p']:.2f}), lambda={info['lambda_inflation']:.2f}, FDR<0.1: {info['n_q_lt_0_1']} "
                  f"[{time.time()-t0:.0f}s]", flush=True)
    pd.DataFrame(infos).to_csv(os.path.join(args.out, "scan_cohort_summary.csv"), index=False, float_format="%.4g")

    for ep, lst in results.items():
        if not lst:
            continue
        if len(lst) < 2:
            print(f"\n{ep}: one cohort only, so no pooled meta-analysis; see scan_<cohort>_{ep}.csv")
            continue
        allr = pd.concat(lst)
        meta_rows = []
        for g, grp in allr.groupby("gene"):
            if len(grp) < 2:
                continue
            mm = meta(grp["logHR_CTLxX"].values, grp["se_CTLxX"].values)
            if "p_random" in mm:
                meta_rows.append(dict(gene=g, k=mm["k"], HR_random=mm["HR_random"], CI_low=mm["CI_low"], CI_high=mm["CI_high"],
                                      p_random=mm["p_random"], I2=mm["I2"], n_cohorts_p_lt_0_05=int((grp["p_CTLxX"] < 0.05).sum()),
                                      n_same_direction=int((np.sign(grp["logHR_CTLxX"]) == np.sign(np.log(mm["HR_random"]))).sum()),
                                      median_expr=float(grp["median_expr"].median())))
        mt = pd.DataFrame(meta_rows).sort_values("p_random")
        mt["q_random"] = bh(mt["p_random"])
        mt["rank_dysfunction"] = np.nan
        up = mt[mt["HR_random"] > 1].sort_values("p_random")
        mt.loc[up.index, "rank_dysfunction"] = np.arange(1, len(up) + 1)
        # annotate hits with partial correlations (adjusted for leukocyte + CAF [+ MSI-H]) in each cohort
        hits = pd.concat([up.head(args.scan_top), mt[mt["HR_random"] < 1].head(args.scan_top // 2)])["gene"].tolist()
        ann = []
        for g in hits:
            r = {"gene": g}
            for s_ in ["CTL", "EXHAUSTION", "TREG", "M2_LIKE", "NEUTROPHIL", "CAF", "ENDOTHELIAL"]:
                vals = []
                for c, df in dfs.items():
                    if g in df:
                        msi = ["MSI_H"] if (args.msi and df["MSI_H"].notna().mean() >= 0.7 and df["MSI_H"].nunique() > 1) else []
                        gf = gene_frame(df, g, ["LEUKOCYTE", "CAF", s_] + msi)
                        if s_ in gf:
                            vals.append(partial_spearman(gf, g, s_, [c_ for c_ in ["LEUKOCYTE", "CAF"] + msi if c_ != s_])[0])
                r[f"pcorr_{s_}_median"] = float(np.nanmedian(vals)) if vals else np.nan
            ann.append(r)
        mt = mt.merge(pd.DataFrame(ann), on="gene", how="left")
        mt.to_csv(os.path.join(args.out, f"scan_meta_{ep}.csv"), index=False, float_format="%.5g")

        print(f"\n===== {ep}: genome-wide CTL x gene scan, pooled across {len(lst)} cohorts ({len(mt)} genes) =====")
        print(f"Calibration: median lambda across cohorts = {np.median([i['lambda_inflation'] for i in infos if i['endpoint']==ep and 'lambda_inflation' in i]):.2f} (1.0 = null)")
        print(f"Genes with pooled FDR<0.1: dysfunction-like (HR>1) {int(((mt.q_random<0.1)&(mt.HR_random>1)).sum())}, "
              f"opposite (HR<1) {int(((mt.q_random<0.1)&(mt.HR_random<1)).sum())}")
        print("\nPositive and negative controls (pooled):")
        for g in SCAN_CONTROLS_POS + SCAN_CONTROLS_NEG:
            x = mt[mt.gene == g]
            if len(x):
                r = x.iloc[0]
                print(f"  {g:9s} {'pos' if g in SCAN_CONTROLS_POS else 'neg'}  HR={r.HR_random:.2f} p={r.p_random:.3g} q={r.q_random:.2g} "
                      f"rank(HR>1)={'' if pd.isna(r.rank_dysfunction) else int(r.rank_dysfunction)}/{len(up)}")
        pos = mt[mt.gene.isin(SCAN_CONTROLS_POS)]["p_random"]
        z = lambda df_: np.log(df_.HR_random) / ((np.log(df_.CI_high) - np.log(df_.CI_low)) / 3.92)
        sc = z(mt).rename("z")
        zpos, zneg = sc[mt.gene.isin(SCAN_CONTROLS_POS)], sc[mt.gene.isin(SCAN_CONTROLS_NEG)]
        if len(zpos) and len(zneg):
            auc = np.mean([(a > b) + 0.5 * (a == b) for a in zpos for b in zneg])
            print(f"Control AUC (positives vs negatives, pooled z): {auc:.2f} (n={len(zpos)} vs {len(zneg)})")
        cols = ["gene", "HR_random", "p_random", "q_random", "k", "n_same_direction", "I2", "pcorr_CTL_median",
                "pcorr_EXHAUSTION_median", "pcorr_CAF_median", "pcorr_M2_LIKE_median"]
        print(f"\nTop dysfunction-like genes (pooled HR>1):")
        print(up.merge(mt[["gene"] + [c for c in cols if c not in up.columns and c != 'gene']], on="gene", how="left")
              .head(30)[[c for c in cols if c in mt.columns or c in up.columns]].round(3).to_string(index=False)
              if False else mt[mt["HR_random"] > 1].head(30)[cols].round(3).to_string(index=False))
        print(f"\nTop opposite genes (HR<1, enhance CTL benefit):")
        print(mt[mt["HR_random"] < 1].head(15)[cols].round(3).to_string(index=False))
    print(f"\nOutputs in {args.out}/ (scan_*.csv)")


# ============================================================================= self-description
SCHEMA = {
    "msi_cohort_selection": {"default": "pooled", "options": ["pooled", "mss", "msi-h"],
        "MSS": "Includes MSI-L unless --msi-low exclude. Unknown labels excluded from both subgroups.",
        "missing": "Pooled retains unknown labels unless --msi-missing exclude/error.",
        "adjustment": "--msi controls adjustment independently; check msi_adjusted and actual model covariates.",
        "audit": "cohort_selection.json; per-model n, events, parameters and status; sparse models skipped"},
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
            "msi_cohort": args.msi_cohort, "msi_low": args.msi_low, "msi_missing": args.msi_missing,
            "cohort_audit_file": "cohort_selection.json",
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
    g.add_argument("--msi-cohort", choices=["pooled", "mss", "msi-h"], default="pooled",
                   help="explicit cohort selection before scores; MSS includes MSI-L by default")
    g.add_argument("--msi-low", choices=["include", "exclude"], default="include",
                   help="whether MSI-L belongs in the MSS-only cohort")
    g.add_argument("--msi-missing", choices=["keep", "exclude", "error"], default="keep",
                   help="unknown MSI: pooled keeps by default; subgroups always exclude")
    g.add_argument("--msi", action="store_true", help="adjust for MSI-H where available")
    g.add_argument("--pcorr-adjust", default="basic",
                   help="'basic', 'full' or a comma list of scores for pcorr_* columns "
                        "(pcorr_full_* is always reported)")
    g.add_argument("--decompose", nargs="*", default=None, metavar="SCORE",
                   help="write decomposition.csv for these scores (flag alone = CTL IFNG MACROPHAGE M2_LIKE NEUTROPHIL)")
    g = ap.add_argument_group("survival")
    g.add_argument("--min-model-samples", type=int, default=40)
    g.add_argument("--min-model-events", type=int, default=20)
    g.add_argument("--events-per-parameter", type=float, default=10)
    g.add_argument("--no-survival", action="store_true", help="correlations only")
    g.add_argument("--endpoints", nargs="+", default=["OS", "PFI"])
    g.add_argument("--interaction-score", default="CTL",
                   help="score whose survival effect is tested for modification by the gene")
    g.add_argument("--adjust-scores", nargs="*", default=[],
                   help="extra Cox covariates, e.g. M2_LIKE TREG (sensitivity; may be mediators)")
    g.add_argument("--score-prognosis", action="store_true",
                   help="also test each signature score as an independent prognostic factor")
    g.add_argument("--km", nargs="*", default=None, metavar="GENE",
                   help="Kaplan-Meier panels for these genes (flag alone = all genes)")
    g = ap.add_argument_group("external scores")
    g.add_argument("--score-file", help="table of per-sample scores (e.g. EcoTyper ecotype abundances); "
                                        "columns become usable in --interaction-score, --adjust-scores, --score-prognosis")
    g.add_argument("--score-transpose", action="store_true", help="score file has scores as rows and samples as columns")
    g = ap.add_argument_group("genome-wide scan")
    g.add_argument("--scan", action="store_true",
                   help="test every expressed gene for the CTL x gene interaction; ignores --genes")
    g.add_argument("--scan-endpoints", nargs="+", default=["OS"], help="endpoints to scan (default OS)")
    g.add_argument("--scan-min-expr", type=float, default=1.0, help="minimum median log2 expression (default 1.0)")
    g.add_argument("--scan-top", type=int, default=100, help="hits annotated with partial correlations (default 100)")
    g.add_argument("--scan-workers", type=int, default=0, help="parallel processes (default: cores-2, max 8)")
    g.add_argument("--describe", action="store_true",
                   help="print the JSON schema of methods, columns, files and interpretation rules, then exit")
    args = ap.parse_args()
    if args.describe:
        describe()
        sys.exit(0)
    if args.msi_cohort != "pooled" and args.source == "tcga" and any(c not in ["COADREAD", "COAD", "READ"] for c in args.cohorts):
        ap.error("MSI subgroup selection is restricted to colorectal TCGA cohorts")
    if args.scan and (args.msi_cohort != "pooled" or args.msi_missing != "keep"):
        ap.error("MSI cohort/missing selection is supported in the targeted pipeline, not --scan")
    if args.min_model_samples < 4 or args.min_model_events < 1 or args.events_per_parameter <= 0:
        ap.error("model guard thresholds must be positive (samples >= 4)")
    if args.no_survival:
        args.endpoints = []
    return args


SCORE_DF = None


def main():
    global SCORE_DF
    args = parse_args()
    if args.score_file:
        SCORE_DF = load_score_file(args.score_file, args.score_transpose)
        print(f"Loaded {SCORE_DF.shape[1]} external scores for {len(SCORE_DF)} patients: {list(SCORE_DF.columns)[:12]}")
    known = set(SIGS) | (set(SCORE_DF.columns) if SCORE_DF is not None else set())
    for name in [args.interaction_score] + list(args.adjust_scores):
        if name not in known:
            sys.exit(f"Unknown score '{name}'. Built-in: {sorted(SIGS)}; supply more with --score-file.")
    if args.scan:
        run_scan(args)
        return
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

    import json
    audit = {}
    for name, df in list(datasets.items()):
        selected, report = select_msi_cohort(df, args.msi_cohort, args.msi_low, args.msi_missing)
        report["survival"] = endpoint_counts(selected, args.endpoints)
        audit[name] = report
        datasets[name] = selected
    # Write even if all subgroup analyses are empty or no model can fit.
    with open(os.path.join(args.out, "cohort_selection.json"), "w") as f:
        json.dump(audit, f, indent=2)
    all_res, prog_rows = [], []
    for name, df in datasets.items():
        df = attach_scores(add_scores(df), SCORE_DF)
        missing = [g for g in genes if g not in df]
        print(f"[{name}] {len(df)} tumors" + (f"; genes not found: {missing}" if missing else ""), flush=True)
        res, has_msi = analyze_cohort(name, df, genes, args, decomp_scores, decomp_rows)
        if not res.empty:
            all_res.append(res)
        if args.score_prognosis:
            prog_rows += score_prognosis(name, df, args.endpoints, has_msi,
                                         list(SCORE_DF.columns) if SCORE_DF is not None else ())
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

