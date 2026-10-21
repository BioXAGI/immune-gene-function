"""
Link EcoTyper cell states to immune / stromal context in TCGA bulk RNA-seq.

AGENT NOTE: run `python state_context.py describe` for a JSON schema of commands, inputs, outputs,
columns and interpretation rules. Every `run` writes run_summary.json. Survival data are not used.

Commands
  build     Collect EcoTyper per-sample state abundances into one table (+ annotation with ecotype
            membership and top marker genes).
              python state_context.py build --ecotyper-dir <.../Carcinoma_Fractions/Cell_States/discovery>
  run       For each query state, partial Spearman correlation with every target (other cell states
            and/or marker-gene signature scores), per TCGA cohort, pooled across cohorts by
            Fisher-z random-effects meta-analysis with 95% confidence intervals, plus an ecotype
            agreement test (do the query's own-ecotype partner states rank highest?).
              python state_context.py run --query-prefix Epi_ --cohorts COADREAD LUAD PAAD BRCA HNSC BLCA --msi
  describe  Print the JSON schema.

Method (run)
  1. Data per cohort: TCGA RNA-seq (UCSC Xena, log2), one primary tumor per patient, signature scores
     (mean z-score of marker genes; SIGS in this file), state abundances matched by the
     12-character patient barcode.
     Confounders adjusted (default): LEUKOCYTE = PTPRC expression (total immune content) and CAF = mean z of
     FAP, COL1A1, COL1A2, COL3A1, PDGFRB, THBS2 (stromal content); MSI-H with --msi. --adjust also accepts
     age, stage, MSI_H, any signature name, or any column of the score table.
  1b. Hypothesis signatures (--signature NAME=G1,G2,... or --signature-file GMT/TSV) and single genes
     (--genes) are scored like the built-in signatures (mean z-score of the member genes present in the
     cohort) and tested as extra targets. Overlap with the query state's marker genes and with the
     confounder genes is reported, because overlap makes a correlation partly true by construction.
  1c. Precomputed signature scores (--signature-scores, e.g. ssGSEA/GSVA/singscore or published per-sample
     scores) are matched to tumors by patient barcode and tested as targets named NAME@pre. When NAME also
     exists as a gene list (--signature or a built-in signature), --signature-method decides:
       compare (default)  test both the precomputed score and the mean z-score, and report how well the two
                          scoring methods agree (Spearman per cohort) in score_method_comparison.csv
       precomputed        use only the precomputed score
       meanz              use only the mean z-score of the gene list
  2. Partial Spearman: rank query, target and covariates; regress query ranks and target ranks on
     covariate ranks; correlate the residuals. Default covariates: LEUKOCYTE (PTPRC), CAF score,
     plus MSI-H with --msi where available. A covariate equal to the target is dropped.
  3. Pooling: z = atanh(r), SE = 1/sqrt(n - k - 3) (k = number of covariates); DerSimonian-Laird
     random effects; pooled r = tanh(z_RE), 95% CI = tanh(z_RE +/- 1.96 SE); two-sided p;
     Benjamini-Hochberg q within each query; I2 heterogeneity.
  4. Ecotype agreement (state targets only): AUC = P(partner state's pooled r > non-partner's);
     Mann-Whitney p; rank of each partner among all state targets.

Interpretation
  - Abundances are shares within a lineage; same-lineage targets are excluded by default because
    shares of one lineage are compositionally linked.
  - Ecotypes were defined from co-occurrence of these same state abundances in TCGA, so ecotype
    agreement here is a consistency check, not independent validation.
  - With thousands of tumors, small correlations become significant; read the size of r and its CI.
  - Correlation is association, not causation; bulk deconvolution cannot localize expression.

Standalone. Requires numpy, pandas, scipy (no other project files). TCGA data are downloaded from
UCSC Xena on first use into --cache (default xena_cache/).
"""
import argparse
import glob
import json
import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats

import gzip
import urllib.parse
import urllib.request

warnings.filterwarnings("ignore")

# ============================================================================= data and statistics (standalone)
XENA = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
TUMOR_CODES = ["01", "03", "06"]           # primary, blood-derived primary, metastatic (preference order)
MSI_COLUMNS = ["CDE_ID_3226963", "MSI_updated_Oct62011", "msi", "MSI", "msi_status"]
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
ALIASES = {"SH2D1B": ["EAT2"], "CMTM2": ["CKLFSF2"], "HAVCR2": ["TIM3"], "TIGIT": ["VSIG9"]}


def fetch(name, cache):
    """Download a UCSC Xena TCGA file once and cache it."""
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


def read_expr(path, needed):
    """Stream a genes x samples matrix keeping only `needed` genes; returns samples x genes."""
    opener = gzip.open if path.endswith(".gz") else open
    rows = {}
    with opener(path, "rt") as f:
        cols = f.readline().rstrip("\n").split("\t")[1:]
        for line in f:
            g, rest = line.split("\t", 1)
            if g in needed:
                vals = rest.rstrip("\n").split("\t")
                try:
                    rows[g] = np.array(vals, dtype=float)
                except ValueError:
                    rows[g] = pd.to_numeric(pd.Series(vals), errors="coerce").values
    return pd.DataFrame(rows, index=cols)


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
    return pd.Series(np.select([st.str.contains("IV"), st.str.contains("III"), st.str.contains("II"),
                                st.str.contains(r"\bI|STAGE I")], [4, 3, 2, 1], np.nan), index=series.index)


def parse_msi(series):
    s = series.astype(str).str.upper()
    out = pd.Series(np.nan, index=series.index)
    out[s.str.contains("MSI-H|MSI_H|MSIH|^H$|HIGH", regex=True)] = 1.0
    out[s.str.contains("MSS|MSI-L|MSI_L|STABLE|^L$|LOW", regex=True)] = 0.0
    return out


def _tnbc(clin):
    neg = lambda col: clin[col].eq("Negative") if col in clin else pd.Series(False, index=clin.index)
    ihc = "lab_proc_her2_neu_immunohistochemistry_receptor_status"
    her2_neg = neg(ihc) | (~clin.get(ihc, pd.Series(index=clin.index)).eq("Positive")
                           & neg("lab_procedure_her2_neu_in_situ_hybrid_outcome_type"))
    return neg("breast_carcinoma_estrogen_receptor_status") & neg("breast_carcinoma_progesterone_receptor_status") & her2_neg


SUBSETS = {"BRCA_TNBC": ("BRCA", _tnbc),
           "BRCA_BASAL": ("BRCA", lambda clin: clin.get("PAM50Call_RNAseq", pd.Series(index=clin.index)).eq("Basal"))}


def load_tcga(cohort, cache, extra_genes=()):
    """Tumors x (signature genes + extra genes + age, stage, MSI_H) for one TCGA cohort or clinical subset."""
    parent, selector = SUBSETS.get(cohort, (cohort, None))
    needed = {g for v in SIGS.values() for g in v} | {a for v in ALIASES.values() for a in v} | set(extra_genes)
    expr = read_expr(fetch(f"TCGA.{parent}.sampleMap/HiSeqV2.gz", cache), needed)
    for g, alts in ALIASES.items():
        if g not in expr:
            hit = next((a for a in alts if a in expr), None)
            if hit:
                expr = expr.rename(columns={hit: g})
    clin = pd.read_csv(fetch(f"TCGA.{parent}.sampleMap/{parent}_clinicalMatrix", cache), sep="\t", index_col=0, low_memory=False)
    expr = expr.loc[pick_tumor_samples(expr.index)]
    if selector is not None:
        keep = set(clin.index[selector(clin).fillna(False).astype(bool)])
        expr = expr.loc[[i for i in expr.index if i in keep]]
    df = expr.copy()
    age_col = next((c for c in ["age_at_initial_pathologic_diagnosis", "age"] if c in clin), None)
    stage_col = next((c for c in ["pathologic_stage", "clinical_stage"] if c in clin), None)
    msi_col = next((c for c in MSI_COLUMNS if c in clin and clin[c].notna().sum() > 20), None)
    df["age"] = pd.to_numeric(clin[age_col], errors="coerce").reindex(df.index) if age_col else np.nan
    df["stage"] = parse_stage(clin[stage_col]).reindex(df.index) if stage_col else np.nan
    df["MSI_H"] = parse_msi(clin[msi_col]).reindex(df.index) if msi_col else np.nan
    return df


def zscore(s):
    return (s - s.mean()) / s.std()


def add_scores(df):
    """Signature score = mean z-score (across tumors of the cohort) of its marker genes."""
    for name, gl in SIGS.items():
        present = [g for g in gl if g in df and df[g].std() > 0]
        df[name] = pd.concat([zscore(df[g]) for g in present], axis=1).mean(axis=1) if present else np.nan
    return df


def _norm_patient(x):
    x = str(x).strip().strip('"').replace(".", "-")
    return x[:12] if x.upper().startswith("TCGA") else x


def load_score_file(path, transpose=False):
    """Per-sample score table; TCGA ids reduced to patient barcodes, primary tumor (code 01) preferred."""
    t = pd.read_csv(path, sep="," if path.endswith(".csv") else "\t", index_col=0)
    if transpose:
        t = t.T
    t.index = t.index.map(lambda v: str(v).strip().strip('"'))
    code = pd.Series([i.replace(".", "-")[13:15] if i.upper().startswith("TCGA") else "01" for i in t.index], index=t.index)
    t = t.assign(_pat=t.index.map(_norm_patient), _pri=(code != "01").astype(int)).sort_values("_pri")
    t = t.drop_duplicates("_pat").set_index("_pat").drop(columns="_pri")
    return t.apply(pd.to_numeric, errors="coerce")


def attach_scores(df, score_df):
    pat = pd.Series(df.index.map(_norm_patient), index=df.index)
    for c in score_df.columns:
        df[c] = pat.map(score_df[c]).values
    return df


def partial_spearman(df, x, y, covars):
    """Rank-based partial correlation: residualize ranks of x and y on ranks of covars, correlate residuals."""
    d = df[[x, y] + covars].dropna()
    if len(d) < 20:
        return np.nan, np.nan
    r = d.rank()
    C = np.column_stack([np.ones(len(r))] + [r[c] for c in covars])
    res = lambda v: v - C @ np.linalg.lstsq(C, v, rcond=None)[0]
    return stats.pearsonr(res(r[x].values), res(r[y].values))


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


LINEAGE_SHORT = {"B.cells": "B", "CD4.T.cells": "CD4", "CD8.T.cells": "CD8", "Dendritic.cells": "DC",
                 "Endothelial.cells": "Endo", "Epithelial.cells": "Epi", "Fibroblasts": "Fibro", "Mast.cells": "Mast",
                 "Monocytes.and.Macrophages": "MoMac", "NK.cells": "NK", "PCs": "PC", "PMNs": "PMN"}
DEFAULT_SIGNATURES = ["CTL", "IFNG", "EXHAUSTION", "TSTEM", "TREG", "NK", "BCELL", "MONOCYTE", "M1_LIKE", "M2_LIKE",
                      "MACROPHAGE", "NEUTROPHIL", "MDSC_LIKE", "CDC", "CAF", "ENDOTHELIAL", "VEGF_ANGIO", "ANTIGEN_PRES"]
DEFAULT_COHORTS = ["COADREAD", "LUAD", "PAAD", "BRCA", "HNSC", "BLCA"]

SCHEMA = {
    "purpose": "Associate cell states (EcoTyper or any per-sample score table) with immune/stromal context "
               "and with other cell states in TCGA tumors, pooled across cohorts with 95% CIs. Survival-free.",
    "commands": {
        "build": "python state_context.py build --ecotyper-dir DIR [--out-prefix cell_states]",
        "run": "python state_context.py run (--query STATE ... | --query-prefix Epi_) [--cohorts ...] "
               "[--targets states signatures] [--adjust LEUKOCYTE CAF] [--msi] [--scores cell_states.tsv] "
               "[--annotation cell_states_annotation.csv] [--include-same-lineage] [--out state_context]",
        "describe": "python state_context.py describe",
    },
    "inputs": {
        "scores": "TSV, first column sample/patient id (TCGA barcodes in dotted or dashed form accepted), "
                  "one column per state, e.g. CD8_S03 (lineage short name + state)",
        "annotation": "CSV with columns state, lineage, ecotype, top_genes, markers (all marker genes; written by build)",
        "cohorts": "TCGA codes on UCSC Xena (COADREAD, LUAD, PAAD, BRCA, HNSC, BLCA, SKCM, ...) or the subsets "
                   "BRCA_TNBC, BRCA_BASAL; downloaded on first use into xena_cache/",
        "confounders": "default LEUKOCYTE (PTPRC) + CAF score; --msi adds MSI-H; --adjust accepts age, stage, MSI_H, "
                       "signature names or score-table columns",
        "hypothesis_signatures": "--signature NAME=G1,G2,... (repeatable), --signature-file (GMT or name/gene table), "
                                 "--genes G1 G2; scored as the mean z-score of member genes present in each cohort",
        "precomputed_scores": "--signature-scores TABLE (e.g. ssGSEA/GSVA/singscore); targets named NAME@pre; "
                              "--signature-method compare|precomputed|meanz when NAME also has a gene list",
    },
    "outputs": {
        "per_cohort.csv": "query, target, target_type (state, signature, custom, gene), cohort, n, k_covariates, r",
        "signature_coverage.csv": "for hypothesis signatures/genes: genes found per cohort and missing genes; for "
                                  "precomputed scores: tumors with a score",
        "score_method_comparison.csv": "(compare mode) per signature and query: pooled r with mean z-score vs "
                                       "precomputed score, CIs, difference, same_sign, cis_overlap, and method agreement "
                                       "(Spearman between the two scores, median/min/max across cohorts)",
        "pooled.csv": "query, target, target_type, target_ecotype, same_ecotype, k_cohorts, r_pooled, ci_low, ci_high, "
                      "p, q, I2, r_median, r_min, r_max, n_positive, marker_overlap (target genes among the query state's "
                      "marker genes), covariate_overlap (target genes among confounder genes)",
        "ecotype_agreement.csv": "query, query_ecotype, n_partners, n_states, auc, mwu_p, partner_ranks, "
                                 "partners_ci_above_0",
        "run_summary.json": "command, cohorts with n, and for each query: top positive/negative targets with CI and "
                            "the ecotype agreement",
    },
    "interpretation_rules": [
        "r_pooled with a 95% CI excluding 0 and consistent sign across cohorts (n_positive = k or 0) is a robust link.",
        "Treat |r| < 0.1 as negligible even when significant; |r| 0.2-0.3 modest; >= 0.4 strong.",
        "High I2 means cohorts disagree; inspect per_cohort.csv.",
        "Ecotype agreement near AUC 1 is expected by construction for EcoTyper states analysed in TCGA.",
        "Signature targets that share genes with the query state's markers are partly correlated by construction "
        "(see marker_overlap).",
        "A hypothesis signature sharing genes with a confounder (covariate_overlap) is partly adjusted away; drop the "
        "overlapping genes or the confounder and compare.",
        "Check signature_coverage.csv: a signature with few genes found in a cohort is unreliable there.",
        "In compare mode, conclusions are robust when the two scoring methods have overlapping CIs (cis_overlap; "
        "same_sign is uninformative when both r are near 0); "
        "low method agreement (Spearman < 0.7) means the result depends on how the signature is scored.",
    ],
}


# ---------------------------------------------------------------------------------------------- build
def build(args):
    rows, ann = [], []
    eco_file = args.ecotypes or os.path.join(args.ecotyper_dir.replace("Cell_States", "Ecotypes"), "ecotypes.txt")
    eco = pd.read_csv(eco_file, sep="\t", index_col=0) if os.path.exists(eco_file) else None
    for lin, sh in LINEAGE_SHORT.items():
        f = glob.glob(os.path.join(args.ecotyper_dir, lin, "*", "state_abundances.txt"))
        if not f:
            print(f"  no state_abundances for {lin}")
            continue
        a = pd.read_csv(f[0], sep="\t", index_col=0).T
        a.columns = [f"{sh}_{c}" for c in a.columns]
        rows.append(a)
        gi_f = glob.glob(os.path.join(args.ecotyper_dir, lin, "*", "gene_info.txt"))
        gi = pd.read_csv(gi_f[0], sep="\t", index_col=0) if gi_f else None
        for s in a.columns:
            st = s.split("_", 1)[1]
            top = gi[gi["State"] == st].sort_values("MaxFC", ascending=False).index[:8].tolist() if gi is not None else []
            ce = eco[(eco["CellType"] == lin) & (eco["State"] == st)]["Ecotype"].tolist() if eco is not None else []
            allm = gi[gi["State"] == st].index.tolist() if gi is not None else []
            ann.append(dict(state=s, lineage=lin, ecotype=ce[0] if ce else "", top_genes=", ".join(top),
                            markers=", ".join(allm)))
    T = pd.concat(rows, axis=1)
    T.index.name = "sample"
    T.to_csv(f"{args.out_prefix}.tsv", sep="\t", float_format="%.5g")
    pd.DataFrame(ann).to_csv(f"{args.out_prefix}_annotation.csv", index=False)
    print(f"Wrote {args.out_prefix}.tsv ({T.shape[0]} samples x {T.shape[1]} states) and {args.out_prefix}_annotation.csv"
          + ("" if eco is not None else " (ecotypes.txt not found: ecotype column empty)"))


# ---------------------------------------------------------------------------------------------- run
def pool(r, n, k):
    r, n, k = np.asarray(r, float), np.asarray(n, float), np.asarray(k, float)
    ok = ~np.isnan(r) & (n - k - 3 > 0)
    r, n, k = r[ok], n[ok], k[ok]
    if len(r) == 0:
        return {}
    z, se = np.arctanh(np.clip(r, -0.999, 0.999)), 1 / np.sqrt(n - k - 3)
    w = 1 / se ** 2
    zf = np.sum(w * z) / w.sum()
    Q = np.sum(w * (z - zf) ** 2)
    kk = len(z)
    tau2 = max(0.0, (Q - (kk - 1)) / (w.sum() - np.sum(w ** 2) / w.sum())) if kk > 1 else 0.0
    wr = 1 / (se ** 2 + tau2)
    zr, ser = np.sum(wr * z) / wr.sum(), np.sqrt(1 / wr.sum())
    return dict(k_cohorts=kk, r_pooled=np.tanh(zr), ci_low=np.tanh(zr - 1.96 * ser), ci_high=np.tanh(zr + 1.96 * ser),
                p=2 * stats.norm.sf(abs(zr / ser)), I2=max(0.0, (Q - (kk - 1)) / Q) if Q > 0 and kk > 1 else 0.0,
                r_median=float(np.median(r)), r_min=float(r.min()), r_max=float(r.max()), n_positive=int((r > 0).sum()))


CONFOUNDER_GENES = {"LEUKOCYTE": set(SIGS["LEUKOCYTE"]), "CAF": set(SIGS["CAF"])}


def parse_signatures(args):
    """Hypothesis signatures from --signature NAME=G1,G2 and --signature-file (GMT: name, description, genes...;
    or two-column TSV/CSV: name, gene). Returns {name: [genes]}."""
    sigs = {}
    for item in args.signature or []:
        if "=" not in item:
            sys.exit(f"--signature expects NAME=GENE1,GENE2,... (got '{item}')")
        name, genes = item.split("=", 1)
        sigs[name.strip()] = [g.strip() for g in genes.replace(";", ",").split(",") if g.strip()]
    if args.signature_file:
        lines = [l.rstrip("\n") for l in open(args.signature_file, encoding="utf-8") if l.strip() and not l.startswith("#")]
        sep = "\t" if "\t" in lines[0] else ","
        fields = [l.split(sep) for l in lines]
        if all(len(f) == 2 for f in fields):
            for name, g in fields:
                if name.lower() in ("name", "signature", "set") and g.lower() in ("gene", "symbol"):
                    continue
                sigs.setdefault(name.strip(), []).append(g.strip())
        else:
            for f in fields:
                sigs[f[0].strip()] = [g.strip() for g in f[2:] if g.strip()]
    return {k: list(dict.fromkeys(v)) for k, v in sigs.items() if v}


def run(args):
    os.makedirs(args.out, exist_ok=True)
    S = load_score_file(args.scores)
    ann = pd.read_csv(args.annotation).fillna("").set_index("state") if args.annotation and os.path.exists(args.annotation) else None
    queries = list(args.query or []) + ([c for c in S.columns if c.startswith(args.query_prefix)] if args.query_prefix else [])
    queries = list(dict.fromkeys(queries))
    missing = [q for q in queries if q not in S.columns and q not in SIGS]
    if not queries or missing:
        sys.exit(f"No usable --query/--query-prefix (missing: {missing}). Available states: {list(S.columns)[:20]}...")
    lineage = lambda s: s.split("_", 1)[0] if s in S.columns else None
    targets = []
    if "states" in args.targets:
        targets += [("state", c) for c in S.columns]
    if "signatures" in args.targets:
        targets += [("signature", s) for s in (args.signatures or DEFAULT_SIGNATURES)]
    custom = parse_signatures(args)
    single = list(dict.fromkeys(args.genes or []))
    clash = [n for n in list(custom) + single if n in SIGS or n in S.columns or n in ("age", "stage", "MSI_H")]
    if clash:
        sys.exit(f"Names clash with built-in scores, state columns or covariates: {clash}; rename them.")
    targets += [("custom", n) for n in custom] + [("gene", g) for g in single]
    extra_genes = {g for v in custom.values() for g in v} | set(single)
    coverage = []
    pre = load_score_file(args.signature_scores, args.signature_scores_transpose) if args.signature_scores else None
    matched = []
    if pre is not None:
        bad = [c for c in pre.columns if c in S.columns]
        if bad:
            sys.exit(f"--signature-scores columns clash with state columns: {bad}")
        listed = set(custom) | {t for tt, t in targets if tt == "signature"} | set(SIGS)
        matched = [c for c in pre.columns if c in listed]
        if args.signature_method == "meanz":
            pre = pre[[c for c in pre.columns if c not in matched]]
            matched = []
        elif args.signature_method == "precomputed":
            targets = [(tt, t) for tt, t in targets if not (tt in ("signature", "custom") and t in matched)]
        targets += [("precomputed", f"{c}@pre") for c in pre.columns]
        print(f"Precomputed scores: {list(pre.columns)} for {len(pre)} patients; with a matching gene list: {matched} "
              f"(method: {args.signature_method})")
    compare_rows = []

    per = []
    cohort_n = {}
    for c in args.cohorts:
        print(f"[{c}] loading", flush=True)
        df = attach_scores(add_scores(load_tcga(c, args.cache, extra_genes)), S)
        for name, genes in custom.items():
            have = [g for g in genes if g in df and df[g].std() > 0]
            coverage.append(dict(signature=name, cohort=c, genes_found=len(have), genes_total=len(genes),
                                 missing=",".join(g for g in genes if g not in have)))
            df[name] = pd.concat([zscore(df[g]) for g in have], axis=1).mean(axis=1) if have else np.nan
        for g in single:
            if g not in df:
                coverage.append(dict(signature=g, cohort=c, genes_found=0, genes_total=1, missing=g))
        if pre is not None:
            df = attach_scores(df, pre.rename(columns=lambda x: f"{x}@pre"))
            for x in pre.columns:
                coverage.append(dict(signature=f"{x}@pre", cohort=c, genes_found=np.nan, genes_total=np.nan,
                                     missing=f"tumors with score: {int(df[f'{x}@pre'].notna().sum())}/{len(df)}"))
            for x in matched:
                if args.signature_method == "compare" and x in df:
                    d2 = df[[x, f"{x}@pre"]].dropna()
                    if len(d2) >= 20:
                        compare_rows.append(dict(signature=x, cohort=c, n=len(d2),
                                                 spearman_meanz_vs_precomputed=stats.spearmanr(d2[x], d2[f"{x}@pre"])[0]))
        msi = ["MSI_H"] if args.msi and df["MSI_H"].notna().mean() >= 0.7 and df["MSI_H"].nunique() > 1 else []
        cohort_n[c] = int(df[queries[0]].notna().sum()) if queries[0] in df else 0
        for q in queries:
            if q not in df or df[q].notna().sum() < args.min_samples:
                continue
            for ttype, t in targets:
                if t == q or t not in df:
                    continue
                if ttype == "state" and not args.include_same_lineage and lineage(t) == lineage(q):
                    continue
                if df[t].notna().sum() < args.min_samples or not df[t].std() > 0:
                    continue
                cov = [x for x in list(args.adjust) + msi if x != t and x != q and x in df]
                d = df[[q, t] + cov].dropna()
                if len(d) < args.min_samples:
                    continue
                r = partial_spearman(d, q, t, cov)[0]
                per.append(dict(query=q, target=t, target_type=ttype, cohort=c, n=len(d), k_covariates=len(cov),
                                covariates="+".join(cov), r=r))
        print(f"  {c}: {cohort_n[c]} tumors with query values", flush=True)
    P = pd.DataFrame(per)
    if P.empty:
        sys.exit("No correlations computed (check that the score table matches the cohorts).")
    P.to_csv(os.path.join(args.out, "per_cohort.csv"), index=False, float_format="%.4g")

    rows = []
    for (q, t, tt), g in P.groupby(["query", "target", "target_type"]):
        res = pool(g["r"], g["n"], g["k_covariates"])
        if res:
            rows.append(dict(query=q, target=t, target_type=tt, **res))
    R = pd.DataFrame(rows)
    R["q"] = np.nan
    for q, idx in R.groupby("query").groups.items():
        R.loc[idx, "q"] = bh(R.loc[idx, "p"].values)
    eco = (lambda s: ann["ecotype"].get(s, "") if ann is not None and s in ann.index else "")
    R["query_ecotype"] = R["query"].map(lambda s: eco(s) or "")
    R["target_ecotype"] = R["target"].map(lambda s: eco(s) or "")
    R["same_ecotype"] = (R["query_ecotype"] != "") & (R["query_ecotype"] == R["target_ecotype"])
    members = {**{k: set(v) for k, v in SIGS.items()}, **{k: set(v) for k, v in custom.items()}, **{g: {g} for g in single}}
    members.update({f"{k}@pre": members[k] for k in matched if k in members})
    mcol = "markers" if ann is not None and "markers" in ann else "top_genes"
    markers = (lambda q: set(x.strip() for x in str(ann[mcol].get(q, "")).split(",") if x.strip())
               if ann is not None and mcol in ann and q in ann.index else set())
    R["marker_overlap"] = [",".join(sorted(members.get(t, set()) & markers(q))) for q, t in zip(R["query"], R["target"])]
    used = lambda t: set().union(*[CONFOUNDER_GENES.get(c, set()) for c in args.adjust if c != t])
    R["covariate_overlap"] = [",".join(sorted(members.get(t, set()) & used(t))) for t in R["target"]]
    if coverage:
        pd.DataFrame(coverage).to_csv(os.path.join(args.out, "signature_coverage.csv"), index=False)
    if compare_rows:
        C = pd.DataFrame(compare_rows)
        side = []
        for x in matched:
            a = R[R["target"] == x].set_index("query")
            b = R[R["target"] == f"{x}@pre"].set_index("query")
            for q in a.index.intersection(b.index):
                side.append(dict(signature=x, query=q, r_meanz=a.loc[q, "r_pooled"], ci_meanz=f"[{a.loc[q, 'ci_low']:.3f}, {a.loc[q, 'ci_high']:.3f}]",
                                 r_precomputed=b.loc[q, "r_pooled"], ci_precomputed=f"[{b.loc[q, 'ci_low']:.3f}, {b.loc[q, 'ci_high']:.3f}]",
                                 difference=b.loc[q, "r_pooled"] - a.loc[q, "r_pooled"],
                                 same_sign=bool(np.sign(a.loc[q, "r_pooled"]) == np.sign(b.loc[q, "r_pooled"])),
                                 cis_overlap=bool(max(a.loc[q, "ci_low"], b.loc[q, "ci_low"]) <= min(a.loc[q, "ci_high"], b.loc[q, "ci_high"]))))
        agree = C.groupby("signature")["spearman_meanz_vs_precomputed"].agg(["median", "min", "max"]).add_prefix("method_agreement_")
        out = pd.DataFrame(side).merge(agree, left_on="signature", right_index=True, how="left") if side else agree.reset_index()
        out.to_csv(os.path.join(args.out, "score_method_comparison.csv"), index=False, float_format="%.4g")
        C.to_csv(os.path.join(args.out, "score_method_agreement_per_cohort.csv"), index=False, float_format="%.4g")
        print("\nScoring-method comparison (mean z-score vs precomputed):")
        print(out.round(3).to_string(index=False))
    R = R.sort_values(["query", "r_pooled"], ascending=[True, False])
    R.to_csv(os.path.join(args.out, "pooled.csv"), index=False, float_format="%.4g")

    agree = []
    for q, g in R[R.target_type == "state"].groupby("query"):
        ce = g["query_ecotype"].iloc[0]
        if not ce:
            continue
        same, other = g[g.same_ecotype], g[~g.same_ecotype]
        if len(same) == 0 or len(other) == 0:
            continue
        auc = np.mean([(a > b) + 0.5 * (a == b) for a in same.r_pooled for b in other.r_pooled])
        ranked = g.sort_values("r_pooled", ascending=False).reset_index(drop=True)
        ranks = sorted(int(ranked.index[ranked.target == t][0]) + 1 for t in same.target)
        agree.append(dict(query=q, query_ecotype=ce, n_partners=len(same), n_states=len(g), auc=auc,
                          mwu_p=stats.mannwhitneyu(same.r_pooled, other.r_pooled, alternative="greater").pvalue,
                          partner_ranks=" ".join(map(str, ranks)), partners_ci_above_0=int((same.ci_low > 0).sum())))
    A = pd.DataFrame(agree, columns=["query", "query_ecotype", "n_partners", "n_states", "auc", "mwu_p",
                                     "partner_ranks", "partners_ci_above_0"])
    A.to_csv(os.path.join(args.out, "ecotype_agreement.csv"), index=False, float_format="%.4g")

    summary = {"command": " ".join(sys.argv), "cohorts": cohort_n, "covariates": list(args.adjust) + (["MSI_H"] if args.msi else []),
               "queries": {}}
    fmt = lambda r: {"target": r.target, "type": r.target_type, "r": round(r.r_pooled, 3),
                     "ci": [round(r.ci_low, 3), round(r.ci_high, 3)], "q": float(r.q), "I2": round(r.I2, 2),
                     "cohorts_positive": f"{r.n_positive}/{r.k_cohorts}"}
    for q, g in R.groupby("query"):
        g = g.sort_values("r_pooled", ascending=False)
        summary["queries"][q] = {"ecotype": g["query_ecotype"].iloc[0], "top_positive": [fmt(r) for r in g.head(args.top).itertuples()],
                                 "top_negative": [fmt(r) for r in g.tail(args.top).iloc[::-1].itertuples()],
                                 "ecotype_agreement": A[A["query"] == q].drop(columns="query").to_dict("records")}
    json.dump(summary, open(os.path.join(args.out, "run_summary.json"), "w"), indent=2, default=float)

    pd.set_option("display.width", 220)
    print(f"\nPooled across {len(cohort_n)} cohorts; covariates {summary['covariates']}.")
    for q, g in R.groupby("query"):
        g = g.sort_values("r_pooled", ascending=False)
        print(f"\n== {q} ({g['query_ecotype'].iloc[0] or 'no ecotype'}) ==")
        show = lambda x: ", ".join(f"{r.target} {r.r_pooled:+.2f} [{r.ci_low:+.2f},{r.ci_high:+.2f}]" for r in x.itertuples())
        print(f"  top +: {show(g.head(args.top))}")
        print(f"  top -: {show(g.tail(args.top).iloc[::-1])}")
        a = A[A["query"] == q]
        if len(a):
            a = a.iloc[0]
            print(f"  ecotype agreement: AUC={a.auc:.2f} (p={a.mwu_p:.2g}); partner ranks {a.partner_ranks} of {a.n_states}")
    print(f"\nOutputs in {args.out}/ (pooled.csv, per_cohort.csv, ecotype_agreement.csv, run_summary.json)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build", help="collect EcoTyper state abundances")
    b.add_argument("--ecotyper-dir", required=True, help=".../Carcinoma_Fractions/Cell_States/discovery")
    b.add_argument("--ecotypes", help="path to ecotypes.txt (auto-detected next to Cell_States)")
    b.add_argument("--out-prefix", default="cell_states")
    r = sub.add_parser("run", help="link query states to immune context and other states")
    r.add_argument("--query", nargs="*", help="state columns, e.g. Epi_S01 CD8_S03")
    r.add_argument("--query-prefix", help="all states starting with this prefix, e.g. Epi_")
    r.add_argument("--cohorts", nargs="+", default=DEFAULT_COHORTS)
    r.add_argument("--targets", nargs="+", default=["states", "signatures"], choices=["states", "signatures"])
    r.add_argument("--signatures", nargs="*", help=f"built-in signature scores to use (default: {' '.join(DEFAULT_SIGNATURES)})")
    r.add_argument("--signature", nargs="*", metavar="NAME=G1,G2",
                   help="hypothesis gene signature(s) to test as targets, e.g. TGFB_RESP=SERPINE1,TGFBI,SMAD7,COL4A1")
    r.add_argument("--signature-file", help="GMT (name, description, genes...) or two-column name/gene table")
    r.add_argument("--genes", nargs="*", help="single genes to test as targets")
    r.add_argument("--signature-scores", help="precomputed per-sample signature scores (rows samples/patients, columns signatures)")
    r.add_argument("--signature-scores-transpose", action="store_true", help="precomputed table has signatures as rows")
    r.add_argument("--signature-method", choices=["compare", "precomputed", "meanz"], default="compare",
                   help="when a precomputed score has a matching gene list: test both and compare (default), "
                        "precomputed only, or mean z-score only")
    r.add_argument("--adjust", nargs="*", default=["LEUKOCYTE", "CAF"],
                   help="covariates (default LEUKOCYTE CAF); also: age stage MSI_H, any signature, any score column")
    r.add_argument("--msi", action="store_true", help="also adjust for MSI-H where available")
    r.add_argument("--include-same-lineage", action="store_true", help="also correlate with states of the same lineage")
    r.add_argument("--scores", default="cell_states.tsv")
    r.add_argument("--annotation", default="cell_states_annotation.csv")
    r.add_argument("--cache", default="xena_cache")
    r.add_argument("--min-samples", type=int, default=50)
    r.add_argument("--top", type=int, default=6)
    r.add_argument("--out", default="state_context")
    sub.add_parser("describe", help="print JSON schema")
    args = ap.parse_args()
    if args.cmd == "describe":
        print(json.dumps(SCHEMA, indent=2))
    elif args.cmd == "build":
        build(args)
    else:
        run(args)


if __name__ == "__main__":
    main()
