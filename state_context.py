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
     (mean z-score of marker genes; see immune_gene_function.SIGS), state abundances matched by the
     12-character patient barcode.
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

Requires: numpy, pandas, scipy; imports immune_gene_function.py from the same folder.
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

import immune_gene_function as igf

warnings.filterwarnings("ignore")

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
        "annotation": "CSV with columns state, lineage, ecotype, top_genes (written by build)",
        "cohorts": "TCGA codes known to immune_gene_function (COADREAD, LUAD, PAAD, BRCA, HNSC, BLCA, SKCM, "
                   "BRCA_TNBC, ...); data are downloaded from UCSC Xena on first use into xena_cache/",
    },
    "outputs": {
        "per_cohort.csv": "query, target, target_type, cohort, n, k_covariates, r",
        "pooled.csv": "query, target, target_type, target_ecotype, same_ecotype, k_cohorts, r_pooled, ci_low, ci_high, "
                      "p, q, I2, r_median, r_min, r_max, n_positive",
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
        "Signature targets that share genes with the query state's markers are partly correlated by construction.",
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
            ann.append(dict(state=s, lineage=lin, ecotype=ce[0] if ce else "", top_genes=", ".join(top)))
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


def run(args):
    os.makedirs(args.out, exist_ok=True)
    S = igf.load_score_file(args.scores)
    ann = pd.read_csv(args.annotation).fillna("").set_index("state") if args.annotation and os.path.exists(args.annotation) else None
    queries = list(args.query or []) + ([c for c in S.columns if c.startswith(args.query_prefix)] if args.query_prefix else [])
    queries = list(dict.fromkeys(queries))
    missing = [q for q in queries if q not in S.columns and q not in igf.SIGS]
    if not queries or missing:
        sys.exit(f"No usable --query/--query-prefix (missing: {missing}). Available states: {list(S.columns)[:20]}...")
    lineage = lambda s: s.split("_", 1)[0] if s in S.columns else None
    targets = []
    if "states" in args.targets:
        targets += [("state", c) for c in S.columns]
    if "signatures" in args.targets:
        targets += [("signature", s) for s in (args.signatures or DEFAULT_SIGNATURES)]

    per = []
    cohort_n = {}
    for c in args.cohorts:
        print(f"[{c}] loading", flush=True)
        df = igf.attach_scores(igf.add_scores(igf.load_tcga(c, ["PTPRC"], args.cache)), S)
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
                r = igf.partial_spearman(d, q, t, cov)[0]
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
        R.loc[idx, "q"] = igf.bh(R.loc[idx, "p"].values)
    eco = (lambda s: ann["ecotype"].get(s, "") if ann is not None and s in ann.index else "")
    R["query_ecotype"] = R["query"].map(lambda s: eco(s) or "")
    R["target_ecotype"] = R["target"].map(lambda s: eco(s) or "")
    R["same_ecotype"] = (R["query_ecotype"] != "") & (R["query_ecotype"] == R["target_ecotype"])
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
    r.add_argument("--signatures", nargs="*", help=f"signature scores to use (default: {' '.join(DEFAULT_SIGNATURES)})")
    r.add_argument("--adjust", nargs="*", default=["LEUKOCYTE", "CAF"], help="covariates (default LEUKOCYTE CAF)")
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
