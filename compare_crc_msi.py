"""Predeclared, paired-bootstrap MSS versus pooled MSI-adjusted comparison.

Run from the repository root; see docs/CRC_MSI_ANALYSIS.md. No TCGA results bundled.
"""
import argparse
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import immune_gene_function as igf


def fdr_families(table, families):
    """BH over all attempted tests in each declared family (failed tests count as p=1)."""
    table = table.copy()
    table["q"] = np.nan
    if not table.empty:
        for _, idx in table.groupby(families, dropna=False).groups.items():
            p = table.loc[idx, "p"]
            q = igf.bh(p.fillna(1))
            table.loc[idx, "q"] = np.where(p.notna(), q, np.nan)
    return table


def bootstrap_difference(frame, estimator, repetitions, rng, minimum=100):
    """Resample the same patients for both fits, preserving the nested cohort overlap.

    Percentile intervals for descriptive differences only: no independent-group p test.
    MSI-stratified resampling preserves observed cohort composition.
    """
    values = []
    groups = [np.flatnonzero(frame.MSI_H.eq(v).to_numpy()) for v in [0, 1]]
    for _ in range(repetitions):
        ix = np.concatenate([rng.choice(g, size=len(g), replace=True) for g in groups if len(g)])
        sample = frame.iloc[ix].reset_index(drop=True)
        pair = estimator(sample)
        if np.isfinite(pair).all():
            values.append(pair[0] - pair[1])  # MSS minus pooled
    out = dict(bootstrap_requested=repetitions, bootstrap_valid=len(values),
               delta_ci_low=np.nan, delta_ci_high=np.nan,
               bootstrap_status="insufficient_valid_replicates")
    if len(values) >= max(minimum, int(np.ceil(0.8*repetitions))):
        out.update(delta_ci_low=float(np.quantile(values, .025)),
                   delta_ci_high=float(np.quantile(values, .975)), bootstrap_status="ok")
    return out


def subgroup(frame, low):
    return frame[frame.MSI_status.isin(["MSS", "MSI-L"] if low == "include" else ["MSS"])]


def correlation_pair(frame, gene, score, covars, low):
    a = igf.partial_spearman_details(subgroup(frame, low), gene, score, covars)
    b = igf.partial_spearman_details(frame, gene, score, covars + ["MSI_H"])
    return a, b


def cox_pair(frame, gene, ep, low, guards):
    fixed = [c for c in ["stage", "age", "CAF"] if c in frame]
    a = igf.cox_interaction(subgroup(frame, low), gene, ep+".time", ep, False, fixed_covariates=fixed, **guards)
    b = igf.cox_interaction(frame, gene, ep+".time", ep, True, fixed_covariates=fixed, **guards)
    return a, b


def compare(raw, config, repetitions=None, no_survival=False):
    genes = config["genes"]
    low = config.get("msi_low", "include")
    df = igf.add_scores(raw.copy())  # identical pooled reference signature scaling
    df, audit = igf.select_msi_cohort(df, "pooled", low, "exclude")
    if not df.index.is_unique:
        raise ValueError("Comparison requires one independent tumor per patient and unique sample IDs")
    if df.MSI_H.nunique() != 2:
        raise ValueError("Pooled MSI-adjusted comparison requires both MSI-H and MSS/MSI-L annotations")
    mss, mss_audit = igf.select_msi_cohort(df, "mss", low)
    _, high_audit = igf.select_msi_cohort(df, "msi-h", low)
    endpoints = [] if no_survival else config["endpoints"]
    audit["survival"] = igf.endpoint_counts(df, config["endpoints"])
    mss_audit["survival"] = igf.endpoint_counts(mss, config["endpoints"])
    high_audit["survival"] = igf.endpoint_counts(df[df.MSI_H.eq(1)], config["endpoints"])
    rng = np.random.default_rng(config["seed"])
    reps = config["bootstrap_repetitions"] if repetitions is None else repetitions
    rows, deltas = [], []
    guards = config["cox_guards"]
    for gene in genes:
        if gene not in df:
            # Retain attempted hypotheses, rather than silently shrinking the FDR family.
            for adjustment in config["adjustments"]:
                for score in config["signatures"]:
                    for cohort in ["mss", "pooled_msi_adjusted"]:
                        rows.append(dict(gene=gene, target=score, family="pcorr:"+adjustment,
                                         cohort=cohort, p=np.nan, estimate=np.nan, status="gene_missing"))
            for ep in endpoints:
                for cohort in ["mss", "pooled_msi_adjusted"]:
                    rows.append(dict(gene=gene, target=ep, family="cox:"+ep, cohort=cohort,
                                     p=np.nan, estimate=np.nan, status="gene_missing"))
            continue
        for adjustment in config["adjustments"]:
            cov = igf.ADJUST_SETS[adjustment]
            for score in config["signatures"]:
                covars = [c for c in cov if c != score]
                gf = igf.gene_frame(df, gene, list(dict.fromkeys([score]+covars+["MSI_H"])))
                gf["MSI_status"] = df.MSI_status
                pair = correlation_pair(gf, gene, score, covars, low)
                for cohort, detail in zip(["mss", "pooled_msi_adjusted"], pair):
                    detail = detail.copy()
                    rows.append(dict(gene=gene, target=score, family="pcorr:"+adjustment,
                                     cohort=cohort, estimate=detail.pop("r"), **detail))
                delta = pair[0]["r"]-pair[1]["r"]
                boot = bootstrap_difference(gf, lambda f: np.array([
                    d["r"] for d in correlation_pair(f, gene, score, covars, low)]), reps, rng) if np.isfinite(delta) else {}
                deltas.append(dict(gene=gene, target=score, family="pcorr:"+adjustment,
                                   delta=delta, units="partial_Spearman_r", **boot))
        for ep in endpoints:
            if not igf.has_endpoint(df, ep):
                for cohort in ["mss", "pooled_msi_adjusted"]:
                    rows.append(dict(gene=gene, target=ep, family="cox:"+ep, cohort=cohort,
                                     p=np.nan, estimate=np.nan, status="endpoint_missing"))
                continue
            # Hold the adjustment set fixed across nested models and bootstrap draws.
            cov = igf.usable_covariates(df, ["stage", "age", "CAF"])
            cf = df[[gene, "CTL", "MSI_H", "MSI_status", ep, ep+".time"]+cov].copy()
            # Restrict both models to common observed covariates; unavailable ones recorded in audit.
            cf = cf.replace([np.inf, -np.inf], np.nan).dropna()
            cf = cf[(cf[ep+".time"] > 0) & cf[ep].isin([0, 1])]
            pair = cox_pair(cf, gene, ep, low, guards)
            # Pipeline beta is per within-model SD product. Convert to a common pooled SD
            # scale before comparing so a shift in subgroup variance does not mimic an effect.
            pooled_scale = cf[gene].std()*cf.CTL.std()
            def common_beta(f, details):
                result = []
                for sub, d in zip([subgroup(f, low), f], details):
                    local = sub[gene].std()*sub.CTL.std()
                    result.append(d.get("logHR_CTLxX", np.nan)*pooled_scale/local if local > 0 else np.nan)
                return np.array(result)
            betas = common_beta(cf, pair)
            for cohort, detail, beta in zip(["mss", "pooled_msi_adjusted"], pair, betas):
                sub = subgroup(cf, low) if cohort == "mss" else cf
                factor = pooled_scale/(sub[gene].std()*sub.CTL.std()) if len(sub) > 1 else np.nan
                se = detail.get("se_CTLxX", np.nan)*factor
                rows.append(dict(gene=gene, target=ep, family="cox:"+ep, cohort=cohort,
                                 estimate=beta, ci_low=beta-1.96*se, ci_high=beta+1.96*se,
                                 HR=np.exp(beta), HR_ci_low=np.exp(beta-1.96*se), HR_ci_high=np.exp(beta+1.96*se),
                                 p=detail.get("p_CTLxX", np.nan), n=detail.get("n", 0),
                                 events=detail.get("events", 0), parameters=detail.get("parameters", 0),
                                 status=detail["status"], covariates=detail.get("covariates", "")))
            delta = betas[0]-betas[1]
            boot = bootstrap_difference(cf, lambda f: common_beta(f, cox_pair(f, gene, ep, low, guards)),
                                        reps, rng) if np.isfinite(delta) else {}
            deltas.append(dict(gene=gene, target=ep, family="cox:"+ep, delta=delta,
                               units="logHR_per_pooled_SD_product", **boot))
    table = fdr_families(pd.DataFrame(rows), ["family", "cohort"])
    audit.update(mss=mss_audit, msi_h=high_audit,
                 cox_covariate_policy="Available at >=70% in annotated pool; common complete cases; no imputation",
                 difference_inference="Paired MSI-stratified bootstrap; descriptive percentile CI, no p or q for delta")
    return table, pd.DataFrame(deltas), audit


def file_hash(path):
    with open(path, "rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="analysis/crc_glycosylation.json")
    ap.add_argument("--cache", default="xena_cache")
    ap.add_argument("--out", default="crc_msi_comparison")
    ap.add_argument("--expr", help="optional local log expression; also requires --clinical")
    ap.add_argument("--clinical")
    ap.add_argument("--bootstrap", type=int)
    ap.add_argument("--no-survival", action="store_true")
    args = ap.parse_args()
    if bool(args.expr) != bool(args.clinical):
        ap.error("--expr and --clinical must be supplied together")
    if args.bootstrap is not None and args.bootstrap < 0:
        ap.error("--bootstrap must be >=0")
    config = json.loads(Path(args.config).read_text())
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    metadata = dict(config=config, command=sys.argv, python=platform.python_version(), status="started")
    metadata["versions"] = {name: importlib.metadata.version(name) if importlib.util.find_spec(name) else None
                            for name in ["numpy", "pandas", "scipy", "lifelines"]}
    try:
        metadata["git_commit"] = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        metadata["git_commit"] = None
    metadata["source_sha256"] = {p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
                                 for p in [__file__, igf.__file__, args.config]}
    metadata["input_sha256"] = {}
    try:
        if not args.no_survival and metadata["versions"]["lifelines"] is None:
            raise RuntimeError("Survival analysis requires lifelines; install requirements.txt or use --no-survival")
        if args.expr:
            for p in [args.expr, args.clinical]:
                metadata["input_sha256"][p] = file_hash(p)
            data_args = argparse.Namespace(expr=args.expr, clinical=args.clinical,
                    endpoints=["OS"], time_col="OS.time", event_col="OS")
            loaded = igf.load_custom(data_args, config["genes"])
            if len(loaded) != 1:
                raise ValueError("Comparison expects one CRC cohort per invocation")
            raw = next(iter(loaded.values()))
        else:
            raw = igf.load_tcga(config["cohort"], config["genes"], args.cache)
            for p in Path(args.cache).iterdir():
                if "COADREAD" in p.name and not p.name.endswith(".part"):
                    with p.open("rb") as f:
                        metadata["input_sha256"][str(p)] = hashlib.file_digest(f, "sha256").hexdigest()
        table, differences, audit = compare(raw, config, args.bootstrap, args.no_survival)
        table.to_csv(out/"estimates.csv", index=False)
        differences.to_csv(out/"mss_minus_pooled.csv", index=False)
        (out/"cohort_selection.json").write_text(json.dumps(audit, indent=2))
        metadata.update(status="completed", bootstrap_effective=args.bootstrap if args.bootstrap is not None else config["bootstrap_repetitions"],
                        no_survival=args.no_survival,
                        cautions=["Exploratory prognosis, not immunotherapy response or causation",
                                  "MSI status does not prove pMMR; bulk RNA cannot identify expressing cells",
                                  "Delta CI is descriptive; overlapping cohorts are not independent",
                                  "Bootstrap holds signature construction fixed; no cell-composition identification"])
    except Exception as e:
        metadata.update(status="failed", error=f"{type(e).__name__}: {e}")
        raise
    finally:
        (out/"run_manifest.json").write_text(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
