import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import immune_gene_function as igf
import compare_crc_msi as comparison


def synthetic(n=240):
    rng = np.random.default_rng(8)
    genes = sorted(igf.genes_needed(["CD276", "C1GALT1", "ACTB"]))
    df = pd.DataFrame(rng.normal(size=(n, len(genes))), columns=genes,
                      index=[f"sample{i}" for i in range(n)])
    labels = pd.Series((["MSS"]*120+["MSI-L"]*20+["MSI-H"]*80+["UNKNOWN"]*20)[:n], index=df.index)
    df["MSI_status"] = labels
    df["MSI_H"] = labels.map({"MSS":0., "MSI-L":0., "MSI-H":1.})
    df["OS.time"] = rng.exponential(500, n)+1
    df["OS"] = (rng.random(n)<.65).astype(int)
    df["age"] = rng.normal(60, 10, n)
    df["stage"] = rng.integers(1, 5, n)
    return df


class MSISelection(unittest.TestCase):
    def setUp(self):
        self.df = synthetic()

    def test_parse_known_and_ambiguous(self):
        x = pd.Series(["MSS", "MSI-L", "msi_h", "stable", "HIGH", None, "not MSS", "MSI-H / MSS", "unknown", "0", "1"])
        result = igf.parse_msi(x)
        np.testing.assert_equal(result[:5].to_numpy(), [0.,0.,1.,0.,1.])
        self.assertTrue(result[5:].isna().all())

    def test_pooled_default_preserves_unknown(self):
        d, r = igf.select_msi_cohort(self.df)
        pd.testing.assert_frame_equal(d, self.df)
        self.assertEqual(r["n_unknown_input"], 20)

    def test_mss_includes_low_and_excludes_unknown(self):
        d, r = igf.select_msi_cohort(self.df, "mss")
        self.assertEqual(len(d), 140)
        self.assertEqual(set(d.MSI_status), {"MSS", "MSI-L"})
        self.assertEqual(r["n_excluded"], 100)

    def test_strict_mss(self):
        d, _ = igf.select_msi_cohort(self.df, "mss", "exclude")
        self.assertEqual(len(d), 120)

    def test_high_only(self):
        d, _ = igf.select_msi_cohort(self.df, "msi-h")
        self.assertTrue(d.MSI_H.eq(1).all())
        self.assertEqual(len(d), 80)

    def test_missing_policy(self):
        d, _ = igf.select_msi_cohort(self.df, missing="exclude")
        self.assertEqual(len(d), 220)
        with self.assertRaisesRegex(ValueError, "Missing"):
            igf.select_msi_cohort(self.df, "mss", missing="error")

    def test_no_annotation_does_not_become_mss(self):
        d, r = igf.select_msi_cohort(self.df.drop(columns=["MSI_status", "MSI_H"]), "mss")
        self.assertTrue(d.empty)
        self.assertEqual(r["n_unknown_input"], 240)

    def test_endpoint_counts_valid_only(self):
        d = pd.DataFrame({"OS":[1,0,1,2,1], "OS.time":[1,2,0,3,np.inf]})
        self.assertEqual(igf.endpoint_counts(d,["OS"])["OS"], {"n":2,"events":1})

    def test_tcga_small_annotation_column_is_kept(self):
        ids = ["TCGA-AB-0001-01A", "TCGA-AB-0002-01A"]
        expr = pd.DataFrame({"CD276":[1,2]},index=ids)
        clin = pd.DataFrame({"msi":["MSS","MSI-H"]},index=ids)
        surv = pd.DataFrame({"OS":[1,0],"OS.time":[100,200]},index=ids)
        with patch.object(igf,"read_expr",return_value=expr), patch.object(igf,"fetch",side_effect=lambda name, cache:name), patch.object(igf.pd,"read_csv",side_effect=[surv,clin]):
            loaded = igf.load_tcga("COADREAD", ["CD276"], "unused")
        self.assertEqual(loaded.MSI_status.tolist(),["MSS","MSI-H"])
        self.assertEqual(loaded.attrs["msi_source"],"msi")


class Statistics(unittest.TestCase):
    def test_partial_degrees_of_freedom(self):
        rng = np.random.default_rng(5)
        d = pd.DataFrame(rng.normal(size=(60, 4)), columns=["x","y","c","d"])
        detail = igf.partial_spearman_details(d,"x","y",["c","d"])
        rho = detail["r"]
        expected = 2*stats.t.sf(abs(rho*np.sqrt(56/(1-rho*rho))), 56)
        self.assertAlmostEqual(detail["p"], expected)
        self.assertLess(detail["ci_low"],rho)
        self.assertGreater(detail["ci_high"],rho)

    def test_partial_singular_constant_and_missing(self):
        d = pd.DataFrame({"x":range(50),"y":range(50),"c":range(50),"d":range(50)})
        self.assertEqual(igf.partial_spearman_details(d,"x","y",["c","d"])["status"],"singular_covariates")
        self.assertEqual(igf.partial_spearman_details(d,"x","y",["c"])["status"],"constant_residual")
        self.assertEqual(igf.partial_spearman_details(d,"x","y",["missing"])["n"],0)

    def test_bh_and_failed_family_members(self):
        np.testing.assert_allclose(igf.bh([.01,.04,.03]),[.03,.04,.04])
        d = pd.DataFrame({"family":["a"]*3,"cohort":["mss"]*3,"p":[.01,.03,np.nan]})
        result = comparison.fdr_families(d,["family","cohort"])
        np.testing.assert_allclose(result.q[:2],[.03,.045])
        self.assertTrue(np.isnan(result.q.iloc[2]))

    def test_nested_bootstrap_uses_same_draw(self):
        d = synthetic()
        d = d[d.MSI_H.notna()]
        result = comparison.bootstrap_difference(d,lambda f:np.array([f.CD276.mean(),f.CD276.mean()]),110,np.random.default_rng(8))
        self.assertEqual(result["delta_ci_low"],0)
        self.assertEqual(result["delta_ci_high"],0)
        self.assertEqual(result["bootstrap_valid"],110)

    def test_bootstrap_unreliable_replicates_suppressed(self):
        result = comparison.bootstrap_difference(synthetic(),lambda f:np.array([np.nan,1]),10,np.random.default_rng(1))
        self.assertEqual(result["bootstrap_valid"],0)
        self.assertTrue(np.isnan(result["delta_ci_low"]))

    def test_comparison_end_to_end_and_seed(self):
        cfg = dict(genes=["CD276","NOT_A_GENE"],msi_low="include",endpoints=[],seed=5,
                   bootstrap_repetitions=2,signatures=["CTL"],adjustments=["basic"],cox_guards={})
        a,b,audit = comparison.compare(synthetic(),cfg,no_survival=True)
        self.assertEqual(len(a),4)
        self.assertEqual(audit["n_selected"],220)
        self.assertEqual(audit["mss"]["n_selected"],140)
        self.assertEqual(a[a.cohort.eq("pooled_msi_adjusted") & a.gene.eq("CD276")].covariates.iloc[0], "LEUKOCYTE+CAF+MSI_H")
        aa,bb,_ = comparison.compare(synthetic(),cfg,no_survival=True)
        pd.testing.assert_frame_equal(a,aa)
        pd.testing.assert_frame_equal(b,bb)

    def test_comparison_requires_two_msi_classes(self):
        d=synthetic();d=d[d.MSI_status.eq("MSS")]
        with self.assertRaisesRegex(ValueError,"requires both"):
            comparison.compare(d,dict(genes=["CD276"],endpoints=[]),no_survival=True)


class CoxGuards(unittest.TestCase):
    def setUp(self):
        self.df=igf.add_scores(synthetic())

    def fit(self,d,**kwargs):
        return igf.cox_interaction(d,"CD276","OS.time","OS",False,**kwargs)

    def test_small_sample(self):
        self.assertEqual(self.fit(self.df.iloc[:12])["status"],"insufficient_samples")

    def test_epv(self):
        d=self.df.copy();d["OS"]=0;d.iloc[:25,d.columns.get_loc("OS")]=1
        out=self.fit(d)
        self.assertEqual(out["status"],"insufficient_events")
        self.assertEqual(out["events"],25)
        self.assertNotIn("HR_CTLxX",out)

    def test_constant_predictor(self):
        d=self.df.copy();d["CD276"]=1
        self.assertEqual(self.fit(d)["status"],"constant_predictor")

    def test_collinear_design(self):
        d=self.df.copy();d["age"]=d.CD276
        self.assertEqual(self.fit(d)["status"],"singular_design")

    def test_requested_msi_not_silently_dropped(self):
        d=self.df.copy();d["MSI_H"]=0
        out=igf.cox_interaction(d,"CD276","OS.time","OS",True)
        self.assertEqual(out["status"],"msi_adjustment_unavailable")

    @unittest.skipUnless(importlib.util.find_spec("lifelines"),"lifelines not installed in this execution environment")
    def test_real_cox_fit_uncertainty(self):
        out=self.fit(self.df)
        self.assertEqual(out["status"],"ok")
        self.assertTrue(out["CI_low_CTLxX"] < out["HR_CTLxX"] < out["CI_high_CTLxX"])

    @unittest.skipUnless(importlib.util.find_spec("lifelines"),"lifelines not installed in this execution environment")
    def test_real_cox_comparison(self):
        cfg=dict(genes=["CD276"],endpoints=["OS"],signatures=[],adjustments=[],seed=9,
                 bootstrap_repetitions=0,cox_guards=dict(min_samples=40,min_events=20,events_per_parameter=10))
        table,delta,_=comparison.compare(synthetic(),cfg)
        self.assertEqual(set(table.status),{"ok"})
        self.assertEqual(len(delta),1)


class CLITests(unittest.TestCase):
    def test_custom_mss_correlation_cli(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);df=synthetic()
            genes=[c for c in df if c in igf.genes_needed(["CD276"])]
            df[genes].T.to_csv(root/"expr.tsv",sep="\t")
            df[["OS","OS.time","age","stage","MSI_status"]].rename(columns={"MSI_status":"msi"}).to_csv(root/"clinical.tsv",sep="\t")
            run=subprocess.run([sys.executable,str(Path(igf.__file__)),"--source","custom","--expr",str(root/"expr.tsv"),"--clinical",str(root/"clinical.tsv"),"--genes","CD276","--msi-cohort","mss","--no-survival","--out",str(root/"out")],capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stderr)
            out=pd.read_csv(root/"out/per_cohort_results.csv")
            self.assertEqual(out.n_samples.iloc[0],140)
            self.assertFalse(out.msi_adjusted.iloc[0])
            audit=json.loads((root/"out/cohort_selection.json").read_text())
            self.assertEqual(audit["CUSTOM"]["n_selected"],140)
            self.assertTrue("pcorr_q_CTL" in out)

class CoxOutputRegression(unittest.TestCase):
    def test_successful_summary_emits_ci_and_counts(self):
        class FakeCox:
            def __init__(self, **kwargs):
                pass
            def fit(self, data, time, event):
                cols = [c for c in data if c not in [time,event]]
                self.summary = pd.DataFrame({"coef":.2,"se(coef)":.1,"exp(coef)":np.exp(.2),"p":.045},index=cols)
                return self
        with patch.dict(sys.modules,{"lifelines":SimpleNamespace(CoxPHFitter=FakeCox)}):
            out=igf.cox_interaction(igf.add_scores(synthetic()),"CD276","OS.time","OS",False)
        self.assertEqual(out["status"],"ok")
        self.assertEqual(out["n"],240)
        self.assertAlmostEqual(out["CI_low_CTLxX"],np.exp(.2-1.96*.1))
        self.assertAlmostEqual(out["CI_high_CTLxX"],np.exp(.2+1.96*.1))

    def test_convergence_warning_suppresses_estimate(self):
        import warnings
        class ConvergenceWarning(UserWarning):
            pass
        class FakeCox:
            def __init__(self, **kwargs):
                pass
            def fit(self, *args):
                warnings.warn("unstable fit",ConvergenceWarning)
                return self
        with patch.dict(sys.modules,{"lifelines":SimpleNamespace(CoxPHFitter=FakeCox)}):
            out=igf.cox_interaction(igf.add_scores(synthetic()),"CD276","OS.time","OS",False)
        self.assertEqual(out["status"],"convergence_warning")
        self.assertNotIn("HR_CTLxX",out)


if __name__=="__main__":
    unittest.main()
