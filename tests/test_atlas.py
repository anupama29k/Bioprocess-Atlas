"""Tests for the Bioprocess Atlas build.   python tests/test_atlas.py"""
import sys
import pandas as pd, tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from atlas.build import build          # noqa: E402
from atlas.schema import VARIABLES     # noqa: E402

res = []
def check(name, cond, detail=""):
    res.append(bool(cond)); print(f"{'PASS' if cond else 'FAIL'}  {name}  {detail}")

with tempfile.TemporaryDirectory() as d:
    r = build(out=Path(d))
    runs, obs, ph, eff = r["runs"], r["observations"], r["phenotypes"], r["effects"]
    check("every dataset is measured, and every run belongs to a listed dataset",
          r["datasets"].measured.all() and set(runs.dataset_id) == set(r["datasets"].dataset_id),
          f"{len(runs)} runs, {len(r['datasets'])} datasets, {runs.organism.nunique()} organisms")
    check("run ids are unique", runs.run_id.is_unique)
    con = pd.read_csv(Path(d) / "contrasts.csv")
    ids = set(runs.run_id)
    check("explicit contrasts: 11 yeast, 8 CHO-K1, and every listed run exists",
          (con.dataset_id == "uhlendorff2026").sum() == 11 and (con.dataset_id == "uhlendorff2026_cho").sum() == 8
          and all(set(s.split(";")) <= ids for s in pd.concat([con.runs, con.reference_runs])))
    check("contaminated CHO-K1 reactor (EXP011 R02) is in no contrast",
          not con.runs.str.contains("EXP011_R02").any() and not con.reference_runs.str.contains("EXP011_R02").any())
    check("every observation uses the controlled vocabulary", set(obs.variable) <= set(VARIABLES))
    check("every run has a DOI and license", runs.source_doi.str.len().gt(0).all() and runs.license.str.len().gt(0).all())
    implicit = runs[runs.condition != runs.reference_condition].campaign.unique()
    check("one reference condition per campaign (runs-table comparisons)",
          runs[runs.is_reference & runs.campaign.isin(implicit)].groupby("campaign").condition.nunique().eq(1).all())
    k = eff[eff.campaign == "CHO-K1 · Uhlendorff"]
    check("contrasts against the O2-failure reference (EXP005) report only early growth",
          set(k[k.perturbation.isin(["glutamine 6 vs 8 mM", "seeding 6e5 vs 4e5 cells/mL"])].phenotype) == {"mu_batch_h"})
    imp = k[(k.perturbation == "lowered impeller (batch)") & (k.phenotype == "vcd_max_1e6_mL")].iloc[0]
    check("CHO-K1: lowered bottom impeller doubles peak VCD in TCX6D batch",
          imp["pct_change"] > 50 and imp.p_welch < 0.05, f"{imp['pct_change']:+.0f} %, p = {imp.p_welch:.3f}")
    m = runs.merge(ph, on="run_id")
    hans = m[m.dataset_id == "hans2020"]
    check("glucose depletion found in the batch phase for every Hans run",
          hans.glucose_depletion_h.between(0, 3.5).all(), f"max {hans.glucose_depletion_h.max():.2f} h")
    sim = m[m.dataset_id == "simensen2022"].iloc[0]
    check("Simensen μ within 10 % of the paper's 0.71 h⁻¹", abs(sim.mu_batch_h / 0.71 - 1) < 0.10,
          f"{sim.mu_batch_h:.3f}")
    aceA = hans[hans.condition == "BW25113 ΔaceA"].acetate_reconsumed_frac.mean()
    par = hans[hans.condition == "BW25113"].acetate_reconsumed_frac.mean()
    check("ΔaceA (no glyoxylate shunt) barely reconsumes acetate vs parent", aceA < 0.2 and par > 0.5,
          f"{aceA:.2f} vs {par:.2f}")
    e = eff[(eff.perturbation == "ΔaceA") & (eff.phenotype == "acetate_reconsumed_frac")].iloc[0]
    check("ΔaceA acetate effect is large", e.effect_size_d < -3, f"d = {e.effect_size_d:.1f}, p = {e.p_welch:.2f}")
    y = eff[eff.campaign == "uhlendorff2026 · W-68"]
    t37 = y[(y.perturbation.str.startswith("37")) & (y.phenotype == "vcd_max_1e6_mL")].iloc[0]
    check("S. cerevisiae W-68 at 37 °C reaches far fewer cells than same-day 22 °C controls",
          t37["pct_change"] < -50 and t37.p_welch < 0.05, f"{t37['pct_change']:.0f} %, p = {t37.p_welch:.3f}")
    l37 = eff[(eff.campaign == "uhlendorff2026 · W-34/70") & eff.perturbation.str.startswith("37")
              & (eff.phenotype == "vcd_max_1e6_mL")].iloc[0]
    check("lager strain W-34/70 keeps its peak cell density at 37 °C", abs(l37["pct_change"]) < 10,
          f"{l37['pct_change']:+.0f} %")
    pre = eff[(eff.perturbation == "3 h vs 24 h pre-culture") & (eff.phenotype == "vcd_max_1e6_mL")]
    check("short 3 h pre-culture lowers peak VCD in both strains",
          (pre["pct_change"] < -30).all() and (pre.p_welch < 0.01).all(), ", ".join(f"{v:.0f} %" for v in pre["pct_change"]))
    c = eff[eff.campaign.str.startswith("schellenberg2021")].set_index(["perturbation", "phenotype"])
    qp = lambda p: c.loc[(p, "qp_pg_cell_day"), "pct_change"]
    check("CHO: glucose overfeeding raises IgG titer and cell-specific productivity (paper's finding)",
          c.loc[("glucose 12 g/L (overfeeding)", "titer_final_gL"), "pct_change"] > 10 and qp("glucose 12 g/L (overfeeding)") > 20,
          f"titer {c.loc[('glucose 12 g/L (overfeeding)', 'titer_final_gL'), 'pct_change']:+.0f} %, qP {qp('glucose 12 g/L (overfeeding)'):+.0f} %")
    check("CHO: CO2 limitation raises cell-specific productivity", qp("CO2 limitation") > 10, f"qP {qp('CO2 limitation'):+.0f} %")
    big = c.loc[("glutamate 2 g/L (low glutamine)", "diameter_end_um"), "pct_change"]
    check("CHO: glutamate arms make bigger, more productive cells", big > 10 and qp("glutamate 2 g/L (low glutamine)") > 50,
          f"diameter {big:+.0f} %, qP {qp('glutamate 2 g/L (low glutamine)'):+.0f} %")
    check("CHO qP is in the usual range for an IgG producer (5–60 pg/cell/day)",
          ph.qp_pg_cell_day.dropna().between(5, 60).all(), f"{ph.qp_pg_cell_day.min():.0f}–{ph.qp_pg_cell_day.max():.0f}")
    check("single-vessel references get no p-value", c.p_welch.isna().all())
    ph_ = eff[eff.campaign == "K. phaffii · hypoxia"]
    act = ph_[ph_.phenotype == "activity_final_AU_mL"]
    check("K. phaffii: hypoxia raises product activity for every product and growth rate (paper's finding)",
          len(act) == 5 and (act["pct_change"] > 50).all(),
          ", ".join(f"{v:+.0f} %" for v in act["pct_change"]))
    check("K. phaffii: growth-phase phenotypes are not compared in fed-batch-only runs",
          not ph_.phenotype.isin(["mu_batch_h", "glucose_depletion_h", "Yxs_batch_gg"]).any())
    check("the misaligned Pichia run (µ = 0.065 normoxia) is not in the atlas",
          not runs.run_id.str.contains("mu0.065_normoxia").any())
    fb = eff[eff.phenotype == "mu_fedbatch_h"]
    nist = eff[eff.campaign == "schapertoens2026 · FB4"].set_index(["perturbation", "phenotype"])
    gal = [p for p in nist.index.get_level_values(0).unique() if "galactose" in p]
    check("NISTCHO: every galactose-supplemented feed raises galactosylation and lowers G0F (paper's finding)",
          len(gal) == 2
          and all(nist.loc[(p, "galactosylation_end_idx"), "pct_change"] > 50 for p in gal)
          and all(nist.loc[(p, "g0f_end_pct"), "pct_change"] < -20 for p in gal)
          and all(nist.loc[(p, "galactosylation_end_idx"), "p_welch"] < 0.01 for p in gal),
          ", ".join(f"{p}: galactosylation {nist.loc[(p, 'galactosylation_end_idx'), 'pct_change']:+.0f} %" for p in gal))
    check("NISTCHO: galactose on the standard feed costs no titer, while low glucose does",
          abs(nist.loc[("standard feed + galactose", "titer_final_gL"), "pct_change"]) < 15
          and nist.loc[("low glucose", "titer_final_gL"), "pct_change"] < -40,
          f"STD+ {nist.loc[('standard feed + galactose', 'titer_final_gL'), 'pct_change']:+.0f} %, "
          f"LoG {nist.loc[('low glucose', 'titer_final_gL'), 'pct_change']:+.0f} %")
    check("the NISTCHO campaign without a standard-feed arm (FB2) produces no effect rows",
          eff[eff.campaign == "schapertoens2026 · FB2"].empty
          and runs[runs.campaign == "schapertoens2026 · FB2"].run_id.nunique() == 9)
    check("feed-dependent phenotypes are flagged, except where the feed is itself the perturbation",
          fb[fb.perturbation_type != "feed"].confounded_by_feed.all()
          and not fb[fb.perturbation_type == "feed"].confounded_by_feed.any())
    bs = eff[eff.campaign == "B. subtilis · surfactin feed rate"].set_index(["perturbation", "phenotype"])
    ti = {p: bs.loc[(p, "titer_final_gL"), "mean"] for p in bs.index.get_level_values(0).unique() if "vs 0.075" in p}
    best = max(ti, key=ti.get)
    check("B. subtilis: surfactin peaks at µ = 0.25 h⁻¹ at about 36 g/L (paper's finding)",
          "0.25" in best and 34 < ti[best] < 38, f"{best} → {ti[best]:.1f} g/L")
    check("B. subtilis: the highest feed rate accumulates glucose (overfeeding)",
          bs.loc[("feed for µ = 0.4 vs 0.075 h⁻¹", "glucose_end_gL"), "mean"] > 20
          and bs.loc[("feed for µ = 0.3 vs 0.075 h⁻¹", "glucose_end_gL"), "mean"] < 5,
          f"{bs.loc[('feed for µ = 0.4 vs 0.075 h⁻¹', 'glucose_end_gL'), 'mean']:.0f} g/L at µ = 0.4")
    check("B. subtilis: batch-phase growth rate does not track the feed setpoint (negative control)",
          abs(bs.xs("mu_batch_h", level=1)["pct_change"]).median() < 10,
          f"median |change| {abs(bs.xs('mu_batch_h', level=1)['pct_change']).median():.1f} %")

print(f"\n{sum(res)}/{len(res)} passed"); sys.exit(0 if all(res) else 1)
