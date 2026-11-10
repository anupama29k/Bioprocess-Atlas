"""
Per-run phenotypes and per-perturbation effects.

Phenotypes are computed only from the phase where runs in a campaign share
the same conditions (the batch phase), unless marked otherwise, so a
difference between conditions is attributable to the perturbation and not to
a different feed. Each phenotype records how it was obtained.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import linregress

from .engine import compute_growth_kinetics, ENGINE_SOURCE

DEPLETED_GL = 0.10      # glucose below this = depleted

PHENOTYPES = {
    # name: (unit, description, confounded_by_feed)
    "mu_batch_h":            ("h⁻¹", "specific growth rate in the batch phase (until glucose depletion; for cell counts, samples below 80 % of peak VCD)", False),
    "mu_batch_se_h":         ("h⁻¹", "standard error of mu_batch_h", False),
    "glucose_depletion_h":   ("h", "first sample with glucose < 0.1 g/L (resolution = sampling interval)", False),
    "acetate_peak_gL":       ("g/L", "highest acetate up to 0.5 h after glucose depletion (overflow)", False),
    "acetate_reconsumed_frac": ("–", "fraction of peak acetate consumed 0.75–2.5 h after glucose depletion", False),
    "Yxs_batch_gg":          ("g/g", "biomass yield on glucose in the batch phase (only if glucose falls monotonically)", False),
    "vcd_max_1e6_mL":        ("10⁶ cells/mL", "peak viable cell density", False),
    "viability_end_pct":     ("%", "viability at the last sample", False),
    "ethanol_end_pct_ww":    ("% w/w", "alcohol at the last sample", False),
    "extract_consumed_frac": ("–", "fraction of the initial real extract consumed by the last sample", False),
    "DO_min_pct":            ("% sat.", "lowest dissolved oxygen (< 5 % = oxygen-limited)", False),
    "lactate_peak_gL":       ("g/L", "highest lactate in the run", False),
    "ivcd_1e9_cells_h_L":    ("10⁹ cells·h/L", "integral of viable cell density over the run (trapezoid)", False),
    "titer_final_gL":        ("g/L", "product titer at the last sample", False),
    "qp_pg_cell_day":        ("pg/cell/day", "cell-specific productivity: product made ÷ integral of viable cells", False),
    "acetate_end_gL":        ("g/L", "acetate at the last sample (late overflow, unlike the batch-phase peak)", False),
    "glucose_end_gL":        ("g/L", "glucose at the last sample (high = fed faster than consumed)", True),
    "product_per_biomass_gg": ("g/g", "product titer ÷ biomass at the last sample", False),
    "galactosylation_end_idx": ("galactoses/glycan", "mean galactose per N-glycan at the last sample (product quality)", False),
    "afucosylation_end_pct": ("%", "afucosylated N-glycans at the last sample (product quality)", False),
    "g0f_end_pct":           ("%", "G0F/G0F share at the last sample (product quality)", False),
    "activity_final_AU_mL":  ("AU/mL", "product activity at the last sample (assay-specific; compare within a product)", False),
    "ethanol_end_gL":        ("g/L", "ethanol at the last sample", False),
    "arabitol_end_gL":       ("g/L", "arabitol at the last sample", False),
    "viability_80_h":        ("h", "first sample after peak VCD with viability below 80 % (culture decline)", False),
    "diameter_end_um":       ("µm", "mean viable cell diameter at the last sample", False),
    "mu_fedbatch_h":         ("h⁻¹", "mean growth rate after glucose depletion (fed-batch); depends on the feed", True),
    "biomass_final_gL":      ("g/L", "last measured biomass", True),
}


def _series(obs, var):
    s = obs[obs.variable == var].sort_values("time_h")
    return s.time_h.to_numpy(float), s.value.to_numpy(float)


def run_phenotypes(run: dict, obs: pd.DataFrame) -> dict:
    """obs = observations of this run only."""
    out = {"run_id": run["run_id"]}
    tx, x = _series(obs, "biomass_gL")
    if len(tx) == 0:
        tx, od = _series(obs, "od600")
        x = od * (run.get("dcw_per_od") or 0.4)
    cell_counts = False
    if len(tx) == 0:   # cell-count based cultures (yeast, mammalian): μ from ln(VCD)
        tx, x = _series(obs, "vcd_cells_mL")
        x = x / 1e6
        cell_counts = True
    tg, g = _series(obs, "glucose_gL")
    ta, a = _series(obs, "acetate_gL")
    tv, v = _series(obs, "vcd_cells_mL")
    if len(tv):
        out["vcd_max_1e6_mL"] = float(v.max() / 1e6)
    for var, key in (("viability_percent", "viability_end_pct"), ("ethanol_pct_ww", "ethanol_end_pct_ww")):
        tt, vv = _series(obs, var)
        if len(tt):
            out[key] = float(vv[-1])
    tvi, vi = _series(obs, "viability_percent")
    if len(tv) >= 3 and len(tvi):
        t_pk = tv[np.argmax(v)]
        late = tvi[(tvi >= t_pk) & (vi < 80)]
        if len(late):
            out["viability_80_h"] = float(late[0])
    if len(ta):
        out["acetate_end_gL"] = float(a[-1])
    tl, lac = _series(obs, "lactate_gL")
    if len(tl):
        out["lactate_peak_gL"] = float(lac.max())
    tdm, dm = _series(obs, "cell_diameter_um")
    if len(tdm):
        out["diameter_end_um"] = float(dm[-1])
    for var, key in (("product_activity_AU_mL", "activity_final_AU_mL"), ("ethanol_gL", "ethanol_end_gL"),
                     ("arabitol_gL", "arabitol_end_gL"), ("galactosylation_index", "galactosylation_end_idx"),
                     ("afucosylation_percent", "afucosylation_end_pct"), ("g0f_percent", "g0f_end_pct")):
        tt, vv = _series(obs, var)
        if len(tt):
            out[key] = float(vv[-1])
    if len(tg):
        out["glucose_end_gL"] = float(g[-1])
    tp, prod = _series(obs, "product_gL")
    if len(tp):
        out["titer_final_gL"] = float(prod[-1])
        tb, bm = _series(obs, "biomass_gL")
        if len(tb) and bm[-1] > 0:
            out["product_per_biomass_gg"] = float(prod[-1] / bm[-1])
    if len(tv) >= 3:
        ivcd = np.trapezoid(v, tv) if hasattr(np, "trapezoid") else np.trapz(v, tv)   # cells·h/mL
        out["ivcd_1e9_cells_h_L"] = float(ivcd * 1e3 / 1e9)
        if len(tp) >= 2 and ivcd > 0:
            # ΔP between first and last product sample over the IVCD of the same window
            w = (tv >= tp[0]) & (tv <= tp[-1])
            iv = np.trapezoid(v[w], tv[w]) if hasattr(np, "trapezoid") else np.trapz(v[w], tv[w])
            if iv > 0:
                # g/L ÷ (cells·h/mL) → pg/cell/day : ×1e12 pg/g ÷ 1e3 mL/L × 24 h/day
                out["qp_pg_cell_day"] = float((prod[-1] - prod[0]) / iv * 1e12 / 1e3 * 24)
    te, ex = _series(obs, "extract_plato")
    if len(te) >= 2 and ex[0] > 0:
        out["extract_consumed_frac"] = float((ex[0] - ex[-1]) / ex[0])
    td, do = _series(obs, "DO_percent")
    if len(td):
        out["DO_min_pct"] = float(np.clip(do, 0, None).min())
    if len(tx) < 3:
        return out

    # batch end = glucose depletion, else the whole run is batch-like
    t_dep = None
    if len(tg):
        # first depleted sample that follows a sample with glucose present (> 0.5 g/L);
        # later fed-batch spikes (overfeeding, flagged outliers) must not reset it.
        # If the first sample is already (nearly) depleted, the first depleted time is an upper bound.
        for i in range(len(g)):
            if g[i] < DEPLETED_GL and (i == 0 or g[:i].max() > 0.5 or g[0] < 0.5):
                t_dep = float(tg[i])
                break
    out["glucose_depletion_h"] = t_dep

    in_batch = tx <= (t_dep if t_dep is not None else tx.max())
    if cell_counts:
        # Sparse daily sampling: fit ln(VCD) over the samples taken while the
        # culture is still growing (< 80 % of its peak), from inoculation on.
        grow = np.arange(len(x)) < max(np.argmax(x >= 0.8 * x.max()), 1)
        if grow.sum() >= 3:
            r = linregress(tx[grow], np.log(np.clip(x[grow], 1e-9, None)))
            out["mu_batch_h"], out["mu_batch_se_h"] = float(r.slope), float(r.stderr)
        in_batch = np.zeros_like(in_batch)   # skip the generic fit below
    if in_batch.sum() >= 3:
        if in_batch.sum() >= 5:
            k = compute_growth_kinetics(tx[in_batch], x[in_batch], dcw_factor=1.0)
            mu, se = k.get("mu"), k.get("mu_se")
        else:   # sparse batch phase: plain log-linear fit over it
            r = linregress(tx[in_batch], np.log(np.clip(x[in_batch], 1e-9, None)))
            mu, se = float(r.slope), float(r.stderr)
        out["mu_batch_h"], out["mu_batch_se_h"] = mu, se

    if len(ta):
        lim = (t_dep + 0.5) if t_dep is not None else ta.max()
        pk = a[ta <= lim]
        if len(pk):
            out["acetate_peak_gL"] = float(pk.max())
        if t_dep is not None:
            late = a[(ta >= t_dep + 0.75) & (ta <= t_dep + 2.5)]
            if len(late) and out.get("acetate_peak_gL", 0) > 0.05:
                out["acetate_reconsumed_frac"] = float(
                    np.clip(1 - late.mean() / out["acetate_peak_gL"], -1, 1))

    if len(tg) and t_dep is not None and not cell_counts:   # g/g yield needs biomass, not cell counts
        m = tg <= t_dep
        gx = np.interp(tg[m], tx, x) if m.sum() else []
        if m.sum() >= 3 and np.all(np.diff(g[m]) < 0):
            r = linregress(g[m], gx)
            out["Yxs_batch_gg"] = float(-r.slope)

    if t_dep is not None and run.get("mode") == "Fed-batch":
        fb = tx > t_dep
        if fb.sum() >= 3:
            out["mu_fedbatch_h"] = float(linregress(tx[fb], np.log(x[fb])).slope)
    if len(_series(obs, "biomass_gL")[0]) or len(_series(obs, "od600")[0]):
        out["biomass_final_gL"] = float(x[-1])
    return out


def _compare(a, b):
    """Stats for perturbed (a) vs reference (b); p and d need ≥ 2 per arm."""
    r = {"mean": a.mean(), "sd": a.std(ddof=1) if len(a) > 1 else np.nan, "n": len(a),
         "ref_mean": b.mean(), "ref_sd": b.std(ddof=1) if len(b) > 1 else np.nan, "ref_n": len(b),
         "difference": a.mean() - b.mean(),
         "pct_change": (a.mean() / b.mean() - 1) * 100 if b.mean() else np.nan,
         "effect_size_d": np.nan, "p_welch": np.nan}
    if len(a) >= 2 and len(b) >= 2:
        sd = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
        r["effect_size_d"] = (a.mean() - b.mean()) / sd if sd > 0 else np.nan
        if sd > 0:
            r["p_welch"] = float(stats.ttest_ind(a, b, equal_var=False).pvalue)
    return r


def effects(runs: pd.DataFrame, pheno: pd.DataFrame, contrasts=None) -> pd.DataFrame:
    """Each perturbed condition vs its reference: from the runs table
    (condition vs reference_condition within a campaign) and from explicit
    contrasts. Phenotypes missing in either arm are skipped."""
    df = runs[["run_id", "campaign", "condition", "reference_condition",
               "perturbation_type", "perturbation"]].merge(pheno, on="run_id")
    names = [p for p in PHENOTYPES if not p.endswith("_se_h")]
    rows = []

    def emit(base, a_df, b_df, only=None, exclude=()):
        for p in names:
            if (only is not None and p not in only) or p in exclude:
                continue
            if p not in a_df or p not in b_df:
                continue
            a, b = a_df[p].dropna().to_numpy(float), b_df[p].dropna().to_numpy(float)
            if len(a) == 0 or len(b) == 0:
                continue
            rows.append({**base, "phenotype": p, "unit": PHENOTYPES[p][0], **_compare(a, b),
                         "confounded_by_feed": PHENOTYPES[p][2] and base["perturbation_type"] != "feed"})

    for (camp, cond), g in df.groupby(["campaign", "condition"]):
        refc = g.reference_condition.iloc[0]
        if cond == refc:
            continue
        ref = df[(df.campaign == camp) & (df.condition == refc)]
        if ref.empty:
            continue
        emit({"campaign": camp, "condition": cond, "reference": refc,
              "perturbation_type": g.perturbation_type.iloc[0],
              "perturbation": g.perturbation.iloc[0] or cond, "dose": np.nan}, g, ref)
    for c in contrasts or []:
        a_df, b_df = df[df.run_id.isin(c["runs"])], df[df.run_id.isin(c["reference_runs"])]
        if a_df.empty or b_df.empty:
            continue
        emit({"campaign": c["campaign"], "condition": a_df.condition.iloc[0],
              "reference": b_df.condition.iloc[0], "perturbation_type": c["perturbation_type"],
              "perturbation": c["perturbation"], "dose": c.get("dose", np.nan)}, a_df, b_df,
             only=c.get("phenotypes"), exclude=c.get("exclude_phenotypes", ()))
    return pd.DataFrame(rows)
