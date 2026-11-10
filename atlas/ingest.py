"""
Dataset ingesters: each turns one published dataset (as stored, unchanged,
under ``data/raw/``) into atlas Runs + observations. Adding a dataset =
adding one function here and registering it in INGESTERS.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .schema import Dataset, Run

ROOT = Path(__file__).parent.parent / "data" / "raw"


def simensen2022() -> Dataset:
    d = pd.read_csv(ROOT / "simensen2022_ecoli_mg1655_batch" / "sampling_with_offgas.csv")
    og = pd.read_csv(ROOT / "simensen2022_ecoli_mg1655_batch" / "offgas_and_control.csv")
    ds = Dataset(
        dataset_id="simensen2022",
        title="E. coli MG1655 aerobic batch on glucose M9, 1 L",
        citation="Simensen V et al. (2022) PLOS ONE 17:e0262450",
        doi="10.1371/journal.pone.0262450", license="CC BY 4.0", measured=True)
    rid = "simensen2022:batch1"
    ds.runs.append(Run(
        run_id=rid, dataset_id=ds.dataset_id, organism="Escherichia coli", strain="K-12 MG1655",
        mode="Batch", vessel="Eppendorf BioFlo 115, 3 L", working_volume_L=1.0,
        medium="modified M9, 2% glucose", temperature_C=37.0, pH_setpoint=7.0, replicate=1,
        campaign="simensen2022", reference_condition="MG1655", condition="MG1655",
        is_reference=True, dcw_per_od=0.5, source_doi=ds.doi, license=ds.license,
        notes="OD not recorded after 5.25 h; harvest at ~7.5 h"))
    cols = {"OD_600nm": "od600", "Glucose_gL": "glucose_gL", "Acetate_gL": "acetate_gL",
            "Lactate_gL": "lactate_gL", "Succinate_gL": "succinate_gL", "Formate_gL": "formate_gL"}
    for c, v in cols.items():
        for t, x in zip(d.Time_h, d[c]):
            if pd.notna(x):
                ds.observations.append((rid, float(t), v, float(x)))
    og = og[og.Time_h <= 7.35]   # after this the culture is being harvested
    for c, v in {"CER_mmol_Lh": "CER_mmol_Lh", "OUR_mmol_Lh": "OUR_mmol_Lh",
                 "DO_percent": "DO_percent", "pH": "pH", "Agitation_rpm": "agitation_rpm"}.items():
        for t, x in zip(og.Time_h, og[c]):
            if pd.notna(x):
                ds.observations.append((rid, float(t), v, float(x)))
    return ds


def hans2020() -> Dataset:
    m = pd.read_csv(ROOT / "hans2020_ecoli_keio_fedbatch" / "measurements_long.csv")
    ds = Dataset(
        dataset_id="hans2020",
        title="8 E. coli K-12 strains (Keio knockouts) in 24 parallel fed-batch mini-bioreactors",
        citation="Hans S et al. (2020) Bioengineering 7:145",
        doi="10.3390/bioengineering7040145", license="CC BY 4.0", measured=True)
    ref = "BW25113"
    for (strain, rep), g in m.groupby(["strain", "reactor_replicate"]):
        ko = strain.replace("BW25113 ", "") if strain.startswith("BW25113 Δ") else ""
        ptype = "gene_knockout" if ko else ("none" if strain == ref else "strain")
        rid = f"hans2020:{strain.replace(' ', '_')}:r{rep}"
        ds.runs.append(Run(
            run_id=rid, dataset_id=ds.dataset_id, organism="Escherichia coli", strain=f"K-12 {strain}",
            mode="Fed-batch", vessel="2mag bioREACTOR 48 on Tecan/Hamilton liquid handlers",
            working_volume_L=None,
            medium="defined medium; batch on 5 g/L glucose (+5 g/L pulse at 1 h), then 400 g/L "
                   "glucose pulses every 5 min",
            temperature_C=37.0, pH_setpoint=None, replicate=int(rep), campaign="hans2020",
            reference_condition=ref, condition=strain, perturbation_type=ptype,
            perturbation=ko or (strain if strain != ref else ""), is_reference=(strain == ref),
            dcw_per_od=None, source_doi=ds.doi, license=ds.license,
            notes="3000 rpm; feed rate adapted per strain by the model (fed-batch phase is "
                  "confounded by feed, batch phase is not); biomass = OD × 2.62 path × 0.33 g/L"))
        # analytical duplicates at the same time point are averaged
        gg = g.groupby(["variable", "time_h"], as_index=False)["value"].mean()
        for v, t, x in gg[["variable", "time_h", "value"]].itertuples(index=False):
            ds.observations.append((rid, float(t), v, float(x)))
    return ds


INGESTERS = {"simensen2022": simensen2022, "hans2020": hans2020}


def uhlendorff2026() -> Dataset:
    """52 yeast propagation runs (S. cerevisiae W-68, S. pastorianus W-34/70)."""
    D = ROOT / "uhlendorff2026_yeast_propagation"
    cond = pd.read_csv(D / "conditions.csv")
    off = pd.read_csv(D / "offline_all_runs.csv")
    inl = pd.read_csv(D / "inline_15min_medians.csv.gz")
    ds = Dataset(
        dataset_id="uhlendorff2026",
        title="Yeast propagation: 2 strains × temperature, aeration, wort, pre-culture, scale (52 runs)",
        citation="Uhlendorff S et al. (2026) Zenodo dataset; data descriptor on Preprints.org 202611.0120",
        doi="10.5281/zenodo.19690248", license="CC BY 4.0", measured=True)

    def label(c):
        vvm = "? vvm" if pd.isna(c.aeration_vvm) else f"{c.aeration_vvm:g} vvm"
        return (f"{c.strain} · {c.preculture_h:g} h pre · {c.wort_plato:g} °P · {vvm} · "
                f"{c.temperature_C:g} °C · {c.scale_L:g} L")

    cond["condition"] = cond.apply(label, axis=1)
    rid = lambda r: f"uhlendorff2026:{r}"
    for c in cond.itertuples():
        same_exp = cond[(cond.experiment == c.experiment)]
        ref_c = same_exp[same_exp.temperature_C == 22].iloc[0] if (same_exp.temperature_C == 22).any() else c
        is_ref = c.temperature_C == 22 or ref_c is c
        pert = "" if is_ref else f"{c.temperature_C:g} °C"
        ds.runs.append(Run(
            run_id=rid(c.run), dataset_id=ds.dataset_id, organism=c.species, strain=c.strain,
            mode="Batch", vessel="stirred-tank bioreactor (4 parallel at 0.7 L; single at 10 L)",
            working_volume_L=float(c.scale_L), medium=f"brewer's wort, {c.wort_plato:g} °P",
            temperature_C=float(c.temperature_C), pH_setpoint=None,
            replicate=int(c.reactor[1:]), campaign=c.experiment,
            # comparisons for this dataset are the explicit one-factor contrasts below
            reference_condition=c.condition, condition=c.condition,
            perturbation_type="none" if is_ref else "temperature", perturbation=pert,
            perturbation_dose=None if is_ref else float(c.temperature_C), is_reference=bool(is_ref),
            source_doi=ds.doi, license=ds.license,
            notes=f"aeration {'unresolved (1 or 0.5 vvm)' if pd.isna(c.aeration_vvm) else str(c.aeration_vvm) + ' vvm'}; "
                  "pH and DO not controlled; 400 rpm"))

    offmap = {"VCD Mean (cells/mL)": "vcd_cells_mL", "Viability Mean (%)": "viability_percent",
              "Alcohol content (% w/w)": "ethanol_pct_ww"}
    for col, var in offmap.items():
        for r, t, v in off[["run", "Time (h)", col]].dropna().itertuples(index=False):
            ds.observations.append((rid(r), float(t), var, float(v)))
    for r, t, sg in off[["run", "Time (h)", "Real extract (SG)"]].dropna().itertuples(index=False):
        ds.observations.append((rid(r), float(t), "extract_plato", float(259 - 259 / sg)))
    inmap = {"Temperature (°C)": "temperature_C", "pH": "pH", "Dissolved O2 (% air sat)": "DO_percent",
             "Dissolved CO2 (%-sat)": "dCO2_percent", "Permittivity (pF/cm)": "permittivity_pF_cm"}
    for col, var in inmap.items():
        for r, t, v in inl[["run", "Time (h)", col]].dropna().itertuples(index=False):
            ds.observations.append((rid(r), float(t), var, float(v)))

    # One-factor contrasts (same-day controls where the design has them)
    def sel(strain, **kw):
        m = cond.strain == strain
        for k, v in kw.items():
            m &= (cond[k].isin(v) if isinstance(v, (list, tuple)) else (cond[k] == v))
        return [rid(r) for r in cond[m].run]

    for strain, e37, e30, eref24, e3h, epilot in (
            ("W-68", "EXP007", "EXP011", ["EXP001", "EXP002"], "EXP005", "EXP013"),
            ("W-34/70", "EXP008", "EXP012", ["EXP003", "EXP004"], "EXP006", None)):
        add = lambda ptype, pert, dose, runs, ref: ds.contrasts.append(dict(
            campaign=f"uhlendorff2026 · {strain}", perturbation_type=ptype, perturbation=pert,
            dose=dose, runs=runs, reference_runs=ref))
        add("temperature", "37 °C vs 22 °C (same day)", 37,
            sel(strain, experiment=e37, temperature_C=37), sel(strain, experiment=e37, temperature_C=22))
        add("temperature", "30 °C vs 22 °C (same day, 24 °P)", 30,
            sel(strain, experiment=e30, temperature_C=30), sel(strain, experiment=e30, temperature_C=22))
        add("preculture", "3 h vs 24 h pre-culture", 3,
            sel(strain, experiment=e3h), sel(strain, experiment=eref24))
        add("aeration", "0.5 vs 1 vvm", 0.5,
            sel(strain, experiment=e37, temperature_C=22), sel(strain, experiment=e3h))
        add("medium", "24 °P vs 12 °P wort", 24,
            sel(strain, experiment=e30, temperature_C=22), sel(strain, experiment=e37, temperature_C=22))
        if epilot:
            add("scale", "10 L vs 0.7 L", 10,
                sel(strain, experiment=epilot), sel(strain, experiment=e37, temperature_C=22))
    return ds


INGESTERS["uhlendorff2026"] = uhlendorff2026


def schellenberg2021() -> Dataset:
    """CHO DG44 (IgG producer), 12-day fed-batch: nutrient, pH and CO2 stress vs standard process."""
    D = ROOT / "schellenberg2021_cho_dg44_stress"
    ds = Dataset(
        dataset_id="schellenberg2021",
        title="CHO DG44 IgG fed-batch: glutamate, glucose, feed, pH and CO2 stress vs standard process (ambr15 and shake flask)",
        citation="Schellenberg J, Solle D et al. (2021) LUH research data; Eng Life Sci 2022, doi:10.1002/elsc.202100062",
        doi="10.25835/0066265", license="CC0 1.0", measured=True)
    # file prefix → (campaign, condition label, perturbation type, dose, vessel, volume, pH set point)
    AMBR, SF = ("ambr15 microbioreactor", 0.015), ("500 mL shake flask", None)
    design = {
        "nu_reference": ("ambr15 · nutrients", "standard process", "none", None, AMBR, None),
        "nu_1.1": ("ambr15 · nutrients", "glutamate 1.2 g/L (low glutamine)", "medium", 1.2, AMBR, None),
        "nu_1.2": ("ambr15 · nutrients", "glutamate 2 g/L (low glutamine)", "medium", 2.0, AMBR, None),
        "nu_2":   ("ambr15 · nutrients", "glucose 12 g/L (overfeeding)", "feed", 12.0, AMBR, None),
        "nu_3":   ("ambr15 · nutrients", "decreased feed", "feed", None, AMBR, None),
        "ph_reference": ("ambr15 · pH", "standard process", "none", None, AMBR, None),
        "ph_0":   ("ambr15 · pH", "pH 7.3 from day 0", "pH", 7.3, AMBR, 7.3),
        "ph_5":   ("ambr15 · pH", "pH 7.3 from day 5", "pH", 7.3, AMBR, 7.3),
        "sf_reference": ("shake flask", "standard process", "none", None, SF, None),
        "sf_nu":  ("shake flask", "nutrient limitation", "medium", None, SF, None),
        "sf_co2": ("shake flask", "CO2 limitation", "other", None, SF, None),
    }
    cols = {"VCC": ("vcd_cells_mL", 1e6), "Viability": ("viability_percent", 1), "Average diameter": ("cell_diameter_um", 1),
            "IgG": ("product_gL", 1), "Glucose": ("glucose_gL", 1), "Lactate": ("lactate_gL", 1),
            "Glutamine": ("glutamine_gL", 1), "Glutamate": ("glutamate_gL", 1)}
    for f in sorted(D.glob("*.csv")):
        stem = f.stem
        key = next(k for k in sorted(design, key=len, reverse=True) if stem.startswith(k + "_") or stem.startswith(k + "-"))
        camp, label, ptype, dose, (vessel, vol), ph = design[key]
        is_ref = ptype == "none"
        rep = 1 if is_ref else int(stem[len(key) + 1])
        rid = f"schellenberg2021:{key}_{rep}"
        ds.runs.append(Run(
            run_id=rid, dataset_id=ds.dataset_id, organism="Cricetulus griseus (CHO)", strain="CHO DG44, IgG producer",
            mode="Fed-batch", vessel=vessel, working_volume_L=vol, medium="chemically defined, fed-batch (12 d)",
            temperature_C=None, pH_setpoint=ph, replicate=rep, campaign=f"schellenberg2021 · {camp}",
            reference_condition="standard process", condition=label, perturbation_type=ptype,
            perturbation="" if is_ref else label, perturbation_dose=dose, is_reference=is_ref,
            source_doi=ds.doi, license=ds.license,
            notes="reference arm has a single vessel: effects get a difference but no p-value"
                  + ("; glutamine near zero at day 0 in these vessels (as measured)" if key.startswith("nu_1") else "")))
        df = pd.read_csv(f, sep=";", skiprows=[1], na_values=["--", ""], encoding="utf-8")
        df = df.loc[:, ~df.columns.str.startswith("Unnamed")].dropna(subset=["ProcessTime"])
        df["t"] = df.ProcessTime.str.replace("Day", "").astype(float) * 24.0
        for col, (var, k) in cols.items():
            for t, v in df[["t", col]].dropna().itertuples(index=False):
                ds.observations.append((rid, float(t), var, float(v) * k))
    return ds


INGESTERS["schellenberg2021"] = schellenberg2021
