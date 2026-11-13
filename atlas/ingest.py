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



# What a replicate number means, per dataset, taken from each source's own description of its
# hardware (see the dataset READMEs). Almost every multi-vessel study in the atlas runs its
# replicates side by side in one parallel system, so they share the seed culture and the run
# day: the atlas's p-values treat them as independent, which is optimistic. Recorded here
# rather than inferred, so that a dataset whose source says nothing stays "unknown".
REPLICATE_SCOPE = {
    # 2mag bioREACTOR 48, 24 vessels run together on one liquid handler.
    "hans2020": "parallel_vessels",
    # ambr15 microbioreactor system; the shake-flask series is 3 flasks from one inoculum.
    "schellenberg2021": "parallel_vessels",
    # DASGIP parallel system, 2 reactors per experiment.
    "uhlendorff2026_cho": "parallel_vessels",
    # 50 mL TubeSpin tubes, 3-4 per strategy within one campaign.
    "schapertoens2026": "parallel_vessels",
    # 2 x 30 L stirred tanks per feed rate. The paper does not say whether the pair ran
    # together or on different days, so this one is left unknown rather than guessed.
    "hiller2024": "unknown",
}


def _tag_replicate_scope(ds: Dataset) -> Dataset:
    """Fill replicate_scope for runs that an ingester left at the default."""
    for r in ds.runs:
        if r.replicate_scope != "unknown":
            continue                                  # the ingester set it per run
        if sum(1 for x in ds.runs if x.condition == r.condition) == 1:
            r.replicate_scope = "single"
        else:
            r.replicate_scope = REPLICATE_SCOPE.get(ds.dataset_id, "unknown")
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
            # the 0.7 L runs are four vessels side by side; the 10 L runs are one vessel each
            replicate_scope="parallel_vessels" if float(c.scale_L) < 5 else "single",
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


def uhlendorff2026_cho() -> Dataset:
    """24 CHO-K1 batch/fed-batch runs in DASGIP reactors with inline turbidity, permittivity, O2 and CO2."""
    D = ROOT / "uhlendorff2026_cho_k1"
    cond = pd.read_csv(D / "conditions.csv").set_index("experiment")
    inl = pd.read_csv(D / "inline_15min_medians.csv.gz")
    ds = Dataset(
        dataset_id="uhlendorff2026_cho",
        title="CHO-K1 batch and fed-batch: agitation/aeration, seeding, glutamine, medium, impeller, feeding (24 runs)",
        citation="Uhlendorff S et al. (2026) Zenodo dataset; data descriptor on Preprints.org 202611.0377",
        doi="10.5281/zenodo.20829178", license="CC BY 4.0", measured=True)
    excluded = {"EXP011_R02": "bacterial contamination from ~51 h", "EXP010_R02": "feed line blocked until ~66 h"}

    def label(c):
        s = (f"{c.mode} · {c.medium} · {c.agitation_rpm} rpm · {c.aeration_vvm:g} vvm · "
             f"seed {c.seeding_cells_mL / 1e5:g}e5 · Gln {c.glutamine_mM:g} mM")
        return s + (" · lowered impeller" if c.lower_impeller.startswith("lowered") else "")

    rid = lambda exp, r: f"uhlendorff2026_cho:{exp}_{r}"
    offmap = {"VCD Mean (cells/mL)": ("vcd_cells_mL", 1), "Viability Mean (%)": ("viability_percent", 1),
              "Avg. diameter Mean (µm)": ("cell_diameter_um", 1), "Glucose (g/L)": ("glucose_gL", 1),
              "Lactate (g/L)": ("lactate_gL", 1), "Ammonium (mM)": ("ammonium_mmolL", 1),
              "Osmolality (mOsmol/kg)": ("osmolality_mOsm_kg", 1)}
    inmap = {"Dissolved CO2 (%-sat)": "dCO2_percent", "Permittivity (pF/cm)": "permittivity_pF_cm",
             "pH": "pH", "Dissolved O2 (% air sat)": "DO_percent", "Temperature (°C)": "temperature_C",
             "Agitation rate (rpm)": "agitation_rpm"}
    for f in sorted((D / "offline").glob("*_offline.csv")):
        exp, mode, r = f.stem.split("_")[:3]
        c = cond.loc[exp]
        run_id = rid(exp, r)
        key = f"{exp}_{r}"
        ds.runs.append(Run(
            run_id=run_id, dataset_id=ds.dataset_id, organism="Cricetulus griseus (CHO)", strain="CHO-K1 (non-producing)",
            mode=c["mode"], vessel="DASGIP parallel stirred-tank bioreactor, 1 L max",
            working_volume_L=0.65 if c["mode"] == "Batch" else 0.55, medium=c.medium, temperature_C=37.0, pH_setpoint=7.0,
            replicate=int(r[1:]), campaign="uhlendorff2026_cho",
            # comparisons for this dataset are the explicit contrasts below
            reference_condition=label(c), condition=label(c), perturbation_type="none", is_reference=True,
            source_doi=ds.doi, license=ds.license,
            notes="; ".join(x for x in [f"experiment {exp}", str(c.events) if pd.notna(c.events) else "",
                                         f"excluded from contrasts: {excluded[key]}" if key in excluded else ""] if x)))
        off = pd.read_csv(f)
        for col, (var, k) in offmap.items():
            for t, v in off[["Time (h)", col]].dropna().itertuples(index=False):
                ds.observations.append((run_id, float(t), var, float(v) * k))
        sub = inl[inl.run == f"{exp}_{mode}_{r}"]
        for col, var in inmap.items():
            for t, v in sub[["Time (h)", col]].dropna().itertuples(index=False):
                ds.observations.append((run_id, float(t), var, float(v)))

    def runs(*exps, only=None):
        out = [rid(e, r) for e in exps for r in ("R01", "R02")]
        out = [x for x in out if x.split(":")[1] not in excluded]
        return [x for x in out if only is None or x.endswith(only)]

    # Runs stopped at different times, so end-of-run viability is not comparable here.
    # EXP005's O2 supply failed from ~104 h to 121 h (peak VCD at 121.5 h): contrasts that use it
    # as reference keep only phenotypes settled before ~100 h.
    EARLY = ["mu_batch_h"]

    def add(ptype, pert, dose, a, b, early_only=False):
        ds.contrasts.append(dict(
            campaign="CHO-K1 · Uhlendorff", perturbation_type=ptype, perturbation=pert, dose=dose, runs=a,
            reference_runs=b, exclude_phenotypes=["viability_end_pct"], **({"phenotypes": EARLY} if early_only else {})))
    add("other", "100 rpm + 0.04 vvm vs 80 rpm + 0.125 vvm", None, runs("EXP005"), runs("EXP001", "EXP002", "EXP003", "EXP004"), early_only=True)
    add("other", "seeding 6e5 vs 4e5 cells/mL", 6e5, runs("EXP006"), runs("EXP005"), early_only=True)
    add("medium", "glutamine 6 vs 8 mM", 6, runs("EXP007"), runs("EXP005"), early_only=True)
    add("medium", "TCX6D vs TC-42 medium (batch)", None, runs("EXP008"), runs("EXP005"), early_only=True)
    add("other", "lowered impeller (batch)", None, runs("EXP009"), runs("EXP008"))
    add("feed", "fed-batch vs batch (TC-42)", None, runs("EXP010"), runs("EXP005"), early_only=True)
    add("medium", "TCX6D vs TC-42 medium (fed-batch)", None, runs("EXP012"), runs("EXP010"))
    add("other", "lowered impeller (fed-batch)", None, runs("EXP011"), runs("EXP012"))
    return ds


INGESTERS["uhlendorff2026_cho"] = uhlendorff2026_cho


def salesvallverdu2024() -> Dataset:
    """K. phaffii P_GAP producers in glucose-limited fed-batch: hypoxia (RQ-controlled) vs normoxia."""
    D = ROOT / "salesvallverdu2024_pichia_hypoxia"
    L = [l.split("\t") for l in (D / "raw_data_microbial_biotech_v1.tab").read_text().split("\n")]
    ds = Dataset(
        dataset_id="salesvallverdu2024",
        title="Komagataella phaffii fed-batch: hypoxia vs normoxia at four growth rates (Crl1) and for two other proteins",
        citation="Sales-Vallverdú A et al. (2024) Microbial Biotechnology 17:e14411; data on CORA",
        doi="10.34810/data1118", license="CC BY 4.0", measured=True)
    # block start row → (product, set µ in h⁻¹)
    blocks = {0: ("Crl1", 0.030), 11: ("Crl1", 0.065), 21: ("Crl1", 0.100), 30: ("Crl1", 0.120),
              39: ("proROL", None), 48: ("CalB", None)}
    starts = sorted(blocks)
    cols = {3: ("biomass_gL", 1), 4: ("product_activity_AU_mL", 1), 5: ("ethanol_gL", 1),
            6: ("arabitol_gL", 1), 8: ("glucose_gL", 1)}
    # Known problems in the source table (see README): one run is misaligned, one column is a copy.
    unusable = {("Crl1", 0.065, "normoxia")}
    drop_cols = {("Crl1", 0.030, "normoxia"): {5}}

    def num(x):
        try:
            return float(x.strip().strip('"'))
        except ValueError:
            return None

    rid_of = {}
    for s in starts:
        e = next((x for x in starts if x > s), len(L))
        product, mu = blocks[s]
        for off in (0, 11):
            label = L[s][off].strip('"').lower()
            o2 = "hypoxia" if "hypoxia" in label else "normoxia"
            if (product, mu, o2) in unusable:
                continue
            tag = f"{product}_{'mu' + format(mu, '.3f') if mu else 'mu-default'}_{o2}"
            rid = f"salesvallverdu2024:{tag}"
            rid_of[(product, mu, o2)] = rid
            ds.runs.append(Run(
                run_id=rid, dataset_id=ds.dataset_id, organism="Komagataella phaffii", strain=f"P_GAP-{product} producer",
                mode="Fed-batch", vessel="stirred-tank bioreactor (bench)", working_volume_L=None,
                medium="defined medium, glucose-limited fed-batch", temperature_C=None, pH_setpoint=None, replicate=1,
                campaign="salesvallverdu2024", reference_condition=tag, condition=tag, perturbation_type="none",
                is_reference=True, source_doi=ds.doi, license=ds.license,
                notes=(f"set growth rate {mu} h⁻¹; " if mu else "") + f"{o2}; time = hours of fed-batch phase"
                      + ("; ethanol column dropped (copy of the hypoxia run)" if (product, mu, o2) in drop_cols else "")))
            for r in L[s + 2:e]:
                if len(r) < off + 9:
                    continue
                v = [num(x) for x in r[off:off + 9]]
                if v[0] is None or v[0] == 0:          # sample 0 is an all-zero placeholder row
                    continue
                t = v[2]
                if t is None:
                    continue
                for j, (var, k) in cols.items():
                    if j in drop_cols.get((product, mu, o2), ()) or v[j] is None:
                        continue
                    ds.observations.append((rid, float(t), var, v[j] * k))

    for product, mu in [("Crl1", 0.030), ("Crl1", 0.100), ("Crl1", 0.120), ("proROL", None), ("CalB", None)]:
        a, b = rid_of.get((product, mu, "hypoxia")), rid_of.get((product, mu, "normoxia"))
        if a and b:
            ds.contrasts.append(dict(
                campaign="K. phaffii · hypoxia", perturbation_type="DO",
                perturbation=f"hypoxia vs normoxia ({product}" + (f", µ = {mu:g} h⁻¹)" if mu else ")"),
                dose=mu, runs=[a], reference_runs=[b],
                exclude_phenotypes=["glucose_depletion_h", "mu_batch_h", "Yxs_batch_gg"]))
    return ds


INGESTERS["salesvallverdu2024"] = salesvallverdu2024


def hiller2024() -> Dataset:
    """B. subtilis surfactin: six growth-rate-controlling feed rates in 30 L fed-batch, plus a batch process."""
    import openpyxl

    D = ROOT / "hiller2024_bsubtilis_surfactin"
    wb = openpyxl.load_workbook(D / "12934_2024_2531_MOESM2_ESM.xlsx", data_only=True)
    ds = Dataset(
        dataset_id="hiller2024",
        title="B. subtilis BMV9 surfactin: six growth-rate-controlling feed rates (0.075–0.4 h⁻¹) in 30 L fed-batch, duplicate bioreactors",
        citation="Hiller E, Off M, Hermann A, Vahidinasab M, Benatto Perino EH, Lilge L, Hausmann R (2024) Microbial Cell Factories 23:260",
        doi="10.1186/s12934-024-02531-w", license="CC BY 4.0", measured=True)
    cols = {"OD600": ("od600", 1), "CDW (g/L)": ("biomass_gL", 1), "Glucose (g/L)": ("glucose_gL", 1),
            "Surfactin (g/L)": ("product_gL", 1), "Acetate (g/L)": ("acetate_gL", 1),
            "Ammonium (g/L)": ("ammonium_gL", 1)}
    REF_MU = 0.075
    rid_of = {}

    def add_run(rid, mu, rep, feed_start):
        is_ref = (mu == REF_MU)
        cond = "batch (no feed)" if mu is None else f"feed for µ = {mu:g} h⁻¹"
        ds.runs.append(Run(
            run_id=rid, dataset_id=ds.dataset_id, organism="Bacillus subtilis",
            strain="BMV9 (spo0A3 trp+ sfp+ ΔmanPA), surfactin producer", mode="Batch" if mu is None else "Fed-batch",
            vessel="30 L stirred-tank bioreactor, 3 Rushton turbines", working_volume_L=12.0,
            medium="chemically defined mineral salt medium, glucose", temperature_C=37.0, pH_setpoint=7.0,
            replicate=rep, campaign="hiller2024",
            # comparisons for this dataset are the explicit one-factor contrasts below
            reference_condition=cond, condition=cond, perturbation_type="none" if is_ref else ("feed" if mu else "none"),
            perturbation="" if is_ref else cond, perturbation_dose=mu, is_reference=is_ref,
            source_doi=ds.doi, license=ds.license,
            notes="batch volume 12 L in a 30 L vessel, final volume ≈ 18.7 L; DO held above 50 %"
                  + (f"; feed started at {feed_start:g} h" if feed_start is not None else "")))

    for ws in wb.worksheets:
        mu = None if ws.title == "Batch" else float(ws.title.split()[0])
        rows = list(ws.iter_rows(values_only=True))
        h_i = next(i for i, r in enumerate(rows) if r and "Time (h)" in [str(c) for c in r])
        header = rows[h_i]
        # the sheets end with a volume-balance table whose numbers must not be read as samples
        end = next((i for i in range(h_i + 1, len(rows))
                    if any("volume" in str(c).lower() for c in rows[i] if c is not None)), len(rows))
        rows = rows[h_i + 1:end]
        # one block per bioreactor; a block starts at each "Time (h)" column
        blocks = [i for i, c in enumerate(header) if str(c) == "Time (h)"]
        for rep, t_col in enumerate(blocks, start=1):
            colmap = {}
            for j in range(t_col, min(t_col + 8, len(header))):
                name = str(header[j])
                if name in cols:
                    colmap[j] = cols[name]
            rid = f"hiller2024:{'batch' if mu is None else 'mu' + format(mu, '.3f')}_R{rep}"
            obs, feed_start = [], None
            for r in rows:
                if len(r) <= t_col or not isinstance(r[t_col], (int, float)):
                    continue
                t = float(r[t_col])
                if str(r[t_col - 1]).upper().startswith("FS"):
                    feed_start = t
                for j, (var, k) in colmap.items():
                    v = r[j] if j < len(r) else None
                    if isinstance(v, (int, float)):
                        obs.append((rid, t, var, float(v) * k))
            if not obs:
                continue
            add_run(rid, mu, rep, feed_start)
            ds.observations += obs
            rid_of.setdefault(mu, []).append(rid)

    for mu in sorted(m for m in rid_of if m not in (None, REF_MU)):
        ds.contrasts.append(dict(
            campaign="B. subtilis · surfactin feed rate", perturbation_type="feed",
            perturbation=f"feed for µ = {mu:g} vs {REF_MU:g} h⁻¹", dose=mu,
            runs=rid_of[mu], reference_runs=rid_of[REF_MU],
            # the feed starts at glucose depletion, so the reconsumption window lies inside the feed phase
            exclude_phenotypes=["acetate_reconsumed_frac"]))
    if None in rid_of:
        ds.contrasts.append(dict(
            campaign="B. subtilis · surfactin feed rate", perturbation_type="feed",
            perturbation=f"fed-batch (µ = {REF_MU:g}) vs batch", dose=None,
            runs=rid_of[REF_MU], reference_runs=rid_of[None],
            exclude_phenotypes=["acetate_reconsumed_frac"]))
    return ds


INGESTERS["hiller2024"] = hiller2024


def schapertoens2026() -> Dataset:
    """NISTCHO fed-batch under seven feeding strategies, with titer and N-glycan quality attributes."""
    D = ROOT / "schapertoens2026_nistcho_feeding"
    ds = Dataset(
        dataset_id="schapertoens2026",
        title="NISTCHO fed-batch: feeding frequency, glucose level and galactose supplementation, with product N-glycan quality",
        citation="Schäpertöns V et al. (2026) npj Systems Biology and Applications 12:67; data on Zenodo",
        doi="10.5281/zenodo.18609042", license="CC BY 4.0", measured=True)
    # letter code in the raw files → (strategy name, campaign/batch)
    NAME = {"A": ("STD", "FB4"), "B": ("STD+", "FB4"), "C": ("LoG+", "FB4"), "G": ("LoG", "FB4"),
            "D": ("HiF", "FB2"), "E": ("HIP", "FB2"), "F": ("HIP+", "FB2")}
    WHAT = {"STD": "standard feed", "STD+": "standard feed + galactose",
            "LoG": "low glucose", "LoG+": "low glucose + galactose",
            "HiF": "high feed frequency", "HIP": "high initial product feed",
            "HIP+": "high initial product feed + galactose"}
    rid = lambda c, r: f"schapertoens2026:{NAME[c][0]}_R{int(r)}"

    growth = pd.read_csv(D / "01_ViCell_growth_data.csv")
    growth["rep"] = growth.Replicate.str.replace("R", "").astype(int)
    titer = pd.read_csv(D / "03_ViCell_titer.csv").rename(columns={"Titer_ug.mL": "titer_ug_mL"})
    lact = pd.read_csv(D / "04_lactate.csv")
    lact.columns = [c.strip() for c in lact.columns]
    lact["rep"] = lact.Replicate.str.replace("R", "").astype(int)
    gluc = pd.read_excel(D / "05_glucose_data.xlsx")

    for c, (name, batch) in NAME.items():
        reps = sorted(growth[growth.Condition == c].rep.unique())
        for r in reps:
            is_ref = (name == "STD")
            ds.runs.append(Run(
                run_id=rid(c, r), dataset_id=ds.dataset_id, organism="Cricetulus griseus (CHO)",
                strain="NISTCHO, cNISTmAb producer", mode="Fed-batch", vessel="50 mL TubeSpin bioreactor tube",
                working_volume_L=0.020, medium=f"chemically defined fed-batch, {WHAT[name]}",
                temperature_C=None, pH_setpoint=None, replicate=int(r),
                campaign=f"schapertoens2026 · {batch}",
                # FB2 was run as a separate campaign with no standard-feed arm, so its runs
                # carry no reference and produce no effect rows (see the dataset README).
                reference_condition=("STD" if batch == "FB4" else name), condition=name,
                perturbation_type="none" if is_ref else "feed", perturbation="" if is_ref else WHAT[name],
                is_reference=is_ref, source_doi=ds.doi, license=ds.license,
                notes=f"campaign {batch}; letter code {c} in the source files"
                      + ("" if batch == "FB4" else "; no same-campaign standard-feed reference")))

    # Cell counts, viability, diameter and titer all come from the strategy-named file,
    # so they do not depend on the letter -> name mapping at all. (The mapping is verified:
    # every run's time grid and VCD series is identical between the two files.)
    BY_NAME = {n: c for c, (n, _) in NAME.items()}
    for row in titer.itertuples():
        run = rid(BY_NAME[row.Condition], row.Replicate)
        for col, var, k in (("Total_VCD", "vcd_cells_mL", 1e6), ("Viability", "viability_percent", 1),
                            ("Average_diameter", "cell_diameter_um", 1),
                            ("titer_ug_mL", "product_gL", 1e-3)):
            v = getattr(row, col)
            if pd.notna(v):
                ds.observations.append((run, float(row.Hours), var, float(v) * k))
    for row in lact.itertuples():
        ds.observations.append((rid(row.Condition, row.rep), float(row.Hours), "lactate_gL",
                                float(row.c_mM) * 0.09008))       # lactate, 90.08 g/mol
    for row in gluc.itertuples():
        ds.observations.append((rid(row.Condition, row.Replicate), float(row.Hour), "glucose_gL",
                                float(getattr(row, "_4"))))

    # Product N-glycans: one file per condition, replicate and sampling hour, as produced by
    # CAFOG. Each row is a glycoform pair across the two heavy chains, with its corrected
    # fractional abundance and the Fuc / Hex / HexNAc counts of the pair.
    #
    # A biantennary glycan at G0 carries 3 hexoses (core mannose), 4 HexNAc and 1 core fucose
    # when fucosylated, and has 2 galactosylation sites. So for a row covering n glycans
    # (n = HexNAc / 4): galactose residues = Hex - 3n, possible sites = 2n, fucosylation
    # sites = n. Rows where one glycan is absent therefore contribute n = 1 rather than
    # being discarded.
    #
    # The galactosylation index follows Schäpertöns et al. 2026, Eq. 7:
    #     100 * sum(A_j * g_j) / sum(A_j * d_j)
    # which reproduces their reported 24 % for the galactose-supplemented strategies.
    for f in sorted((D / "glyco").glob("*.csv")):
        c, rep, hour = f.stem.split("_")
        g = pd.read_csv(f)
        ab = g.corr_abundance.fillna(g.abundance).astype(float)
        ok = g.HexNAc.notna() & g.Hex.notna()
        if not ok.any() or ab[ok].sum() <= 0:
            continue
        n_gly = g.HexNAc / 4.0                       # glycans represented by this row
        gal_res = g.Hex - 3 * n_gly                  # galactose residues
        gal_sites = 2 * n_gly                        # possible galactosylation sites
        fuc_sites = n_gly                            # one core fucose site per glycan
        run, t = rid(c, rep), float(hour)
        ds.observations.append((run, t, "galactosylation_percent",
                                float(100 * (ab[ok] * gal_res[ok]).sum() / (ab[ok] * gal_sites[ok]).sum())))
        ds.observations.append((run, t, "afucosylation_percent",
                                float(100 * (ab[ok] * (fuc_sites[ok] - g.Fuc[ok].fillna(0))).sum()
                                      / (ab[ok] * fuc_sites[ok]).sum())))
        g0f = ab[g.glycoform == "A2G0F/A2G0F"].sum()
        ds.observations.append((run, t, "g0f_percent", float(100 * g0f / ab.sum())))
    return ds


INGESTERS["schapertoens2026"] = schapertoens2026


def calzadilla2023() -> Dataset:
    """One 2 L HEK293 chemostat held at five dilution-rate set points in sequence."""
    import openpyxl

    D = ROOT / "calzadilla2023_hek293_chemostat"
    ws = openpyxl.load_workbook(D / "Continuous_culture_data_Multiplicity_HEK293.xlsx",
                                data_only=True)["Sheet1"]
    rows = [list(r) for r in ws.iter_rows(values_only=True)]
    ds = Dataset(
        dataset_id="calzadilla2023",
        title="HEK293 producing the Her1 extracellular domain: one 2 L chemostat through five dilution-rate set points",
        citation="Calzadilla N, Boggiano T, Vazquez E, Mulet A (2023) Data in Brief 50:109604",
        doi="10.17632/t9rcjv5362.4", license="CC BY 4.0", measured=True)
    rid = "calzadilla2023:chemostat"
    ds.runs.append(Run(
        run_id=rid, dataset_id=ds.dataset_id, organism="Homo sapiens (HEK293)",
        strain="HEK293 secreting the Her1 extracellular domain", mode="Continuous",
        vessel="2 L stirred-tank bioreactor", working_volume_L=2.0,
        medium="chemically defined, 24 mmol/L glucose and 5.5 mmol/L glutamine in the feed",
        temperature_C=None, pH_setpoint=None, replicate=1, campaign="calzadilla2023",
        # A single vessel: the five steady states are sequential periods of one run, not
        # independent reactors, so no comparison is asserted between them.
        reference_condition="chemostat", condition="chemostat",
        perturbation_type="none", is_reference=True, source_doi=ds.doi, license=ds.license,
        notes="single vessel; dilution rate set to 0.45, 0.40, 0.35, 0.40 then 0.45 1/d "
              "(steady states SS1-SS5), recorded as the dilution_rate_h time series; "
              "the study reports different steady states at the same dilution rate"))
    # column index -> (variable, factor to the atlas unit)
    COLS = {4: ("vcd_cells_mL", 1.0), 12: ("dilution_rate_h", 1.0),
            13: ("product_gL", 1e-3),                 # µg/mL -> g/L
            26: ("glucose_gL", 0.18016),              # mmol/L -> g/L
            30: ("lactate_gL", 0.09008),
            34: ("glutamine_gL", 0.14615)}
    for r in rows[4:]:
        if len(r) < 35 or not isinstance(r[3], (int, float)):
            continue
        t_h = float(r[3])
        for j, (var, k) in COLS.items():
            v = r[j]
            if isinstance(v, (int, float)):
                ds.observations.append((rid, t_h, var, float(v) * k))
    return ds


INGESTERS["calzadilla2023"] = calzadilla2023
