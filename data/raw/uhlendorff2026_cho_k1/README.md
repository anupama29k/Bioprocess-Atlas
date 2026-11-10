# CHO-K1 batch and fed-batch cultivations with inline turbidity, permittivity, O2 and CO2 (measured)

**Source:** Uhlendorff S, Fulek R, Eimler J, Pein-Hackelbusch M, Frahm B (2026). *Dataset Based on Chinese Hamster Ovary (CHO) Cultivations including Turbidity, Permittivity, O2 and CO2 Measurements.* Zenodo, https://doi.org/10.5281/zenodo.20829178 (v1). Data descriptor: Preprints.org, https://doi.org/10.20944/preprints202611.0377.v1 (not peer reviewed).
**License:** CC BY 4.0. Please cite the authors.

## What is here

- `offline/`: the authors' 24 time-aligned offline tables, unchanged. Sampled 1–2 times a day: VCD, TCD, viability, diameter, compactness and aggregation (Countstar, technical triplicates; mean and SEM), plus glucose, lactate, ammonium and osmolality.
- `inline_15min_medians.csv.gz`: the authors' inline tables (logged every 30–60 s) reduced to 15-minute medians for selected channels: dissolved CO2, transmission (turbidity), permittivity, pH, DO, temperature, agitation, cumulative feed and liquid volume. The multifrequency capacitance spectra and gas-flow channels are left out. They are in the original Zenodo files.
- `conditions.csv`: the experiment design, transcribed from the data descriptor's table.

## Design

CHO-K1 (Bielefeld University) in a DASGIP parallel system with two reactors (R01, R02) per experiment. All runs are at 37 °C, pH 7.0 and 30 % DO set point. Batch runs start at 650 mL; fed-batch runs start at 550 mL and are kept near that volume.

The four standard batches (EXP001–004: TC-42, 80 rpm, 0.125 vvm) are repeats of the same process on different days. The later experiments each change one or two factors from a 100 rpm, 0.04 vvm base. The atlas compares:

| contrast | perturbed | reference |
|---|---|---|
| 100 rpm + 0.04 vvm vs 80 rpm + 0.125 vvm | EXP005 | EXP001–004 |
| seeding 6 × 10⁵ vs 4 × 10⁵ cells/mL | EXP006 | EXP005 |
| glutamine 6 vs 8 mM | EXP007 | EXP005 |
| TCX6D vs TC-42 medium (batch) | EXP008 | EXP005 |
| lowered vs standard bottom impeller (batch) | EXP009 | EXP008 |
| fed-batch vs batch (TC-42) | EXP010 | EXP005 |
| TCX6D vs TC-42 medium (fed-batch) | EXP012 | EXP010 R01 |
| lowered vs standard bottom impeller (fed-batch) | EXP011 R01 | EXP012 |

## Caveats

- Agitation and aeration change together between the standard batches and the 100 rpm base.
- **EXP005 is in five contrasts (the 100 rpm base).** Its O2 supply failed: DO falls below 10 % at about 104 h and supply was restored at 121 h. Peak VCD in both EXP005 reactors falls at 121.5 h, so peak VCD, IVCD and anything later are capped by the failure. These contrasts therefore report only the growth rate, which is fitted from earlier samples (`phenotypes: mu_batch_h` in `contrasts.csv`).
- The four standard batches (80 rpm, 0.125 vvm) also fall below 10 % DO from about 110–118 h, before their peak VCD, despite the 30 % set point. They are probably oxygen-limited at high density. From EXP006 on, DO stays above about 12 %.
- The runs were stopped at different times (161–328 h), so viability at the last sample is not compared. The atlas uses `viability_80_h` (the first sample after peak with viability below 80 %) instead.
- **EXP010 R02:** the feed line was blocked, so feeding effectively started at about 66 h. It is kept in the data but left out of the fed-batch medium contrast.
- **EXP011 R02** was contaminated with bacteria from about 51 h. It is kept in the data but excluded from every contrast.
- No product titer: CHO-K1 here is a non-producing host.
