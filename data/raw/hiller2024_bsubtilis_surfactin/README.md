# B. subtilis BMV9 surfactin: six growth-rate-controlling feed rates in 30 L fed-batch (measured)

**Source:** Hiller E, Off M, Hermann A, Vahidinasab M, Benatto Perino EH, Lilge L, Hausmann R (2024). *The influence of growth rate-controlling feeding strategy on the surfactin production in Bacillus subtilis bioreactor processes.* Microbial Cell Factories 23:260. https://doi.org/10.1186/s12934-024-02531-w
**File:** `12934_2024_2531_MOESM2_ESM.xlsx`, the article's Supplementary Material 2, reproduced unchanged.
**License:** CC BY 4.0 (the article and its supplementary material).

## Design

*B. subtilis* BMV9 (*spo0A3*, *trp*⁺, *sfp*⁺, Δ*manPA*), a surfactin producer, in a 30 L stirred-tank bioreactor with three Rushton turbines: 12 L of chemically defined mineral salt medium at 37 °C and pH 7.0, with dissolved oxygen held above 50 %. The final volume is about 18.7 L.

Every run shares the same batch phase. When the batch glucose is exhausted, an exponential glucose feed starts, set to hold a chosen growth rate:

| condition | set µ (h⁻¹) | bioreactors |
|---|---|---|
| batch, no feed | – | 1 |
| feed for µ = 0.075 (atlas reference) | 0.075 | 2 |
| feed for µ = 0.15 | 0.15 | 2 |
| feed for µ = 0.2 | 0.2 | 2 |
| feed for µ = 0.25 | 0.25 | 2 |
| feed for µ = 0.3 | 0.3 | 2 |
| feed for µ = 0.4 | 0.4 | 2 |

Measured per sample: OD600, cell dry weight, glucose, surfactin, acetate and ammonium (ammonium is not in the batch sheet). The feed start (FS) and feed end (FE) samples are labelled in the source; the feed start is recorded in each run's notes.

This is the atlas's first graded **process** dose series with a product titer, and its largest vessel.

## What the atlas finds, against the paper

| paper | atlas |
|---|---|
| highest surfactin ≈ 36 g/L at µ = 0.25 h⁻¹ | 36.0 g/L mean at µ = 0.25, the maximum of the series |
| acetate up to ≈ 3 g/L at the higher feed rates | acetate at the last sample reaches 2.9 g/L (µ = 0.3) and 3.0 g/L (µ = 0.4), against 0.03 g/L at µ = 0.075 |
| feeding above 0.25 h⁻¹ accumulates glucose | glucose at the last sample: 0.0–0.3 g/L up to µ = 0.3, then 44–53 g/L at µ = 0.4 |

## Caveats

- **`product_per_biomass_gg` runs 7–17 % above the paper's Y_P/X, and the reason is their smoothing, not the definition.** The paper defines Y_P/X as P/X at the time of maximum surfactin (their Eq. 3). Surfactin peaks at the last sample in every run here, so that definition, "last over last", "P_max over X_max" and the increment form all give the same number — checked, they agree to three decimals. The remaining gap comes from the paper's Methods: the series were fitted with polynomial functions and the production parameters read off those fits, whereas the atlas computes from the measured points. A polynomial fit pulls in the final point, which is the maximum, so their values come out lower. The atlas value is therefore transparent and reproducible but not identical to the published figure, and the phenotype is deliberately not called Y_P/X. Y_P/S is not computed at all: the feed's glucose concentration is not in the file.

- **The parsing is validated against the paper's own batch-phase figures.** The paper reports that the initial 25 g/L glucose was exhausted between 11.5 and 16 h at a mean CDW of 7.22 ± 0.82 g/L. From the ingested data the atlas finds 24.0 g/L initial glucose and 6.98 ± 1.04 g/L at exhaustion. One run (µ = 0.075 R1) falls outside their time range at 20 h, because the atlas requires glucose below 0.1 g/L while that run still measured 0.44 g/L at 16 h.
- **The feed is the perturbation here, so feed-dependent phenotypes are the result rather than a confounder.** The atlas marks `confounded_by_feed` false for any contrast whose perturbation type is `feed`; `mu_fedbatch_h` is the direct check on whether the set growth rate was reached (it is not: the actual rate falls short of the set rate, increasingly so above 0.2 h⁻¹).
- **Batch-phase phenotypes act as a negative control.** The feed starts after the batch phase, so `mu_batch_h` should not differ between arms, and in the atlas it does not, apart from a single arm where the sparse batch sampling makes the fit unreliable.
- Two runs (µ = 0.2 R2, µ = 0.4 R2) have only two samples before the feed starts, which is too few for a batch growth-rate fit, so that phenotype is missing for them.
- At **µ = 0.4 R1** glucose never falls below 0.1 g/L, so no glucose-depletion time and no fed-batch growth rate are reported for that run. The feed outran consumption from the start.
- A few glucose values in the µ = 0.15 sheet are slightly negative (down to −0.0009 g/L), which is measurement noise around zero. They are kept as reported.
- `acetate_reconsumed_frac` is not compared for this dataset. It measures acetate taken back up 0.75–2.5 h after glucose runs out, but here the feed starts at exactly that moment, so the window lies inside the feed phase. The batch-phase acetate peak is still compared, and acts as another negative control.
- The volume-balance table at the foot of each sheet is not part of the time course and is not ingested.
