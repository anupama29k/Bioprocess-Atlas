# Changelog

## v0.5.0 — 2026-11-13

9 measured datasets · 172 runs · 7 organisms · 48 one-factor comparisons · 222,437 observations

**New dataset**

- `calzadilla2023` — one 2 L HEK293 chemostat producing the Her1 extracellular domain, held at five dilution-rate set points over 85 days. CC BY 4.0. The atlas's first continuous cultivation, and its first human cell line.

**Schema 0.5**

- New variable: dilution rate.
- Phenotypes that assume a batch or fed-batch accumulation (specific growth rate, substrate depletion time, biomass yield, byproduct peak and reconsumption, integral viable cell density, cell-specific productivity) are no longer computed for continuous runs, where they are meaningless.
- The chemostat is a single vessel, so its five steady-state periods are stored as one run with the dilution rate as a time series, and no comparison is asserted between them.

**Statistics, after a verification pass against the source papers**

- **FDR-corrected q-values.** `effects.csv` gains `q_bh`, a Benjamini-Hochberg false-discovery rate computed within each campaign. Of 211 tests, 78 reach p < 0.05 but only 51 survive at q < 0.05, and about 11 of the 78 are expected by chance. Two effects previously quoted on p alone do not survive and are now stated as effect sizes: ΔaceA acetate reconsumption (d = −5.6, q = 0.20) and W-68 peak VCD at 37 °C (d = −20.5, q = 0.062).
- **`replicate_scope` on every run**, recording what its replicate number means, taken from each source's own description of its hardware. 44 of the 51 surviving effects rest on vessels run side by side in one parallel system, 7 on replicates of unstated kind, and none on replication a source calls independent. The Welch test treats replicates as independent, so its p-values are anticonservative.

**Corrections**

- **Galactosylation is now the source study's own metric.** `galactosylation_index` (galactoses per glycan) is replaced by `galactosylation_percent`, implemented to Schäpertöns et al.'s Eq. 7. The atlas returns 24.3 % for the galactose arms against the paper's reported 24 %. Incomplete glycoform pairs are now weighted rather than dropped.
- **NISTCHO titer was being silently dropped** after a refactor, because `itertuples` renames a column containing a dot. Cell counts and titer are now read from the named file, which also removed the atlas's dependence on an inferred condition-letter mapping.
- **B. subtilis volume-balance rows were being parsed as data**, which produced a glucose depletion time of 0.39 h. Row parsing is now bounded at the volume table.
- `Yxs_batch_gg` is no longer computed from cell counts, where it has no meaning.
- Ammonium in mM and in g/L are kept as separate variables and never converted, since no source states whether a g/L figure is on an NH4⁺, elemental-N or salt basis.

**External validation**

- Galactosylation 24.3 % vs the paper's 24 %; B. subtilis biomass at glucose exhaustion 6.98 ± 1.04 g/L vs the paper's 7.22 ± 0.82; NISTCHO condition mapping confirmed across all 25 runs. The remaining 7–17 % gap in B. subtilis Y_P/X is explained: the paper fits polynomials and reads parameters off the fits.

41 tests.

## v0.4.0 — 2026-11-10

8 measured datasets · 171 runs · 6 organisms · 48 one-factor comparisons · 222,062 observations

**New dataset**

- `schapertoens2026` — NISTCHO fed-batch under seven feeding strategies (feed frequency, glucose level, galactose supplementation), 25 vessels, with cNISTmAb titer **and N-glycan product quality**. CC BY 4.0. The atlas's first product-quality data: galactose supplementation roughly doubles galactosylation and halves the G0F share, with no titer penalty on the standard feed.

**Schema 0.4**

- New variables: galactosylation index, afucosylation percent, G0F percent.
- New phenotypes: galactosylation, afucosylation and G0F at the last sample.

31 tests.

## v0.3.0 — 2026-11-10

7 measured datasets · 146 runs · 6 organisms · 45 one-factor comparisons · 220,453 observations

**New datasets**

- `schellenberg2021` — CHO DG44 IgG producer, 12-day fed-batch: glutamate, glucose overfeeding, reduced feed, pH and CO2 stress against a standard process (ambr15 and shake flask), 21 vessels. CC0.
- `uhlendorff2026_cho` — CHO-K1 batch and fed-batch in DASGIP reactors with inline turbidity, permittivity, DO and CO2: agitation/aeration, seeding density, glutamine, medium, impeller position and feeding, 24 runs in duplicate. CC BY 4.0.
- `salesvallverdu2024` — *K. phaffii* glucose-limited fed-batch, hypoxia against matched normoxia at four growth rates and for three secreted lipases, 11 runs. CC BY 4.0.
- `hiller2024` — *B. subtilis* BMV9 surfactin in a 30 L bioreactor: six growth-rate-controlling feed rates (0.075–0.4 h⁻¹) in duplicate plus a batch process, 13 runs. CC BY 4.0.

**Schema 0.3**

- New variables: glutamine, glutamate, ammonium (both mass and molar), osmolality, cell diameter, arabitol, product activity.
- New phenotypes: final titer, cell-specific productivity (qP), integral viable cell density, peak lactate, end cell diameter, viability-decline time, end acetate, end glucose, product per biomass, end ethanol and arabitol, final product activity.
- A contrast can now declare which phenotypes it is valid for (`phenotypes` / `exclude_phenotypes`), for designs where a reference run is compromised or a phenotype does not apply.
- Feed-dependent phenotypes are no longer flagged as confounded when the perturbation is itself the feed, because there they are the result.
- Biomass yield is no longer computed from cell counts.

**Checks**

28 tests, including each source study's own headline finding: the ΔaceA glyoxylate-shunt defect, the ale strain's collapse at 37 °C, CHO productivity under glucose overfeeding and CO2 limitation, the *Pichia* hypoxia effect, and the *B. subtilis* surfactin optimum at µ = 0.25 h⁻¹.

## v0.2 — 2026-11-10

Seed release: 3 measured datasets (*E. coli* batch with off-gas, 24 Keio-knockout mini-bioreactors, 52 brewing-yeast propagations), 77 runs, schema with explicit contrasts, phenotype engine vendored from the Bioprocess Analyzer, 13 tests.
