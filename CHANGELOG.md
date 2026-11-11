# Changelog

## v0.5.0 — 2026-11-10

9 measured datasets · 172 runs · 7 organisms · 48 one-factor comparisons · 222,437 observations

**New dataset**

- `calzadilla2023` — one 2 L HEK293 chemostat producing the Her1 extracellular domain, held at five dilution-rate set points over 85 days. CC BY 4.0. The atlas's first continuous cultivation, and its first human cell line.

**Schema 0.5**

- New variable: dilution rate.
- Phenotypes that assume a batch or fed-batch accumulation (specific growth rate, substrate depletion time, biomass yield, byproduct peak and reconsumption, integral viable cell density, cell-specific productivity) are no longer computed for continuous runs, where they are meaningless.
- The chemostat is a single vessel, so its five steady-state periods are stored as one run with the dilution rate as a time series, and no comparison is asserted between them.

32 tests.

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
