# Changelog

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
