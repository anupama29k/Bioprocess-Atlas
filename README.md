# Bioprocess Atlas

**Real fermentation runs from published studies, changed one factor at a time and reduced to the same set of phenotypes, so that the effect of a gene, a temperature or a feed can be compared across labs, organisms and scales.**

**[Interactive explainer →](https://anupama29k.github.io/bioprocess-atlas/)** (architecture diagram you can drag around, one comparison followed through the pipeline, effects explorer). The source is `docs/index.html`.

| release | datasets | bioreactor runs | organisms | one-factor comparisons | observations |
|---|---|---|---|---|---|
| v0.2 | 3 | 77 | 3 | 18 | 100,239 |

## Why

Open upstream bioprocess data is scarce. The largest public benchmark is simulated, and real industrial datasets usually keep only harvest-day values. Single-cell biology fixed this with perturbation atlases. The Bioprocess Atlas applies the same idea to bioreactors:

- **Compare.** Every perturbation is measured against its own reference, ideally reactors run on the same day, and every number carries replicate statistics (Welch t-test, effect size d). Arms with a single run report the difference only.
- **Reproduce.** Every table is rebuilt from the unchanged published files with one command. Each release is stamped with a checksum and the version of the analysis engine.
- **Model.** Set points become doses (°C, vvm, °P, feed rate). This is the bioprocess version of graded genetic perturbations, and it is the training data a virtual bioreactor needs.

## How it is built

```
published studies ──► data/raw/<dataset>/       unchanged files + README (citation, license, caveats)
                          │
                          ▼
                     atlas/ingest.py            one function per dataset → runs, observations, contrasts
                          │
                          ▼
                     atlas/schema.py            controlled vocabulary; the unit is part of the variable name
                          │
                          ▼
                     atlas/phenotypes.py        phenotypes per run (growth engine from Bioprocess Analyzer)
                          │                     + effects per contrast (Welch t, effect size d)
                          ▼
                     release/                   CSV tables + manifest.json (schema, engine, checksum)
                          │
          ┌───────────────┼──────────────────┐
   Bioprocess Analyzer  notebooks/papers   virtual bioreactor (planned)
```

### Release tables

| file | one row per | key columns |
|---|---|---|
| `runs.csv` | bioreactor run | run_id, dataset_id, organism, strain, condition, perturbation_type, dose, reference_condition, vessel, working_volume_L |
| `observations.csv.gz` | measurement | run_id, time_h, variable, value |
| `phenotypes.csv` | run × phenotype | run_id, phenotype, value, confounded_by_feed |
| `effects.csv` | perturbation × phenotype | mean, ref_mean, pct_change, effect_size_d, p_welch, n, ref_n |
| `contrasts.csv` | designed comparison | dataset_id, label, perturbed condition, reference condition |
| `manifest.json` | release | schema_version, counts, vocabulary, phenotype definitions, checksum, engine |

Phenotypes are taken only from phases in which both arms had the same conditions (for example, the batch phase before feeding). Phenotypes that depend on the feed are flagged `confounded_by_feed`.

## Datasets

| id | what | runs | source | license |
|---|---|---|---|---|
| `simensen2022` | *E. coli* MG1655, aerobic batch on glucose M9, 1 L, with off-gas | 1 | Simensen V et al. (2022) PLOS ONE 17:e0262450. [doi:10.1371/journal.pone.0262450](https://doi.org/10.1371/journal.pone.0262450) | CC BY 4.0 |
| `hans2020` | 8 *E. coli* K-12 strains (Keio knockouts and parents) in 24 parallel fed-batch mini-bioreactors | 24 | Hans S et al. (2020) Bioengineering 7:145. [doi:10.3390/bioengineering7040145](https://doi.org/10.3390/bioengineering7040145) | CC BY 4.0 |
| `uhlendorff2026` | Brewing yeast propagation: 2 strains × temperature, aeration, wort gravity, pre-culture, scale | 52 | Uhlendorff S et al. (2026) Zenodo. [doi:10.5281/zenodo.19690248](https://doi.org/10.5281/zenodo.19690248) | CC BY 4.0 |

Each `data/raw/<dataset>/README.md` lists the conditions and the known caveats.

### Example findings in v0.2

- The *E. coli* ΔaceA strain reconsumes 5% of its acetate after glucose runs out, against 78% for the parent (p = 0.015, d = −5.6). This is the expected loss of the glyoxylate shunt.
- The ale yeast W-68 reaches a 68% lower peak viable cell density at 37 °C than at 22 °C (p = 0.027), and its viability falls to 26% against 90%. The lager strain W-34/70 is unaffected (+1%).
- A 3 h pre-culture lowers peak viable cell density by 37–41% in both strains, compared with a 24 h pre-culture.
- All Keio knockouts show a similar −9 to −14% shift in peak acetate. A shift this uniform probably comes from the shared background or the plate rather than from the individual genes; treat it as a caveat.

## Quick start

```bash
pip install -r requirements.txt
python -m atlas.build          # rebuilds release/ from data/raw/
python tests/test_atlas.py     # 13 checks, including known biology
```

```python
import pandas as pd
eff = pd.read_csv("release/effects.csv")
eff[eff.p_welch < 0.05].sort_values("effect_size_d")
```

The [Bioprocess Analyzer](https://github.com/anupama29k/bioprocessanalyzer) has an Atlas page that browses a release, and it can analyse your own runs with the same engine.

## Adding a dataset

See [CONTRIBUTING.md](CONTRIBUTING.md). In short: unchanged files plus a README go in `data/raw/<name>/`, you write one ingester function, and you rebuild. The build rejects any variable that is not in the vocabulary.

## Roadmap

| phase | goal | gate |
|---|---|---|
| 0 · Seed (now) | schema, build, tests, 3 measured datasets | engine reproduces published values |
| 1 · Aggregate | 20–50 open datasets, 3+ organisms, field names mapped to the PREFER ontology | 200+ runs |
| 2 · Generate | designed campaign on 24–48 parallel mini-bioreactors, with the parent strain in every plate column | batch effect smaller than the biological effects |
| 3 · Scale and model | partner labs, graded process doses, a virtual bioreactor tested on held-out perturbations | |

## Repository layout

```
atlas/        schema, ingesters, phenotype engine, build
data/raw/     source files as published, one folder per dataset
release/      built tables (committed so they can be used without building)
docs/         interactive explainer (GitHub Pages)
tests/        test_atlas.py
```

## License and citation

Created and maintained by **Anupama Kozhiyalam** ([@anupama29k](https://github.com/anupama29k)).


- **Code:** MIT, © 2026 Anupama Kozhiyalam (see `LICENSE`).
- **Data:** each dataset keeps its original license (all CC BY 4.0 so far). If you use the release, cite the original studies listed above and this repository (see `CITATION.cff`).

The growth-rate engine in `atlas/engine.py` is vendored from Bioprocess Analyzer at commit 2439bf9.
