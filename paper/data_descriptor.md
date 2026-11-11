# The Bioprocess Atlas: comparable perturbation effects from published bioreactor cultivations

**Anupama Kozhiyalam**¹

¹ [affiliation to fill in] · ORCID [to fill in] · correspondence: [email to fill in]

*Draft data descriptor, prepared for submission to* Scientific Data *(or* Data in Brief*). Numbers refer to release v0.5.0 (v0.4.0 is archived at doi:10.5281/zenodo.24125548). Every figure in the text is generated from that release and should be regenerated before submission if the release changes.*

---

## Abstract

Upstream bioprocess development still begins, for most new strains and products, from a blank sheet. The experiments that would inform it have been run many times, in many laboratories, but their results remain locked in figures and in formats that cannot be compared across studies, and the largest openly available benchmark for the field is simulated rather than measured. Here we describe the Bioprocess Atlas, an open collection of measured bioreactor cultivations drawn from published studies in which a single factor was varied against a reference condition. Release v0.5.0 holds 9 datasets, 172 bioreactor runs and 222,437 measurements spanning seven organisms, from microbial fermentations to mammalian cell culture, in batch, fed-batch and continuous mode, and including both product titer and N-glycan product-quality attributes. Every run is reduced to the same set of phenotypes by a single versioned analysis engine, and every perturbation is reported against its own reference with replicate statistics. All tables are rebuilt from the unchanged published source files by one command, and each release carries a checksum and the engine version. The atlas is intended as a reference set for process development, for benchmarking, and as training data for predictive models of bioreactor behaviour.

---

## Background & Summary

Process development for a new strain or product is largely empirical. Teams run screening campaigns across temperature, pH, feed rate, medium and dissolved oxygen, discover which factors matter for their particular system, and then keep that knowledge inside the organisation. The same lessons are relearned repeatedly, at considerable cost, because the results that would generalise are not available in a form that allows comparison.

This is not for want of experiments. The published literature contains a large number of well-controlled bioreactor studies in which one factor was varied against a control. What it lacks is any convention for reporting them such that the effect of a temperature shift in one laboratory can be placed alongside the effect of a gene deletion in another. Values are reported in figures, in idiosyncratic spreadsheet layouts, in units that differ between papers, and with phenotypes defined differently each time. Where numeric data are deposited at all, they are deposited per study, with no shared schema.

The consequences are visible in the resources the field does have. The most widely used open benchmark for bioprocess modelling is a simulation. Industrial datasets, when released, typically retain only harvest-day values, which removes the dynamics that process models need. Public repositories contain scattered deposits of real time courses, but finding them requires reading each paper, and using them requires rebuilding the analysis each time.

Single-cell biology faced a structurally similar problem and addressed it by building perturbation atlases: large, uniformly processed collections in which each perturbation is reported against a matched control, enabling comparison across experiments that were never designed together. The Bioprocess Atlas applies that approach to bioreactors.

The design rests on four decisions.

**Measured data only.** Simulated runs are excluded. The atlas is small as a result, and deliberately so: its value depends on every value in it having been measured in a vessel.

**One factor against a reference.** Each comparison names a perturbed condition, a reference condition, and the single factor that differs. Where a study's design does not support such a comparison — most often because the reference was run in a separate campaign — the runs are included but no effect is computed, and the reason is recorded.

**One analysis engine for every run.** Phenotypes such as specific growth rate are computed by a single implementation, vendored into the repository at a fixed commit and stamped on every release. A change to the engine therefore produces a new release rather than silently altering published values.

**Rebuild from source, not from intermediates.** The published files are stored unchanged. Every table in a release is regenerated from them by one command, and a checksum over all observations is published with each release.

Release v0.5.0 contains 9 datasets and 172 runs across *Escherichia coli*, *Bacillus subtilis*, *Saccharomyces cerevisiae*, *Saccharomyces pastorianus*, *Komagataella phaffii*, Chinese hamster ovary (CHO) and HEK293 cells; 72 batch, 99 fed-batch and one continuous cultivation; working volumes from 15 mL to 12 L; and 48 one-factor comparisons, of which 211 phenotype-level effects carry replicate statistics and 77 reach p < 0.05. The perturbations include gene deletions, strain background, temperature, pH, dissolved oxygen, medium composition, inoculum age, aeration, agitation, impeller geometry, vessel scale and six graded feed rates.

---

## Methods

### Dataset selection

Candidate studies were required to meet four criteria: (i) measurements from real cultivations, with multiple sampling points per run rather than end-point values only; (ii) at least one factor varied against a reference condition, preferably with replicate vessels in each arm; (iii) numeric values available for download, in supplementary files or a repository, since values readable only from figures cannot be ingested without re-digitisation; and (iv) a licence permitting redistribution. Seven of the eight datasets are CC BY 4.0 and one is CC0 1.0.

Studies were found by searching bibliographic databases and data repositories (Zenodo, figshare, Mendeley Data, Dryad, 4TU, Dataverse instances and institutional repositories), and by searching open-access journals that routinely attach numeric supplementary files. The search is ongoing; candidate datasets that failed one criterion are recorded so that they are not re-examined.

### Ingestion

Each dataset occupies a directory under `data/raw/` holding the published files unchanged, together with a README recording the citation, the digital object identifier, the licence, the organism, strain, medium, vessel and working volume, the experimental design, and the known caveats. Reshaping from a publication layout into long format is performed in code at build time rather than by editing the stored files, so that the stored files remain comparable with the publisher's.

One ingester function per dataset converts its files into runs, observations and explicit contrasts. Adding a dataset requires adding one such function and no changes elsewhere.

### Schema and controlled vocabulary

Four tables are keyed by `run_id`. A **run** is one vessel: it records organism, strain, mode, vessel, working volume, medium, temperature and pH set points, replicate number, campaign, condition, reference condition, perturbation type, perturbation and dose. An **observation** is a single measurement: run, time in hours, variable, value. A **phenotype** is a derived quantity for one run. An **effect** is one perturbation against its reference for one phenotype.

Units are carried by the variable name — `glucose_gL`, `DO_percent`, `vcd_cells_mL` — rather than in a separate column, so values are directly comparable across datasets without unit reconciliation at analysis time. The build rejects any variable outside the vocabulary, which prevents silent divergence as datasets accumulate. Release v0.5.0 defines 34 variables, of which 33 are populated, from standard process measurements (biomass, substrates, byproducts, dissolved oxygen, pH, off-gas) through inline sensor channels (permittivity, dissolved carbon dioxide, turbidity) to product measurements (titer, enzyme activity) and product-quality attributes (galactosylation index, afucosylation, G0F share).

### Phenotypes

Twenty-eight phenotypes are computed per run. Specific growth rate is obtained with a changepoint-segmented, sliding-window log-linear fit; for cultures quantified by cell count rather than dry weight, the fit uses samples below 80% of peak viable cell density. Other phenotypes include substrate depletion time, peak and final byproduct concentrations, byproduct reconsumption, biomass yield, peak viable cell density, integral viable cell density, viability decline, final titer, cell-specific productivity and product-quality attributes at harvest.

Phenotypes are computed only from phases in which the compared runs experienced the same conditions. In fed-batch cultivations the batch phase precedes the perturbation in several designs, which makes batch-phase phenotypes a useful internal control.

### Effects and replicate statistics

For each comparison, the perturbed and reference arms are summarised by mean, standard deviation and n; the difference, the percentage change, the standardised effect size and a Welch two-sample t-test are reported, together with a Benjamini-Hochberg false-discovery rate (`q_bh`) computed **within each campaign**, that being the family of tests a reader actually scans when asking which effects in a given screen are real. Effect size and p-value require at least two runs per arm; where an arm has a single vessel, the difference is reported and the inferential fields are left empty rather than filled with a value the design cannot support. Of 365 effect rows, 211 carry a p-value. Seventy-eight of these fall below p = 0.05, but only fifty-one survive at q = 0.05, and roughly eleven of the seventy-eight would be expected by chance alone. This is reported prominently rather than buried: with two to four vessels per arm, the p-value is a weak instrument, and several effects that are large, reproducible and consistent with established biology nonetheless fail correction. The ΔaceA acetate-reconsumption result is the clearest case: the effect size is d = −5.6 and the underlying biochemistry is not in doubt, yet with three vessels per arm it reaches only p = 0.015, which becomes q = 0.20 among the fifty-nine tests of that screen. Users are directed to read effect sizes alongside q, and not to filter on p alone.

Two features of the effect model deserve note.

**Scoped contrasts.** A contrast may declare which phenotypes it is valid for. This was introduced for a CHO-K1 dataset in which the reference experiment suffered an oxygen-supply failure between roughly 104 and 121 hours, with peak viable cell density in both reference vessels falling at 121.5 hours — that is, the peak was set by the failure rather than by biology. Comparisons against that reference are restricted to the specific growth rate, which is fitted from earlier samples. Recording such a restriction in the contrast, rather than omitting the dataset or silently reporting a compromised number, keeps the usable part of a flawed experiment usable.

**Feed as perturbation.** Phenotypes that depend on the feed are flagged as confounded, so that a difference in final biomass between two strains in fed-batch is not read as a strain effect. Where the perturbation is itself the feed, however, those phenotypes are the result rather than a confounder, and the flag is not applied. This distinction matters for the graded feed-rate series described below, in which the fed-batch growth rate is the direct measure of whether the intended growth rate was achieved.

### Reproducibility

`python -m atlas.build` regenerates every table from `data/raw/`. Each release carries a manifest recording the schema version, the engine version, counts, the full variable and phenotype definitions, and a SHA-256 checksum over all observations. The analysis engine is vendored at a fixed commit of an independent repository and identified in the manifest.

---

## Data Records

The release is archived at Zenodo (doi:10.5281/zenodo.24125547, which always resolves to the latest version) and developed at https://github.com/anupama29k/Bioprocess-Atlas.

| File | One row per | Key fields |
|---|---|---|
| `runs.csv` | bioreactor run | run_id, dataset_id, organism, strain, mode, vessel, working_volume_L, medium, temperature_C, pH_setpoint, replicate, campaign, condition, reference_condition, perturbation_type, perturbation, perturbation_dose, source_doi, license, notes |
| `observations.csv.gz` | measurement | run_id, time_h, variable, value |
| `phenotypes.csv` | run | run_id and 28 phenotype columns |
| `effects.csv` | perturbation × phenotype | campaign, condition, reference, perturbation_type, perturbation, dose, phenotype, unit, mean, sd, n, ref_mean, ref_sd, ref_n, difference, pct_change, effect_size_d, p_welch, q_bh, confounded_by_feed |
| `contrasts.csv` | designed comparison | campaign, perturbation_type, perturbation, dose, dataset_id, runs, reference_runs, phenotype scoping |
| `manifest.json` | release | schema_version, counts, variable and phenotype definitions, organisms, engine, observations checksum |

### Source datasets in v0.4.0

| id | organism | design | runs | source | licence |
|---|---|---|---|---|---|
| `simensen2022` | *E. coli* MG1655 | aerobic batch on glucose minimal medium with off-gas, 1 L | 1 | PLOS ONE 17:e0262450 | CC BY 4.0 |
| `hans2020` | *E. coli* K-12 | eight strains including six Keio single-gene knockouts, 24 parallel adaptive fed-batch mini-bioreactors | 24 | Bioengineering 7:145 | CC BY 4.0 |
| `uhlendorff2026` | *S. cerevisiae*, *S. pastorianus* | brewing-yeast propagation: two strains × temperature, aeration, wort gravity, pre-culture age and scale (0.7 L and 10 L) | 52 | Zenodo 10.5281/zenodo.19690248 | CC BY 4.0 |
| `schellenberg2021` | CHO DG44, IgG producer | 12-day fed-batch: glutamate, glucose overfeeding, reduced feed, pH shift and CO₂ limitation against a standard process, in 15 mL bioreactors and shake flasks | 21 | LUH 10.25835/0066265 | CC0 1.0 |
| `uhlendorff2026_cho` | CHO-K1 | batch and fed-batch with inline turbidity, permittivity, dissolved oxygen and CO₂: agitation and aeration, seeding density, glutamine, medium, impeller position, feeding | 24 | Zenodo 10.5281/zenodo.20829178 | CC BY 4.0 |
| `salesvallverdu2024` | *K. phaffii* | glucose-limited fed-batch, oxygen limitation against matched normoxia at four growth rates and for three secreted lipases | 11 | Microb Biotechnol 17:e14411 | CC BY 4.0 |
| `hiller2024` | *B. subtilis* BMV9 | surfactin production in a 30 L vessel: six growth-rate-controlling feed rates in duplicate, plus a batch process | 13 | Microb Cell Fact 23:260 | CC BY 4.0 |
| `calzadilla2023` | HEK293 | one 2 L chemostat held at five dilution-rate set points over 85 days; single vessel, so no comparison is asserted | 1 | Data in Brief 50:109604 | CC BY 4.0 |
| `schapertoens2026` | NISTCHO | seven feeding strategies (frequency, glucose level, galactose) with titer and N-glycan quality attributes | 25 | npj Syst Biol Appl 12:67 | CC BY 4.0 |

---

## Technical Validation

The atlas is validated by requiring that it reproduce, from the raw files and through the common pipeline, the principal finding that each source study reported independently. Thirty-three automated checks run on every build; the build fails if any does not pass.

| Source study | Finding reported by the authors | Recovered by the atlas |
|---|---|---|
| Hans 2020 | loss of isocitrate lyase (Δ*aceA*) removes the glyoxylate shunt, so acetate cannot be reconsumed | 5% of peak acetate reconsumed against 78% for the parent (d = −5.6; p = 0.015, q = 0.20) |
| Uhlendorff 2026 | the ale strain fails at elevated temperature while the lager strain does not | W-68 peak viable cell density 68% lower at 37 °C than same-day 22 °C controls (d = −20.5; p = 0.027, q = 0.062), viability 26% against 90%; W-34/70 unaffected (+1%) |
| Schellenberg 2021 | glucose overfeeding and CO₂ limitation raise cell-specific antibody productivity | qP +45% and titer +20% under overfeeding; qP +21% under CO₂ limitation |
| Sales-Vallverdú 2024 | oxygen limitation improves recombinant protein production | final lipase activity 75–153% higher than matched normoxic runs, in all five comparisons across three products and four growth rates |
| Hiller 2024 | surfactin titer is maximal at a feed set for µ = 0.25 h⁻¹; higher feed rates accumulate glucose and acetate | maximum of the series at µ = 0.25 h⁻¹, 36.0 g/L; 49 g/L residual glucose at µ = 0.4 h⁻¹; end acetate 2.9–3.0 g/L at the two highest rates against 0.03 g/L at the lowest |
| Schäpertöns 2026 | galactose supplementation shifts N-glycosylation | galactosylation index 0.22 → 0.46–0.51 and G0F share 56% → ~30% in all three galactose arms (p < 0.001), with no titer penalty on the standard feed (−4%, p = 0.48) |

Validation additionally uses negative controls that the designs supply. In the *B. subtilis* feed-rate series the feed begins only after a batch phase common to all arms, so batch-phase growth rate should not track the feed set point; across the five comparisons the median absolute change is 4.5%, against changes of 80–385% in the fed-batch growth rate. Independent published values are reproduced where available: for the *E. coli* batch cultivation the engine returns µ = 0.768 h⁻¹ against 0.71 h⁻¹ reported, and 0.74 h⁻¹ when computed from the carbon dioxide evolution rate. Cell-specific productivities fall between 13 and 36 pg/cell/day, within the range expected of an industrial antibody-producing line.

The pipeline also surfaces artefacts that single-study analysis tends to hide. In the Keio knockout screen, every knockout shows a 9–14% shift in peak acetate in the same direction, a uniformity across unrelated genes that points to a shared cause — the resistance cassette, plate position or assay batch — rather than to gene function; this is recorded as a caveat on the dataset. The oxygen-supply failure in the CHO-K1 reference experiment, described above, was identified by comparing the time of peak viable cell density with the dissolved-oxygen trace, and is handled by scoping rather than by discarding the data. In the NISTCHO dataset, three feeding strategies were run in a campaign that contained no standard-feed arm; those nine runs are included with their phenotypes but produce no effect rows, because comparing them with the reference would confound strategy with campaign.

---

## Usage Notes

Release tables are plain CSV and load without special tooling:

```python
import pandas as pd
eff = pd.read_csv("release/effects.csv")
eff[eff.p_welch < 0.05].sort_values("effect_size_d")
```

Several conventions should be observed when using the data.

**Respect the campaign boundary.** Runs carry a campaign identifier. Comparisons are constructed within a campaign; constructing one across campaigns reintroduces precisely the batch confounding the atlas is built to avoid.

**Check `confounded_by_feed` before interpreting an effect**, and note that it is deliberately not set where the feed is the perturbation.

**Do not compare product measurements across products.** Enzyme activity is reported in assay units that differ per enzyme, under the variable `product_activity_AU_mL`, and is comparable only within one product.

**Filter on `q_bh`, not on `p_welch`, and read the effect size.** The release reports 211 tests; filtering on p < 0.05 alone collects an expected eleven false positives. Conversely, a large and mechanistically expected effect may fail correction simply because the design has three vessels per arm, so a non-significant q is not evidence of absence.

**Prefer effect sizes to percentage changes when a reference value is near zero.** Several residual-substrate comparisons produce very large percentages that are arithmetically correct but not informative.

**Read the dataset README before relying on a dataset.** Each records the caveats found during ingestion, including unresolved ambiguities in the source files.

The principal limitation of v0.5.0 is size and coverage. Nine datasets cannot characterise the space of bioprocess perturbations, and the distribution reflects what is openly available rather than what is industrially important: product-quality attributes appear in one dataset, scale is varied in one, and several industrially central hosts are absent. A systematic search for datasets in which one process is measured at two or more scales found essentially nothing in the public domain, which is itself a finding: the comparability question on which scale-up depends cannot presently be studied with open data.

The atlas is designed to grow by contribution. Adding a dataset requires its unchanged files, a README, and one ingester function; the build rejects variables outside the vocabulary, so the schema cannot drift silently as datasets accumulate.

---

## Code Availability

All code is at https://github.com/anupama29k/Bioprocess-Atlas under the MIT licence; the data retain the licences of their source publications. The release archived at doi:10.5281/zenodo.24125547 contains the code, the unchanged source files and the built tables. The phenotype engine is vendored from the Bioprocess Analyzer (https://github.com/anupama29k/bioprocessanalyzer) at commit 2439bf9 and identified in every release manifest. An interactive explainer is published at https://anupama29k.github.io/Bioprocess-Atlas/.

---

## References

To be formatted for the target journal. The eight source studies are cited in the Data Records table and in `release/datasets.csv`; the perturbation-atlas precedent in single-cell biology and the simulated bioprocess benchmark should be cited in the Background section.

---

## Notes for the author (remove before submission)

- Fill in affiliation, ORCID and correspondence address.
- Decide the venue. *Scientific Data* expects a dataset of broad reuse value and has an article-processing charge; *Data in Brief* is cheaper and faster but lower profile; a bioprocess-engineering journal would reach the intended readers most directly. Check current charges before committing.
- Consider contacting the authors of the source datasets before submission. It is a courtesy, it sometimes yields corrections or additional data, and some may wish to be acknowledged.
- The Background section should cite the single-cell perturbation-atlas literature and the simulated bioprocess benchmark explicitly; the draft refers to both without citations.
- A figure would strengthen it: the pipeline schematic from the explainer, and a panel showing one comparison from raw time course through to effect.
- Regenerate every number in this draft from the release at submission time.
