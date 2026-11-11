# HEK293 chemostat producing the Her1 extracellular domain (measured)

**Source:** Calzadilla N, Boggiano T, Vazquez E, Mulet A (2023). *Continuous culture of HEK293 producing the ECD-Her1: existence of steady-state multiplicity at the same external conditions. Cell culture and metabolic datasets.* Mendeley Data, https://doi.org/10.17632/t9rcjv5362.4 (v4). Related article: Data in Brief 50:109604.
**License:** CC BY 4.0.
**File:** `Continuous_culture_data_Multiplicity_HEK293.xlsx`, the dataset's own `Continuous culture data_Multiplicity_HEK293` workbook, unchanged. The companion LC-MS metabolite workbooks (76 metabolites) are left at the source.

## What this is, and what it is not

One 2 L HEK293 chemostat, sampled daily for about 85 days (64 samples), producing the extracellular domain of Her1. The dilution rate was changed between five held periods, labelled SS1–SS5 in the source:

| period | dilution rate (1/d) |
|---|---|
| SS1 | 0.45 |
| SS2 | 0.40 |
| SS3 | 0.35 |
| SS4 | 0.40 |
| SS5 | 0.45 |

**This is a single vessel.** The five periods are sequential states of one reactor, not independent reactors, so samples within a period are repeated measures rather than biological replicates. The atlas therefore stores this as **one run** and asserts **no comparison** between the periods: there are no effect rows for this dataset. The dilution rate is recorded as a time series (`dilution_rate_h`), so anyone who wants to look at the response can slice the run themselves with the design in view.

The study's own point is that the same dilution rate did not always give the same steady state — SS1 and SS5 are both at 0.45 1/d — which is a further reason not to treat the periods as interchangeable arms.

## Measurements ingested

Viable cell density, product titer, glucose, lactate, glutamine and the dilution rate. Concentrations reported in mmol/L are converted to g/L (glucose 180.16, lactate 90.08, glutamine 146.15 g/mol) and titer from µg/mL to g/L.

## Caveats

- **Continuous mode.** Batch and fed-batch phenotypes do not apply here: in a chemostat the growth rate is set by the dilution rate and product washes out rather than accumulating. The atlas does not compute specific growth rate, substrate depletion time, biomass yield, integral viable cell density or cell-specific productivity for continuous runs, and this is the dataset that prompted that rule.
- The source workbook also carries per-period means, standard deviations and coefficients of variation. Only the measured values are ingested; the atlas recomputes any summary it needs.
- Glucose falls to about 0.05 g/L at the end. In a chemostat that is a steady-state residual, not depletion.
- The 76-metabolite LC-MS tables in the source could be added later if the atlas gains a vocabulary for intracellular metabolites.
