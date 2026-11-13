# K. phaffii (Pichia pastoris) fed-batch: hypoxia vs normoxia (measured)

**Source:** Sales-Vallverdú A, Gasset A, Requena-Moreno G, Valero F, Montesinos-Seguí JL, Garcia-Ortega X (2024). *Synergic kinetic and physiological control to improve the efficiency of Komagataella phaffii recombinant protein production bioprocesses.* Microbial Biotechnology 17(2):e14411, https://doi.org/10.1111/1751-7915.14411. Data: CORA Repositori de Dades de Recerca, https://doi.org/10.34810/data1118
**License:** CC BY 4.0. Files reproduced unchanged (`raw_data_microbial_biotech_v1.tab`, `Readme_original.txt`); the atlas ingester reshapes the side-by-side blocks into long format.

## Design

Glucose-limited fed-batch cultivations of constitutive (P_GAP) K. phaffii producer strains, sampled through the fed-batch phase. Each condition is one bioreactor; every hypoxic run has a matched normoxic run:

| product | comparison | set growth rate |
|---|---|---|
| Crl1 (*Candida rugosa* lipase 1) | hypoxia vs normoxia | 0.030, 0.065, 0.100, 0.120 h⁻¹ |
| proROL (*Rhizopus oryzae* lipase) | hypoxia vs normoxia | the study's reference rate |
| CalB (*Candida antarctica* lipase B) | hypoxia vs normoxia | the study's reference rate |

Hypoxia is held by controlling the respiratory quotient (the paper targets RQ ≈ 1.4), rather than by a DO set point. Measured per sample: biomass (g/L), extracellular lipolytic activity (AU/mL), ethanol, arabitol, and HPLC by-products.

## Caveats

- **One run per arm**, so the atlas reports differences without p-values.
- **Product activity is in assay units (AU/mL) and differs per enzyme.** Compare only within one product, never between Crl1, proROL and CalB. The atlas keeps it as `product_activity_AU_mL` for this reason.
- **µ = 0.065 h⁻¹ normoxia is left out.** In the source table its columns are shifted: "FB time", "X" and "Act" repeat the same pair of numbers, so time and biomass cannot be recovered. The hypoxic run at that rate is kept, but it has no reference.
- **µ = 0.030 h⁻¹ normoxia: the ethanol column is dropped.** It is identical to the hypoxic run's ethanol column, which is a copy rather than a measurement. Arabitol in that block is zero throughout and is kept as reported.
- **Some cells are written in scientific notation** (`5.8477922E-2`), mixed in among plain decimals in the same column. A reader that drops the exponent turns 0.058 g/L of arabitol into 5.85 g/L, which looks like a plausible measurement rather than an error. The ingester parses the exponent; a test asserts that arabitol rises monotonically in every hypoxic run, which is what fails if that ever regresses.
- Time is hours into the fed-batch phase, not from inoculation. The first row of each block is an all-zero placeholder and is skipped. Growth-phase phenotypes (batch µ, glucose depletion, biomass yield) do not apply: the atlas excludes them from these comparisons and reports the fed-batch growth rate instead.
