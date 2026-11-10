# CHO DG44 IgG producer: nutrient, pH and CO2 stress vs standard process (measured)

**Source data:** Solle D, Schellenberg J (2021). *CHO stress induced product production.* Research data repository of Leibniz Universität Hannover. https://doi.org/10.25835/0066265
**Related paper:** Schellenberg J et al. (2022). *Stress-induced increase of monoclonal antibody production in CHO cells.* Engineering in Life Sciences 22(5):427–436. https://doi.org/10.1002/elsc.202100062
**License:** CC0 1.0 (public domain dedication). The 21 CSV files are reproduced unchanged except that the µ sign in the unit row was restored (it was mis-encoded in the download).

## Design

12-day fed-batch cultures of an IgG-producing CHO DG44 line, sampled daily. Three series, each with its own standard-process reference:

| series (atlas campaign) | vessel | arms vs "standard process" | vessels per arm |
|---|---|---|---|
| ambr15 · nutrients | ambr15, 15 mL | glutamate 1.2 g/L; glutamate 2 g/L; glucose 12 g/L (overfeeding); decreased feed | 2 (reference 1) |
| ambr15 · pH | ambr15, 15 mL | pH 7.3 from day 0; pH 7.3 from day 5 | 2 (reference 1) |
| shake flask | 500 mL shake flask | nutrient limitation; CO2 limitation | 3 (reference 1) |

Measured each day: viable cell concentration, viability, mean cell diameter, IgG, glucose, lactate, glutamine and glutamate.

## Caveats

- **Each reference is a single vessel**, so the atlas reports differences for these contrasts but no p-values or effect sizes.
- In the glutamate arms, glutamine is close to zero at day 0 (0.04–0.18 g/L, against 0.72 g/L in the reference). Glutamate therefore appears to replace glutamine rather than being added on top. The condition labels say "low glutamine" for this reason. This is read from the data, not from the paper's methods.
- The nutrient series has no samples on days 1–4. The shake-flask nutrient arms are missing days 1–3, and the shake-flask reference is missing days 1, 2 and 4.
- IgG at day 0 in the nutrient series is about 0.25 g/L, which is probably carry-over from the inoculum. Specific productivity (qP) uses the product made after the first sample, so this offset cancels out.
- The standard process's pH set point and the exact nature of the "nutrient limitation" and "CO2 limitation" arms are described in the paper, not in the data files.
