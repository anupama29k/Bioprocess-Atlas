# NISTCHO fed-batch: seven feeding strategies, with product N-glycan quality (measured)

**Source:** Schäpertöns V et al. (2026). *Effects of feeding strategies on culture performance and product quality in NISTCHO.* npj Systems Biology and Applications 12:67. Data: Zenodo, https://doi.org/10.5281/zenodo.18609042 (v3; the concept DOI 10.5281/zenodo.17046013 always resolves to the latest version).
**License:** CC BY 4.0.

## What is here

The Zenodo record is a single 7.3 GB archive that is mostly raw mass-spectrometry files (`data/mzml/`, 327 mzML files). Only the process and glycan result tables are reproduced here, unchanged:

- `01_ViCell_growth_data.csv`, `02_ViCell_growth_data_averaged.csv` — viable cell density, viability and mean cell diameter.
- `03_titer.csv`, `03_ViCell_titer.csv` — cNISTmAb titer (µg/mL) with the matching cell counts.
- `04_lactate.csv` — lactate (mM).
- `05_glucose_data.xlsx` — glucose (g/L).
- `glyco/<condition>_<replicate>_<hour>.csv` — the authors' CAFOG output: relative abundance of each glycoform pair, with Fuc, Hex and HexNAc counts, at two sampling times per run.

The mzML files, the intermediate R objects and the per-run CAFOG inputs are not copied; fetch the Zenodo archive for those.

## Design

NISTCHO (the open NIST reference cell line) producing cNISTmAb, in 50 mL TubeSpin tubes at 20 mL working volume, sampled roughly daily for about 11 days. Two separate fed-batch campaigns:

| campaign | strategy | what changes | vessels |
|---|---|---|---|
| FB4 | **STD** (atlas reference) | standard feed | 4 |
| FB4 | STD+ | standard feed + galactose | 4 |
| FB4 | LoG | low glucose | 4 |
| FB4 | LoG+ | low glucose + galactose | 4 |
| FB2 | HiF | high feed frequency | 3 |
| FB2 | HIP | high initial product feed | 3 |
| FB2 | HIP+ | high initial product feed + galactose | 3 |

In the source files the strategies are letter-coded: A = STD, B = STD+, C = LoG+, G = LoG, D = HiF, E = HIP, F = HIP+.

## Product quality

This is the atlas's first dataset with product **quality** attributes, not just titer. The glycoform table reports pairs of glycans (both heavy chains), so the atlas derives three standard attributes from the complete pairs (HexNAc = 8), renormalised:

- `galactosylation_index` — mean galactose residues per N-glycan, from the Hex count above the G0/G0 baseline of 6, halved to put it per glycan (0–2).
- `afucosylation_percent` — share of glycans lacking core fucose, from the Fuc count.
- `g0f_percent` — share of the product carrying G0F on both heavy chains.

These are stored as observations at each sampling time, so each run's phenotype is the value at harvest.

## What the atlas finds

Galactose supplementation roughly doubles galactosylation (index 0.22 → 0.46–0.51) and cuts the G0F share from 56 % to about 30 %, in all three galactose arms, with very tight replicates. Adding galactose to the standard feed does this without a titer penalty (−4 %, p = 0.48), whereas lowering glucose costs 58 % of the titer and 31 % of the integral viable cell density.

## Caveats

- **The FB2 campaign has no standard-feed arm.** HiF, HIP and HIP+ were run separately from the STD reference, so comparing them against STD would confound the strategy with the campaign. Those nine runs are in the atlas with their phenotypes, but they carry no reference and produce no effect rows. Use them for cross-run comparison at your own risk, or to compare the three FB2 strategies with each other.
- **Harvest times differ.** The LoG and LoG+ runs end at about 245 h, the others at about 263–265 h, and the final glycan sample is at 240 h rather than 264 h for those arms. The titer comparison for LoG and LoG+ is therefore slightly against them, since they had about 20 h less to accumulate product.
- Quality attributes exist at only two time points per run (about 120 h and harvest), so they are not a time course in the way the cell counts are.
- Lactate is converted from mM to g/L at 90.08 g/mol. Titer is converted from µg/mL to g/L.
- Incomplete glycoform pairs (HexNAc = 4, one glycan missing) are left out of the three derived attributes and the rest renormalised.
- These are TubeSpin tubes, not instrumented bioreactors: there is no controlled pH or dissolved oxygen, and no inline data.
