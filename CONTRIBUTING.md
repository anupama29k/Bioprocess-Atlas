# Contributing a dataset

The atlas grows one dataset at a time. A dataset qualifies when it has:

- **time courses** from real bioreactor or shake-flask runs (not harvest-day values only),
- **at least one factor varied against a reference**, ideally with two or more runs per arm,
- **an open license** that allows redistribution (CC BY 4.0 or CC0 preferred).

## Steps

1. **Raw files.** Create `data/raw/<firstauthor><year>_<short>/`. Put the published files there unchanged, or reshaped only if they cannot be read otherwise. Add a `README.md` with:
   - the citation and DOI, and the license;
   - organism, strain(s), medium, vessel and working volume;
   - the conditions and which run is the reference for each;
   - known caveats (missing samples, feed starts, unit conversions you did).
2. **Ingester.** Add a function to `atlas/ingest.py` that returns a `Dataset` (see `atlas/schema.py`):
   - one `Run` per bioreactor, with `perturbation_type`, `dose` and `reference_condition`;
   - observations named only with variables from `VARIABLES` (the unit is part of the name, so convert values to it);
   - explicit `contrasts` if the reference is not obvious, for example same-day controls.

   Register the function in `INGESTERS`.
3. **Build and test.**
   ```bash
   python -m atlas.build
   python tests/test_atlas.py
   ```
   Add a test that checks one result the original paper reports (for example a growth rate or a final titer).
4. **Pull request.** Open a PR with the new folder, the ingester, the test and the rebuilt `release/`.

## New variables or phenotypes

If you need a variable that is not in the vocabulary, add it to `VARIABLES` in `atlas/schema.py`, with its unit and a one-line description, in the same PR. Phenotypes live in `atlas/phenotypes.py`. If a phenotype depends on the feed or another factor that differs between arms, set `confounded_by_feed`.
