"""Regenerate docs/demo_data.json (and the data embedded in docs/index.html) from release/.

    python -m atlas.build && python scripts/make_demo_data.py

The two worked examples in the stepper (ΔaceA and W-68 at 37 °C) are kept from the
existing demo_data.json; everything else is read from the release tables.
"""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
REL, DOCS = ROOT / "release", ROOT / "docs"
CAMPAIGN_NAMES = {
    "hans2020": "E. coli Keio screen",
    "uhlendorff2026 · W-68": "yeast · W-68",
    "uhlendorff2026 · W-34/70": "yeast · W-34/70",
    "CHO-K1 · Uhlendorff": "CHO-K1 · process changes",
    "schellenberg2021 · ambr15 · nutrients": "CHO DG44 · nutrients (ambr15)",
    "schellenberg2021 · ambr15 · pH": "CHO DG44 · pH (ambr15)",
    "schellenberg2021 · shake flask": "CHO DG44 · stress (shake flask)",
}


def num(x):
    return None if x is None or (isinstance(x, float) and np.isnan(x)) else round(float(x), 4)


def main():
    man = json.loads((REL / "manifest.json").read_text())
    runs = pd.read_csv(REL / "runs.csv")
    pheno = pd.read_csv(REL / "phenotypes.csv").set_index("run_id")
    eff = pd.read_csv(REL / "effects.csv")
    con = pd.read_csv(REL / "contrasts.csv")
    dsets = pd.read_csv(REL / "datasets.csv")
    old = json.loads((DOCS / "demo_data.json").read_text())

    def arm_runs(row):
        c = con[(con.campaign == row.campaign) & (con.perturbation == row.perturbation)]
        if len(c):
            return c.iloc[0].runs.split(";"), c.iloc[0].reference_runs.split(";")
        g = runs[runs.campaign == row.campaign]
        return (g[g.condition == row.condition].run_id.tolist(), g[g.condition == row.reference].run_id.tolist())

    out = []
    for row in eff[~eff.confounded_by_feed].itertuples():
        a, b = arm_runs(row)
        vals = lambda ids: [[r, num(pheno.at[r, row.phenotype])] for r in ids
                            if r in pheno.index and pd.notna(pheno.at[r, row.phenotype])]
        out.append(dict(campaign=CAMPAIGN_NAMES.get(row.campaign, row.campaign), perturbation=row.perturbation,
                        type=row.perturbation_type, phenotype=row.phenotype, unit=row.unit, mean=num(row.mean),
                        ref=num(row.ref_mean), pct=num(row.pct_change), d=num(row.effect_size_d), p=num(row.p_welch),
                        n=int(row.n), rn=int(row.ref_n), a=vals(a), b=vals(b)))
    data = dict(
        manifest=dict(version=man["schema_version"], datasets=man["datasets"], runs=man["runs"],
                      conditions=man["conditions"], contrasts=man["contrasts"], organisms=man["organisms"],
                      observations=man["observations"], engine=man["engine"], checksum=man["observations_sha256_16"],
                      variables=len(man["variables"]), phenotypes=man["phenotypes"]),
        datasets=json.loads(dsets.to_json(orient="records")),
        effects=out, examples=old["examples"])
    (DOCS / "demo_data.json").write_text(json.dumps(data, ensure_ascii=False))
    html = (DOCS / "index.html").read_text()
    blob = json.dumps(data, ensure_ascii=False).replace("</", "<\\/")
    html = re.sub(r'(<script type="application/json" id="data">).*?(</script>)',
                  lambda m: m.group(1) + blob + m.group(2), html, count=1, flags=re.S)
    (DOCS / "index.html").write_text(html)
    print(f"{len(out)} effects, {man['runs']} runs, {man['datasets']} datasets")


if __name__ == "__main__":
    main()
