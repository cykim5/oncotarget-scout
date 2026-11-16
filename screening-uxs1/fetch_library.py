"""Fetch an approved-drug repurposing library from the ChEMBL API.

Source: ChEMBL web services, molecule endpoint, filtered to
max_phase=4 (approved) and molecule_type=Small molecule.
Filters applied here: has a canonical SMILES, has a preferred name,
150 <= full molecular weight <= 600 Da. From the filtered pool a
deterministic random sample of N_LIBRARY compounds is drawn
(seed fixed, so the library is reproducible).

Output: library.csv with columns chembl_id, pref_name,
canonical_smiles, full_mwt. No compound identities are invented;
every row comes from the ChEMBL record.
"""
import csv
import json
import os
import random
import time
import urllib.parse
import urllib.request

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "library.csv")

N_LIBRARY = 100
SEED = 42
MW_MIN = 150.0
MW_MAX = 600.0
PAGE_SIZE = 1000

API = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"


def fetch_page(offset):
    params = {
        "max_phase": 4,
        "molecule_type": "Small molecule",
        "limit": PAGE_SIZE,
        "offset": offset,
    }
    url = API + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except Exception as e:
            print(f"  fetch offset={offset} attempt {attempt + 1} failed: {e}", flush=True)
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"ChEMBL fetch failed at offset {offset}")


pool = []
offset = 0
while True:
    data = fetch_page(offset)
    mols = data.get("molecules", [])
    if not mols:
        break
    for m in mols:
        struct = m.get("molecule_structures") or {}
        smiles = struct.get("canonical_smiles")
        props = m.get("molecule_properties") or {}
        try:
            mw = float(props.get("full_mwt"))
        except (TypeError, ValueError):
            continue
        name = m.get("pref_name")
        if not smiles or not name:
            continue
        if not (MW_MIN <= mw <= MW_MAX):
            continue
        pool.append({
            "chembl_id": m.get("molecule_chembl_id"),
            "pref_name": name,
            "canonical_smiles": smiles,
            "full_mwt": mw,
        })
    total = data.get("page_meta", {}).get("total_count")
    print(f"  fetched offset {offset}, pool so far {len(pool)} (total {total})", flush=True)
    offset += len(mols)
    if total is not None and offset >= total:
        break

print(f"Filtered pool: {len(pool)} approved small molecules", flush=True)
random.seed(SEED)
random.shuffle(pool)
library = pool[:N_LIBRARY]
library.sort(key=lambda r: r["chembl_id"])

with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["chembl_id", "pref_name", "canonical_smiles", "full_mwt"])
    w.writeheader()
    for row in library:
        w.writerow(row)

print(f"Wrote {len(library)} compounds to {OUT}", flush=True)
