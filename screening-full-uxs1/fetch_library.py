"""Fetch every approved small molecule in ChEMBL, not a sample of 100.

Same source and filters as screening-uxs1/fetch_library.py, with the
sampling step removed, so this is the whole approved drug space inside
the molecular weight window.
"""
import csv
import json
import os
import time
import urllib.parse
import urllib.request

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "library.csv")

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
    req = urllib.request.Request(url, headers={"Accept": "application/json",
                                               "User-Agent": "research-script"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                return json.loads(resp.read().decode())
        except Exception as err:
            if attempt == 3:
                raise
            print(f"  retry {attempt + 1} after {err}", flush=True)
            time.sleep(5)


rows = []
seen = set()
offset = 0
total = None
while True:
    page = fetch_page(offset)
    if total is None:
        total = page["page_meta"]["total_count"]
        print(f"approved small molecules in ChEMBL: {total}", flush=True)
    molecules = page.get("molecules", [])
    if len(molecules) == 0:
        break
    for m in molecules:
        structures = m.get("molecule_structures") or {}
        smiles = structures.get("canonical_smiles")
        name = m.get("pref_name")
        props = m.get("molecule_properties") or {}
        mwt = props.get("full_mwt")
        cid = m.get("molecule_chembl_id")
        if smiles is None or name is None or mwt is None:
            continue
        if cid in seen:
            continue
        mwt = float(mwt)
        if mwt < MW_MIN or mwt > MW_MAX:
            continue
        seen.add(cid)
        rows.append({
            "chembl_id": cid,
            "pref_name": name,
            "canonical_smiles": smiles,
            "full_mwt": mwt,
        })
    offset = offset + PAGE_SIZE
    print(f"  fetched {offset} of {total}, kept {len(rows)}", flush=True)
    if offset >= total:
        break
    time.sleep(0.2)

rows.sort(key=lambda r: r["chembl_id"])
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["chembl_id", "pref_name",
                                      "canonical_smiles", "full_mwt"])
    w.writeheader()
    for r in rows:
        w.writerow(r)
print(f"wrote {len(rows)} compounds to {OUT}", flush=True)
