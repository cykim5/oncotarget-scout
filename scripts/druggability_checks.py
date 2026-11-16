"""Pull Open Targets tractability and ChEMBL compound counts for the top 5.

Two gaps the druggability note flagged:
  1. no Open Targets tractability score was retrieved
  2. "no known inhibitor found" was a search result, not a database answer
"""

import json
import time
import urllib.request

GENES = ["INTS6", "NFE2L2", "UXS1", "NCAPG2", "SPRR1B"]

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
AGENT = "research-script"


def get_json(url, data=None):
    headers = {"User-Agent": AGENT}
    body = None
    if data is not None:
        headers["Content-Type"] = "application/json"
        body = json.dumps(data).encode()
    req = urllib.request.Request(url, data=body, headers=headers)
    with urllib.request.urlopen(req, timeout=40) as r:
        return json.loads(r.read().decode())


def open_targets(symbol):
    """Ensembl id, UniProt accessions and the tractability flags."""
    found = get_json(OT_URL, {
        "query": '{ search(queryString:"%s", entityNames:["target"]) '
                 '{ hits { id name } } }' % symbol
    })
    ensembl = None
    for hit in found["data"]["search"]["hits"]:
        if hit["name"] == symbol:
            ensembl = hit["id"]
            break
    if ensembl is None:
        return None

    detail = get_json(OT_URL, {
        "query": '{ target(ensemblId:"%s") { approvedSymbol biotype '
                 'proteinIds { id source } '
                 'tractability { label modality value } } }' % ensembl
    })
    t = detail["data"]["target"]
    t["ensemblId"] = ensembl
    return t


def chembl_targets(uniprot):
    """ChEMBL targets for a UniProt accession, with compound counts."""
    out = []
    data = get_json(f"{CHEMBL}/target.json?target_components__accession={uniprot}&limit=20")
    for t in data.get("targets", []):
        cid = t["target_chembl_id"]
        acts = get_json(f"{CHEMBL}/activity.json?target_chembl_id={cid}&limit=1")
        mechs = get_json(f"{CHEMBL}/mechanism.json?target_chembl_id={cid}&limit=20")
        out.append({
            "chembl_id": cid,
            "pref_name": t.get("pref_name"),
            "target_type": t.get("target_type"),
            "activity_count": acts["page_meta"]["total_count"],
            "mechanisms": [
                {
                    "molecule_chembl_id": m.get("molecule_chembl_id"),
                    "action_type": m.get("action_type"),
                    "mechanism_of_action": m.get("mechanism_of_action"),
                    "max_phase": m.get("max_phase"),
                }
                for m in mechs.get("mechanisms", [])
            ],
        })
        time.sleep(0.3)
    return out


results = {}
for gene in GENES:
    print("...", gene, flush=True)
    row = {"symbol": gene}
    t = open_targets(gene)
    if t is None:
        row["error"] = "no Open Targets match"
        results[gene] = row
        continue

    row["ensembl_id"] = t["ensemblId"]
    row["biotype"] = t.get("biotype")
    row["uniprot"] = [p["id"] for p in t.get("proteinIds", [])
                      if p["source"] == "uniprot_swissprot"]
    # Keep only the flags that came back true, which is what the field means
    # by "tractable".
    row["tractability_true"] = [
        {"modality": x["modality"], "label": x["label"]}
        for x in t.get("tractability", []) if x["value"] == True
    ]
    row["tractability_all_count"] = len(t.get("tractability", []))

    row["chembl"] = []
    for acc in row["uniprot"]:
        row["chembl"].extend(chembl_targets(acc))
    results[gene] = row
    time.sleep(0.3)

out_path = "results/tables/druggability_open_targets_chembl.json"
with open(out_path, "w") as f:
    json.dump(results, f, indent=2)
print("written to", out_path)
