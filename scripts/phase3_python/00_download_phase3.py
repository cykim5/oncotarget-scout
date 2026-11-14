"""Phase 3 downloads for OncoTarget Scout.

Fetches, with verification where the source offers it:
  1. A lung tumor single cell dataset from CZ CELLxGENE Discover:
     dataset 4ca422ae-e4ad-4884-9073-7c3799cd9143 (35,954 cells,
     Chromium 10x, lung adenocarcinoma and squamous cell lung
     carcinoma, 4 donors), from the collection "Transcriptome
     Analysis of Archived Tumor Tissues by Visium, GeoMx DSP, and
     Chromium Methods Reveals Inter- and Intra-Patient
     Heterogeneity" (bioRxiv 2024.11.01.621259). Chosen because its
     author cell type labels include an explicit malignant cell
     compartment, which is what the single cell layer needs; the
     other LUAD datasets found in the CELLxGENE scan either label
     tumor cells only as lung epithelial cells or are multi-million
     cell atlases too large for this environment. The CELLxGENE
     download URL publishes no checksum, so the file is verified
     by size (502,582,278 bytes) and by opening it with anndata
     in the processing script. Recorded in results/run_log.md.
  2. DepMap Public 20Q2 release files from the Broad's figshare
     deposit (article 12280541, published 2020-06-04), for the
     time split check: Achilles_gene_effect.csv (CERES gene
     effect scores; Chronos did not exist yet in 20Q2),
     CCLE_mutations.csv (cell line mutation calls),
     sample_info.csv (cell line annotations), and
     Achilles_common_essentials.csv (common essential gene list).
     File md5s are verified against the figshare-supplied values.

Outputs (data/raw is gitignored; nothing here is committed):
  data/raw/cellxgene/luad_sc_4ca422ae.h5ad
  data/raw/depmap20q2/{Achilles_gene_effect,CCLE_mutations,
                       sample_info,Achilles_common_essentials}.csv

Run from the repository root with the project venv:
  .venv/bin/python scripts/phase3_python/00_download_phase3.py
"""

import hashlib
import json
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CELLXGENE = ROOT / "data" / "raw" / "cellxgene"
DEPMAP20Q2 = ROOT / "data" / "raw" / "depmap20q2"
for d in (CELLXGENE, DEPMAP20Q2):
    d.mkdir(parents=True, exist_ok=True)

CX_DATASET = "4ca422ae-e4ad-4884-9073-7c3799cd9143"
CX_URL = f"https://datasets.cellxgene.cziscience.com/{CX_DATASET}.h5ad"
CX_SIZE = 502582278

FIGSHARE_ARTICLE = 12280541  # DepMap Public 20Q2
DEPMAP_FILES = ["Achilles_gene_effect.csv", "CCLE_mutations.csv",
                "sample_info.csv", "Achilles_common_essentials.csv"]


def fetch(url, dest, expected_md5=None, expected_size=None):
    """Download url to dest (via a temp file, verified, then renamed).

    Size and md5 are checked on the temp file BEFORE it replaces the
    destination, so a truncated download never lands under the final
    name (the convention from scripts/phase2_python/01_download_phase2.py).
    """
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    req = urllib.request.Request(url, headers={"User-Agent": "oncotarget-scout/phase3"})
    with urllib.request.urlopen(req, timeout=900) as resp, open(tmp, "wb") as fh:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
    if expected_size is not None and tmp.stat().st_size != expected_size:
        raise SystemExit(f"size mismatch for {dest.name}: got "
                         f"{tmp.stat().st_size}, expected {expected_size}")
    if expected_md5 is not None and md5_of(tmp) != expected_md5:
        raise SystemExit(f"md5 mismatch for {dest.name}")
    tmp.replace(dest)
    return dest.stat().st_size


def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


print("[1/2] CELLxGENE single cell dataset (LUAD + LUSC, 35,954 cells) ...")
size = fetch(CX_URL, CELLXGENE / "luad_sc_4ca422ae.h5ad", expected_size=CX_SIZE)
print(f"      luad_sc_4ca422ae.h5ad: {size} bytes, size verified "
      f"(no checksum published by the source)")

print("[2/2] DepMap Public 20Q2 via figshare ...")
with urllib.request.urlopen(
        f"https://api.figshare.com/v2/articles/{FIGSHARE_ARTICLE}/files?page_size=100",
        timeout=120) as resp:
    inventory = {f["name"]: f for f in json.load(resp)}
for name in DEPMAP_FILES:
    entry = inventory[name]
    dest = DEPMAP20Q2 / name
    size = fetch(entry["download_url"], dest,
                 expected_md5=entry["supplied_md5"],
                 expected_size=entry["size"])
    print(f"      {name}: {size} bytes, size and md5 verified")

print("Phase 3 downloads complete.")
