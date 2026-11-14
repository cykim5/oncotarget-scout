"""Phase 2 downloads for OncoTarget Scout.

Fetches, with verification where the source offers it:
  1. TCGA-LUAD somatic mutations for STK11, KEAP1, KRAS from the
     cBioPortal public API (PanCancer Atlas study
     luad_tcga_pan_can_atlas_2018). The Xena GDC hub mutation file
     (TCGA-LUAD.mutect2_snv.tsv.gz) returns HTTP 403 and the legacy
     Xena hub mutation files do too, so the cBioPortal API is the
     documented substitute; it serves the same underlying TCGA
     mutation calls. Recorded as a deviation in results/run_log.md.
  2. GTEx v8 median gene TPM by tissue (safety layer).
  3. DepMap Public 24Q4 release files from the Broad's figshare
     deposit (article 27993248): CRISPRGeneEffect.csv (Chronos),
     Model.csv (cell line annotations), OmicsSomaticMutations.csv
     (cell line mutation calls), AchillesCommonEssentialControls.csv
     (common essential gene list, for the pan-essential penalty).
     The DepMap portal download API returns an HTML human
     verification page to scripted requests, so the figshare API is
     the documented keyless route; file md5s are verified against
     the figshare-supplied values. Figshare hosting stopped after
     24Q4, so 24Q4 is the newest release available this way.

Outputs (data/raw is gitignored; nothing here is committed):
  data/raw/cbioportal/mutations_STK11.json (etc.), sample_list.json
  data/raw/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
  data/raw/depmap24q4/{CRISPRGeneEffect,Model,OmicsSomaticMutations,
                       AchillesCommonEssentialControls}.csv

Run from the repository root with the project venv:
  .venv/bin/python scripts/phase2_python/01_download_phase2.py
"""

import hashlib
import json
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CBIO = ROOT / "data" / "raw" / "cbioportal"
GTEX = ROOT / "data" / "raw" / "gtex"
DEPMAP = ROOT / "data" / "raw" / "depmap24q4"
for d in (CBIO, GTEX, DEPMAP):
    d.mkdir(parents=True, exist_ok=True)

GENES = {"STK11": 6794, "KEAP1": 9817, "KRAS": 3845}
STUDY = "luad_tcga_pan_can_atlas_2018"
PROFILE = f"{STUDY}_mutations"
SAMPLE_LIST = f"{STUDY}_all"

GTEX_URL = ("https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/"
            "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz")

FIGSHARE_ARTICLE = 27993248  # DepMap Public 24Q4
DEPMAP_FILES = ["CRISPRGeneEffect.csv", "Model.csv",
                "OmicsSomaticMutations.csv",
                "AchillesCommonEssentialControls.csv"]


def fetch(url, dest, expect_json=False, expected_md5=None, expected_size=None):
    """Download url to dest (via a temp file, verified, then renamed).

    The size and md5 are checked on the temp file BEFORE it replaces
    the destination, so a truncated download never lands under the
    final name. (The first run of this script caught a truncated
    OmicsSomaticMutations.csv exactly this way; see run_log.md.)
    """
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    req = urllib.request.Request(url, headers={"User-Agent": "oncotarget-scout/phase2"})
    with urllib.request.urlopen(req, timeout=600) as resp, open(tmp, "wb") as fh:
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
    if expect_json:
        json.loads(dest.read_text())  # raises if the body is not JSON
    return dest.stat().st_size


def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


print("[1/3] cBioPortal mutations (TCGA-LUAD PanCancer Atlas) ...")
for symbol, entrez in GENES.items():
    url = (f"https://www.cbioportal.org/api/molecular-profiles/{PROFILE}/mutations"
           f"?sampleListId={SAMPLE_LIST}&entrezGeneId={entrez}&projection=SUMMARY")
    dest = CBIO / f"mutations_{symbol}.json"
    size = fetch(url, dest, expect_json=True)
    n = len(json.loads(dest.read_text()))
    print(f"      {symbol}: {n} mutation records, {size} bytes")
url = f"https://www.cbioportal.org/api/sample-lists/{SAMPLE_LIST}/sample-ids"
size = fetch(url, CBIO / "sample_list.json", expect_json=True)
print(f"      sample list: {len(json.loads((CBIO / 'sample_list.json').read_text()))} samples")

print("[2/3] GTEx v8 median TPM ...")
size = fetch(GTEX_URL, GTEX / "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz")
print(f"      {size} bytes")

print("[3/3] DepMap Public 24Q4 via figshare ...")
with urllib.request.urlopen(
        f"https://api.figshare.com/v2/articles/{FIGSHARE_ARTICLE}/files?page_size=100",
        timeout=120) as resp:
    inventory = {f["name"]: f for f in json.load(resp)}
for name in DEPMAP_FILES:
    entry = inventory[name]
    dest = DEPMAP / name
    size = fetch(entry["download_url"], dest,
                 expected_md5=entry["supplied_md5"],
                 expected_size=entry["size"])
    print(f"      {name}: {size} bytes, size and md5 verified")

print("Phase 2 downloads complete.")
