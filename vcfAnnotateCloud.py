#!/usr/bin/env python3
import sys, json, time, requests
from typing import List, Tuple, Iterable, Optional

# -----------------------------
# VEP API CONFIG
# -----------------------------
VEP_URL = "https://rest.ensembl.org/vep/human/region"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}

# -----------------------------
# CancerVar API CONFIG
# -----------------------------
CANCERVAR_URL = "http://cancervar.wglab.org/api_new.php"


# ============================================================
# VCF PARSER
# ============================================================
def parse_vcf(path: str) -> Iterable[Tuple[str, int, str, str]]:
    """
    Yields (chrom, pos, ref, alt) for each ALT allele in the VCF.
    Splits multi-allelic lines into separate records.
    """
    with open(path, "rt") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            chrom, pos, _, ref, alts = cols[:5]
            pos = int(pos)
            for alt in alts.split(","):
                if alt == ".":
                    continue
                yield (chrom.replace("chr", ""), pos, ref, alt)


# ============================================================
# VEP REGION STRING BUILDER
# ============================================================
def to_vep_region_strings(batch: List[Tuple[str, int, str, str]]) -> List[str]:
    regions = []
    for chrom, pos, ref, alt in batch:
        regions.append(f"{chrom} {pos} {pos} {ref}/{alt}")
    return regions


# ============================================================
# VEP POST REQUEST
# ============================================================
def vep_post(regions: List[str],
             assembly: str = "GRCh38",
             include_mane: bool = True,
             pick_best: bool = True,
             retries: int = 3) -> Optional[list]:

    payload = {
        "variants": regions,
        "assembly": assembly,
        "hgvs": 1,
        "mane": 1 if include_mane else 0,
        "pick": 1 if pick_best else 0,
        "canonical": 1
    }

    for attempt in range(retries):
        try:
            r = requests.post(VEP_URL, headers=HEADERS, data=json.dumps(payload), timeout=60)
            if r.status_code == 429:
                time.sleep(1 + attempt)
                continue
            r.raise_for_status()
            return r.json()
        except Exception as e:
            if attempt == retries - 1:
                sys.stderr.write(f"[VEP] Failed after {retries} attempts: {e}\n")
                return None
            time.sleep(2)
    return None


# ============================================================
# PICK BEST TRANSCRIPT CONSEQUENCE
# ============================================================
def pick_transcript(tc_list: list) -> dict:
    if not tc_list:
        return {}
    for t in tc_list:
        if t.get("pick") == 1:
            return t
    return tc_list[0]


# ============================================================
# CANCERVAR QUERY
# ============================================================
def query_cancervar(chrom: str, pos: int, ref: str, alt: str, build="hg38") -> dict:
    url = (
        f"{CANCERVAR_URL}?queryType=position&build={build}"
        f"&chr={chrom}&pos={pos}&ref={ref}&alt={alt}"
    )

    headers = {
        "User-Agent": "Mozilla/5.0",
        "Accept": "application/json"
    }

    try:
        r = requests.get(url, headers=headers, timeout=20)
        r.raise_for_status()

        # If HTML returned, no data
        if "text/html" in r.headers.get("Content-Type", ""):
            return {"cancervar": None, "opai": None}

        data = r.json()

        # Your JSON has top-level keys:
        # "Cancervar" and "OPAI"
        cv_raw = data.get("Cancervar")
        opai_raw = data.get("OPAI")

        def clean(x):
            if x in ("NA", "", None):
                return None
            return x

        return {
            "cancervar": clean(cv_raw),
            "opai": clean(opai_raw)
        }

    except Exception:
        return {"cancervar": None, "opai": None}
# ============================================================
# PROCESS BATCH → DICT RECORDS
# ============================================================
def process_batch_to_dict(batch: List[Tuple[str,int,str,str]]) -> List[dict]:
    region_inputs = to_vep_region_strings(batch)
    res = vep_post(region_inputs)

    records = []

    if res is None:
        # VEP failed → return empty annotations
        for chrom, pos, ref, alt in batch:
            cv = query_cancervar(chrom, pos, ref, alt)
            records.append({
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "gene_symbol": None,
                "transcript_id": None,
                "hgvsc": None,
                "hgvsp": None,
                "consequence": None,
                "cancervar": cv["cancervar"],
                "opai": cv["opai"]
            })
        return records

    # VEP returned results
    for (chrom, pos, ref, alt), rec in zip(batch, res):
        tc_list = rec.get("transcript_consequences", []) or []
        tc = pick_transcript(tc_list)

        gene_symbol = tc.get("gene_symbol")
        transcript_id = tc.get("transcript_id")
        hgvsc = tc.get("hgvsc")
        hgvsp = tc.get("hgvsp")
        consequence = ",".join(tc.get("consequence_terms", []))

        # CancerVar annotation
        cv = query_cancervar(chrom, pos, ref, alt)

        records.append({
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "gene_symbol": gene_symbol,
            "transcript_id": transcript_id,
            "hgvsc": hgvsc,
            "hgvsp": hgvsp,
            "consequence": consequence,
            "cancervar": cv["cancervar"],
            "opai": cv["opai"]
        })

    return records


# ============================================================
# MAIN ANNOTATION FUNCTION
# ============================================================
def annotate_vcf_to_json(in_vcf: str, out_json: str, batch_size: int = 200):
    all_records = []
    batch: List[Tuple[str,int,str,str]] = []

    for var in parse_vcf(in_vcf):
        batch.append(var)
        if len(batch) >= batch_size:
            all_records.extend(process_batch_to_dict(batch))
            batch.clear()

    if batch:
        all_records.extend(process_batch_to_dict(batch))

    with open(out_json, "w") as f:
        json.dump(all_records, f, indent=2)

    print(f"\nAnnotation complete → {len(all_records)} variants")
    print(f"Output written to: {out_json}")


# ============================================================
# CLI ENTRY POINT
# ============================================================
if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: python annotate_simple_cloud.py input.vcf output.json\n")
        sys.exit(1)

    annotate_vcf_to_json(sys.argv[1], sys.argv[2])