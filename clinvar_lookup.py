import gzip

# Try to import pytabix-lite
try:
    import tabixlite as tbx
    TABIX_AVAILABLE = True
except ImportError:
    TABIX_AVAILABLE = False


# ------------------------------------------------------------
# Parse INFO field
# ------------------------------------------------------------
def parse_info_field(info_str):
    info = {}
    for item in info_str.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True
    return info


# ------------------------------------------------------------
# FAST: pytabix-lite lookup (works on Windows + Streamlit Cloud)
# ------------------------------------------------------------
def get_clinsig_tabixlite(chrom, pos, ref, alt, vcf_path="clinvar.vcf.gz"):
    if not TABIX_AVAILABLE:
        return None

    chrom = chrom.replace("chr", "")
    pos = int(pos)

    try:
        tb = tbx.TabixFile(vcf_path)
        records = tb.query(chrom, pos, pos)
    except Exception:
        return None

    for fields in records:
        v_chrom = fields[0]
        v_pos = int(fields[1])
        v_ref = fields[3]
        v_alt = fields[4]
        info_str = fields[7]

        if v_pos != pos or v_ref != ref:
            continue

        alts = v_alt.split(",")
        if alt not in alts:
            continue

        idx = alts.index(alt)
        info = parse_info_field(info_str)

        clnsig = info.get("CLNSIG", "").split("|")
        clnrev = info.get("CLNREVSTAT", "").split("|")
        clnvid = info.get("CLNVID", "").split("|")
        allele_id = info.get("ALLELEID", "").split("|")

        return {
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "clinical_significance": clnsig[idx] if idx < len(clnsig) else None,
            "review_status": clnrev[idx] if idx < len(clnrev) else None,
            "variation_id": clnvid[idx] if idx < len(clnvid) else None,
            "allele_id": allele_id[idx] if idx < len(allele_id) else None
        }

    return None


# ------------------------------------------------------------
# SLOW: Pure-Python fallback (linear scan)
# ------------------------------------------------------------
def get_clinsig_pure_python(chrom, pos, ref, alt, vcf_path="clinvar.vcf.gz"):
    chrom = chrom.replace("chr", "")
    pos = int(pos)

    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            v_chrom, v_pos, _, v_ref, v_alt, _, _, info_str = fields[:8]

            if v_chrom != chrom:
                continue
            if int(v_pos) != pos:
                continue
            if v_ref != ref:
                continue

            alts = v_alt.split(",")
            if alt not in alts:
                continue

            idx = alts.index(alt)
            info = parse_info_field(info_str)

            clnsig = info.get("CLNSIG", "").split("|")
            clnrev = info.get("CLNREVSTAT", "").split("|")
            clnvid = info.get("CLNVID", "").split("|")
            allele_id = info.get("ALLELEID", "").split("|")

            return {
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "clinical_significance": clnsig[idx] if idx < len(clnsig) else None,
                "review_status": clnrev[idx] if idx < len(clnrev) else None,
                "variation_id": clnvid[idx] if idx < len(clnvid) else None,
                "allele_id": allele_id[idx] if idx < len(allele_id) else None
            }

    return None


# ------------------------------------------------------------
# Unified helper: try tabix-lite → fallback to pure Python
# ------------------------------------------------------------
def get_clinsig(chrom, pos, ref, alt, vcf_path="clinvar.vcf.gz"):
    """
    Recommended entry point.
    Uses fast pytabix-lite lookup if available,
    otherwise falls back to pure Python linear scan.
    """
    res = get_clinsig_tabixlite(chrom, pos, ref, alt, vcf_path=vcf_path)
    if res:
        return res
    return get_clinsig_pure_python(chrom, pos, ref, alt, vcf_path=vcf_path)


# ------------------------------------------------------------
# CLI test
# ------------------------------------------------------------
if __name__ == "__main__":
    result = get_clinsig("1", 45331556, "C", "T", vcf_path="clinvar.vcf.gz")
    print(result)
