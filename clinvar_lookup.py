import gzip

def parse_info_field(info_str):
    info = {}
    for item in info_str.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True
    return info

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

if __name__ == "__main__":
    result = get_clinsig_pure_python(
        chrom="1",
        pos=45331556,
        ref="C",
        alt="T",
        vcf_path="clinvar.vcf.gz"
    )
    print(result)
