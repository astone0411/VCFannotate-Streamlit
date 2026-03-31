import streamlit as st
import tempfile
import json
import sqlite3
from pathlib import Path

# Import your annotators
from vcfAnnotateCloud import annotate_vcf_to_json
from clinvar_lookup import get_clinsig_pure_python

# Path to your SQLite DB (can live in your GitHub repo)
DB_PATH = Path("lims.db")

# ---------------------------------------
# Initialize SQLite database
# ---------------------------------------
def init_db():
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    cur.execute("""
        CREATE TABLE IF NOT EXISTS variants (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            chrom TEXT,
            pos INTEGER,
            ref TEXT,
            alt TEXT,
            gene_symbol TEXT,
            transcript_id TEXT,
            hgvsc TEXT,
            hgvsp TEXT,
            consequence TEXT,
            cancervar TEXT,
            opai TEXT,
            clinvar_clinsig TEXT,
            clinvar_review TEXT,
            clinvar_variation_id TEXT,
            clinvar_allele_id TEXT
        )
    """)
    conn.commit()
    conn.close()

init_db()

# ---------------------------------------
# Insert annotated variants into SQLite
# ---------------------------------------
def insert_variants(records):
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()

    for r in records:
        cur.execute("""
            INSERT INTO variants (
                chrom, pos, ref, alt,
                gene_symbol, transcript_id,
                hgvsc, hgvsp, consequence,
                cancervar, opai,
                clinvar_clinsig, clinvar_review,
                clinvar_variation_id, clinvar_allele_id
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            r["chrom"], r["pos"], r["ref"], r["alt"],
            r["gene_symbol"], r["transcript_id"],
            r["hgvsc"], r["hgvsp"], r["consequence"],
            r["cancervar"], r["opai"],
            r["clinvar_clinsig"], r["clinvar_review"],
            r["clinvar_variation_id"], r["clinvar_allele_id"]
        ))

    conn.commit()
    conn.close()

# ---------------------------------------
# STREAMLIT UI
# ---------------------------------------
st.title("VCF Annotation Portal (VEP + CancerVar + ClinVar)")

uploaded_file = st.file_uploader("Upload a .vcf file", type=["vcf"])

if uploaded_file:
    st.write(f"Uploaded: {uploaded_file.name}")

    if st.button("Run Annotation"):
        with st.spinner("Annotating variants... this may take a moment"):
            # Save uploaded VCF to temp file
            with tempfile.NamedTemporaryFile(delete=False, suffix=".vcf") as tmp:
                tmp.write(uploaded_file.read())
                tmp_path = tmp.name

            # Output JSON path
            out_json = tmp_path + ".json"

            # Run your VEP + CancerVar annotator
            annotate_vcf_to_json(tmp_path, out_json)

            # Load JSON
            with open(out_json) as f:
                records = json.load(f)

            # Add ClinVar annotation to each record
            for r in records:
                clin = get_clinsig_pure_python(
                    r["chrom"], r["pos"], r["ref"], r["alt"]
                )

                r["clinvar_clinsig"] = clin.get("clinical_significance") if clin else None
                r["clinvar_review"] = clin.get("review_status") if clin else None
                r["clinvar_variation_id"] = clin.get("variation_id") if clin else None
                r["clinvar_allele_id"] = clin.get("allele_id") if clin else None

            # Insert into SQLite
            insert_variants(records)

        st.success(f"Annotated {len(records)} variants and saved to database")

        st.subheader("Preview")
        st.dataframe(records)
