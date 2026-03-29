import streamlit as st
import tempfile
import json
import sqlite3
from pathlib import Path

# Import your annotator
from vcfAnnotateCloud import annotate_vcf_to_json

# Path to your SQLite DB (can be inside your GitHub repo)
DB_PATH = Path("lims.db")

# ---------------------------------------
# Create table if not exists
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
            opai TEXT
        )
    """)
    conn.commit()
    conn.close()

init_db()

# ---------------------------------------
# Insert annotated variants
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
                cancervar, opai
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            r["chrom"], r["pos"], r["ref"], r["alt"],
            r["gene_symbol"], r["transcript_id"],
            r["hgvsc"], r["hgvsp"], r["consequence"],
            r["cancervar"], r["opai"]
        ))

    conn.commit()
    conn.close()

# ---------------------------------------
# STREAMLIT UI
# ---------------------------------------
st.title("VCF Annotation Portal (VEP + CancerVar)")

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

            # Run your annotator
            annotate_vcf_to_json(tmp_path, out_json)

            # Load JSON
            with open(out_json) as f:
                records = json.load(f)

            # Insert into SQLite
            insert_variants(records)

        st.success(f"Annotated {len(records)} variants and saved to database")

        st.subheader("Preview")
        st.dataframe(records)
