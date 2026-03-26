#!/usr/bin/env python3

"""
Extract sequences and CSV rows matching protein IDs.

Inputs:
- snakemake.input.fasta
- snakemake.input.protein_file
- snakemake.input.protein_ids

Outputs:
- snakemake.output.outfasta
- snakemake.output.protein_csv
"""

from Bio import SeqIO
import pandas as pd
import sys


# -------------------------------
# Snakemake I/O
# -------------------------------

FASTA_FILE = snakemake.input.fasta
CSV_FILE = snakemake.input.protein_file
IDS_FILE = snakemake.input.protein_ids
GENOME_FILE = snakemake.input.genome_file
PROTEIN = snakemake.params.protein_name

OUT_FASTA = snakemake.output.outfasta
OUT_CSV = snakemake.output.protein_csv
REMOVE_HYPOTHETICALS = snakemake.params.remove_hypotheticals

# -------------------------------
# Load Protein IDs
# -------------------------------

print("📥 Loading protein IDs...")

with open(IDS_FILE) as f:
    protein_ids = set(line.strip() for line in f if line.strip())

if not protein_ids:
    sys.exit("❌ No protein IDs found. Stopping.")

print(f"🔢 {len(protein_ids)} protein IDs loaded.")


# -------------------------------
# Filter FASTA
# -------------------------------

print("🧬 Filtering FASTA sequences...")

matched_records = []
skipped_hypotheticals = 0

for record in SeqIO.parse(FASTA_FILE, "fasta"):
    if record.id not in protein_ids:
        continue

    header = record.description.lower()

    if REMOVE_HYPOTHETICALS:
        # Remove the ID from header
        # record.description looks like:
        # "ABC_001 hypothetical protein"
        header_without_id = header.replace(record.id.lower(), "", 1).strip()

        # Check if hypothetical
        if "hypothetical" in header_without_id:
            skipped_hypotheticals += 1
            continue

    matched_records.append(record)

if not matched_records:
    print("⚠️ No matching FASTA records found.")

SeqIO.write(matched_records, OUT_FASTA, "fasta")

print(f"✅ FASTA written → {OUT_FASTA}")
if REMOVE_HYPOTHETICALS:
    print(f"🚫 Skipped {skipped_hypotheticals} hypothetical proteins.")


# -------------------------------
# Filter CSV
# -------------------------------

print("📊 Filtering CSV file...")

df = pd.read_csv(CSV_FILE)
gf = pd.read_csv(GENOME_FILE)

merged = pd.merge(gf, df, on="genome_file", how="inner")
# Assumes column name is 'protein'
if "locus_tag" not in df.columns:
    sys.exit("❌ CSV does not contain a 'locus_tag' column.")

filtered_df = merged[merged["locus_tag"].isin(protein_ids)]

if filtered_df.empty:
    print("⚠️ No matching rows in CSV.")

filtered_df['Manual_annotation'] = PROTEIN



filtered_df.to_csv(OUT_CSV, index=False)

print(f"✅ CSV written → {OUT_CSV}")

print("🎉 Extraction complete.")