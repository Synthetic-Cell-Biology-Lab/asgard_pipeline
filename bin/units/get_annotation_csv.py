#!/usr/bin/env python3

import sys
import pandas as pd
from Bio import SeqIO

############################################
# Inputs
############################################

fasta_file = sys.argv[1]
protein_csv = sys.argv[2]
cluster_csv = sys.argv[3]
output_csv = sys.argv[4]

############################################
# Extract protein IDs from FASTA
############################################

print("[INFO] Reading FASTA...")

protein_ids = {record.id for record in SeqIO.parse(fasta_file, "fasta")}

print(f"[INFO] Proteins in FASTA: {len(protein_ids)}")

############################################
# Load protein CSV
############################################

print("[INFO] Loading protein annotation table...")

df = pd.read_csv(protein_csv)

print(f"[INFO] Proteins in CSV: {len(df)}")

############################################
# Filter rows by FASTA proteins
############################################

filtered = df[df["locus_tag"].isin(protein_ids)]

print(f"[INFO] Proteins retained: {len(filtered)}")

############################################
# Keep only taxonomy columns
############################################

taxonomy_columns = [
    "locus_tag",
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]

missing_cols = [c for c in taxonomy_columns if c not in filtered.columns]

if missing_cols:
    raise ValueError(f"Missing required columns: {missing_cols}")

filtered = filtered[taxonomy_columns]

print("[INFO] Loading cluster table...")

cluster_df = pd.read_csv(cluster_csv)

print(f"[INFO] Cluster entries: {len(cluster_df)}")

merged = pd.merge(filtered, cluster_df, left_on='locus_tag', right_on="id", how="left")

############################################
# Write output
############################################

merged.to_csv(output_csv, index=False)

print(f"[INFO] Annotation table written → {output_csv}")