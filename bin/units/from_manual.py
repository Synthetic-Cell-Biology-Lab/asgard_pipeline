#!/usr/bin/env python3

import sys
import pandas as pd
from Bio import SeqIO
import ast

############################################
# INPUTS
############################################

protein_file = sys.argv[1]
genome_file = sys.argv[2]
fasta_file = sys.argv[3]

ids_out = sys.argv[4]
fasta_out = sys.argv[5]
csv_out = sys.argv[6]

tax_level = sys.argv[7]
tax_filter = sys.argv[8]
annotation_filter = sys.argv[9:]

############################################
# Normalize parameters
############################################

if tax_level == "None":
    tax_level = None
if tax_filter == "None":
    tax_filter = None
if annotation_filter == "None":
    annotation_filter = None

############################################
# Load tables
############################################

protein_df = pd.read_csv(protein_file)
genome_df = pd.read_csv(genome_file)

str_cols = protein_df.select_dtypes(include="object").columns
protein_df[str_cols] = protein_df[str_cols].apply(lambda col: col.str.strip())

str_cols2 = genome_df.select_dtypes(include="object").columns
genome_df[str_cols2] = genome_df[str_cols2].apply(lambda col: col.str.strip())

print(f"[DEBUG] Protein table loaded: {protein_df.shape[0]} rows", file=sys.stderr)
print(f"[DEBUG] Genome table loaded: {genome_df.shape[0]} rows", file=sys.stderr)

############################################
# Taxonomic filtering
############################################

if tax_level and tax_filter:

    print(
        f"[DEBUG] Applying taxonomic filter: {tax_level} == {tax_filter}",
        file=sys.stderr
    )

    before = genome_df.shape[0]

    genome_df = genome_df[genome_df[tax_level] == tax_filter]

    after = genome_df.shape[0]

    print(
        f"[DEBUG] Genomes before filter: {before} | after filter: {after}",
        file=sys.stderr
    )

############################################
# Merge protein and genome tables
############################################

before_merge = protein_df.shape[0]

merged = protein_df.merge(genome_df, on="genome_file", how="inner")

print(
    f"[DEBUG] Proteins before merge: {before_merge} | after merge: {merged.shape[0]}",
    file=sys.stderr
)

############################################
# Annotation distribution before filtering
############################################

if "Manual_annotation" in merged.columns:

    print("[DEBUG] Manual annotation counts BEFORE filtering:", file=sys.stderr)

    counts = merged["Manual_annotation"].fillna("unknown").value_counts()

    for ann, count in counts.items():
        print(f"    {ann}: {count}", file=sys.stderr)
############################################
# Annotation filtering
############################################

if annotation_filter:

    annotations = [x.strip() for x in annotation_filter]

    print(f"[DEBUG] Applying annotation filter: {annotations}", file=sys.stderr)

    before = merged.shape[0]

    merged = merged[merged["Manual_annotation"].isin(annotations)]

    after = merged.shape[0]

    print(
        f"[DEBUG] Proteins before annotation filter: {before} | after: {after}",
        file=sys.stderr
    )
############################################
# Annotation distribution after filtering
############################################

if "Manual_annotation" in merged.columns:

    print("[DEBUG] Manual annotation counts AFTER filtering:", file=sys.stderr)

    counts = merged["Manual_annotation"].fillna("unknown").value_counts()

    for ann, count in counts.items():
        print(f"    {ann}: {count}", file=sys.stderr)

############################################
# Save filtered protein table
############################################

merged.to_csv(csv_out, index=False)

print(f"[DEBUG] Filtered protein CSV written: {csv_out}", file=sys.stderr)

############################################
# Extract protein IDs
############################################

protein_ids = set(merged["locus_tag"])

print(f"[DEBUG] Total protein IDs retained: {len(protein_ids)}", file=sys.stderr)

with open(ids_out, "w") as f:
    for pid in protein_ids:
        f.write(f"{pid}\n")

print(f"[DEBUG] Protein ID list written: {ids_out}", file=sys.stderr)

############################################
# Map locus_tag → manual annotation
############################################

annotation_map = dict(
    zip(merged["locus_tag"], merged["Manual_annotation"].fillna("unknown"))
)

############################################
# Extract FASTA and rewrite headers
############################################

records = []

for record in SeqIO.parse(fasta_file, "fasta"):

    locus = record.id

    if locus in protein_ids:

        annotation = annotation_map.get(locus, "unknown")

        record.id = locus
        record.description = f"{locus} {annotation}"

        records.append(record)

SeqIO.write(records, fasta_out, "fasta")

print(f"[DEBUG] FASTA sequences written: {len(records)}", file=sys.stderr)
print(f"[DEBUG] Output FASTA: {fasta_out}", file=sys.stderr)