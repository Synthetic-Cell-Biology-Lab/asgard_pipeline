#!/usr/bin/env python3

import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ---------------------------
# Files
# ---------------------------
dali_file = (
    "/home/anirudh/asgard_pipeline/database/structures/SepF/comparison/s001A.dali.tsv"
)

reference_fasta = "/home/anirudh/asgard_pipeline/database/collated/Version1/filtered/85comp10con/fasta/v1_cp85_con10.fasta"

subjects_faa = "/home/anirudh/asgard_pipeline/database/structures/SepF/comparison/dali_subjects.faa"

blast_db = "v1_cp85_con10"

blast_out = "/home/anirudh/asgard_pipeline/database/structures/SepF/comparison/dali_vs_reference.tsv"
filtered_blast = "/home/anirudh/asgard_pipeline/database/structures/SepF/comparison/dali_vs_reference.filtered.tsv"

mapping_out = "/home/anirudh/asgard_pipeline/database/structures/SepF/comparison/subject_accession_map.tsv"

# ---------------------------
# User thresholds
# ---------------------------
MIN_PIDENT = 100.0
MIN_QCOV = 70.0
MAX_EVALUE = 1e-5

# ---------------------------
# Read DALI output
# ---------------------------
print("Reading DALI results...")

df = pd.read_csv(dali_file, sep="\t", comment="#")

subjects = df[["sbjct", "sbjct-sequence"]].drop_duplicates(subset="sbjct")

records = [
    SeqRecord(Seq(row["sbjct-sequence"]), id=str(row["sbjct"]), description="")
    for _, row in subjects.iterrows()
]

SeqIO.write(records, subjects_faa, "fasta")

print(f"Wrote {len(records)} subject sequences")

# ---------------------------
# Build BLAST database
# ---------------------------
print("Building BLAST database...")

subprocess.run(
    [
        "makeblastdb",
        "-in",
        reference_fasta,
        "-dbtype",
        "prot",
        "-out",
        blast_db,
    ],
    check=True,
)

# ---------------------------
# Run BLASTP
# ---------------------------
print("Running BLASTP...")

subprocess.run(
    [
        "blastp",
        "-query",
        subjects_faa,
        "-db",
        blast_db,
        "-out",
        blast_out,
        "-outfmt",
        "6 qseqid sseqid pident length qlen slen bitscore evalue",
        "-evalue",
        str(MAX_EVALUE),
        "-max_hsps",
        "1",
    ],
    check=True,
)

# ---------------------------
# Read BLAST results
# ---------------------------
print("Reading BLAST output...")

hits = pd.read_csv(
    blast_out,
    sep="\t",
    header=None,
    names=[
        "subject",
        "accession",
        "pident",
        "length",
        "qlen",
        "slen",
        "bitscore",
        "evalue",
    ],
)

if hits.empty:
    print("No BLAST hits found.")
    exit()

# ---------------------------
# Calculate coverage
# ---------------------------
hits["qcov"] = 100 * hits["length"] / hits["qlen"]

# ---------------------------
# Apply filters
# ---------------------------
filtered = hits[
    (hits["pident"] >= MIN_PIDENT)
    & (hits["qcov"] >= MIN_QCOV)
    & (hits["evalue"] <= MAX_EVALUE)
].copy()

print(f"Hits before filtering: {len(hits)}")
print(f"Hits after filtering: {len(filtered)}")

if filtered.empty:
    print("No hits passed filtering.")
    exit()

# ---------------------------
# Keep best hit per subject
# ---------------------------
filtered = filtered.sort_values(["subject", "bitscore"], ascending=[True, False])

filtered.to_csv(filtered_blast)

best_hits = filtered.drop_duplicates(subset="subject", keep="first")

# ---------------------------
# Save mapping
# ---------------------------
best_hits[["subject", "accession"]].to_csv(
    mapping_out,
    sep="\t",
    index=False,
)

print(f"Saved {len(best_hits)} mappings to:")
print(mapping_out)

print("\nTop mappings:")
print(best_hits[["subject", "accession"]].head())
