#!/usr/bin/env python3

from Bio import SeqIO
import pandas as pd
import sys
import argparse
from distutils.util import strtobool

# -------------------------------
# Input Handling
# -------------------------------


def get_args():
    """Handle both Snakemake and CLI usage."""

    if "snakemake" in globals():
        return {
            "FASTA_FILE": snakemake.input.fasta,
            "CSV_FILE": snakemake.input.protein_file,
            "IDS_FILE": snakemake.input.protein_ids,
            "GENOME_FILE": snakemake.input.genome_file,
            "OUT_FASTA": snakemake.output.outfasta,
            "OUT_CSV": snakemake.output.protein_csv,
            "PROTEIN": snakemake.params.protein_name,
            "REMOVE_HYPOTHETICALS": snakemake.params.remove_hypotheticals,
            "TAXON_LEVEL": snakemake.params.taxon_level,
            "TAXON_FILTER": snakemake.params.taxon_filter,
        }

    # CLI mode
    parser = argparse.ArgumentParser(
        description="Extract sequences and CSV rows matching protein IDs"
    )

    parser.add_argument("--fasta", required=True)
    parser.add_argument("--csv", required=True)
    parser.add_argument("--ids", required=True)
    parser.add_argument("--genome_file", required=True)
    parser.add_argument("--outfasta", required=True)
    parser.add_argument("--outcsv", required=True)
    parser.add_argument("--protein_name", required=True)
    parser.add_argument(
        "--remove_hypotheticals",
        help="Remove hypothetical proteins",
        type=lambda x: bool(strtobool(x)),
        default=False,
    )
    parser.add_argument("--taxon_level")
    parser.add_argument("--taxon_filter")

    args = parser.parse_args()

    return {
        "FASTA_FILE": args.fasta,
        "CSV_FILE": args.csv,
        "IDS_FILE": args.ids,
        "GENOME_FILE": args.genome_file,
        "OUT_FASTA": args.outfasta,
        "OUT_CSV": args.outcsv,
        "PROTEIN": args.protein_name,
        "REMOVE_HYPOTHETICALS": args.remove_hypotheticals,
        "TAXON_LEVEL": args.taxon_level,
        "TAXON_FILTER": args.taxon_filter,
    }


# Load args
cfg = get_args()

FASTA_FILE = cfg["FASTA_FILE"]
CSV_FILE = cfg["CSV_FILE"]
IDS_FILE = cfg["IDS_FILE"]
GENOME_FILE = cfg["GENOME_FILE"]
OUT_FASTA = cfg["OUT_FASTA"]
OUT_CSV = cfg["OUT_CSV"]
PROTEIN = cfg["PROTEIN"]
REMOVE_HYPOTHETICALS = cfg["REMOVE_HYPOTHETICALS"]
TAXON_LEVEL = cfg["TAXON_LEVEL"]
TAXON_FILTER = cfg["TAXON_FILTER"]


if bool(TAXON_LEVEL) ^ bool(TAXON_FILTER):
    raise ValueError(
        "Taxon filtering requires both --taxon_level and --taxon_filter inputs."
    )


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
        header_without_id = header.replace(record.id.lower(), "", 1).strip()

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

if "locus_tag" not in df.columns:
    sys.exit("❌ CSV does not contain a 'locus_tag' column.")

filtered_df = merged[merged["locus_tag"].isin(protein_ids)]

if filtered_df.empty:
    print("⚠️ No matching rows in CSV.")


filtered_df["Manual_annotation"] = PROTEIN

if TAXON_FILTER and TAXON_LEVEL:
    print("Before Taxon filtering: ", filtered_df.shape)

    tax_filtered = filtered_df[filtered_df[f"{TAXON_LEVEL}"] == f"{TAXON_FILTER}"]

    print("After Taxon filtering: ", tax_filtered.shape)

    tax_filtered.to_csv(OUT_CSV, index=False)
else:
    filtered_df.to_csv(OUT_CSV, index=False)


print(f"✅ CSV written → {OUT_CSV}")
print("🎉 Extraction complete.")
