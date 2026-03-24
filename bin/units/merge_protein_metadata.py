#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def parse_fasta(fasta_path):
    """
    Extract locus_tag + sequence info from FASTA
    """
    records = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        header = record.description

        # --- Extract locus_tag (FIRST token safest) ---
        locus_tag = header.split()[0]

        records.append({
            "locus_tag": locus_tag,
            "sequence": str(record.seq),
            "seq_length": len(record.seq)
        })

    return pd.DataFrame(records)


def main():
    parser = argparse.ArgumentParser(description="Merge FASTA + annotation + metadata")

    parser.add_argument("--fasta", required=True)
    parser.add_argument("--annotation", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    # ---------------- LOAD ----------------
    fasta_df = parse_fasta(args.fasta)
    ann_df = pd.read_csv(args.annotation)
    meta_df = pd.read_csv(args.metadata)

    # ---------------- DEBUG SAFETY ----------------
    if "locus_tag" not in ann_df.columns:
        raise ValueError("Annotation CSV must contain 'locus_tag'")

    if "genome_file" not in ann_df.columns:
        raise ValueError("Annotation CSV must contain 'genome_file'")

    if "genome_file" not in meta_df.columns:
        raise ValueError("Metadata CSV must contain 'genome_file'")

    # ---------------- MERGE ----------------
    merged = fasta_df.merge(ann_df, on="locus_tag", how="left")

    merged = merged.merge(
        meta_df,
        on="genome_file",
        how="left",
        suffixes=("", "_meta")
    )

    # ---------------- OPTIONAL CLEANUP ----------------
    # Remove completely empty columns
    merged = merged.dropna(axis=1, how="all")

    # ---------------- SAVE ----------------
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()