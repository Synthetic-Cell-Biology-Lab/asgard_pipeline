#!/usr/bin/env python3
"""
trim_codon_aln.py — Trim a codon alignment using column indices from trimAl.

Workflow:
    1. Run MAFFT on protein sequences -> aligned.faa
    2. Run trimAl on aligned.faa with -colnumbering -> trimmed.faa + cols.txt
    3. Run PAL2NAL on aligned.faa (untrimmed) + nucleotides -> codon.fna
    4. Run this script: trim_codon_aln.py codon.fna cols.txt trimmed_codon.fna

trimAl -colnumbering command:
    trimal -in aligned.faa -out trimmed.faa -automated1 -colnumbering > cols.txt

The -colnumbering output is a single line of 0-based column indices, e.g.:
    #ColumnsMap 0, 1, 2, 5, 6, 7, 10, 11 ...
or just:
    0, 1, 2, 5, 6, 7, 10, 11 ...
Both formats are handled.
"""

import sys
import argparse
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def parse_args():
    p = argparse.ArgumentParser(
        description="Trim a codon alignment using trimAl -colnumbering indices."
    )
    p.add_argument("codon_aln",  help="PAL2NAL codon alignment (FASTA, nucleotide)")
    p.add_argument("cols_file",  help="trimAl -colnumbering output (kept protein column indices)")
    p.add_argument("output",     help="Output trimmed codon alignment (FASTA)")
    p.add_argument(
        "--codon-table", type=int, default=1,
        help="NCBI codon table number for optional translation check (default: 1)"
    )
    p.add_argument(
        "--no-check", action="store_true",
        help="Skip post-trim sanity checks"
    )
    return p.parse_args()


def parse_colnumbering(cols_file: str) -> list[int]:
    """
    Parse trimAl -colnumbering output.
    Handles both plain integers and the '#ColumnsMap ...' header format.
    Returns a sorted list of 0-based protein column indices.
    """
    cols = []
    with open(cols_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Strip optional '#ColumnsMap' prefix
            if line.startswith("#"):
                line = line.split(None, 1)[1] if " " in line else ""
            # Parse comma- or whitespace-separated integers
            for token in line.replace(",", " ").split():
                try:
                    cols.append(int(token))
                except ValueError:
                    continue
    if not cols:
        raise ValueError(f"No column indices parsed from '{cols_file}'. "
                         "Make sure trimAl was run with -colnumbering.")
    return sorted(cols)


def protein_cols_to_nuc_cols(protein_cols: list[int]) -> list[int]:
    """Convert 0-based protein column indices to nucleotide column indices (×3)."""
    nuc_cols = []
    for c in protein_cols:
        nuc_cols += [c * 3, c * 3 + 1, c * 3 + 2]
    return nuc_cols


def trim_codon_alignment(
    codon_aln: str,
    nuc_cols: list[int],
    output: str,
    codon_table: int = 1,
    check: bool = True,
) -> None:
    records = list(SeqIO.parse(codon_aln, "fasta"))
    if not records:
        raise ValueError(f"No sequences found in '{codon_aln}'.")

    seq_len = len(records[0].seq)

    # Validate all sequences are the same length
    bad = [r.id for r in records if len(r.seq) != seq_len]
    if bad:
        raise ValueError(
            f"Sequences are not the same length in '{codon_aln}'. "
            f"Offending IDs: {bad[:5]}"
        )

    # Validate codon alignment is a multiple of 3
    if seq_len % 3 != 0:
        raise ValueError(
            f"Codon alignment length ({seq_len}) is not a multiple of 3. "
            "Check your PAL2NAL output."
        )

    # Validate all requested columns exist
    max_col = max(nuc_cols)
    if max_col >= seq_len:
        raise IndexError(
            f"Column index {max_col} is out of range for alignment of length {seq_len}. "
            "Are the column indices from the correct trimAl run?"
        )

    # Validate output will be a multiple of 3
    if len(nuc_cols) % 3 != 0:
        raise ValueError(
            f"Number of nucleotide columns to keep ({len(nuc_cols)}) is not a multiple of 3. "
            "Something went wrong with the column index parsing."
        )

    # Extract columns
    trimmed_records = []
    for rec in records:
        seq = str(rec.seq)
        trimmed_seq = "".join(seq[i] for i in nuc_cols)
        trimmed_records.append(
            SeqRecord(Seq(trimmed_seq), id=rec.id, description="")
        )

    SeqIO.write(trimmed_records, output, "fasta")

    # Sanity checks
    if check:
        n_codons_in  = seq_len // 3
        n_codons_out = len(nuc_cols) // 3
        pct_kept = 100 * n_codons_out / n_codons_in

        print(f"  Sequences:       {len(records)}")
        print(f"  Input codons:    {n_codons_in}")
        print(f"  Kept codons:     {n_codons_out} ({pct_kept:.1f}%)")
        print(f"  Output length:   {len(nuc_cols)} nt")

        # Spot-check: translate first sequence and confirm no internal stops
        first_seq = str(trimmed_records[0].seq).replace("-", "N")
        # Pad to multiple of 3 if gaps cause partial codon at end
        pad = (-len(first_seq)) % 3
        first_seq += "N" * pad
        try:
            protein = Seq(first_seq).translate(table=codon_table, to_stop=False)
            internal_stops = str(protein)[:-1].count("*")
            if internal_stops > 0:
                print(
                    f"  WARNING: {internal_stops} internal stop codon(s) detected in "
                    f"'{trimmed_records[0].id}' after trimming. "
                    "This may indicate a frame shift or a genuine pseudogene.",
                    file=sys.stderr,
                )
            else:
                print(f"  Translation check: OK (no internal stops in first sequence)")
        except Exception as e:
            print(f"  WARNING: translation check failed: {e}", file=sys.stderr)

        print(f"  Written to:      {output}")


def main():
    args = parse_args()

    print("Parsing trimAl column indices...")
    protein_cols = parse_colnumbering(args.cols_file)
    print(f"  Protein columns to keep: {len(protein_cols)}")

    nuc_cols = protein_cols_to_nuc_cols(protein_cols)
    print(f"  Nucleotide columns to keep: {len(nuc_cols)}")

    print("Trimming codon alignment...")
    trim_codon_alignment(
        codon_aln=args.codon_aln,
        nuc_cols=nuc_cols,
        output=args.output,
        codon_table=args.codon_table,
        check=not args.no_check,
    )

    print("Done.")


if __name__ == "__main__":
    main()