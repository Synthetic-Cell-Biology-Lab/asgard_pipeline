#!/usr/bin/env python3

import argparse
from Bio import SeqIO


def filter_by_length(input_fasta, output_fasta, min_len, max_len):
    """
    Filter protein sequences based on sequence length.

    Parameters
    ----------
    input_fasta : str
        Input FASTA file.
    output_fasta : str
        Output FASTA file.
    min_len : int
        Minimum sequence length (inclusive).
    max_len : int
        Maximum sequence length (inclusive).
    """

    count_in = 0
    count_out = 0

    with open(output_fasta, "w") as out_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            count_in += 1
            seq_len = len(record.seq)

            if min_len <= seq_len <= max_len:
                SeqIO.write(record, out_handle, "fasta")
                count_out += 1

    print(f"Total sequences read      : {count_in}")
    print(f"Sequences retained        : {count_out}")
    print(f"Sequences removed         : {count_in - count_out}")
    print(f"Output written to         : {output_fasta}")


def main():
    parser = argparse.ArgumentParser(
        description="Filter protein sequences based on sequence length."
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input protein FASTA file",
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output FASTA file",
    )

    parser.add_argument(
        "--min_len",
        type=int,
        required=True,
        help="Minimum allowed sequence length",
    )

    parser.add_argument(
        "--max_len",
        type=int,
        required=True,
        help="Maximum allowed sequence length",
    )

    args = parser.parse_args()

    if args.min_len > args.max_len:
        parser.error("min_len cannot be greater than max_len.")

    filter_by_length(
        args.input,
        args.output,
        args.min_len,
        args.max_len,
    )


if __name__ == "__main__":
    main()