#!/usr/bin/env python3

import sys

def read_fasta(fasta_file):
    sequences = {}
    current_seq = None

    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_seq = line[1:].split()[0]
                sequences[current_seq] = ""
            else:
                sequences[current_seq] += line

    return sequences


def check_alignment(sequences):
    lengths = set(len(seq) for seq in sequences.values())
    if len(lengths) != 1:
        raise ValueError(f"Alignment error: sequences have different lengths: {lengths}")
    return lengths.pop()


def write_itol_fasta(sequences, output_file, label="MSA"):
    with open(output_file, "w") as out:
        out.write("DATASET_ALIGNMENT\n")
        out.write("SEPARATOR COMMA\n")
        out.write(f"DATASET_LABEL,{label}\n")
        out.write("COLOR,#ff0000\n")
        out.write("COLOR_SCHEME,clustal\n")
        out.write("DISPLAY_CONSENSUS,1\n")
        out.write("DISPLAY_CONSERVATION,1\n\n")
        out.write("DATA\n\n")

        for name, seq in sequences.items():
            out.write(f">{name}\n")
            out.write(f"{seq}\n\n")


def main():
    if len(sys.argv) < 3:
        print("Usage: python fasta_to_itol_fasta.py input.fasta output.txt")
        sys.exit(1)

    fasta_file = sys.argv[1]
    output_file = sys.argv[2]

    sequences = read_fasta(fasta_file)
    aln_length = check_alignment(sequences)

    print(f"Loaded {len(sequences)} sequences")
    print(f"Alignment length: {aln_length}")

    write_itol_fasta(sequences, output_file)

    print(f"iTOL FASTA-style dataset written to: {output_file}")


if __name__ == "__main__":
    main()