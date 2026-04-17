#!/usr/bin/env python3

import os
import glob
import numpy as np
import pandas as pd
import logomaker
import matplotlib.pyplot as plt

# ==============================
# CONFIG
# ==============================
INPUT_FOLDER = "/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz2_run2/clade_wise"   # folder with aligned FASTA files
OUTPUT_FOLDER = "/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz2_run2/clade_wise"
WINDOW_SIZE = 30                # max positions per logo (for long alignments)
STEP_SIZE = 30                  # slide window (can overlap if < WINDOW_SIZE)
CENTER_ON_NON_GAP = False        # trims gap-heavy columns

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

# ==============================
# UTIL FUNCTIONS
# ==============================

def read_fasta(filepath):
    sequences = []
    with open(filepath) as f:
        seq = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    sequences.append(seq)
                    seq = ""
            else:
                seq += line
        if seq:
            sequences.append(seq)
    return sequences


def filter_alignment(seqs):
    """Remove columns that are mostly gaps"""
    arr = np.array([list(s) for s in seqs])
    keep_cols = []

    for i in range(arr.shape[1]):
        col = arr[:, i]
        gap_fraction = np.sum(col == "-") / len(col)
        if gap_fraction < 0.8:  # keep if <80% gaps
            keep_cols.append(i)

    return ["".join(arr[:, i] for i in keep_cols)]


def alignment_to_matrix(seqs):
    """Convert alignment to count matrix"""
    seq_len = len(seqs[0])
    matrix = []

    for i in range(seq_len):
        col = [s[i] for s in seqs]
        counts = {aa: col.count(aa) for aa in AMINO_ACIDS}
        matrix.append(counts)

    return pd.DataFrame(matrix)


def plot_logo(df, outpath, title):
    plt.figure(figsize=(len(df)*2, 3))

    logo = logomaker.Logo(df, shade_below=.5, fade_below=.5)

    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)

    logo.ax.set_ylabel("Count")
    logo.ax.set_title(title, fontsize=10)

    plt.tight_layout()
    plt.savefig(outpath, format="svg")
    plt.close()


# ==============================
# MAIN
# ==============================

def process_file(filepath):
    filename = os.path.basename(filepath).replace(".fasta", "")
    print(f"Processing {filename}...")

    seqs = read_fasta(filepath)

    # sanity check
    lengths = set(len(s) for s in seqs)
    if len(lengths) != 1:
        print(f"Skipping {filename} (unaligned sequences)")
        return

    seq_len = len(seqs[0])

    # sliding window for long alignments
    for start in range(0, seq_len, STEP_SIZE):
        end = start + WINDOW_SIZE
        window_seqs = [s[start:end] for s in seqs]

        # skip empty windows
        if all(set(s) == {"-"} for s in window_seqs):
            continue

        df = alignment_to_matrix(window_seqs)

        outname = f"{filename}_pos{start}-{end}.svg"
        outpath = os.path.join(OUTPUT_FOLDER, outname)

        plot_logo(df, outpath, f"{filename} ({start}-{end})")


def main():
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    fasta_files = glob.glob(os.path.join(INPUT_FOLDER, "*.fasta"))

    for f in fasta_files:
        process_file(f)

    print("All logos generated!")


if __name__ == "__main__":
    main()