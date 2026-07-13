#!/usr/bin/env python3

"""
Create conserved-block alignment figure from an MSA.

Inputs
------
msa.fasta
metadata.csv

metadata.csv columns:
locus_tag,class,order

locus_tag must match FASTA header

Example:
locus_tag,class,order
seq1,Lokiarchaeia,Lokiarchaeales
seq2,Lokiarchaeia,Lokiarchaeales

Usage
-----
python conserved_blocks.py \
    --msa msa.fasta \
    --metadata metadata.csv \
    --output conserved_blocks.png
"""

from Bio import AlignIO
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from collections import Counter
import argparse

# =====================================================
# USER SETTINGS
# =====================================================
# #ftsz1
# BLOCKS = [
#     (0, 0),
#     (178, 196),
#     (209, 224),
#     (272, 300),
#     (318, 347),
#     (355, 371),
#     (380, 390),
#     (602, 640),
# ]

# ftsz2
# BLOCKS = [
#     (0, 0),
#     (78, 96),
#     (104, 119),
#     (157, 190),
#     (210, 240),
#     (272, 282),
#     (295, 340),
#     (555, 562),
#     (576, 576),
# ]


# actin_selected


# ftsz1-2
# BLOCKS = [(223, 242), (255, 284), (290, 310), (322, 349), (410, 444)]

# odin ftsz

BLOCKS = [(0, 0), (360, 420), (625, 625)]

FONT_SIZE = 8
LABEL_WIDTH = 30
CELL_W = 1
CELL_H = 1

# Fixed number of columns reserved for every inter-block spacer number.
# 5 works for up to 3-digit counts with one blank cell on each side.
SPACER_WIDTH = 5

# Minimum fraction of non-gap residues that must share the same
# amino-acid identity for a column to receive background colour.
CONSERVATION_THRESHOLD = 0.90

# =====================================================
# COLOR SCHEME
# =====================================================

GREEN_SET = set("ACGNPSTV")
YELLOW_SET = set("ACFILMVWY")
RED_SET = set("DEHKNQR")

GREEN_BG = "#66ff66"
YELLOW_BG = "#ffff66"


# =====================================================
# HELPERS
# =====================================================


def compute_column_conservation(sequences, col_index):
    """
    Return the most-common non-gap residue and its frequency
    (as a fraction of non-gap positions) for a given MSA column.

    Returns (dominant_aa, frequency).
    If the column is all gaps, returns ('-', 0.0).
    """
    residues = [seq[col_index] for seq in sequences if seq[col_index] != "-"]

    if not residues:
        return "-", 0.0

    most_common_aa, count = Counter(residues).most_common(1)[0]
    return most_common_aa, count / len(residues)


def build_conservation_map(sequences, blocks):
    """
    Pre-compute per-column conservation for every column
    that appears inside any block.

    Returns a dict: { msa_column_index -> (dominant_aa, frequency) }
    """
    conservation = {}
    for start, end in blocks:
        for col in range(start, end + 1):
            conservation[col] = compute_column_conservation(sequences, col)
    return conservation


def residue_style(aa, col_index, conservation_map):
    """
    Returns (background_color, text_color).

    Colour is applied only when the dominant residue at this column
    meets CONSERVATION_THRESHOLD *and* the current residue matches
    that dominant residue.
    """
    if aa == "-":
        return "white", "black"

    dominant_aa, freq = conservation_map.get(col_index, ("-", 0.0))

    if freq >= CONSERVATION_THRESHOLD and aa == dominant_aa:
        if aa in GREEN_SET:
            return GREEN_BG, "black"
        if aa in YELLOW_SET:
            return YELLOW_BG, "black"
        if aa in RED_SET:
            return "white", "red"

    return "white", "black"


def extract_blocks(sequence, blocks):
    """
    Extract conserved blocks and insert the number of
    non-gap residues between consecutive blocks.

    Returns a list of tuples:
        ("block",  <substring>, msa_start_col)
        ("spacer", <count_string>, None)
    """
    pieces = []

    for i, (start, end) in enumerate(blocks):

        pieces.append(("block", sequence[start : end + 1], start))

        if i < len(blocks) - 1:
            next_start = blocks[i + 1][0]
            omitted_region = sequence[end + 1 : next_start]
            residue_count = sum(aa != "-" for aa in omitted_region)
            pieces.append(("spacer", str(residue_count), None))

    return pieces


# =====================================================
# MAIN
# =====================================================


def main(msa_file, metadata_file, output):

    aln = AlignIO.read(msa_file, "fasta")

    metadata = pd.read_csv(metadata_file)

    meta_dict = metadata.set_index("locus_tag").to_dict("index")

    rows = []

    for record in aln:

        seq_id = record.id.split("/")[0]

        if seq_id not in meta_dict:
            print(f"Skipping {seq_id}: not found in metadata")
            continue

        tax = meta_dict[seq_id]

        rows.append(
            {
                "class": tax["class"],
                "order": tax["order"],
                "locus_tag": seq_id,
                "sequence": str(record.seq),
            }
        )

    rows = sorted(
        rows,
        key=lambda x: (
            str(x["class"]) if pd.notna(x["class"]) else "",
            str(x["order"]) if pd.notna(x["order"]) else "",
        ),
    )

    # -------------------------------------------------
    # Pre-compute conservation across ALL sequences in
    # the MSA (not just the filtered subset) so the
    # threshold reflects the full alignment.
    # -------------------------------------------------

    all_sequences = [str(record.seq) for record in aln]
    conservation_map = build_conservation_map(all_sequences, BLOCKS)

    # -------------------------------------------------
    # Pre-render rows and determine figure size
    # -------------------------------------------------

    total_rows = len(rows)
    longest_line = 0
    rendered_rows = []

    for row in rows:

        blocks = extract_blocks(row["sequence"], BLOCKS)
        label = (
            f"{row['order']}_{row['locus_tag']}"
            if pd.notna(row["order"]) and row["order"] != ""
            else row["locus_tag"]
        )
        rendered_rows.append((row, label, blocks))

        width = sum(
            SPACER_WIDTH if piece_type == "spacer" else len(piece)
            for piece_type, piece, *_ in blocks
        )

        longest_line = max(longest_line, width)

    fig_w = max(25, longest_line * 0.15)
    fig_h = max(8, total_rows * 0.25)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=300)
    ax.axis("off")

    y = 0
    current_class = None

    # -------------------------------------------------
    # Draw rows
    # -------------------------------------------------

    for row, label, blocks in rendered_rows:

        # Class header
        if row["class"] != current_class:
            y += 1
            ax.text(
                0,
                -y,
                row["class"],
                fontsize=FONT_SIZE + 2,
                fontweight="bold",
                family="monospace",
            )
            y += 1
            current_class = row["class"]

        # Row label
        ax.text(0, -y, label, fontsize=FONT_SIZE, family="monospace", va="center")

        x = LABEL_WIDTH

        for piece_type, piece, msa_start in blocks:

            # ---- inter-block spacer number ----
            if piece_type == "spacer":
                ax.text(
                    x + SPACER_WIDTH / 2,
                    -y,
                    piece,
                    fontsize=FONT_SIZE,
                    family="monospace",
                    va="center",
                    ha="center",
                )
                x += SPACER_WIDTH
                continue

            # ---- residue block ----
            for offset, aa in enumerate(piece):

                col_index = msa_start + offset
                bg, fg = residue_style(aa, col_index, conservation_map)

                rect = Rectangle(
                    (x, -y - 0.5), CELL_W, CELL_H, facecolor=bg, edgecolor="none"
                )
                ax.add_patch(rect)

                ax.text(
                    x + 0.5,
                    -y,
                    aa,
                    ha="center",
                    va="center",
                    fontsize=FONT_SIZE,
                    color=fg,
                    family="monospace",
                )

                x += 1

        y += 1

    ax.set_xlim(0, LABEL_WIDTH + longest_line + 10)
    ax.set_ylim(-y - 1, 2)

    plt.tight_layout()
    plt.savefig(output, bbox_inches="tight")
    print(f"Saved to {output}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--msa", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    main(args.msa, args.metadata, args.output)
