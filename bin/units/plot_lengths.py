import argparse
import logging
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

BIN_LABELS = ["short", "mid_short", "mid_long", "long"]


def parse_args():
    p = argparse.ArgumentParser(
        description="Bin FASTA sequences by length, split into files, and plot length distribution."
    )
    p.add_argument("fasta", help="Input FASTA file")
    p.add_argument("output", help="Output plot image (e.g. summary.png)")
    p.add_argument("split_dir", help="Directory to write length-binned FASTA files")
    p.add_argument(
        "--max-swarm",
        type=int,
        default=2000,
        help="Skip swarmplot and use stripplot above this many sequences (default: 2000)",
    )
    return p.parse_args()


def load_records(fasta_path):
    if not os.path.isfile(fasta_path):
        log.error(f"FASTA file not found: {fasta_path}")
        sys.exit(1)

    try:
        records = list(SeqIO.parse(fasta_path, "fasta"))
    except Exception as e:
        log.error(f"Failed to parse FASTA file '{fasta_path}': {e}")
        sys.exit(1)

    if not records:
        log.error(f"No sequences found in '{fasta_path}'. Is it a valid FASTA file?")
        sys.exit(1)

    return records


def assign_bins(df):
    """Assign length bins, falling back gracefully for small/degenerate data."""
    n_unique = df["Length"].nunique()

    if n_unique == 1:
        log.warning("All sequences have the same length; assigning a single bin.")
        df["quantile_bin"] = "all"
        return df

    if n_unique < 10:
        n_bins = min(4, n_unique)
        log.warning(
            f"Only {n_unique} unique lengths; using {n_bins} equal-width bins instead of quartiles."
        )
        codes = pd.cut(df["Length"], bins=n_bins, labels=False, duplicates="drop")
        # give consistent, human-readable names regardless of how many bins we got
        n_actual = codes.nunique()
        names = (
            BIN_LABELS[:n_actual]
            if n_actual <= len(BIN_LABELS)
            else [f"bin_{i}" for i in range(n_actual)]
        )
        df["quantile_bin"] = codes.map(dict(enumerate(names)))
    else:
        # duplicates="drop" avoids the ValueError when many sequences share
        # a length and quartile edges collide; may yield <4 bins.
        codes = pd.qcut(df["Length"], q=4, labels=False, duplicates="drop")
        n_actual = codes.nunique()
        if n_actual < 4:
            log.warning(
                f"Quartile binning collapsed to {n_actual} bins due to repeated lengths."
            )
        names = (
            BIN_LABELS[:n_actual]
            if n_actual <= len(BIN_LABELS)
            else [f"bin_{i}" for i in range(n_actual)]
        )
        df["quantile_bin"] = codes.map(dict(enumerate(names)))

    return df


def main():
    args = parse_args()

    records = load_records(args.fasta)
    lengths = [len(r.seq) for r in records]
    df = pd.DataFrame({"Length": lengths, "Record": records})

    df = assign_bins(df)

    os.makedirs(args.split_dir, exist_ok=True)
    out_dir_for_plot = os.path.dirname(args.output)
    if out_dir_for_plot:
        os.makedirs(out_dir_for_plot, exist_ok=True)

    for q in df["quantile_bin"].dropna().unique():
        subset = df[df["quantile_bin"] == q]
        out_path = os.path.join(args.split_dir, f"{q}.fasta")
        SeqIO.write(subset["Record"].tolist(), out_path, "fasta")
        log.info(f"{q}: {len(subset)} sequences → {out_path}")

    # -----------------------------
    # Style
    # -----------------------------
    sns.set_theme(style="whitegrid", context="talk")
    plt.rcParams["figure.dpi"] = 120

    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.25)

    # 1. Histogram
    ax1 = fig.add_subplot(gs[0, 0])
    sns.histplot(
        df["Length"],
        bins=min(40, df["Length"].nunique()),
        kde=(df["Length"].nunique() > 1),
        ax=ax1,
    )
    ax1.set_title("Histogram + KDE")
    ax1.set_xlabel("Sequence Length")
    ax1.set_ylabel("Frequency")
    for q in df["Length"].quantile([0.25, 0.5, 0.75]):
        ax1.axvline(q, linestyle="--", alpha=0.7)

    # 2. Violin plot
    ax2 = fig.add_subplot(gs[0, 1])
    sns.violinplot(y=df["Length"], inner="quartile", ax=ax2)
    ax2.set_title("Violin Plot")
    ax2.set_ylabel("Sequence Length")

    # 3. Bar plot (binned counts)
    ax3 = fig.add_subplot(gs[1, 0])
    n_bins = min(15, max(2, df["Length"].nunique()))
    bins = np.linspace(df["Length"].min(), df["Length"].max(), n_bins)
    bin_col = pd.cut(df["Length"], bins=bins, duplicates="drop")
    bin_counts = bin_col.value_counts().sort_index()
    bin_labels = [f"{int(b.left)}-{int(b.right)}" for b in bin_counts.index]
    ax3.bar(bin_labels, bin_counts.values)
    ax3.set_title("Binned Counts")
    ax3.set_xlabel("Length Range")
    ax3.set_ylabel("Count")
    ax3.tick_params(axis="x", rotation=45)

    # 4. Beeswarm / strip plot
    ax4 = fig.add_subplot(gs[1, 1])
    if len(df) > args.max_swarm:
        log.warning(
            f"{len(df)} sequences exceeds --max-swarm ({args.max_swarm}); using stripplot instead of swarmplot."
        )
        sns.stripplot(y=df["Length"], size=3, alpha=0.5, ax=ax4)
        ax4.set_title("Strip Plot")
    else:
        sns.swarmplot(y=df["Length"], size=3, ax=ax4)
        ax4.set_title("Beeswarm Plot")
    ax4.set_ylabel("Sequence Length")

    fig.suptitle("Sequence Length Distribution Overview", fontsize=18)

    plt.savefig(args.output, bbox_inches="tight")
    plt.close()

    df[["Length", "quantile_bin"]].to_csv(
        os.path.join(args.split_dir, "length_bins.csv"), index=False
    )
    log.info(f"Saved plot → {args.output}")


if __name__ == "__main__":
    main()
