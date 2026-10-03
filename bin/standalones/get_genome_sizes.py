# script to plot the genome sizes comparing at a given taxonomic level
#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


TAXONOMIC_LEVELS = [
    "domain",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]


def plot_genome_sizes(
    input_file,
    taxonomic_level,
    output_file=None,
    genome_size_col="genome_size",
    figsize=None,
    dpi=300,
    max_groups=None,
):
    """
    Plot comparative genome-size distributions across taxonomic groups.

    Parameters
    ----------
    input_file : str
        Path to input TSV file.
    taxonomic_level : str
        Taxonomic level to compare (e.g. phylum, genus, species).
    output_file : str, optional
        Output figure path. Defaults to a PNG in the current directory.
    genome_size_col : str
        Column containing genome sizes.
    figsize : tuple, optional
        Figure dimensions (width, height).
    dpi : int
        Resolution of saved figure.
    max_groups : int, optional
        Display only the N groups with the most observations.
        Useful for highly resolved taxonomic levels.
    """

    taxonomic_level = taxonomic_level.lower()

    if taxonomic_level not in TAXONOMIC_LEVELS:
        raise ValueError(
            f"Invalid taxonomic level: {taxonomic_level}\n"
            f"Choose from: {', '.join(TAXONOMIC_LEVELS)}"
        )

    # --------------------------------------------------
    # 1. Read input
    # --------------------------------------------------

    df = pd.read_csv(
        input_file,
        sep="\t",
        low_memory=False,
    )

    # Normalize header names
    df.columns = df.columns.str.strip().str.lower()

    required_cols = [taxonomic_level, genome_size_col.lower()]
    missing_cols = [col for col in required_cols if col not in df.columns]

    if missing_cols:
        raise ValueError(
            f"Missing required columns: {missing_cols}\n"
            f"Available columns: {list(df.columns)}"
        )

    # --------------------------------------------------
    # 2. Clean data
    # --------------------------------------------------

    df = df[required_cols].copy()

    df.columns = ["taxon", "genome_size"]

    df["taxon"] = df["taxon"].astype("string").str.strip()

    # Convert genome sizes to numeric
    df["genome_size"] = pd.to_numeric(
        df["genome_size"],
        errors="coerce",
    )

    # Remove missing, empty, and invalid observations
    df = df.dropna(subset=["taxon", "genome_size"])

    df = df[
        (df["taxon"] != "")
        & (df["genome_size"] > 0)
        & (~df["taxon"].str.lower().isin(
            ["nan", "none", "na", "unknown", "unclassified"]
        ))
    ].copy()

    if df.empty:
        raise ValueError("No valid genome-size observations remain.")

    # --------------------------------------------------
    # 3. Optional group filtering
    # --------------------------------------------------

    group_counts = df["taxon"].value_counts()

    if max_groups is not None:
        selected_groups = group_counts.head(max_groups).index
        df = df[df["taxon"].isin(selected_groups)].copy()

    # Order groups by median genome size
    group_order = (
        df.groupby("taxon")["genome_size"]
        .median()
        .sort_values()
        .index
        .tolist()
    )

    n_groups = len(group_order)

    # Automatically scale figure width
    if figsize is None:
        width = max(8, min(0.65 * n_groups + 3, 24))
        figsize = (width, 6)

    # --------------------------------------------------
    # 4. Scientific plotting style
    # --------------------------------------------------

    sns.set_theme(
        context="paper",
        style="ticks",
        font="sans-serif",
        font_scale=1.2,
        rc={
            "axes.linewidth": 1.0,
            "axes.labelweight": "medium",
            "axes.titleweight": "bold",
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "savefig.bbox": "tight",
        },
    )

    fig, ax = plt.subplots(figsize=figsize)

    palette = sns.color_palette(
        "bright",
        n_colors=n_groups,
    )

    # --------------------------------------------------
    # 5. Box plot
    # --------------------------------------------------

    sns.boxplot(
        data=df,
        x="taxon",
        y="genome_size",
        order=group_order,
        hue="taxon",
        hue_order=group_order,
        palette=palette,
        dodge=False,
        legend=False,
        width=0.65,
        linewidth=1.1,
        fliersize=0,
        medianprops={
            "color": "black",
            "linewidth": 1.5,
        },
        whiskerprops={
            "linewidth": 1.1,
        },
        capprops={
            "linewidth": 1.1,
        },
        boxprops={
            "edgecolor": "black",
            "alpha": 0.85,
        },
        ax=ax,
    )

    # --------------------------------------------------
    # 6. Overlay individual observations
    # --------------------------------------------------

    sns.stripplot(
        data=df,
        x="taxon",
        y="genome_size",
        order=group_order,
        color="black",
        alpha=0.35,
        size=2.5,
        jitter=0.22,
        zorder=2,
        ax=ax,
    )

    # --------------------------------------------------
    # 7. Labels and formatting
    # --------------------------------------------------

    ax.set_xlabel(
        taxonomic_level.capitalize(),
        fontsize=13,
        labelpad=10,
    )

    ax.set_ylabel(
        "Genome size (bp)",
        fontsize=13,
        labelpad=10,
    )

    ax.set_title(
        f"Genome size across {taxonomic_level} groups",
        fontsize=15,
        pad=15,
    )

    ax.tick_params(
        axis="x",
        labelrotation=60,
        labelsize=9,
    )

    ax.tick_params(
        axis="y",
        labelsize=10,
    )

    # Scientific notation for large genome sizes
    ax.ticklabel_format(
        axis="y",
        style="sci",
        scilimits=(6, 6),
        useMathText=True,
    )

    ax.yaxis.get_offset_text().set_fontsize(10)

    # Remove top and right spines
    sns.despine(ax=ax)

    # Light horizontal gridlines
    ax.set_axisbelow(True)

    ax.yaxis.grid(
        True,
        linestyle="--",
        linewidth=0.5,
        alpha=0.35,
    )

    ax.xaxis.grid(False)

    plt.tight_layout()

    # --------------------------------------------------
    # 8. Save figure
    # --------------------------------------------------

    if output_file is None:
        output_file = (
            f"genome_sizes_by_{taxonomic_level}.png"
        )

    output_file = Path(output_file)

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fig.savefig(
        output_file,
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
    )

    print(f"Taxonomic level : {taxonomic_level}")
    print(f"Valid genomes   : {len(df)}")
    print(f"Taxonomic groups: {n_groups}")
    print(f"Saved figure    : {output_file.resolve()}")

    plt.close(fig)


# ------------------------------------------------------
# Command-line interface
# ------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=(
            "Plot genome-size distributions across "
            "taxonomic groups from a TSV file."
        )
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input TSV file",
    )

    parser.add_argument(
        "-t", "--taxon",
        required=True,
        choices=TAXONOMIC_LEVELS,
        help="Taxonomic level to compare",
    )

    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output figure filename",
    )

    parser.add_argument(
        "--genome-size-col",
        default="genome_size",
        help="Genome-size column name (default: genome_size)",
    )

    parser.add_argument(
        "--max-groups",
        type=int,
        default=None,
        help="Plot only the N most abundant taxonomic groups",
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Figure resolution (default: 300)",
    )

    args = parser.parse_args()

    plot_genome_sizes(
        input_file=args.input,
        taxonomic_level=args.taxon,
        output_file=args.output,
        genome_size_col=args.genome_size_col,
        max_groups=args.max_groups,
        dpi=args.dpi,
    )