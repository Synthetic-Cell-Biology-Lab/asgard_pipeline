#!/usr/bin/env python3

import argparse
import os
import re

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker
import numpy as np


############################################################
# Arguments
############################################################

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",       required=True)
    parser.add_argument("--outdir",      required=True)
    parser.add_argument("--protein",     required=True)
    parser.add_argument("--genome_file", required=True,
                        help="CSV with columns: genome_file + taxonomy ranks")
    parser.add_argument("--tax_level",   required=True,
                        help="Column name to group/split plots by "
                             "(e.g. phylum, class, order)")
    parser.add_argument("--top_n",       type=int, default=25)
    parser.add_argument("--min_genomes", type=int, default=1)
    parser.add_argument("--pfam_map", required=True)
    return parser.parse_args()


############################################################
# Color palette
############################################################

DISTINCT_20 = [
    "#E63946", "#F4A261", "#2A9D8F", "#457B9D", "#6A0572",
    "#F1C453", "#264653", "#A8DADC", "#E76F51", "#8338EC",
    "#06D6A0", "#FFB703", "#219EBC", "#FB8500", "#023047",
    "#8ECAE6", "#C77DFF", "#D62828", "#52B788", "#B5838D",
]

def make_palette(names):
    n = len(names)
    if n <= len(DISTINCT_20):
        colors = DISTINCT_20[:n]
    else:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            "distinct", DISTINCT_20
        )
        colors = [
            matplotlib.colors.to_hex(cmap(i / (n - 1)))
            for i in range(n)
        ]
    return dict(zip(names, colors))


############################################################
# Helpers
############################################################

def safe_name(s):
    return re.sub(r"[^A-Za-z0-9._-]", "_", str(s))


def explode_pfam(df):
    df = df[
        df["Pfam_acc"].notna() &
        (df["Pfam_acc"].astype(str).str.strip() != "") &
        (df["Pfam_acc"].astype(str).str.strip() != "nan")
    ].copy()
    df["Pfam_acc"] = df["Pfam_acc"].astype(str).str.split("; ")
    df = df.explode("Pfam_acc")
    df["Pfam_acc"] = df["Pfam_acc"].str.strip()
    return df[df["Pfam_acc"] != ""]


############################################################
# Plot
############################################################

def build_plot(df_slice, color_col, color_labels, palette,
               protein_name, tax_level, taxon_label, top_n, min_genomes):

    total_genomes = df_slice["genome"].nunique()

    # Genome count per (domain x color group)
    counts = (
        df_slice
        .groupby(["Pfam_label", color_col])["genome"]
        .nunique()
        .reset_index(name="genome_count")
    )

    totals = (
        counts.groupby("Pfam_label")["genome_count"]
        .sum()
        .reset_index(name="total_genomes")
        .query("total_genomes >= @min_genomes")
        .nlargest(top_n, "total_genomes")
    )

    if totals.empty:
        print(f"  Skipping '{taxon_label}': no domains after filtering")
        return None

    top_domains = totals["Pfam_label"].tolist()
    counts      = counts[counts["Pfam_label"].isin(top_domains)]

    pivot = (
        counts
        .pivot(index="Pfam_label", columns=color_col, values="genome_count")
        .reindex(index=top_domains, columns=color_labels)
        .fillna(0)
        .iloc[::-1]   # highest count at top
    )

    n_domains  = len(pivot)
    fig_height = max(6, n_domains * 0.55 + 3)
    fig, axes  = plt.subplots(
        1, 2, figsize=(22, fig_height),
        gridspec_kw={"width_ratios": [2.5, 1]}
    )

    # --- Panel 1: stacked genome count ---
    ax1   = axes[0]
    lefts = np.zeros(n_domains)

    for label in color_labels:
        vals = pivot[label].values
        ax1.barh(range(n_domains), vals, left=lefts,
                 color=palette.get(label, "#cccccc"), label=label,
                 height=0.7, edgecolor="white", linewidth=0.4)
        lefts += vals

    ax1_top = ax1.twiny()
    ax1_top.set_xlim(0, ax1.get_xlim()[1] / total_genomes)
    ax1_top.xaxis.set_major_formatter(
        matplotlib.ticker.PercentFormatter(xmax=1, decimals=0)
    )
    ax1_top.set_xlabel("Genome Proportion", fontsize=12)

    ax1.set_yticks(range(n_domains))
    ax1.set_yticklabels(pivot.index, fontsize=11)
    ax1.set_xlabel("Genome Count", fontsize=12)
    ax1.set_ylabel("Pfam Domain", fontsize=12)
    ax1.set_title("Genome Count / Proportion", fontsize=13, fontweight="bold")
    ax1.spines[["top", "right"]].set_visible(False)


    # --- Legend ---
    active = set(counts[color_col])
    handles = [
        mpatches.Patch(color=palette.get(t, "#cccccc"), label=t)
        for t in color_labels if t in active
    ]
    fig.legend(handles=handles, title=color_col.capitalize(),
               bbox_to_anchor=(1.01, 0.5), loc="center left",
               fontsize=10, title_fontsize=11, frameon=True)

    fig.suptitle(
        f"Top {top_n} Pfam domains in neighborhood of {protein_name}"
        f" | {tax_level}: {taxon_label}",
        fontsize=16, fontweight="bold", y=1.02
    )
    plt.tight_layout()

    return fig


############################################################
# Main
############################################################

def main():

    args = parse_args()

    print("📂 Loading data")
    df = pd.read_csv(args.input)

    tax_df = pd.read_csv(args.genome_file)

    if "genome_file" not in tax_df.columns:
        raise ValueError("genome_file CSV must have a 'genome_file' column")
    if args.tax_level not in tax_df.columns:
        raise ValueError(f"genome_file CSV missing column '{args.tax_level}'")

    # Color one rank below the split level
    tax_hierarchy = ["species", "genus", "family", "order", "class", "phylum", "domain"]
    tax_idx       = tax_hierarchy.index(args.tax_level)
    color_col     = tax_hierarchy[tax_idx - 1]

    if color_col not in tax_df.columns:
        raise ValueError(f"genome_file CSV missing color column '{color_col}'")

    df = df.merge(
        tax_df[["genome_file", args.tax_level, color_col]]
              .rename(columns={"genome_file": "genome"}),
        on="genome", how="left"
    )

    df[args.tax_level] = df[args.tax_level].fillna("Unknown").astype(str).str.strip()
    df[color_col]      = df[color_col].fillna("Unknown").astype(str).str.strip()

    print("🔎 Exploding Pfam accessions")
    df = explode_pfam(df)
    

    if df.empty:
        print("⚠ No Pfam annotations found")
        return
    
    print("📖 Loading Pfam descriptions")

    pfam_map = pd.read_csv(
        args.pfam_map,
        sep="\t",
        names=["Pfam_acc", "Pfam_desc"]
    )

    # assumes columns: Pfam_acc, Pfam_desc
    pfam_dict = dict(
        zip(
            pfam_map["Pfam_acc"],
            pfam_map["Pfam_desc"]
        )
    )

    # Replace accession with description if available
    df["Pfam_label"] = df["Pfam_acc"].map(pfam_dict)

    # fallback to accession if no description found
    df["Pfam_label"] = df["Pfam_label"].fillna(df["Pfam_acc"])

    # Global palette built before any split so colors are consistent
    all_color_labels = sorted(df[color_col].unique())
    palette          = make_palette(all_color_labels)

    os.makedirs(args.outdir, exist_ok=True)

    groups = (
        df.groupby(args.tax_level)["genome"]
        .nunique()
        .reset_index(name="n")
        .sort_values("n", ascending=False)
    )

    print(f"\nFound {len(groups)} groups at '{args.tax_level}' level\n")

    for taxon_label in groups[args.tax_level]:

        df_slice = df[df[args.tax_level] == taxon_label]
        n        = df_slice["genome"].nunique()
        print(f"  '{taxon_label}': {n} genomes")

        fig = build_plot(
            df_slice, color_col, all_color_labels, palette,
            args.protein, args.tax_level, taxon_label,
            args.top_n, args.min_genomes
        )

        if fig is None:
            continue

        out = os.path.join(args.outdir, f"{safe_name(taxon_label)}.pfam_frequency.svg")
        fig.savefig(out, format="svg", bbox_inches="tight")
        plt.close(fig)
        print(f"  ✅ {out}")

    print("\nDone.")


if __name__ == "__main__":
    main()