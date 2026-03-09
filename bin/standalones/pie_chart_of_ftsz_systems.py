#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import argparse


############################################################
# HARD CODED CONFIG
############################################################

USE_HARDCODED = True

INPUT_CSV = "/home/anirudh/asgard_pipeline/database/collated/Version1/filtered/85comp10con/protein_file/jan2026_85comp10con_pf.csv"

OUTPUT_PREFIX = "/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsZ_system_summary"


############################################################
# Load and process data
############################################################

def process_csv(input_csv):

    df = pd.read_csv(input_csv)

    # keep only FtsZ1 and FtsZ2
    df = df[df["Manual_annotation"].isin(["FtsZ1", "FtsZ2"])]

    # group per genome
    grouped = df.groupby("genome_file")["Manual_annotation"].agg(list).reset_index()

    grouped["has_FtsZ1"] = grouped["Manual_annotation"].apply(lambda x: "FtsZ1" in x)
    grouped["has_FtsZ2"] = grouped["Manual_annotation"].apply(lambda x: "FtsZ2" in x)

    grouped["has_both"] = grouped["has_FtsZ1"] & grouped["has_FtsZ2"]

    grouped["system_type"] = grouped.apply(
        lambda r: "2 FtsZ system" if r["has_both"] else "1 FtsZ system",
        axis=1
    )

    grouped = grouped.drop(columns=["Manual_annotation"])

    return grouped


############################################################
# Plot pie chart
############################################################
def plot_pie(grouped, output_png):

    import matplotlib.pyplot as plt

    counts = grouped["system_type"].value_counts()

    labels = counts.index.tolist()
    sizes = counts.values
    total = sum(sizes)

    # Okabe–Ito colorblind-safe palette
    colors = [
        "#7fcdbb",  # pastel teal
        "#f4a6a6"   # pastel coral
    ]


    fig, ax = plt.subplots(figsize=(6,6), facecolor="white")

    wedges, texts, autotexts = ax.pie(
        sizes,
        autopct=lambda pct: f"{pct:.1f}%\n(n={int(round(pct*total/100))})",
        startangle=90,
        colors=colors,
        wedgeprops=dict(edgecolor="white", linewidth=2)
    )

    ax.set_title(
        "Distribution of FtsZ Systems Across Genomes",
        fontsize=16,
        pad=20
    )

    ax.legend(
        wedges,
        labels,
        title="System Type",
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        frameon=False
    )

    plt.tight_layout()

    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close()

############################################################
# MAIN
############################################################

def main():

    if USE_HARDCODED:

        input_csv = INPUT_CSV
        output_prefix = OUTPUT_PREFIX

    else:

        parser = argparse.ArgumentParser()

        parser.add_argument(
            "input_csv",
            help="Protein metadata CSV"
        )

        parser.add_argument(
            "--output_prefix",
            default="ftsZ_systems"
        )

        args = parser.parse_args()

        input_csv = args.input_csv
        output_prefix = args.output_prefix

    grouped = process_csv(input_csv)

    summary_file = f"{output_prefix}_genome_summary.csv"
    pie_file = f"{output_prefix}_piechart.png"

    grouped.to_csv(summary_file, index=False)

    plot_pie(grouped, pie_file)

    print("Genome summary written to:", summary_file)
    print("Pie chart written to:", pie_file)


if __name__ == "__main__":
    main()