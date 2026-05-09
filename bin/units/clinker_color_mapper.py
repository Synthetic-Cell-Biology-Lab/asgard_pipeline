#!/usr/bin/env python3

import argparse
import pandas as pd
import colorcet as cc
import matplotlib.colors as mcolors

############################################################
# MAIN
############################################################

def main():

    parser = argparse.ArgumentParser(
        description="Generate clinker color map"
    )

    parser.add_argument(
        "--classified",
        required=True,
        help="Classified dataframe CSV"
    )

    parser.add_argument(
        "--synteny",
        required=True,
        help="Synteny CSV containing Manual_annotation"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output color mapping CSV"
    )

    args = parser.parse_args()

    ########################################################
    # Load files
    ########################################################

    classified_df = pd.read_csv(args.classified)

    synteny_df = pd.read_csv(args.synteny)

    ########################################################
    # Clean strings
    ########################################################

    synteny_df["locus_tag"] = (
        synteny_df["locus_tag"]
        .astype(str)
        .str.strip()
    )

    classified_df["neighbor_locus"] = (
        classified_df["neighbor_locus"]
        .astype(str)
        .str.strip()
    )

    ########################################################
    # Manual annotations take priority
    ########################################################

    manual_df = synteny_df[
        ["locus_tag", "Manual_annotation"]
    ].dropna()

    manual_df = manual_df.rename(
        columns={
            "Manual_annotation": "function"
        }
    )

    ########################################################
    # Track manually annotated proteins
    ########################################################

    manual_locus_tags = set(
        manual_df["locus_tag"]
    )

    ########################################################
    # Add classifier annotations
    ########################################################

    classified_subset = classified_df[
        ~classified_df["neighbor_locus"].isin(
            manual_locus_tags
        )
    ][["neighbor_locus", "functional_bin"]]

    classified_subset = classified_subset.rename(
        columns={
            "neighbor_locus": "locus_tag",
            "functional_bin": "function"
        }
    )

    ########################################################
    # Combine
    ########################################################

    combined = pd.concat(
        [
            manual_df,
            classified_subset
        ],
        ignore_index=True
    )

    ########################################################
    # Clean labels
    ########################################################

    combined["function"] = (
        combined["function"]
        .astype(str)
        .str.strip()
        .str.lower()
    )

    ########################################################
    # Remove duplicates
    ########################################################

    combined = combined.drop_duplicates(
        subset=["locus_tag"],
        keep="first"
    )

    ########################################################
    # Generate Glasbey colors
    ########################################################

    unique_functions = sorted(
        combined["function"].unique()
    )

    glasbey = cc.glasbey_bw

    color_map = {
        func: mcolors.to_hex(
            glasbey[i % len(glasbey)]
        )
        for i, func in enumerate(unique_functions)
    }

    combined["color"] = (
        combined["function"]
        .map(color_map)
    )

    ########################################################
    # Export only locus_tag + color
    ########################################################

    combined.to_csv(
        args.output,
        index=False
    )

    ########################################################
    # Summary
    ########################################################

    print("\nGenerated color groups:\n")

    for func, color in color_map.items():
        print(f"{func:30s} {color}")

    print(f"\nSaved:")
    print(args.output)


if __name__ == "__main__":
    main()