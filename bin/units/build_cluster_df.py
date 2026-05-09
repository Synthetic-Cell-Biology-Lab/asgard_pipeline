#!/usr/bin/env python3

import argparse
import pandas as pd


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        required=True,
        help="Input clustered neighborhood CSV"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output cluster summary CSV"
    )

    args = parser.parse_args()

    df = pd.read_csv(args.input)

    # remove rows without cluster labels
    df = df.dropna(subset=["cluster_label"])

    # remove empty labels
    df = df[
        df["cluster_label"].astype(str).str.strip().ne("")
    ]

    grouped = (
        df.groupby(["annotation_cluster", "cluster_label"])
          ["neighbor_locus"]
          .apply(list)
          .reset_index()
    )

    grouped["protein_count"] = grouped["neighbor_locus"].apply(len)

    grouped = grouped.rename(columns={
        "cluster_label": "cluster_label",
        "neighbor_locus": "protein_ids",
        "annotation_cluster": "cluster_id",
        
    })

    grouped = grouped[
        ["cluster_label", "protein_ids", "cluster_id", "protein_count"]
    ]

    # sort descending
    grouped = grouped.sort_values(
        by="protein_count",
        ascending=False
    )

    grouped.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()