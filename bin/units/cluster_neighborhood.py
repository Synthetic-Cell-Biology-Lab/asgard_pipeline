#!/usr/bin/env python3

import argparse
import pandas as pd
import networkx as nx



threshold = 0.4
############################################################
# Cluster based on shared IPS IDs
############################################################

def cluster_ips(ids_series):

    ids_series = (
        ids_series.fillna("")
        .astype(str)
    )

    # convert to sets of domains
    domain_sets = {
    idx: {
        x_clean
        for x in val.split(";")
        if (x_clean := x.strip().replace("-", "")) != ""
    }
    for idx, val in ids_series.items()
}


    ########################################################
    # Build graph
    ########################################################

    G = nx.Graph()

    for i, dom_i in domain_sets.items():
        G.add_node(i)

        for j, dom_j in domain_sets.items():

            if i >= j:
                continue

            intersection = len(dom_i & dom_j)
            union = len(dom_i | dom_j)

            if union == 0:
                continue

            jaccard = intersection / union

            if jaccard >= threshold:
                G.add_edge(i, j)

    ########################################################
    # Connected components = clusters
    ########################################################

    cluster_map = {}

    for cid, comp in enumerate(nx.connected_components(G)):
        for node in comp:
            cluster_map[node] = cid

    return cluster_map


############################################################
# Generate cluster labels
############################################################

def assign_cluster_labels(df, column, cluster_col):

    label_map = (
        df.groupby(cluster_col)[column]
        .first()
        .to_dict()
    )

    return df[cluster_col].map(label_map)


############################################################
# MAIN
############################################################

def main():

    parser = argparse.ArgumentParser(
        description="Cluster proteins by shared InterPro IDs"
    )

    parser.add_argument("--input", required=True)
    parser.add_argument("--column", required=True, help="Column containing IPS IDs")
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    print("📂 Loading dataframe")

    df = pd.read_csv(args.input)

    if args.column not in df.columns:
        raise ValueError(f"{args.column} not found in dataframe")

    ########################################################
    # Cluster IPS IDs
    ########################################################

    print("🔎 Clustering based on IPS IDs")

    cluster_map = cluster_ips(df[args.column])

    df["annotation_cluster"] = df.index.map(cluster_map)

    ########################################################
    # Assign labels
    ########################################################

    df["cluster_label"] = assign_cluster_labels(
        df,
        args.column,
        "annotation_cluster"
    )

    ########################################################
    # Save
    ########################################################

    df.to_csv(args.output, index=False)

    print(f"✅ Saved clustered annotations → {args.output}")


if __name__ == "__main__":
    main()