#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import duckdb
from Bio import SeqIO


############################################################
# Load protein annotations from parquet
############################################################

def load_ips_annotations(proteins, parquet_file):

    if not proteins:
        return {}, {}

    con = duckdb.connect()

    protein_df = pd.DataFrame({"protein": list(proteins)})
    con.register("protein_list", protein_df)

    result = con.execute(f"""
        SELECT
            p.protein,
            COALESCE(s.ipr_acc, s.sig_acc) AS acc,
            COALESCE(s.ipr_desc, s.sig_desc) AS desc
        FROM protein_list p
        LEFT JOIN read_parquet('{parquet_file}') s
        ON p.protein = s.protein
    """).fetchdf()

    acc_lookup = (
        result
        .dropna(subset=["acc"])
        .groupby("protein")["acc"]
        .apply(lambda x: "; ".join(sorted(set(x))))
        .to_dict()
    )

    desc_lookup = (
        result
        .dropna(subset=["desc"])
        .groupby("protein")["desc"]
        .apply(lambda x: "; ".join(sorted(set(x))))
        .to_dict()
    )

    return acc_lookup, desc_lookup


############################################################
# Extract neighborhood
############################################################

def extract_neighborhood(filtered_df, protein_name, protein_col, gbk_path, window):

    rows = []

    grouped = filtered_df.groupby("genome_file")

    for genome, gdf in grouped:

        gbk_file = os.path.join(gbk_path, genome, f"{genome}.gbk")

        if not os.path.exists(gbk_file):
            print(f"⚠ Missing {gbk_file}")
            continue

        records = list(SeqIO.parse(gbk_file, "genbank"))

        target_loci = set(gdf["locus_tag"])

        for rec in records:

            cds_features = [
                f for f in rec.features if f.type == "CDS"
            ]

            for idx, feat in enumerate(cds_features):

                locus = feat.qualifiers.get("locus_tag", [None])[0]

                if locus not in target_loci:
                    continue

                center_start = int(feat.location.start)
                center_end = int(feat.location.end)
                center_mid = (center_start + center_end) // 2

                region_start = max(0, center_start - window)
                region_end = center_end + window

                for n_idx, neighbor in enumerate(cds_features):

                    n_start = int(neighbor.location.start)
                    n_end = int(neighbor.location.end)

                    if n_end < region_start or n_start > region_end:
                        continue

                    n_locus = neighbor.qualifiers.get("locus_tag", [""])[0]

                    # skip center gene itself
                    if n_locus == locus:
                        continue

                    product = neighbor.qualifiers.get("product", [""])[0]

                    n_mid = (n_start + n_end) // 2
                    distance_bp = n_mid - center_mid

                    gene_offset = n_idx - idx

                    rows.append({
                        "genome": genome,
                        "center_protein": protein_name,
                        "center_locus": locus,

                        "neighbor_locus": n_locus,
                        "neighbor_product": product,

                        "start": n_start,
                        "end": n_end,
                        "strand": neighbor.location.strand,

                        "center_start": center_start,
                        "center_end": center_end,
                        "center_mid": center_mid,
                        "neighbor_mid": n_mid,

                        "distance_bp": distance_bp,
                        "gene_offset": gene_offset
                    })

    return pd.DataFrame(rows)


############################################################
# MAIN
############################################################

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--filtered_proteins", required=True)
    parser.add_argument("--protein_col", required=True)
    parser.add_argument("--ips", required=True)
    parser.add_argument("--protein", required=True)
    parser.add_argument("--window", type=int, default=5000)
    parser.add_argument("--gbk_path", required=True)
    parser.add_argument("--out", required=True)

    args = parser.parse_args()

    print("📂 Loading filtered protein CSV")

    df = pd.read_csv(args.filtered_proteins)

    ########################################################
    # Filter for the protein class
    ########################################################

    df = df[df[args.protein_col] == args.protein]

    if df.empty:
        print(f"⚠ No proteins found for {args.protein}")
        pd.DataFrame().to_csv(args.out, index=False)
        return

    ########################################################
    # Extract neighborhood
    ########################################################

    print(f"🔎 Extracting neighborhoods for {args.protein}")

    neighborhood_df = extract_neighborhood(
        df,
        args.protein,
        args.protein_col,
        args.gbk_path,
        args.window
    )

    if neighborhood_df.empty:
        print("⚠ No neighborhoods found")
        neighborhood_df.to_csv(args.out, index=False)
        return

    ########################################################
    # Annotate with InterProScan
    ########################################################

    print("🔎 Annotating domains from InterProScan")

    proteins = set(neighborhood_df["neighbor_locus"])

    acc_lookup, desc_lookup = load_ips_annotations(proteins, args.ips)

    neighborhood_df["IPS_acc"] = neighborhood_df["neighbor_locus"].map(
        lambda x: acc_lookup.get(x, "")
    )

    neighborhood_df["IPS_desc"] = neighborhood_df["neighbor_locus"].map(
        lambda x: desc_lookup.get(x, "")
    )

    ########################################################
    # Write output
    ########################################################

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    neighborhood_df.to_csv(args.out, index=False)

    print(f"✅ Saved neighborhood table → {args.out}")


if __name__ == "__main__":
    main()