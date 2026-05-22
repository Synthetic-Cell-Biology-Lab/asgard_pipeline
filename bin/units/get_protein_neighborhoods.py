#!/usr/bin/env python3

import argparse
import os

import duckdb
import pandas as pd
from Bio import SeqIO


############################################################
# Load protein annotations from parquet
############################################################
def load_ips_annotations(proteins, parquet_file):

    if not proteins:
        return {}, {}, {}

    con = duckdb.connect()

    protein_df = pd.DataFrame({
        "protein": list(proteins)
    })

    con.register(
        "protein_list",
        protein_df
    )

    result = con.execute(f"""
        SELECT
            p.protein,
            COALESCE(s.ipr_acc, s.sig_acc)                           AS acc,
            COALESCE(s.ipr_desc, s.sig_desc)                         AS desc,
            CASE 
                WHEN LOWER(s.analysis) = 'pfam' 
                THEN s.sig_acc 
            END AS pfam_acc,

            CASE 
                WHEN LOWER(s.analysis) = 'pfam' 
                THEN s.sig_desc 
            END AS pfam_desc            
        FROM protein_list p
        LEFT JOIN read_parquet('{parquet_file}') s
        ON p.protein = s.protein
    """).fetchdf()

    ########################################################
    # Aggregate
    ########################################################

    acc_lookup = (
        result
        .dropna(subset=["acc"])
        .groupby("protein")["acc"]
        .apply(lambda x:
            "; ".join(sorted(set(x)))
        )
        .to_dict()
    )

    desc_lookup = (
        result
        .dropna(subset=["desc"])
        .groupby("protein")["desc"]
        .apply(lambda x:
            "; ".join(sorted(set(x)))
        )
        .to_dict()
    )

    pfam_lookup = (
        result
        .dropna(subset=["pfam_acc"])
        .groupby("protein")["pfam_acc"]
        .apply(lambda x:
            "; ".join(sorted(set(x)))
        )
        .to_dict()
    )

    pfam_desc_lookup = (
        result
        .dropna(subset=["pfam_acc", "pfam_desc"])
        .drop_duplicates(subset=["pfam_acc"])
        .set_index("pfam_acc")["pfam_desc"]
        .to_dict()
    )
    return acc_lookup, desc_lookup, pfam_lookup, pfam_desc_lookup


############################################################
# Build GBK index
############################################################

def build_gbk_index(base_path, max_depth=2):

    gbk_index = {}

    for root, dirs, files in os.walk(base_path):

        depth = root[len(base_path):].count(
            os.sep
        )

        if depth > max_depth:
            continue

        for f in files:

            if f.endswith((".gbk", ".gb")):

                genome_name = (
                    os.path.splitext(f)[0]
                    .replace(".domains", "")
                    .strip()
                )

                gbk_index[genome_name] = (
                    os.path.join(root, f)
                )

    print("\nIndexed GBKs:")
    print(sorted(gbk_index.keys()))

    return gbk_index


############################################################
# Extract neighborhood
############################################################

def extract_neighborhood(
    filtered_df,
    protein_name,
    gbk_index,
    window,
    subset_gbk=None
):

    rows = []

    ########################################################
    # Normalize genome names
    ########################################################

    filtered_df["genome_file"] = (
        filtered_df["genome_file"]
        .astype(str)
        .str.replace(
            ".domains",
            "",
            regex=False
        )
        .str.strip()
    )

    ########################################################
    # Optional subset
    ########################################################

    if subset_gbk:

        subset_list = [
            x.replace(".domains", "").strip()
            for x in subset_gbk.split(",")
        ]

        filtered_df = filtered_df[
            filtered_df["genome_file"]
            .isin(subset_list)
        ]

    ########################################################
    # Group genomes
    ########################################################

    grouped = filtered_df.groupby(
        "genome_file"
    )

    print("\nGenomes requested:")
    print(sorted(grouped.groups.keys()))

    ########################################################
    # Iterate genomes
    ########################################################

    for genome, gdf in grouped:

        gbk_file = gbk_index.get(genome)

        ####################################################
        # Missing GBK
        ####################################################

        if not gbk_file:

            print(
                f"⚠ Missing GBK for genome: {genome}"
            )

            continue

        print(
            f"\nProcessing genome:"
            f" {genome}"
        )

        ####################################################
        # Parse GBK
        ####################################################

        records = list(
            SeqIO.parse(
                gbk_file,
                "genbank"
            )
        )

        ####################################################
        # Normalize target loci
        ####################################################

        target_loci = set(
            gdf["locus_tag"]
            .astype(str)
            .str.strip()
        )

        print(
            f"Target loci count:"
            f" {len(target_loci)}"
        )

        ####################################################
        # Iterate records
        ####################################################

        for rec in records:

            cds_features = [
                f
                for f in rec.features
                if f.type == "CDS"
            ]

            ################################################
            # Iterate target genes
            ################################################

            for idx, feat in enumerate(
                cds_features
            ):

                locus = (
                    feat.qualifiers
                    .get("locus_tag", [None])[0]
                )

                if locus:

                    locus = str(locus).strip()

                if locus not in target_loci:
                    continue

                ################################################
                # Center gene coordinates
                ################################################

                center_strand = (
                    feat.location.strand
                )

                center_start = int(
                    feat.location.start
                )

                center_end = int(
                    feat.location.end
                )

                center_mid = (
                    center_start + center_end
                ) // 2

                ################################################
                # Neighborhood window
                ################################################

                region_start = max(
                    0,
                    center_start - window
                )

                region_end = (
                    center_end + window
                )

                ################################################
                # Iterate neighbors
                ################################################

                for n_idx, neighbor in enumerate(
                    cds_features
                ):

                    n_start = int(
                        neighbor.location.start
                    )

                    n_end = int(
                        neighbor.location.end
                    )

                    ################################################
                    # Window filter
                    ################################################

                    if (
                        n_end < region_start
                        or
                        n_start > region_end
                    ):
                        continue

                    n_locus = (
                        neighbor.qualifiers
                        .get("locus_tag", [""])[0]
                    )

                    n_locus = str(n_locus).strip()

                    ################################################
                    # Skip self
                    ################################################

                    if n_locus == locus:
                        continue

                    ################################################
                    # Metadata
                    ################################################

                    product = (
                        neighbor.qualifiers
                        .get("product", [""])[0]
                    )

                    n_strand = (
                        neighbor.location.strand
                    )

                    n_mid = (
                        n_start + n_end
                    ) // 2

                    ################################################
                    # Relative distances
                    ################################################

                    distance_bp = (
                        n_mid - center_mid
                    )

                    gene_offset = (
                        n_idx - idx
                    )

                    ################################################
                    # Strand normalization
                    ################################################

                    if center_strand == -1:

                        distance_bp = -distance_bp

                        gene_offset = -gene_offset

                    ################################################
                    # Store row
                    ################################################

                    row = {

                        "genome":
                            genome,

                        "center_protein":
                            protein_name,

                        "center_locus":
                            locus,

                        "neighbor_locus":
                            n_locus,

                        "neighbor_product":
                            product,

                        "start":
                            n_start,

                        "end":
                            n_end,

                        "center_start":
                            center_start,

                        "center_end":
                            center_end,

                        "center_mid":
                            center_mid,

                        "neighbor_mid":
                            n_mid,

                        "center_strand":
                            center_strand,

                        "neighbor_strand":
                            n_strand,

                        "same_strand":
                            center_strand == n_strand,

                        "distance_bp":
                            distance_bp,

                        "gene_offset":
                            gene_offset,
                    }

                    rows.append(row)

    return pd.DataFrame(rows)


############################################################
# MAIN
############################################################

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--filtered_proteins",
        required=True
    )

    parser.add_argument(
        "--protein_col",
        required=True
    )

    parser.add_argument(
        "--ips",
        required=True
    )

    parser.add_argument(
        "--protein",
        required=True
    )

    parser.add_argument(
        "--window",
        type=int,
        default=5000
    )

    parser.add_argument(
        "--gbk_path",
        required=True
    )

    parser.add_argument(
        "--out",
        required=True
    )

    parser.add_argument(
        "--subset_gbk",
        default=None
    )

    parser.add_argument(
        "--pfam_out",
        required=True
    )

    args = parser.parse_args()

    ########################################################
    # Load protein CSV
    ########################################################

    print(
        "📂 Loading filtered protein CSV"
    )

    df = pd.read_csv(
        args.filtered_proteins
    )

    ########################################################
    # Normalize protein labels
    ########################################################

    df[args.protein_col] = (
        df[args.protein_col]
        .astype(str)
        .str.strip()
        .str.lower()
    )

    protein_query = (
        str(args.protein)
        .strip()
        .lower()
    )

    ########################################################
    # Filter target protein class
    ########################################################

    df = df[
        df[args.protein_col]
        == protein_query
    ]

    print(
        f"\nProteins matched:"
        f" {len(df)}"
    )

    ########################################################
    # Empty case
    ########################################################

    if df.empty:

        print(
            f"⚠ No proteins found "
            f"for {args.protein}"
        )

        pd.DataFrame().to_csv(
            args.out,
            index=False
        )

        return

    ########################################################
    # Build GBK index
    ########################################################

    gbk_index = build_gbk_index(
        args.gbk_path
    )

    ########################################################
    # Extract neighborhoods
    ########################################################

    print(
        f"\n🔎 Extracting neighborhoods "
        f"for {args.protein}"
    )

    neighborhood_df = extract_neighborhood(
        df,
        args.protein,
        gbk_index,
        args.window,
        args.subset_gbk
    )

    ########################################################
    # Empty neighborhood
    ########################################################

    if neighborhood_df.empty:

        print(
            "⚠ No neighborhoods found"
        )

        neighborhood_df.to_csv(
            args.out,
            index=False
        )

        return

    ########################################################
    # Annotate InterPro
    ########################################################

    print(
        "\n🔎 Annotating domains "
        "from InterProScan"
    )

    proteins = set(
        neighborhood_df[
            "neighbor_locus"
        ]
    )

    acc_lookup, desc_lookup, pfam_lookup, pfam_desc_lookup = (
        load_ips_annotations(
            proteins,
            args.ips
        )
    )

    neighborhood_df["IPS_acc"] = (
        neighborhood_df["neighbor_locus"]
        .map(lambda x: acc_lookup.get(x, ""))
    )

    neighborhood_df["IPS_desc"] = (
        neighborhood_df["neighbor_locus"]
        .map(lambda x: desc_lookup.get(x, ""))
    )

    neighborhood_df["Pfam_acc"] = (
        neighborhood_df["neighbor_locus"]
        .map(lambda x: pfam_lookup.get(x, ""))
    )

    with open(args.pfam_out, 'w') as out:
        for k, v in pfam_desc_lookup.items():
            out.write(f"{k}\t{v}\n")

    ########################################################
    # Save
    ########################################################

    os.makedirs(
        os.path.dirname(args.out),
        exist_ok=True
    )

    neighborhood_df.to_csv(
        args.out,
        index=False
    )

    print(
        f"\n✅ Saved neighborhood table "
        f"→ {args.out}"
    )


if __name__ == "__main__":
    main()