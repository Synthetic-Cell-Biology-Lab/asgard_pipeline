#!/usr/bin/env python3

import argparse
import duckdb
import pandas as pd
from Bio import SeqIO


############################################################
# ANNOTATE ONE GBK FILE WITH DOMAIN DATA
############################################################

def annotate_gbk(gbk_input, gbk_output, parquet_file):

    print("📂 Connecting DuckDB...")
    con = duckdb.connect()

    print(f"🧬 Reading {gbk_input}")

    records = list(SeqIO.parse(gbk_input, "genbank"))

    proteins = []

    ############################################################
    # Extract CDS locus tags
    ############################################################

    for rec in records:
        for feat in rec.features:

            if feat.type != "CDS":
                continue

            locus = feat.qualifiers.get("locus_tag", [None])[0]

            if locus:
                proteins.append(locus)

    if not proteins:

        print("⚠️ No CDS proteins found")
        SeqIO.write(records, gbk_output, "genbank")
        return

    print(f"   Found {len(proteins)} CDS proteins")

    ############################################################
    # Create dataframe of proteins
    ############################################################

    protein_df = pd.DataFrame({"protein": proteins})

    con.register("protein_list", protein_df)

    ############################################################
    # Query parquet directly
    ############################################################

    print("🔎 Querying domain annotations...")

    result = con.execute(f"""
        SELECT p.protein, s.domains
        FROM protein_list p
        LEFT JOIN read_parquet('{parquet_file}') s
        ON p.protein = s.protein
    """).fetchall()

    domain_lookup = {
        protein: domains
        for protein, domains in result
        if domains is not None
    }

    print(f"   Matched {len(domain_lookup)} proteins with domains")

    ############################################################
    # Annotate CDS features
    ############################################################

    for rec in records:
        for feat in rec.features:

            if feat.type != "CDS":
                continue

            locus = feat.qualifiers.get("locus_tag", [None])[0]

            if locus in domain_lookup:

                feat.qualifiers["domains"] = [str(domain_lookup[locus])]

    ############################################################
    # Write output GBK
    ############################################################

    SeqIO.write(records, gbk_output, "genbank")

    print(f"✅ Saved annotated file → {gbk_output}")


############################################################
# MAIN
############################################################

def main():

    parser = argparse.ArgumentParser(
        description="Annotate GenBank CDS features with InterProScan domain data"
    )

    parser.add_argument(
        "--input_gbk",
        required=True,
        help="Input GenBank file"
    )

    parser.add_argument(
        "--output_gbk",
        required=True,
        help="Output annotated GenBank file"
    )

    parser.add_argument(
        "--parquet",
        required=True,
        help="Protein summary parquet file containing domain annotations"
    )

    args = parser.parse_args()

    annotate_gbk(
        args.input_gbk,
        args.output_gbk,
        args.parquet
    )


if __name__ == "__main__":
    main()