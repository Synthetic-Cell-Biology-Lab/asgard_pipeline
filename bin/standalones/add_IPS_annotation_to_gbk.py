#!/usr/bin/env python3

import os
import argparse
import duckdb
from Bio import SeqIO


USE_HARDCODED = True

GBK_FOLDER = "/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/temp_synteny"
PARQUET_FILE = "/home/anirudh/asgard_pipeline/database/collated/Version1/filtered/85comp10con/IPS/protein_summary.parquet"


def annotate_gbk_folder(gbk_folder, parquet_file):

    print("📂 Connecting DuckDB...")
    con = duckdb.connect()

    print("📂 Registering parquet...")
    con.execute(f"""
    CREATE TABLE protein_summary AS
    SELECT * FROM '{parquet_file}'
    """)

    gbk_files = [
        os.path.join(gbk_folder, f)
        for f in os.listdir(gbk_folder)
        if f.endswith(".gbk")
    ]

    print(f"📁 Found {len(gbk_files)} GBK files")

    for gbk_path in gbk_files:

        print(f"🧬 Processing {os.path.basename(gbk_path)}")

        records = list(SeqIO.parse(gbk_path, "genbank"))

        proteins = []

        for rec in records:
            for feat in rec.features:

                if feat.type != "CDS":
                    continue

                locus = feat.qualifiers.get("locus_tag", [None])[0]

                if locus:
                    proteins.append(locus)

        if not proteins:
            continue

        print(f"   Found {len(proteins)} CDS proteins")

        ################################
        # DuckDB query
        ################################

        protein_list = ",".join([f"'{p}'" for p in proteins])

        result = con.execute(f"""
            SELECT protein, domains
            FROM protein_summary
            WHERE protein IN ({protein_list})
        """).fetchall()

        domain_lookup = {p: d for p, d in result if d is not None}

        print(f"   Matched {len(domain_lookup)} proteins with domains")

        ################################
        # Annotate CDS
        ################################

        for rec in records:
            for feat in rec.features:

                if feat.type != "CDS":
                    continue

                locus = feat.qualifiers.get("locus_tag", [None])[0]

                if locus in domain_lookup:

                    feat.qualifiers["domains"] = [str(domain_lookup[locus])]

        ################################
        # Write GBK
        ################################

        SeqIO.write(records, gbk_path, "genbank")

    print("✅ Domain annotation complete")


def main():

    if USE_HARDCODED:

        gbk_folder = GBK_FOLDER
        parquet_file = PARQUET_FILE

    else:

        parser = argparse.ArgumentParser()

        parser.add_argument("--gbk_folder", required=True)
        parser.add_argument("--parquet", required=True)

        args = parser.parse_args()

        gbk_folder = args.gbk_folder
        parquet_file = args.parquet

    annotate_gbk_folder(gbk_folder, parquet_file)


if __name__ == "__main__":
    main()