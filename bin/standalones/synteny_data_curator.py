#!/usr/bin/env python3

import os
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq

############################################################
# DEFAULT CONFIG (used if CLI not provided)
############################################################

GENOME_IDS = [
    "GCA_041410205.1_ASM4141020v1_genomic",
    "GCA_008000775.2_ASM800077v2_genomic",
    "GCA_001940665.2_ASM194066v2_genomic"
]

TARGET_CATEGORIES = [
    "Tubulin",
    "FtsZ1",
    "FtsZ2",
    "CetZ",
    "SepF"
]

WINDOW = 5000

GENOME_METADATA = "/home/anirudh/asgard_pipeline/database/collated/Version1/filtered/85comp10con/genome_file/jan2026_85comp10con_gf.csv"
PROTEIN_METADATA = "/home/anirudh/asgard_pipeline/database/collated/Version1/filtered/85comp10con/protein_file/jan2026_85comp10con_pf.csv"

DATABASE_PATH = "/home/anirudh/asgard_pipeline/database/cds_genomes"
OUTPUT_DIR = "/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/temp_synteny"

############################################################
# COLOR MAP
############################################################

COLOR_MAP = {
    "FtsZ1": "255 0 0",
    "FtsZ2": "255 120 0",
    "Tubulin": "0 120 255",
    "CetZ": "0 200 100",
    "SepF": "200 0 200"
}

############################################################
# INTERVAL MERGE
############################################################

def merge_intervals(intervals):

    if not intervals:
        return []

    intervals = sorted(intervals)
    merged = [intervals[0]]

    for start, end in intervals[1:]:

        last_start, last_end = merged[-1]

        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))

    return merged


############################################################
# LOAD METADATA
############################################################

def load_metadata(genome_ids, genome_meta, protein_meta, categories):

    genome_df = pd.read_csv(genome_meta)
    protein_df = pd.read_csv(protein_meta)

    genome_df = genome_df[genome_df["genome_file"].isin(genome_ids)]
    protein_df = protein_df[protein_df["genome_file"].isin(genome_ids)]

    protein_df = protein_df[
        protein_df["Manual_annotation"].isin(categories)
    ]

    merged = pd.merge(protein_df, genome_df, on="genome_file")

    locus_to_annotation = dict(
        zip(merged["locus_tag"], merged["Manual_annotation"])
    )

    return merged, locus_to_annotation


############################################################
# EXTRACT SYNTENY
############################################################

def process_genome(genome_id, locus_to_annotation, combined_rows, db_path, window):

    gbk_path = os.path.join(
        db_path,
        genome_id,
        f"{genome_id}.gbk"
    )

    if not os.path.exists(gbk_path):
        print("Missing:", gbk_path)
        return

    records = list(SeqIO.parse(gbk_path, "genbank"))

    target_loci = set(locus_to_annotation.keys())

    seq_chunks = []
    new_features = []

    genome_pos = 0

    for rec in records:

        intervals = []

        for feat in rec.features:

            if feat.type != "CDS":
                continue

            locus = feat.qualifiers.get("locus_tag", [None])[0]

            if locus in target_loci:

                start = int(feat.location.start)
                end = int(feat.location.end)

                region_start = max(0, start - window)
                region_end = min(len(rec.seq), end + window)

                intervals.append((region_start, region_end))

        intervals = merge_intervals(intervals)

        for region_start, region_end in intervals:

            region_seq = rec.seq[region_start:region_end]
            seq_chunks.append(str(region_seq))

            offset = genome_pos - region_start

            for feat in rec.features:

                if feat.type != "CDS":
                    continue

                f_start = int(feat.location.start)
                f_end = int(feat.location.end)

                if f_end < region_start or f_start > region_end:
                    continue

                new_start = max(f_start, region_start) + offset
                new_end = min(f_end, region_end) + offset

                locus = feat.qualifiers.get("locus_tag", [""])[0]

                qualifiers = dict(feat.qualifiers)

                annotation = locus_to_annotation.get(locus)

                if annotation:

                    qualifiers["final_annotation"] = [annotation]
                    qualifiers["gene"] = [annotation]
                    qualifiers["product"] = [annotation]

                    if annotation in COLOR_MAP:
                        qualifiers["colour"] = [COLOR_MAP[annotation]]

                new_feat = SeqFeature(
                    FeatureLocation(
                        new_start,
                        new_end,
                        strand=feat.location.strand
                    ),
                    type="CDS",
                    qualifiers=qualifiers
                )

                new_features.append(new_feat)

                combined_rows.append({
                    "genome": genome_id,
                    "locus_tag": locus,
                    "start": f_start,
                    "end": f_end,
                    "strand": feat.location.strand,
                    "annotation": annotation if annotation else "",
                    "product": qualifiers.get("product", [""])[0]
                })

            genome_pos += len(region_seq)

    new_record = SeqRecord(
        Seq("".join(seq_chunks)),
        id=genome_id,
        name=genome_id,
        description="FtsZ system synteny regions"
    )

    new_record.features = new_features
    new_record.annotations["molecule_type"] = "DNA"
    new_record.annotations["topology"] = "linear"

    return new_record


############################################################
# MAIN
############################################################

def main():

    parser = argparse.ArgumentParser(description="Extract synteny regions around FtsZ systems")

    parser.add_argument("--genomes", nargs="+", default=GENOME_IDS)
    parser.add_argument("--categories", nargs="+", default=TARGET_CATEGORIES)
    parser.add_argument("--genome_metadata", default=GENOME_METADATA)
    parser.add_argument("--protein_metadata", default=PROTEIN_METADATA)
    parser.add_argument("--database_path", default=DATABASE_PATH)
    parser.add_argument("--window", type=int, default=WINDOW)
    parser.add_argument("--output_dir", default=OUTPUT_DIR)

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    merged, locus_to_annotation = load_metadata(
        args.genomes,
        args.genome_metadata,
        args.protein_metadata,
        args.categories
    )

    combined_rows = []

    for genome in args.genomes:

        record = process_genome(
            genome,
            locus_to_annotation,
            combined_rows,
            args.database_path,
            args.window
        )

        if record:

            out_path = os.path.join(
                args.output_dir,
                f"{genome}.ftsz_synteny.gbk"
            )

            SeqIO.write(record, out_path, "genbank")

    combined_df = pd.DataFrame(combined_rows)

    combined_df.to_csv(
        os.path.join(args.output_dir, "all_synteny_proteins.csv"),
        index=False
    )

    print("Finished.")


if __name__ == "__main__":
    main()