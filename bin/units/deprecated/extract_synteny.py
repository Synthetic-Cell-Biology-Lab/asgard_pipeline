#!/usr/bin/env python3

import os
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq


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
# MERGE INTERVALS
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
# LOAD PROTEINS FOR A GENOME
############################################################

def load_locus_annotations(csv_file, genome):

    df = pd.read_csv(csv_file)

    df = df[df["genome_file"] == genome]

    if df.empty:
        return {}

    return dict(zip(df["locus_tag"], df["Manual_annotation"]))


############################################################
# PROCESS GENOME
############################################################

def process_genome(genome_id, locus_to_annotation, db_path, window):

    gbk_path = os.path.join(
        db_path,
        genome_id,
        f"{genome_id}.gbk"
    )

    if not os.path.exists(gbk_path):

        print(f"WARNING: Missing genome file {gbk_path}")
        return None

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

                    qualifiers["gene"] = [annotation]
                    qualifiers["product"] = [annotation]
                    qualifiers["final_annotation"] = [annotation]

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

            genome_pos += len(region_seq)

    if not seq_chunks:
        return None

    new_record = SeqRecord(
        Seq("".join(seq_chunks)),
        id=genome_id,
        name=genome_id,
        description="Synteny regions"
    )

    new_record.features = new_features
    new_record.annotations["molecule_type"] = "DNA"
    new_record.annotations["topology"] = "linear"

    return new_record


############################################################
# MAIN
############################################################

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--input_csv", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--database_path", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--window", type=int, default=5000)

    args = parser.parse_args()

    locus_to_annotation = load_locus_annotations(
        args.input_csv,
        args.genome
    )

    if not locus_to_annotation:

        print(f"No loci found for genome {args.genome}")
        return

    record = process_genome(
        args.genome,
        locus_to_annotation,
        args.database_path,
        args.window
    )

    if record:

        os.makedirs(os.path.dirname(args.output), exist_ok=True)

        SeqIO.write(record, args.output, "genbank")

        print(f"Saved: {args.output}")

    else:

        print(f"No synteny regions found for {args.genome}")


if __name__ == "__main__":
    main()