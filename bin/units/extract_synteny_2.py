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
# BUILD ONE SEQRECORD FROM A GROUP OF INTERVALS
############################################################

def build_record(genome_id, contig_idx, intervals, rec_map, locus_to_annotation):
    """
    Given a list of (region_start, region_end, source_record) tuples,
    concatenate their sequences and remap all CDS features into the
    new coordinate space.  Returns a SeqRecord.
    """
    seq_chunks = []
    new_features = []
    genome_pos = 0

    for region_start, region_end, rec in intervals:
        region_seq = rec.seq[region_start:region_end]
        seq_chunks.append(str(region_seq))

        offset = genome_pos - region_start

        for feat in rec.features:
            if feat.type != "CDS":
                continue

            f_start = int(feat.location.start)
            f_end   = int(feat.location.end)

            if f_end < region_start or f_start > region_end:
                continue

            new_start = max(f_start, region_start) + offset
            new_end   = min(f_end,   region_end)   + offset

            locus      = feat.qualifiers.get("locus_tag", [""])[0]
            qualifiers = dict(feat.qualifiers)
            annotation = locus_to_annotation.get(locus)

            if annotation:
                qualifiers["gene"]             = [annotation]
                qualifiers["product"]          = [annotation]
                qualifiers["final_annotation"] = [annotation]
                if annotation in COLOR_MAP:
                    qualifiers["colour"] = [COLOR_MAP[annotation]]

            new_features.append(SeqFeature(
                FeatureLocation(new_start, new_end, strand=feat.location.strand),
                type="CDS",
                qualifiers=qualifiers,
            ))

        genome_pos += len(region_seq)

    record_id = f"{genome_id}_ctg{contig_idx + 1}"
    new_record = SeqRecord(
        Seq("".join(seq_chunks)),
        id=record_id,
        name=record_id,
        description="Synteny region",
    )
    new_record.features = new_features
    new_record.annotations["molecule_type"] = "DNA"
    new_record.annotations["topology"]      = "linear"
    return new_record


############################################################
# PROCESS GENOME  →  returns a list of SeqRecords
############################################################

def process_genome(genome_id, locus_to_annotation, db_path, window, gap_threshold):
    """
    Extracts synteny windows, merges overlapping ones, then splits into
    separate contigs wherever the gap between consecutive intervals
    exceeds gap_threshold.  Returns one SeqRecord per contig.
    """
    gbk_path = os.path.join(db_path, genome_id, f"{genome_id}.gbk")

    if not os.path.exists(gbk_path):
        print(f"WARNING: Missing genome file {gbk_path}")
        return []

    records     = list(SeqIO.parse(gbk_path, "genbank"))
    target_loci = set(locus_to_annotation.keys())

    # --- collect (region_start, region_end, source_record) per chromosome ---
    raw_intervals = []   # list of (start, end, rec)

    for rec in records:
        windows = []
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            locus = feat.qualifiers.get("locus_tag", [None])[0]
            if locus in target_loci:
                s = int(feat.location.start)
                e = int(feat.location.end)
                windows.append((max(0, s - window), min(len(rec.seq), e + window)))

        for start, end in merge_intervals(windows):
            raw_intervals.append((start, end, rec))

    if not raw_intervals:
        return []

    # --- split into contig groups based on gap_threshold ---
    contig_groups   = []   # list of lists of (start, end, rec)
    current_group   = [raw_intervals[0]]

    for i in range(1, len(raw_intervals)):
        prev_start, prev_end, prev_rec = raw_intervals[i - 1]
        cur_start,  cur_end,  cur_rec  = raw_intervals[i]

        # Gap is defined between the end of the previous interval and the
        # start of the current one.  Intervals on different chromosomes are
        # always treated as separate contigs.
        gap = cur_start - prev_end
        same_chrom = (prev_rec.id == cur_rec.id)

        if same_chrom and gap <= gap_threshold:
            current_group.append(raw_intervals[i])
        else:
            contig_groups.append(current_group)
            current_group = [raw_intervals[i]]

    contig_groups.append(current_group)   # flush last group

    # --- build one SeqRecord per contig group ---
    seqrecords = []
    for idx, group in enumerate(contig_groups):
        rec = build_record(genome_id, idx, group, None, locus_to_annotation)
        seqrecords.append(rec)

    return seqrecords


############################################################
# MAIN
############################################################

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_csv",      required=True)
    parser.add_argument("--genome",         required=True)
    parser.add_argument("--database_path",  required=True)
    parser.add_argument("--output",         required=True)
    parser.add_argument("--window",         type=int, default=5000)
    parser.add_argument("--gap_threshold",  type=int, default=0,
                        help="Maximum gap (bp) between intervals before "
                             "splitting into a new contig (default: 0 = "
                             "always split on any gap)")
    args = parser.parse_args()

    locus_to_annotation = load_locus_annotations(args.input_csv, args.genome)

    if not locus_to_annotation:
        print(f"No loci found for genome {args.genome}")
        return

    records = process_genome(
        args.genome,
        locus_to_annotation,
        args.database_path,
        args.window,
        args.gap_threshold,
    )

    if records:
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        SeqIO.write(records, args.output, "genbank")
        print(f"Saved {len(records)} contig(s) → {args.output}")
    else:
        print(f"No synteny regions found for {args.genome}")


if __name__ == "__main__":
    main()