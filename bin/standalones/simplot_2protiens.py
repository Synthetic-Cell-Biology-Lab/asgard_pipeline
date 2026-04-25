#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO, AlignIO
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import tempfile
import argparse
import os

############################################################
# Arguments
############################################################

parser = argparse.ArgumentParser(description="Protein1 vs Protein2 similarity profile")

parser.add_argument("--metadata", required=True)
parser.add_argument("--fasta", required=True)
parser.add_argument("--protein_col", default="Manual_annotation")
parser.add_argument("--genome_col", default="genome_file")
parser.add_argument("--seqid_col", default="locus_tag")

parser.add_argument("--protein1", required=True)
parser.add_argument("--protein2", required=True)

parser.add_argument("--nterm_len", type=int, default=150)
parser.add_argument("--use_half_protein2", action="store_true",
                    help="Use first half of protein2 (e.g., for MreB-like proteins)")

parser.add_argument("--output_prefix", default="protein_compare")

args = parser.parse_args()

############################################################
# Load metadata
############################################################

df = pd.read_csv(args.metadata)

############################################################
# Filter genomes containing BOTH proteins
############################################################

g1 = set(df[df[args.protein_col] == args.protein1][args.genome_col])
g2 = set(df[df[args.protein_col] == args.protein2][args.genome_col])

common_genomes = g1.intersection(g2)

print(f"Genomes with both proteins: {len(common_genomes)}")

df = df[df[args.genome_col].isin(common_genomes)]

############################################################
# Extract IDs
############################################################

p1_ids = set(
    df[df[args.protein_col] == args.protein1][args.seqid_col]
)

p2_ids = set(
    df[df[args.protein_col] == args.protein2][args.seqid_col]
)

print(f"{args.protein1} sequences: {len(p1_ids)}")
print(f"{args.protein2} sequences: {len(p2_ids)}")

############################################################
# Load FASTA
############################################################

seq_dict = {rec.id: rec for rec in SeqIO.parse(args.fasta, "fasta")}

############################################################
# Extract regions
############################################################

def get_nterm(seq, length):
    return seq[:length] if len(seq) >= length else seq

p1_seqs = []
p2_seqs = []

for sid in p1_ids:
    if sid in seq_dict:
        p1_seqs.append((sid, get_nterm(seq_dict[sid].seq, args.nterm_len)))

for sid in p2_ids:
    if sid in seq_dict:
        seq = seq_dict[sid].seq

        if args.use_half_protein2:
            seq = seq[:len(seq)//2]

        p2_seqs.append((sid, get_nterm(seq, args.nterm_len)))

print(f"Total sequences retained: {len(p1_seqs) + len(p2_seqs)}")

############################################################
# Write temp FASTA
############################################################

tmp_fasta = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta").name

with open(tmp_fasta, "w") as f:
    for sid, seq in p1_seqs + p2_seqs:
        f.write(f">{sid}\n{seq}\n")

############################################################
# Run MAFFT
############################################################

aligned_file = tmp_fasta + ".aln.fasta"

print("Running MAFFT...")

with open(aligned_file, "w") as out:
    subprocess.run(
        ["mafft", "--auto", tmp_fasta],
        stdout=out,
        stderr=subprocess.DEVNULL
    )

############################################################
# Load alignment
############################################################

alignment = AlignIO.read(aligned_file, "fasta")

############################################################
# Define groups
############################################################

p1_ids = set([sid for sid, _ in p1_seqs])
p2_ids = set([sid for sid, _ in p2_seqs])

############################################################
# Similarity calculation
############################################################

scores = []

for col in range(alignment.get_alignment_length()):

    p1_res = []
    p2_res = []

    for rec in alignment:
        if rec.id in p1_ids:
            p1_res.append(rec.seq[col])
        elif rec.id in p2_ids:
            p2_res.append(rec.seq[col])

    p1_res = [r for r in p1_res if r != "-"]
    p2_res = [r for r in p2_res if r != "-"]

    if not p1_res or not p2_res:
        scores.append(np.nan)
        continue

    matches = 0
    total = 0

    for a in p1_res:
        for b in p2_res:
            total += 1
            if a == b:
                matches += 1

    scores.append(matches / total)

############################################################
# Smooth
############################################################

smooth_scores = pd.Series(scores).rolling(10, min_periods=1).mean()

############################################################
# Plot
############################################################

plt.figure(figsize=(12,4))
plt.plot(smooth_scores, linewidth=2)

plt.xlabel("Alignment Position")
plt.ylabel("Mean Identity")
plt.title(f"{args.protein1} vs {args.protein2} similarity (N-term)")

plt.ylim(0,1)
plt.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(f"{args.output_prefix}_similarity_profile.png", dpi=300)

############################################################
# Save alignment
############################################################

os.rename(aligned_file, f"{args.output_prefix}_aligned.fasta")

print("Done.")