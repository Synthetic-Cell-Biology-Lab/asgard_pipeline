#!/usr/bin/env bash
set -euo pipefail

FAA=$1
FFN=$2

NR_FAA=$3
NR_FFN=$4

ALIGNED=$5
TRIMMED=$6
CODON=$7

TREE_PREFIX=$8
THREADS=$9

########################################
# 1. Remove identical sequences
########################################

echo "Running CD-HIT"

cd-hit \
    -i "$FAA" \
    -o "$NR_FAA" \
    -c 0.99 \
    -n 5 \
    -T "$THREADS" \
    -M 0

########################################
# 2. Extract nucleotide sequences
########################################

echo "Extract nucleotide sequences"

seqkit grep \
    -f <(grep "^>" "$NR_FAA" | sed 's/>//') \
    "$FFN" \
    > "$NR_FFN"

########################################
# 3. MAFFT alignment (accurate mode)
########################################

echo "Running MAFFT"

mafft --localpair --maxiterate 1000 --thread "$THREADS" \
    "$NR_FAA" \
    > "$ALIGNED"

########################################
# 4. Trim alignment
########################################

echo "Running trimAl"

trimal \
    -in "$ALIGNED" \
    -out "$TRIMMED" \
    -automated1

########################################
# 5. Codon alignment
########################################

echo "Running PAL2NAL"

pal2nal.pl \
    "$TRIMMED" \
    "$NR_FFN" \
    -output fasta \
    > "$CODON"

########################################
# 6. Phylogenetic tree
########################################

echo "Running IQ-TREE"

iqtree \
    -s "$TRIMMED" \
    -m MFP \
    -bb 1000 \
    -nt "$THREADS" \
    -pre "$TREE_PREFIX"

echo "Selection preparation complete"