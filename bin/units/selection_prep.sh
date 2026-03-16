#!/usr/bin/env bash
set -euo pipefail

FAA=$1
FFN=$2

NR_FAA=$3
NR_FFN=$4

ALIGNED=$5
TRIMMED=$6
# TRIMMED_COD=$7
CODON=$7

TREE_PREFIX=$8
THREADS=$9
COL_NUMBERING=${10}

# Directory of this script — used to locate trim_codon_aln.py
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

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
# 2. Extract nucleotide sequences,
#    strip terminal stop codon,
#    and clean headers to locus tag only
########################################

echo "Extracting and cleaning nucleotide sequences"

seqkit seq -ni "$NR_FAA" \
    | seqkit grep -f - "$FFN" \
    | seqkit subseq --region 1:-4 \
    | seqkit replace -p " .*" -r "" \
    > "$NR_FFN"

########################################
# 3. Translation consistency filter
########################################

echo "Checking FAA/FFN consistency"

TRANSLATED=$(mktemp --suffix=".faa")
GOOD_IDS=$(mktemp --suffix=".txt")
CLEAN_FAA=$(mktemp --suffix=".faa")
CLEAN_FFN=$(mktemp --suffix=".ffn")
trap 'rm -f "$TRANSLATED" "$GOOD_IDS" "$CLEAN_FAA" "$CLEAN_FFN"' EXIT

seqkit translate "$NR_FFN" \
    | seqkit replace -p " .*" -r "" \
    > "$TRANSLATED"

python3 - "$NR_FAA" "$TRANSLATED" "$GOOD_IDS" << 'EOF'
import sys
from Bio import SeqIO

faa   = {r.id: str(r.seq).rstrip("*") for r in SeqIO.parse(sys.argv[1], "fasta")}
trans = {r.id: str(r.seq).rstrip("*") for r in SeqIO.parse(sys.argv[2], "fasta")}

good = []
bad  = []
for seq_id, prot in faa.items():
    if trans.get(seq_id) == prot:
        good.append(seq_id)
    else:
        bad.append(seq_id)

with open(sys.argv[3], "w") as f:
    f.write("\n".join(good) + "\n")

print(f"  Consistent:   {len(good)}")
print(f"  Inconsistent: {len(bad)} (removed)")
if bad:
    print("  Removed IDs:")
    for i in bad:
        print(f"    {i}")
EOF

N_GOOD=$(wc -l < "$GOOD_IDS")
if [ "$N_GOOD" -lt 4 ]; then
    echo "ERROR: fewer than 4 consistent sequences remain — aborting" >&2
    exit 1
fi

seqkit grep -f "$GOOD_IDS" "$NR_FAA" > "$CLEAN_FAA"
seqkit grep -f "$GOOD_IDS" "$NR_FFN" > "$CLEAN_FFN"

########################################
# 4. MAFFT alignment (accurate mode)
########################################

echo "Running MAFFT"

mafft --localpair --maxiterate 1000 --thread -1 \
    "$CLEAN_FAA" \
    > "$ALIGNED"

########################################
# 5. Trim protein alignment,
#    saving kept column indices
########################################

echo "Running trimAl"

trimal \
    -in "$ALIGNED" \
    -out "$TRIMMED" \
    -automated1 \
    -colnumbering > "$COL_NUMBERING"

########################################
# 6. Codon alignment using UNTRIMMED
#    protein alignment
########################################

echo "Running PAL2NAL"

pal2nal.pl \
    "$ALIGNED" \
    "$CLEAN_FFN" \
    -output fasta \
    > "$CODON"

########################################
# 7. Trim codon alignment using
#    protein column indices
########################################

# echo "Trimming codon alignment"

# python3 "$SCRIPT_DIR/trim_codon_aln.py" \
#     "$CODON" \
#     "${TRIMMED}.cols" \
#     "$TRIMMED_COD"

########################################
# 8. Sanity check sequence counts
########################################

echo "Sanity check"

N_ALIGNED=$(grep -c "^>" "$ALIGNED")
N_TRIMMED=$(grep -c "^>" "$TRIMMED")
N_CODON=$(grep -c "^>" "$CODON")
# N_TRIMMED_COD=$(grep -c "^>" "$TRIMMED_COD")

echo "  Input (consistent):  $N_GOOD"
echo "  Protein alignment:   $N_ALIGNED"
echo "  Trimmed protein:     $N_TRIMMED"
echo "  Codon alignment:     $N_CODON"
# echo "  Trimmed codon:       $N_TRIMMED_COD"

if [ "$N_ALIGNED" -ne "$N_TRIMMED" ] || \
   [ "$N_ALIGNED" -ne "$N_CODON" ]; then
    echo "ERROR: sequence count mismatch across pipeline stages" >&2
    exit 1
fi

# Verify trimmed codon length is a multiple of 3
# COD_LEN=$(awk '/^>/{next} {gsub(/-/,""); printf $0} END{print ""}' "$TRIMMED_COD" \
#     | head -1 | tr -d '\n' | wc -c)
# if [ $(( COD_LEN % 3 )) -ne 0 ]; then
#     echo "ERROR: trimmed codon alignment length ($COD_LEN) is not a multiple of 3" >&2
#     exit 1
# fi
# echo "  Codon length check:  OK ($COD_LEN nt, $(( COD_LEN / 3 )) codons)"

########################################
# 9. Phylogenetic tree
########################################

echo "Running IQ-TREE"

iqtree \
    -s "$TRIMMED" \
    -m MFP \
    -bb 1000 \
    -T AUTO \
    -pre "$TREE_PREFIX"

echo "Selection preparation complete"