#!/usr/bin/env bash
set -euo pipefail

FAA=$1
FFN=$2

NR_FAA=$3
NR_FFN=$4

ALIGNED=$5
TRIMMED=$6
TRIMMED_COD=$7
CODON=$8

TREE_PREFIX=$9
THREADS=$10

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

# Build a lookup of translated sequences keyed by ID
# Keep only IDs where protein sequences match exactly
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
# 6. Trim codon alignment
########################################

echo "Running trimAl"

# trimal \
#     -in "$CODON" \
#     -out "$TRIMMED_COD" \
#     -automated1

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
    "$CLEAN_FFN" \
    -output fasta \
    > "$CODON"

########################################
# 7. Sanity check sequence counts
########################################

echo "Sanity check"

N_ALIGNED=$(grep -c "^>" "$ALIGNED")
N_CODON=$(grep -c "^>" "$CODON")
# N_TRIMMED_COD=$(grep -c "^>" "$TRIMMED_COD")
N_TRIMMED=$(grep -c "^>" "$TRIMMED")

echo "  Input (consistent): $N_GOOD"
echo "  Protein alignment:  $N_ALIGNED"
echo "  Codon alignment:    $N_CODON"
# echo "  Trimmed codon:      $N_TRIMMED_COD"
echo "  Trimmed:            $N_TRIMMED"

if [ "$N_ALIGNED" -ne "$N_CODON" ] || [ "$N_ALIGNED" -ne "$N_TRIMMED" ]; then
    echo "ERROR: sequence count mismatch across pipeline stages" >&2
    exit 1
fi

########################################
# 8. Phylogenetic tree
########################################

echo "Running IQ-TREE"

iqtree \
    -s "$TRIMMED" \
    -m MFP \
    -bb 1000 \
    -T AUTO \
    -pre "$TREE_PREFIX"

echo "Selection preparation complete"