#!/usr/bin/env bash

TARGET_DIR="/home/anirudh/asgard_pipeline/database/structures/SepF/comparison/Target"
REF_DIR="/home/anirudh/asgard_pipeline/database/structures/SepF/comparison/Reference"
OUT_DIR="/home/anirudh/asgard_pipeline/database/structures/SepF/comparison/Output"

mkdir -p "$OUT_DIR"

# Output CSV for parsed values
CSV_FILE="$OUT_DIR/usalign_summary.csv"
echo "target,reference,TMscore1,TMscore2,RMSD,Aligned_length,Seq_ID" > "$CSV_FILE"

for target in "$TARGET_DIR"/*; do
    tbase=$(basename "$target")
    tname="${tbase%.*}"

    for ref in "$REF_DIR"/*; do
        rbase=$(basename "$ref")
        rname="${rbase%.*}"

        prefix="$OUT_DIR/${tname}_vs_${rname}"
        outfile="${prefix}_alignment.txt"

        # Run USalign
        USalign "$target" "$ref" -o "$prefix" > "$outfile"

        # --- Extract values ---
        TM1=$(grep -m1 "TM-score=" "$outfile" | head -n1 | awk '{print $2}')
        TM2=$(grep -m2 "TM-score=" "$outfile" | tail -n1 | awk '{print $2}')

        RMSD=$(grep "RMSD=" "$outfile" | head -n1 | sed -E 's/.*RMSD= *([0-9.]+).*/\1/')
        ALN_LEN=$(grep "Aligned length=" "$outfile" | sed -E 's/Aligned length= *([0-9]+).*/\1/')
        SEQ_ID=$(grep "Seq_ID=" "$outfile" | sed -E 's/.*Seq_ID= *([0-9.]+).*/\1/')

        # Append to CSV
        echo "$tname,$rname,$TM1,$TM2,$RMSD,$ALN_LEN,$SEQ_ID" >> "$CSV_FILE"

    done
done

echo "Done. Results in $OUT_DIR"