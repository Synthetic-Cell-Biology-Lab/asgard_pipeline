#!/usr/bin/env python3

import pandas as pd

# ============================================================
# HARDCODED PATHS
# ============================================================

INPUT_TABLE = "/home/anirudh/asgard_pipeline/database/collated/Version1/filtered/85comp10con/protein_file/jan2026_85comp10con_pf.csv"


# dictionary: annotation_name -> file containing protein ids
annotation_files = {
    "Dynein_like" : "/home/anirudh/asgard_pipeline/database/protein_sets/protein_ids/Dynein.ids"
}

# ============================================================
# LOAD TABLE
# ============================================================

df = pd.read_csv(INPUT_TABLE)

# ensure column exists
if "Manual_annotation" not in df.columns:
    df["Manual_annotation"] = ""

# ============================================================
# APPLY ANNOTATIONS
# ============================================================

for annotation, file in annotation_files.items():

    with open(file) as f:
        protein_ids = set(line.strip().replace(" ", "_") for line in f if line.strip())

    df.loc[df["locus_tag"].isin(protein_ids), "Manual_annotation"] = annotation

# ============================================================
# EXPORT
# ============================================================

df.to_csv(INPUT_TABLE, index=False)

print("Annotated table written to:", INPUT_TABLE)