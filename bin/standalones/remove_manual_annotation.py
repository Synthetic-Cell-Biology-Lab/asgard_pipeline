#!/usr/bin/env python3

import pandas as pd

# ============================================================
# HARDCODED PATHS
# ============================================================

INPUT_TABLE = "/home/anirudh/asgard_pipeline/database/collated/Version1/filtered/85comp10con/protein_file/jan2026_85comp10con_pf.csv"


# dictionary: annotation_name -> file containing protein ids
tags_to_remove = ["DipA"]

# ============================================================
# LOAD TABLE
# ============================================================

df = pd.read_csv(INPUT_TABLE)


# ensure column exists
if "Manual_annotation" not in df.columns:
    df["Manual_annotation"] = ""

# ============================================================
# REMOVE ANNOTATIONS
# ============================================================
for tag in tags_to_remove:
    df["Manual_annotation"] = df["Manual_annotation"].replace(tag, None)


# ============================================================
# EXPORT
# ============================================================

df.to_csv(INPUT_TABLE, index=False)
entries = df["Manual_annotation"].unique()
print(entries)
print("Annotated table written to:", INPUT_TABLE)
