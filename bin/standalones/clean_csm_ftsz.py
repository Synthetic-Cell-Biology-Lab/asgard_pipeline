import pandas as pd
import re

INPUT_CSV = "/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz_fin_v1/extraction_exploration/ftsz_csm.rev.csv"
OUTPUT_CSV = "/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz_fin_v1/extraction_exploration/ftsz_csm_clean.rev.csv"

def clean_text(s):
    if pd.isna(s):
        return s
    s = re.sub(r"_+", " ", str(s))   # replace multiple underscores
    s = re.sub(r"\s+", " ", s).strip()  # clean extra spaces
    s = re.sub(r" ", "_", s)
    return s

# Load CSV
df = pd.read_csv(INPUT_CSV, on_bad_lines="skip")

# Apply only to first column
first_col = df.columns[0]
df[first_col] = df[first_col].apply(clean_text)

# Save
df.to_csv(OUTPUT_CSV, index=False)

print(f"✅ Cleaned first column and saved → {OUTPUT_CSV}")