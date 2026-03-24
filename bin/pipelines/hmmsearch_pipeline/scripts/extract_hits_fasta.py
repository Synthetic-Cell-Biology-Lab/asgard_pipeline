import sys
import pandas as pd
from Bio import SeqIO

input_csv = sys.argv[1]
input_fasta = sys.argv[2]
output_fasta = sys.argv[3]

df = pd.read_csv(input_csv)

# ---- Column detection (robust) ----
if "protein_id" in df.columns:
    id_col = "protein_id"
elif "target_name" in df.columns:
    id_col = "target_name"
else:
    raise ValueError(f"No valid ID column found. Columns: {list(df.columns)}")

hit_ids = set(df[id_col])

print(f"[INFO] Using column: {id_col}")
print(f"[INFO] Unique hit IDs: {len(hit_ids)}")

# ---- FASTA filtering ----
records = []
for record in SeqIO.parse(input_fasta, "fasta"):
    record_id = record.id.split()[0]  # safer
    if record_id in hit_ids:
        records.append(record)

SeqIO.write(records, output_fasta, "fasta")

print(f"[INFO] Extracted {len(records)} sequences")