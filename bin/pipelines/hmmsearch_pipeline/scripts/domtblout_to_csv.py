import sys
import pandas as pd

input_file = sys.argv[1]
output_file = sys.argv[2]

rows = []

with open(input_file) as f:
    for line in f:
        if line.startswith("#"):
            continue

        parts = line.strip().split(maxsplit=22)

        # --- Core fields ---
        target_name = parts[0]
        target_acc = parts[1]
        tlen = int(parts[2])

        query_name = parts[3]
        query_acc = parts[4]
        qlen = int(parts[5])

        full_evalue = float(parts[6])
        full_score = float(parts[7])
        full_bias = float(parts[8])

        dom_idx = int(parts[9])
        dom_count = int(parts[10])

        c_evalue = float(parts[11])
        i_evalue = float(parts[12])
        dom_score = float(parts[13])
        dom_bias = float(parts[14])

        hmm_from = int(parts[15])
        hmm_to = int(parts[16])

        ali_from = int(parts[17])
        ali_to = int(parts[18])

        env_from = int(parts[19])
        env_to = int(parts[20])

        acc = float(parts[21])
        description = parts[22] if len(parts) > 22 else ""

        # --- Derived metric ---
        coverage = (ali_to - ali_from + 1) / qlen if qlen > 0 else 0

        rows.append({
            "protein_id": target_name,
            "gene": target_name,
            "ascog_id": query_name,

            "coverage": coverage,
            "dom_score": dom_score,
            "c_evalue": c_evalue,
            "i_evalue": i_evalue,

            "tacc": target_acc,
            "tlen": tlen,
            "qacc": query_acc,
            "qlen": qlen,

            "full_evalue": full_evalue,
            "full_score": full_score,
            "full_bias": full_bias,

            "dom_idx": dom_idx,
            "dom_count": dom_count,
            "dom_bias": dom_bias,

            "hmm_from": hmm_from,
            "hmm_to": hmm_to,
            "ali_from": ali_from,
            "ali_to": ali_to,
            "env_from": env_from,
            "env_to": env_to,

            "acc": acc,
            "hmm_description": description,
            "source_file": input_file,

            "arcog_id": "",
            "category": "",
            "description": description
        })

df = pd.DataFrame(rows)

# Optional: sort by best hits
df = df.sort_values(by=["protein_id", "i_evalue"])

df.to_csv(output_file, index=False)

print(f"Saved {len(df)} domain hits → {output_file}")