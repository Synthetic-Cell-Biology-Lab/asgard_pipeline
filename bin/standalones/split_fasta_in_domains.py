import re
import time
from collections import defaultdict

import pandas as pd
from Bio import SeqIO, Entrez
from tqdm import tqdm

# =========================
# CONFIG
# =========================
Entrez.email = "anirudhbantwalbaliga@gmail.com"  # REQUIRED by NCBI
FASTA_FILE = "/home/anirudh/asgard_pipeline/database/protein_sets/SepF/SepF_1/SepF.rev.fasta"
ASGARD_CSV = "/home/anirudh/asgard_pipeline/database/protein_sets/SepF/SepF_1/extraction_exploration/SepF.rev.csv"
OUTPUT_CSV = "/home/anirudh/asgard_pipeline/database/protein_sets/SepF/SepF_1/extraction_exploration/SepF.full.rev.csv"

# =========================
# FUNCTIONS
# =========================

def extract_ox(header):
    match = re.search(r"OX=(\d+)", header)
    return match.group(1) if match else None


def fetch_taxonomy(taxid):
    """Fetch taxonomy lineage from NCBI"""
    try:
        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        record = Entrez.read(handle)[0]
        lineage = record["LineageEx"]
       
        ranks = {
            "domain": None,  # domain
            "kingdom": None,
            "phylum": None,
            "class": None,
            "order": None
        }

        for item in lineage:
            if item["Rank"] in ranks:
                ranks[item["Rank"]] = item["ScientificName"]

        return {
            "domain": ranks["domain"],
            "kingdom": ranks["kingdom"],
            "phylum": ranks["phylum"],
            "class": ranks["class"],
            "order": ranks["order"]
        }

    except Exception:
        return {
            "domain": None,
            "kingdom": None,
            "phylum": None,
            "class": None,
            "order": None
        }


# =========================
# STEP 1: Parse FASTA
# =========================

records = list(SeqIO.parse(FASTA_FILE, "fasta"))

data = []
tax_cache = {}

print("Extracting taxonomy...")

for record in tqdm(records):
    header = record.description
    ox = extract_ox(header)

    if ox:
        if ox not in tax_cache:
            tax_cache[ox] = fetch_taxonomy(ox)
            time.sleep(0.34)  # NCBI rate limit (~3 req/sec)

        tax_info = tax_cache[ox]
    else:
        tax_info = {
            "domain": None,
            "kingdom": None,
            "phylum": None,
            "class": None,
            "order": None
        }

    data.append({
        "id": record.id,
        "description": header,
        "OX": ox,
        **tax_info
    })

df = pd.DataFrame(data)

# =========================
# STEP 2: Append Asgard CSV
# =========================

if ASGARD_CSV:
    asgard_df = pd.read_csv(ASGARD_CSV)

    # Ensure same columns
    for col in df.columns:
        if col not in asgard_df.columns:
            asgard_df[col] = None

    df = pd.concat([df, asgard_df], ignore_index=True)

# =========================
# STEP 3: Save CSV
# =========================

df.to_csv(OUTPUT_CSV, index=False)
print(f"Saved taxonomy table → {OUTPUT_CSV}")

# =========================
# STEP 4: Split FASTA by domain
# =========================

domain_groups = defaultdict(list)

# Map ID → domain
id_to_domain = dict(zip(df["id"], df["domain"]))

for record in records:
    domain = id_to_domain.get(record.id, "Unknown")
    if domain is None:
        domain = "Unknown"

    domain_groups[domain].append(record)

# Write files
for domain, recs in domain_groups.items():
    safe_name = domain.replace(" ", "_")
    out_file = f"/home/anirudh/asgard_pipeline/database/protein_sets/SepF/SepF_1/sepf_{safe_name}.rev.fasta"
    SeqIO.write(recs, out_file, "fasta")
    print(f"Wrote {len(recs)} sequences → {out_file}")