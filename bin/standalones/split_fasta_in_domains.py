import re
import time
from collections import defaultdict

import pandas as pd
from Bio import SeqIO, Entrez
from tqdm import tqdm

# =========================
# CONFIG
# =========================
Entrez.email = "anirudhbantwalbaliga@gmail.com"
Entrez.api_key = "eed94abceca160fc601883bc7e508f604608"  # get free at https://www.ncbi.nlm.nih.gov/account/
RATE_LIMIT = 0.11  # ~10 req/sec with API key, use 0.34 without

FASTA_FILE = "/home/anirudh/asgard_pipeline/database/protein_sets/SepF/SepF_1/SepF.rev.fasta"
ASGARD_CSV = "/home/anirudh/asgard_pipeline/database/protein_sets/SepF/SepF_1/extraction_exploration/SepF.rev.csv"
OUTPUT_CSV = "/home/anirudh/asgard_pipeline/database/protein_sets/SepF/SepF_1/extraction_exploration/SepF.full.rev.csv"
OUT_DIR = "/home/anirudh/asgard_pipeline/database/protein_sets/SepF/SepF_1"

BATCH_SIZE = 500  # NCBI max per request

# =========================
# FUNCTIONS
# =========================

def extract_ox(header):
    match = re.search(r"OX=(\d+)", header)
    return match.group(1) if match else None


def fetch_taxonomy_batch(taxids):
    """
    Fetch taxonomy for a batch of taxids in one NCBI request.
    Returns dict: taxid -> taxonomy dict
    """
    results = {}
    try:
        handle = Entrez.efetch(
            db="taxonomy",
            id=",".join(taxids),
            retmode="xml"
        )
        records = Entrez.read(handle)

        for record in records:
            taxid = record["TaxId"]
            lineage = record.get("LineageEx", [])
            print(lineage)
            # NCBI uses "superkingdom" not "domain"
            rank_map = {
                "domain": "domain",
                "kingdom": "kingdom",
                "phylum": "phylum",
                "class": "class",
                "order": "order",
            }
            ranks = {v: None for v in rank_map.values()}

            for item in lineage:
                ncbi_rank = item["Rank"]
                if ncbi_rank in rank_map:
                    ranks[rank_map[ncbi_rank]] = item["ScientificName"]

            results[str(taxid)] = ranks

    except Exception as e:
        print(f"  ⚠️ Batch fetch failed: {e}")
        # Return empty results for all taxids in this batch
        empty = {"domain": None, "kingdom": None, "phylum": None,
                 "class": None, "order": None}
        for tid in taxids:
            results[tid] = empty

    return results


# =========================
# STEP 1: Parse FASTA
# =========================

records = list(SeqIO.parse(FASTA_FILE, "fasta"))
print(f"Loaded {len(records)} sequences")

# Extract OX for each record
record_data = []
for record in records:
    ox = extract_ox(record.description)
    record_data.append({
        "id": record.id,
        "description": record.description,
        "OX": ox,
    })

# =========================
# STEP 2: Batch fetch taxonomy
# =========================

# Collect unique taxids that need fetching
all_oxids = list({r["OX"] for r in record_data if r["OX"] is not None})
print(f"Unique taxids to fetch: {len(all_oxids)}")

tax_cache = {}
batches = [all_oxids[i:i + BATCH_SIZE] for i in range(0, len(all_oxids), BATCH_SIZE)]

for batch in tqdm(batches, desc="Fetching taxonomy"):
    batch_results = fetch_taxonomy_batch(batch)
    tax_cache.update(batch_results)
    time.sleep(RATE_LIMIT)

# =========================
# STEP 3: Build dataframe
# =========================

empty_tax = {"domain": None, "kingdom": None, "phylum": None,
             "class": None, "order": None}

data = []
for row in record_data:
    ox = row["OX"]
    tax_info = tax_cache.get(ox, empty_tax) if ox else empty_tax
    data.append({**row, **tax_info})

df = pd.DataFrame(data)

# =========================
# STEP 4: Append Asgard CSV
# =========================
# =========================
# STEP 4: Replace with Asgard CSV (based on locus_tag)
# =========================

if ASGARD_CSV:
    asgard_df = pd.read_csv(ASGARD_CSV)
    print(f"Asgard CSV columns: {list(asgard_df.columns)}")

    # Ensure both have same columns
    for col in df.columns:
        if col not in asgard_df.columns:
            asgard_df[col] = None

    for col in asgard_df.columns:
        if col not in df.columns:
            df[col] = None

    # Set index to locus_tag for matching
    if "locus_tag" not in df.columns or "locus_tag" not in asgard_df.columns:
        raise ValueError("locus_tag column missing in one of the dataframes")

    df.set_index("locus_tag", inplace=True)
    asgard_df.set_index("locus_tag", inplace=True)

    # Replace rows in df where locus_tag exists in asgard_df
    overlap = df.index.intersection(asgard_df.index)
    print(f"Replacing {len(overlap)} rows from Asgard CSV")

    df.update(asgard_df)

    # (Optional) If you also want to ADD new rows from Asgard not in df:
    missing = asgard_df.index.difference(df.index)
    print(f"Adding {len(missing)} new rows from Asgard CSV")

    df = pd.concat([df, asgard_df.loc[missing]])

    # Reset index back
    df.reset_index(inplace=True)

# =========================
# STEP 5: Save CSV
# =========================

df.to_csv(OUTPUT_CSV, index=False)
print(f"Saved taxonomy table → {OUTPUT_CSV}")
print(df["domain"].value_counts(dropna=False))

# =========================
# STEP 6: Split FASTA by domain
# =========================

id_to_domain = dict(zip(df["id"], df["domain"]))

domain_groups = defaultdict(list)
for record in records:
    domain = id_to_domain.get(record.id) or "Unknown"
    domain_groups[domain].append(record)

for domain, recs in domain_groups.items():
    safe_name = re.sub(r'[^\w]', '_', domain)
    out_file = f"{OUT_DIR}/sepf_{safe_name}.rev.fasta"
    SeqIO.write(recs, out_file, "fasta")
    print(f"Wrote {len(recs)} sequences → {out_file}")