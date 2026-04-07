import re
import time
import pandas as pd
import requests
from Bio import SeqIO, Entrez
from tqdm import tqdm
import json

# =========================
# CONFIG
# =========================
Entrez.email = "anirudhbantwalbaliga@gmail.com"
Entrez.api_key = "eed94abceca160fc601883bc7e508f604608"

FASTA_FILE = "/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz_fin_v1/with_csm_seq/ftsz.95.rev.fasta"
ASGARD_CSV = "/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz_fin_v1/extraction_exploration/ftsz.rev.csv"
OUTPUT_CSV = "/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz_fin_v1/extraction_exploration/ftsz_csm.rev.csv"

NCBI_BATCH_SIZE = 100
TAX_BATCH_SIZE = 200

# =========================
# HELPERS
# =========================
def normalize_accession(acc):
    return acc.strip().replace("-", "_")

def is_valid_ncbi(acc):
    return bool(re.match(r"^[A-Z]{3}\d+\.\d+$|^(WP_|XP_|NP_|YP_|AP_|ZP_)", acc))

def is_mag_record(rec):
    keywords = rec.get("GBSeq_keywords", [])
    return any(k in ["MAG", "ENV"] for k in keywords)

# =========================
# PARSE HEADERS
# =========================
def classify_and_parse(record):
    rid = record.id
    desc = record.description

    if rid.startswith("sp|") or rid.startswith("tr|"):
        return {"seq_id": rid, "source": "uniprot", "accession": rid.split("|")[1]}

    if re.match(r'^[A-Z]{8}_\d+$', rid):
        return {"seq_id": rid, "source": "asgard", "accession": rid}

    accession = rid.split("|")[0]
    return {"seq_id": rid, "source": "ncbi", "accession": accession}

# =========================
# UNIPROT FETCH
# =========================
def fetch_uniprot_single(accession):
    empty = {"product": None, "taxid": None}

    try:
        r = requests.get(f"https://rest.uniprot.org/uniprotkb/{accession}", timeout=20)
        if r.status_code == 404:
            return empty

        data = json.loads(r.text)

        prot = data.get("proteinDescription", {})
        rec = prot.get("recommendedName", {})

        if rec:
            product = rec.get("fullName", {}).get("value")
        else:
            sub = prot.get("submissionNames", [])
            product = sub[0].get("fullName", {}).get("value") if sub else None

        taxid = data.get("organism", {}).get("taxonId")

        return {"product": product, "taxid": str(taxid) if taxid else None}

    except Exception as e:
        print(f"⚠️ UniProt error ({accession}): {e}")
        return empty

# =========================
# NCBI FETCH (MAG-AWARE)
# =========================
def fetch_ncbi_batch(accessions):
    results = {}
    empty = {"product": None, "taxid": None, "is_mag": False}

    accessions = [normalize_accession(a) for a in accessions]
    accessions = [a for a in accessions if is_valid_ncbi(a)]

    if not accessions:
        return results

    try:
        handle = Entrez.efetch(
            db="protein",
            id=",".join(accessions),
            retmode="xml"
        )
        records = Entrez.read(handle)
        handle.close()

    except Exception as e:
        print(f"⚠️ NCBI fetch failed: {e}")
        return {a: empty for a in accessions}

    for rec in records:
        acc = rec.get("GBSeq_accession-version", "") or rec.get("GBSeq_primary-accession", "")
        bare = acc.split(".")[0]

        product = None
        taxid = None
        is_mag = is_mag_record(rec)

        for feat in rec.get("GBSeq_feature-table", []):

            key = feat.get("GBFeature_key")

            if key == "Protein":
                for qual in feat.get("GBFeature_quals", []):
                    if qual["GBQualifier_name"] == "product":
                        product = qual["GBQualifier_value"]

            if key == "source":
                for qual in feat.get("GBFeature_quals", []):
                    if qual["GBQualifier_name"] == "db_xref":
                        val = qual["GBQualifier_value"]
                        if val.startswith("taxon:"):
                            taxid = val.replace("taxon:", "")

        data = {
            "product": product,
            "taxid": taxid,
            "is_mag": is_mag
        }

        results[acc] = data
        results[bare] = data

    time.sleep(0.12)
    return results

# =========================
# TAXONOMY FETCH
# =========================
def fetch_taxonomy_batch(taxids):
    results = {}

    try:
        handle = Entrez.efetch(db="taxonomy", id=",".join(taxids), retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        for rec in records:
            tid = str(rec["TaxId"])

            tax = {
                "domain": None,
                "phylum": None,
                "class": None,
                "order": None,
                "family": None,
                "genus": None,
                "species": rec.get("ScientificName")
            }

            for item in rec.get("LineageEx", []):
                rank = item["Rank"]
                name = item["ScientificName"]

                if rank == "superkingdom":
                    tax["domain"] = name
                elif rank in tax:
                    tax[rank] = name

            results[tid] = tax

    except Exception as e:
        print(f"⚠️ Taxonomy fetch failed: {e}")

    return results

# =========================
# MAIN
# =========================
records = list(SeqIO.parse(FASTA_FILE, "fasta"))
parsed = [classify_and_parse(r) for r in records]

uniprot_rows = [p for p in parsed if p["source"] == "uniprot"]
ncbi_rows = [p for p in parsed if p["source"] == "ncbi"]

print(f"UniProt: {len(uniprot_rows)} | NCBI: {len(ncbi_rows)}")

# UniProt
uniprot_cache = {}
for r in tqdm(uniprot_rows, desc="UniProt"):
    uniprot_cache[r["accession"]] = fetch_uniprot_single(r["accession"])
    time.sleep(0.2)

# NCBI
ncbi_cache = {}
if ncbi_rows:
    accs = [normalize_accession(r["accession"]) for r in ncbi_rows]

    batches = [accs[i:i+NCBI_BATCH_SIZE] for i in range(0, len(accs), NCBI_BATCH_SIZE)]

    for batch in tqdm(batches, desc="NCBI"):
        ncbi_cache.update(fetch_ncbi_batch(batch))

# Collect taxids
all_taxids = set()

for v in uniprot_cache.values():
    if v.get("taxid"):
        all_taxids.add(v["taxid"])

for v in ncbi_cache.values():
    if v.get("taxid"):
        all_taxids.add(v["taxid"])

print(f"Unique taxids: {len(all_taxids)}")

# Taxonomy
tax_cache = {}
taxids = list(all_taxids)

batches = [taxids[i:i+TAX_BATCH_SIZE] for i in range(0, len(taxids), TAX_BATCH_SIZE)]

for batch in tqdm(batches, desc="Taxonomy"):
    tax_cache.update(fetch_taxonomy_batch(batch))
    time.sleep(0.1)

# Final table
rows = []

for p, rec in zip(parsed, records):
    row = {
        "seq_id": p["seq_id"],
        "accession": p["accession"],
        "source": p["source"],
        "length": len(rec.seq),
        "is_mag": False
    }

    if p["source"] == "uniprot":
        meta = uniprot_cache.get(p["accession"], {})
        tax = tax_cache.get(meta.get("taxid"), {})
        row.update(meta)
        row.update(tax)

    elif p["source"] == "ncbi":
        meta = ncbi_cache.get(normalize_accession(p["accession"]), {})
        tax = tax_cache.get(meta.get("taxid"), {})
        row.update(meta)
        row.update(tax)

    rows.append(row)

df = pd.DataFrame(rows)
df.to_csv(OUTPUT_CSV, index=False)

print(f"\n✅ Saved → {OUTPUT_CSV}")