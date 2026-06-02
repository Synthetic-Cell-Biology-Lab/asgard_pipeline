import re
import time
import json
import pandas as pd
import requests
from Bio import SeqIO, Entrez
from tqdm import tqdm

# =========================
# CONFIG (injected by Snakemake)
# =========================
FASTA_FILE = snakemake.input.fasta
PIPELINE_CSV = snakemake.input.pipeline_csv
OUTPUT_CSV = snakemake.output.final_csv

Entrez.email = snakemake.config["entrez"]["email"]
Entrez.api_key = snakemake.config["entrez"]["api_key"]

NCBI_BATCH_SIZE = 10
TAX_BATCH_SIZE = 200
TAX_COLS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


# =========================
# HELPERS
# =========================
def normalize_accession(acc):
    return acc.strip().replace("-", "_")


def is_valid_ncbi(acc):
    return bool(re.match(r"^[A-Z]{3}\d+\.\d+$|^(WP_|XP_|NP_|YP_|AP_|ZP_)", acc))


def is_mag_record(rec):
    return any(k in ["MAG", "ENV"] for k in rec.get("GBSeq_keywords", []))


def tax_is_complete(row_dict, tax_cols):
    for col in tax_cols:
        val = row_dict.get(col)
        if val is None or (isinstance(val, float) and pd.isna(val)):
            return False
        if str(val).strip() in ("", "nan", "None", "unknown", "Unknown"):
            return False
    return True


def classify_and_parse(record):
    rid = record.id
    if rid.startswith("sp|") or rid.startswith("tr|"):
        return {"seq_id": rid, "source": "uniprot", "accession": rid.split("|")[1]}
    if re.match(r"^[A-Z]{8}_\d+$", rid):
        return {"seq_id": rid, "source": "pipeline", "accession": rid}
    accession = rid.split("|")[0]
    return {"seq_id": rid, "source": "ncbi", "accession": accession}


# =========================
# FETCH FUNCTIONS
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


def fetch_ncbi_batch(accessions, max_retries=4):
    results = {}
    empty = {"product": None, "taxid": None, "is_mag": False}
    accessions = [normalize_accession(a) for a in accessions]
    accessions = [a for a in accessions if is_valid_ncbi(a)]
    if not accessions:
        return results

    wait = 1.0
    for attempt in range(1, max_retries + 1):
        try:
            handle = Entrez.efetch(db="protein", id=",".join(accessions), retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            break
        except Exception as e:
            print(f"⚠️ NCBI fetch failed (attempt {attempt}/{max_retries}): {e}")
            if attempt == max_retries:
                if len(accessions) > 1:
                    mid = len(accessions) // 2
                    results.update(fetch_ncbi_batch(accessions[:mid], max_retries))
                    results.update(fetch_ncbi_batch(accessions[mid:], max_retries))
                else:
                    results[accessions[0]] = empty
                return results
            time.sleep(wait)
            wait *= 2

    for rec in records:
        acc = rec.get("GBSeq_accession-version", "") or rec.get(
            "GBSeq_primary-accession", ""
        )
        bare = acc.split(".")[0]
        product, taxid = None, None
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

        data = {"product": product, "taxid": taxid, "is_mag": is_mag}
        results[acc] = data
        results[bare] = data

    time.sleep(0.3)
    return results


def fetch_taxonomy_batch(taxids):
    results = {}
    try:
        handle = Entrez.efetch(db="taxonomy", id=",".join(taxids), retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        for rec in records:
            tid = str(rec["TaxId"])
            tax = {col: None for col in TAX_COLS}
            tax["species"] = rec.get("ScientificName")
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
# LOAD PIPELINE CSV
# =========================
print("Loading pipeline CSV …")
pip_df = pd.read_csv(PIPELINE_CSV, dtype=str, low_memory=False)
pip_df.columns = pip_df.columns.str.strip()

CARRY_COLS = ["product"] + TAX_COLS
CARRY_COLS = [c for c in CARRY_COLS if c in pip_df.columns]

id_col = "locus_tag" if "locus_tag" in pip_df.columns else "seq_id"

pip_lookup = (
    pip_df[[id_col] + CARRY_COLS]
    .drop_duplicates(subset=id_col)
    .set_index(id_col)
    .to_dict(orient="index")
)
print(f"  Pipeline CSV rows loaded: {len(pip_lookup)}")

# =========================
# PARSE FASTA + CLASSIFY
# =========================
fasta_records = list(SeqIO.parse(FASTA_FILE, "fasta"))
parsed = [classify_and_parse(r) for r in fasta_records]

cached_rows = []  # in pipeline CSV + taxonomy complete
partial_rows = []  # in pipeline CSV + taxonomy incomplete
fetch_rows = []  # not in pipeline CSV → full remote fetch

for p, rec in zip(parsed, fasta_records):
    cached = pip_lookup.get(p["seq_id"])

    if p["source"] == "pipeline" or cached is not None:
        if cached and tax_is_complete(cached, TAX_COLS[:3]):
            cached_rows.append((p, rec, cached))
        else:
            partial_rows.append((p, rec, cached))
    else:
        fetch_rows.append((p, rec, None))

print(f"\nFully cached  : {len(cached_rows)}")
print(f"Partial cache : {len(partial_rows)}")
print(f"Remote fetch  : {len(fetch_rows)}")

# =========================
# REMOTE FETCH
# =========================
needs_fetch = partial_rows + fetch_rows
fetch_parsed = [p for p, _, _ in needs_fetch]
fetch_records = [rec for _, rec, _ in needs_fetch]
fetch_cached = [c for _, _, c in needs_fetch]

# Filter to only non-pipeline sources
uniprot_rows = [p for p in fetch_parsed if p["source"] == "uniprot"]
ncbi_rows = [p for p in fetch_parsed if p["source"] == "ncbi"]

print(f"\nRemote — UniProt: {len(uniprot_rows)} | NCBI: {len(ncbi_rows)}")

ncbi_cache = {}
if ncbi_rows:
    accs = [normalize_accession(r["accession"]) for r in ncbi_rows]
    batches = [
        accs[i : i + NCBI_BATCH_SIZE] for i in range(0, len(accs), NCBI_BATCH_SIZE)
    ]
    for batch in tqdm(batches, desc="NCBI"):
        ncbi_cache.update(fetch_ncbi_batch(batch))

uniprot_cache = {}
for r in tqdm(uniprot_rows, desc="UniProt"):
    uniprot_cache[r["accession"]] = fetch_uniprot_single(r["accession"])
    time.sleep(0.2)

all_taxids = set()
for v in uniprot_cache.values():
    if v.get("taxid"):
        all_taxids.add(v["taxid"])
for v in ncbi_cache.values():
    if v.get("taxid"):
        all_taxids.add(v["taxid"])

print(f"\nUnique taxids to resolve: {len(all_taxids)}")

tax_cache = {}
taxids = list(all_taxids)
batches = [
    taxids[i : i + TAX_BATCH_SIZE] for i in range(0, len(taxids), TAX_BATCH_SIZE)
]
for batch in tqdm(batches, desc="Taxonomy"):
    tax_cache.update(fetch_taxonomy_batch(batch))
    time.sleep(0.1)

# =========================
# ASSEMBLE FINAL TABLE
# =========================
output_rows = []

for p, rec, cached in cached_rows:
    row = {
        "locus_tag": p["seq_id"],
        "accession": p["accession"],
        "source": p["source"],
        "length": len(rec.seq),
        "is_mag": False,
        "product": cached.get("product"),
    }
    for tc in TAX_COLS:
        row[tc] = cached.get(tc)
    output_rows.append(row)

for p, rec, prior_cache in zip(fetch_parsed, fetch_records, fetch_cached):
    row = {
        "locus_tag": p["seq_id"],
        "accession": p["accession"],
        "source": p["source"],
        "length": len(rec.seq),
        "is_mag": False,
    }

    if p["source"] == "uniprot":
        meta = uniprot_cache.get(p["accession"], {})
    elif p["source"] == "ncbi":
        meta = ncbi_cache.get(normalize_accession(p["accession"]), {})
    else:
        meta = {}

    tax = tax_cache.get(meta.get("taxid"), {})

    if prior_cache:
        row.update({k: v for k, v in prior_cache.items() if v is not None})

    row.update(meta)
    row.update(tax)

    output_rows.append(row)

df = pd.DataFrame(output_rows)
df.to_csv(OUTPUT_CSV, index=False)

print(f"\n✅ Saved → {OUTPUT_CSV}")
print(f"   Total rows      : {len(df)}")
print(f"   Fully cached    : {len(cached_rows)}")
print(f"   Partial re-fetch: {len(partial_rows)}")
print(f"   Full remote     : {len(fetch_rows)}")
