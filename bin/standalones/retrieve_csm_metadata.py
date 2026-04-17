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
Entrez.email   = "anirudhbantwalbaliga@gmail.com"
Entrez.api_key = "eed94abceca160fc601883bc7e508f604608"

FASTA_FILE  = "/home/anirudh/asgard_pipeline/database/protein_sets/SepF/SepF_1/SepF.rev.fasta"
ASGARD_CSV  = "/home/anirudh/asgard_pipeline/database/protein_sets/SepF/SepF_1/extraction_exploration/SepF.rev.csv"
OUTPUT_CSV  = "/home/anirudh/asgard_pipeline/database/protein_sets/SepF/SepF_1/extraction_exploration/SepF_meta.rev.csv"

NCBI_BATCH_SIZE = 10
TAX_BATCH_SIZE  = 200

TAX_COLS = ["domain", "phylum", "class", "order", "family", "genus", "species"]

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

def tax_is_complete(row_dict, tax_cols):
    """Return True only if every taxonomy column is present and non-empty."""
    for col in tax_cols:
        val = row_dict.get(col)
        if val is None or (isinstance(val, float) and pd.isna(val)):
            return False
        if str(val).strip() in ("", "nan", "None", "unknown", "Unknown"):
            return False
    return True

# =========================
# PARSE HEADERS
# =========================
def classify_and_parse(record):
    rid = record.id

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
        rec  = prot.get("recommendedName", {})
        if rec:
            product = rec.get("fullName", {}).get("value")
        else:
            sub     = prot.get("submissionNames", [])
            product = sub[0].get("fullName", {}).get("value") if sub else None
        taxid = data.get("organism", {}).get("taxonId")
        return {"product": product, "taxid": str(taxid) if taxid else None}
    except Exception as e:
        print(f"⚠️ UniProt error ({accession}): {e}")
        return empty

# =========================
# NCBI FETCH (with retry + exponential backoff)
# =========================
def fetch_ncbi_batch(accessions, max_retries=4):
    results = {}
    empty   = {"product": None, "taxid": None, "is_mag": False}

    accessions = [normalize_accession(a) for a in accessions]
    accessions = [a for a in accessions if is_valid_ncbi(a)]
    if not accessions:
        return results

    wait = 1.0  # initial backoff in seconds

    for attempt in range(1, max_retries + 1):
        try:
            handle  = Entrez.efetch(db="protein", id=",".join(accessions), retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            break   # success — exit retry loop

        except Exception as e:
            print(f"⚠️ NCBI fetch failed (attempt {attempt}/{max_retries}): {e}")

            if attempt == max_retries:
                # Final attempt failed — if batch is splittable, recurse on halves
                if len(accessions) > 1:
                    print(f"  ↳ Splitting batch of {len(accessions)} and retrying halves …")
                    mid = len(accessions) // 2
                    results.update(fetch_ncbi_batch(accessions[:mid], max_retries))
                    results.update(fetch_ncbi_batch(accessions[mid:], max_retries))
                else:
                    print(f"  ↳ Single accession {accessions[0]} failed — marking empty.")
                    results[accessions[0]] = empty
                return results

            time.sleep(wait)
            wait *= 2   # exponential backoff: 1 → 2 → 4 → 8 s

    for rec in records:
        acc  = rec.get("GBSeq_accession-version", "") or rec.get("GBSeq_primary-accession", "")
        bare = acc.split(".")[0]

        product = None
        taxid   = None
        is_mag  = is_mag_record(rec)

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
        results[acc]  = data
        results[bare] = data

    time.sleep(0.3)
    return results

# =========================
# TAXONOMY FETCH
# =========================
def fetch_taxonomy_batch(taxids):
    results = {}
    try:
        handle  = Entrez.efetch(db="taxonomy", id=",".join(taxids), retmode="xml")
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
# LOAD ASGARD_CSV LOOKUP
# =========================
print("Loading ASGARD_CSV …")
rev_df = pd.read_csv(ASGARD_CSV, dtype=str, low_memory=False)

CARRY_COLS = ["product"] + TAX_COLS
CARRY_COLS = [c for c in CARRY_COLS if c in rev_df.columns]

# Index on seq_id (the column name used in the previous output)
id_col = "seq_id" if "seq_id" in rev_df.columns else "locus_tag"

rev_lookup = (
    rev_df[[id_col] + CARRY_COLS]
    .drop_duplicates(subset=id_col)
    .set_index(id_col)
    .to_dict(orient="index")
)
print(f"  ASGARD_CSV rows loaded: {len(rev_lookup)}")

# =========================
# PARSE FASTA + SPLIT
# =========================
records = list(SeqIO.parse(FASTA_FILE, "fasta"))
parsed  = [classify_and_parse(r) for r in records]

cached_rows = []   # seq_id found AND all taxonomy complete
partial_rows = []  # seq_id found BUT taxonomy incomplete → re-fetch taxonomy only
fetch_rows  = []   # seq_id not found → full remote fetch

threshold_tax_cols = TAX_COLS[:3]

for p, rec in zip(parsed, records):
    cached = rev_lookup.get(p["seq_id"])

    if cached is None:
        # Not in ASGARD_CSV at all
        fetch_rows.append((p, rec, None))

    elif tax_is_complete(cached, [c for c in threshold_tax_cols if c in CARRY_COLS]):
        # Found and taxonomy complete — use as-is
        cached_rows.append((p, rec, cached))

    else:
        # Found but taxonomy incomplete — keep product, re-fetch taxonomy
        partial_rows.append((p, rec, cached))

print(f"\nFully cached (seq_id + complete taxonomy) : {len(cached_rows)}")
print(f"Partial cache (seq_id found, tax missing)  : {len(partial_rows)}")
print(f"No cache (full remote fetch needed)        : {len(fetch_rows)}")

# =========================
# REMOTE FETCH
# =========================
# For partial rows we only need taxonomy — but we need a taxid first,
# so we still hit NCBI/UniProt for the taxid then fetch taxonomy.
needs_fetch = partial_rows + fetch_rows   # combined list needing some network call

fetch_parsed  = [p   for p, _, _ in needs_fetch]
fetch_records = [rec for _, rec, _ in needs_fetch]
fetch_cached  = [c   for _, _, c  in needs_fetch]   # may be None or partial dict

uniprot_rows = [p for p in fetch_parsed if p["source"] == "uniprot"]
ncbi_rows    = [p for p in fetch_parsed if p["source"] == "ncbi"]

print(f"\nRemote — UniProt: {len(uniprot_rows)} | NCBI: {len(ncbi_rows)}")

# NCBI
ncbi_cache = {}
if ncbi_rows:
    accs    = [normalize_accession(r["accession"]) for r in ncbi_rows]
    batches = [accs[i:i+NCBI_BATCH_SIZE] for i in range(0, len(accs), NCBI_BATCH_SIZE)]
    for batch in tqdm(batches, desc="NCBI"):
        ncbi_cache.update(fetch_ncbi_batch(batch))

# UniProt
uniprot_cache = {}
for r in tqdm(uniprot_rows, desc="UniProt"):
    uniprot_cache[r["accession"]] = fetch_uniprot_single(r["accession"])
    time.sleep(0.2)

# Collect taxids
all_taxids = set()
for v in uniprot_cache.values():
    if v.get("taxid"):
        all_taxids.add(v["taxid"])
for v in ncbi_cache.values():
    if v.get("taxid"):
        all_taxids.add(v["taxid"])

print(f"\nUnique taxids to resolve: {len(all_taxids)}")

# Taxonomy
tax_cache = {}
taxids    = list(all_taxids)
batches   = [taxids[i:i+TAX_BATCH_SIZE] for i in range(0, len(taxids), TAX_BATCH_SIZE)]
for batch in tqdm(batches, desc="Taxonomy"):
    tax_cache.update(fetch_taxonomy_batch(batch))
    time.sleep(0.1)

# =========================
# ASSEMBLE FINAL TABLE
# =========================
output_rows = []

# ── Fully cached rows ──
for p, rec, cached in cached_rows:
    row = {
        "seq_id":    p["seq_id"],
        "accession": p["accession"],
        "source":    p["source"],
        "length":    len(rec.seq),
        "is_mag":    False,
        "product":   cached.get("product"),
    }
    for tc in TAX_COLS:
        row[tc] = cached.get(tc)
    output_rows.append(row)

# ── Partially cached + fully remote rows ──
for p, rec, prior_cache in zip(fetch_parsed, fetch_records, fetch_cached):
    row = {
        "seq_id":    p["seq_id"],
        "accession": p["accession"],
        "source":    p["source"],
        "length":    len(rec.seq),
        "is_mag":    False,
    }

    if p["source"] == "uniprot":
        meta = uniprot_cache.get(p["accession"], {})
    elif p["source"] == "ncbi":
        meta = ncbi_cache.get(normalize_accession(p["accession"]), {})
    else:
        meta = {}

    tax = tax_cache.get(meta.get("taxid"), {})

    # Start from prior cache if available (preserves product if already known)
    if prior_cache:
        row.update({k: v for k, v in prior_cache.items() if v is not None})

    row.update(meta)   # remote meta wins for product / taxid / is_mag
    row.update(tax)    # freshly fetched taxonomy always wins

    output_rows.append(row)

df = pd.DataFrame(output_rows)
df.to_csv(OUTPUT_CSV, index=False)

print(f"\n✅ Saved → {OUTPUT_CSV}")
print(f"   Total rows      : {len(df)}")
print(f"   Fully cached    : {len(cached_rows)}")
print(f"   Partial re-fetch: {len(partial_rows)}")
print(f"   Full remote     : {len(fetch_rows)}")