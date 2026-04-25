#!/usr/bin/env python3

from Bio import SeqIO
from datetime import datetime
# ==============================
# 🔧 HARDCODED INPUTS
# ==============================
FASTA_PATH = "/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz_fin_v1/with_csm_seq/ftsz.95.rev.fasta"
IDS_PATH = "/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz_fin_v1/with_csm_seq/v3/ftsz.csm.ids"
OUTPUT_PATH = "/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz_fin_v1/with_csm_seq/v3/ftsz.csm.fasta"
LOG_PATH = "/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz_fin_v1/phylogeny_pipeline.log"

MIN_LEN = 0
MAX_LEN = 10000000


# ==============================
# 🎨 COLOR DEFINITIONS (terminal only)
# ==============================
class Color:
    GREEN = "\033[92m"
    RED = "\033[91m"
    YELLOW = "\033[93m"
    CYAN = "\033[96m"
    RESET = "\033[0m"

def log(message, color=Color.RESET):
    """Print colored message + append plain text to log file"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Terminal output (colored)
    print(f"{color}[{timestamp}] {message}{Color.RESET}")
    
    # File output (no color)
    with open(LOG_PATH, "a") as log_f:
        log_f.write(f"[{timestamp}] {message}\n")

# ==============================
# 🚀 START
# ==============================
log("Starting protein filtering run", Color.CYAN)

# ==============================
# 📂 LOAD IDS
# ==============================
with open(IDS_PATH) as f:
    ids = set(line.strip().replace(" ", "_") for line in f if line.strip())

log(f"Loaded {len(ids)} IDs", Color.CYAN)

# ==============================
# 🔍 FILTER FASTA
# ==============================
kept = 0
total = 0

with open(OUTPUT_PATH, "w") as out_f:
    for record in SeqIO.parse(FASTA_PATH, "fasta"):
        total += 1
        
        record_id = record.id
        seq_len = len(record.seq)
        
        if record_id in ids and MIN_LEN <= seq_len <= MAX_LEN:
            SeqIO.write(record, out_f, "fasta")
            kept += 1
            log(f"KEPT: {record_id} | length={seq_len}", Color.GREEN)
        else:
            log(f"SKIPPED: {record_id} | length={seq_len}", Color.RED)

# ==============================
# 📊 SUMMARY
# ==============================
log(f"Total sequences scanned: {total}", Color.YELLOW)
log(f"Sequences kept: {kept}", Color.GREEN)
log(f"Output written to: {OUTPUT_PATH}", Color.CYAN)
log("Run complete\n", Color.CYAN)