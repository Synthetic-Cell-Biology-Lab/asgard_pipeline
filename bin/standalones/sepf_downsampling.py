
import csv
import subprocess
import shutil
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO

# ---------------------------
# INPUT FILES
# ---------------------------
csv_file = "/home/anirudh/asgard_pipeline/database/protein_sets/SepF/SepF_1/extraction_exploration/SepF.rev.csv"
fasta_file = "/home/anirudh/asgard_pipeline/database/protein_sets/SepF/SepF_1/SepF.rev.fasta"

output_fasta = "/home/anirudh/asgard_pipeline/database/protein_sets/SepF/SepF_1/SepF.rep.rev.fasta"
log_file = "/home/anirudh/asgard_pipeline/database/protein_sets/SepF/SepF_1/SepF_downsampling_log.txt"
output_csv = "/home/anirudh/asgard_pipeline/database/protein_sets/SepF/SepF_1/extraction_exploration/SepF.rep.rev.csv"

work_dir = Path("/home/anirudh/asgard_pipeline/database/protein_sets/SepF/SepF_1/mmseqs_work")

# ---------------------------
# LOAD CSV
# ---------------------------
with open(csv_file) as f:
    reader = csv.DictReader(f)
    records = list(reader)

# ---------------------------
# GROUP BY GENOME
# ---------------------------
genome_groups = defaultdict(list)
for r in records:
    genome_groups[r["genome_file"]].append(r)

selected = set()
remaining = []
log_lines = []

# ---------------------------
# RULE 1: KEEP GENOMES WITH 2 ENTRIES
# ---------------------------
for genome, group in genome_groups.items():
    if len(group) == 2:
        for r in group:
            selected.add(r["locus_tag"])
        log_lines.append(f"[PAIR] {genome}: kept both ({len(group)})")
    else:
        remaining.extend(group)

# ---------------------------
# LOAD FASTA
# ---------------------------
seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

# ---------------------------
# GROUP BY GENUS
# ---------------------------
genus_groups = defaultdict(list)
for r in remaining:
    genus_groups[r["genus"]].append(r)

# ---------------------------
# FUNCTION: CLEAN MMSEQS RUN
# ---------------------------
def pick_medoid(genus, group):
    genus_dir = work_dir / genus.replace(" ", "_")
    genus_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = genus_dir / "input.fasta"
    db = genus_dir / "db"
    cluster = genus_dir / "cluster"
    tmp = genus_dir / "tmp"
    rep = genus_dir / "rep"
    rep_fasta = genus_dir / "rep.fasta"

    # write genus fasta
    with open(fasta_path, "w") as out:
        for r in group:
            tag = r["locus_tag"]
            if tag in seq_dict:
                SeqIO.write(seq_dict[tag], out, "fasta")

    try:
        # create db
        subprocess.run(["mmseqs", "createdb", fasta_path, db], check=True)

        # cluster
        subprocess.run([
            "mmseqs", "cluster",
            db, cluster, tmp,
            "--min-seq-id", "0.3",
            "-c", "0.8"
        ], check=True)

        # representative sequences
        subprocess.run([
            "mmseqs", "createseqfiledb",
            db, cluster, rep
        ], check=True)

        subprocess.run([
            "mmseqs", "convert2fasta",
            rep, rep_fasta
        ], check=True)

        reps = list(SeqIO.parse(rep_fasta, "fasta"))
        if reps:
            return reps[0].id

    except subprocess.CalledProcessError:
        return None

    finally:
        # 🔥 CLEANUP genus-specific files
        shutil.rmtree(genus_dir, ignore_errors=True)

    return None

# ---------------------------
# RULE 2: GENUS REPRESENTATIVE
# ---------------------------
work_dir.mkdir(exist_ok=True)

for genus, group in genus_groups.items():
    if len(group) == 1:
        tag = group[0]["locus_tag"]
        selected.add(tag)
        log_lines.append(f"[GENUS-single] {genus}: {tag}")
    else:
        medoid = pick_medoid(genus, group)
        if medoid:
            selected.add(medoid)
            log_lines.append(f"[GENUS-medoid] {genus}: {medoid}")
        else:
            fallback = group[0]["locus_tag"]
            selected.add(fallback)
            log_lines.append(f"[GENUS-fallback] {genus}: {fallback}")

# ---------------------------
# RULE 3: ENSURE ONE PER ORDER
# ---------------------------
order_map = defaultdict(list)
for r in records:
    order_map[r["order"]].append(r)

selected_orders = set(
    r["order"] for r in records if r["locus_tag"] in selected
)

for order, group in order_map.items():
    if order not in selected_orders:
        tag = group[0]["locus_tag"]
        selected.add(tag)
        log_lines.append(f"[ORDER-fill] {order}: {tag}")

# ---------------------------
# WRITE FASTA
# ---------------------------
with open(output_fasta, "w") as out:
    for tag in selected:
        if tag in seq_dict:
            SeqIO.write(seq_dict[tag], out, "fasta")

# ---------------------------
# WRITE FILTERED CSV
# ---------------------------
selected_records = [r for r in records if r["locus_tag"] in selected]

with open(output_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=records[0].keys())
    writer.writeheader()
    writer.writerows(selected_records)

# ---------------------------
# WRITE LOG
# ---------------------------
with open(log_file, "w") as f:
    f.write("\n".join(log_lines))

# ---------------------------
# FINAL CLEANUP (just in case)
# ---------------------------
shutil.rmtree(work_dir, ignore_errors=True)

print(f"Done.")
print(f"Sequences selected: {len(selected)}")
print(f"FASTA: {output_fasta}")
print(f"CSV: {output_csv}")
print(f"LOG: {log_file}")