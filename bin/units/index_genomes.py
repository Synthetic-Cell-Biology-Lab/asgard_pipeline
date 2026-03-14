from pathlib import Path
from Bio import SeqIO

cds_dir = Path(snakemake.input.cds_dir)

for fasta in cds_dir.rglob("*"):

    if fasta.suffix not in {".faa", ".ffn"}:
        continue

    idx = fasta.with_suffix(fasta.suffix + ".idx")

    if idx.exists():
        print(f"Index exists: {idx}")
        continue

    print(f"Building index for {fasta}")

    SeqIO.index_db(str(idx), str(fasta), "fasta")