from Bio import SeqIO

input_fastas = [snakemake.input.rev_fasta] + list(snakemake.input.additional_fastas)
output_fasta = snakemake.output.full_fasta

records = []
seen = set()

for fasta in input_fastas:
    for record in SeqIO.parse(fasta, "fasta"):
        clean_id = (
            record.id.split()[0]
            .replace("|", "_")
            .replace("/", "_")
            .replace(":", "_")
            .replace(" ", "_")
        )
        if clean_id in seen:
            continue
        seen.add(clean_id)
        record.id = clean_id
        record.name = clean_id
        record.description = ""
        records.append(record)

SeqIO.write(records, output_fasta, "fasta")
print(f"Merged and normalized {len(records)} sequences from {len(input_fastas)} files")
