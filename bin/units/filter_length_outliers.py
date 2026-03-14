from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

faa_in = snakemake.input.faa
ffn_in = snakemake.input.ffn

faa_out = snakemake.output.faa
ffn_out = snakemake.output.ffn
plotfile = snakemake.output.plot

proteins = list(SeqIO.parse(faa_in, "fasta"))
lengths = np.array([len(r.seq) for r in proteins])

# statistics
median = np.median(lengths)
q1 = np.percentile(lengths, 25)
q3 = np.percentile(lengths, 75)
iqr = q3 - q1

lower_iqr = q1 - 1.5 * iqr
upper_iqr = q3 + 1.5 * iqr

lower_median = median * 0.5

lower = max(lower_iqr, lower_median)
upper = upper_iqr

# filter proteins
filtered_proteins = [r for r in proteins if lower <= len(r.seq) <= upper]

SeqIO.write(filtered_proteins, faa_out, "fasta")

valid_ids = {r.id for r in filtered_proteins}

# filter nucleotides using same IDs
ffn_records = (
    r for r in SeqIO.parse(ffn_in, "fasta")
    if r.id in valid_ids
)

SeqIO.write(ffn_records, ffn_out, "fasta")

# plot
plt.hist(lengths, bins=30)
plt.axvline(lower)
plt.axvline(upper)
plt.axvline(median)
plt.savefig(plotfile)