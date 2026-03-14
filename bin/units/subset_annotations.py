import pandas as pd

protein = snakemake.wildcards.protein
annotation_col = snakemake.params.annotation_col

df = pd.read_csv(
    snakemake.input[0],
    usecols=["locus_tag","genome_file",annotation_col]
)

subset = df[df[annotation_col] == protein]

subset.to_csv(snakemake.output[0], index=False)