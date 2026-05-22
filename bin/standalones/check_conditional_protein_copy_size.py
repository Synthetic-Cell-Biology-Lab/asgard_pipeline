import pandas as pd




csv = "/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz_fin_v1/extraction_exploration/ftsz.rev.csv"
csv_df = pd.read_csv(csv)

genomes_with_tubulin = csv_df[csv_df['Manual_annotation'] == "Tubulin"]['genome_file'].to_list()

print(genomes_with_tubulin)


genome_subset = csv_df[csv_df['genome_file'].isin(genomes_with_tubulin)]

print(genome_subset)

heim_genomes = genome_subset[genome_subset['class'] == "Heimdallarchaeia"]
print(heim_genomes)