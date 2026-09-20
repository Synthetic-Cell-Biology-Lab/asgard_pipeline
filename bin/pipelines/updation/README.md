# Database Updation Pipeline

`bin/pipelines/updation/Snakefile` builds either an isolated genome/protein set or a cumulative update of the local database. It prepares genome download lists, downloads/extracts/link genomes, runs quality/annotation tools, and collates genome/protein summary tables.

## Configuration

Set `download.source` to `gtotree` to obtain GTDB-verified accessions and
metadata with `gtt-get-accessions-from-GTDB -t <download.gtt_taxon>` (the
default taxon is `asgardarchaeota`), or to `ncbi` to retain the NCBI
taxon-search workflow. `gtt` remains accepted as a backwards-compatible alias
for `gtotree`. Downloads are still retrieved with NCBI Datasets after the
accession list is selected, so both sources produce the same normalized FASTA
layout.

`run.update: false` makes a fresh isolated set: prior store and annotation
contents are ignored. Set it to `true` to include previously selected genomes
for the same organism alongside the new selection; the store membership
manifest prevents genomes from another configured organism from being added.

## Logs

The long-running acquisition, quality-control, classification, and annotation
rules create individual log files beneath
`<run.log_dir>/updation/<run.id>/`. Logs are further grouped by rule and, for
Bakta and InterProScan, by organism and genome. This keeps diagnostics from a
given run separate from previous runs while preserving the output produced by
each external tool.

With `download.source: gtotree`, the workflow does not invoke CheckM2 or
GTDB-Tk. It reuses the retrieved `checkm2_*` and GTDB taxonomy fields, writes
a normalized CheckM2-compatible quality report, and writes the selected
metadata to `<run>/gtdbtk/<org>/gtdbtk.ar53.summary.tsv`. The latter is the
same phylogeny metadata handoff expected from GTDB-Tk, so downstream placement
and annotation steps use the retrieved metadata directly. GToTree metadata
must include an accession column, `checkm2_completeness`,
`checkm2_contamination`, and `gtdb_taxonomy` (or `classification`).

For non-GToTree runs, use `deduplication.enabled: false` (the default) to
skip skder and pass every CheckM2-passing genome onward. Set it to `true` to
enable ANI dereplication.

To use a cached GTDB metadata release with a non-GToTree source, set
`gtdb.run_gtdbtk: false`, set `gtdb.metadata_url` to the archaeal GTDB
metadata release URL, and choose the local cache path with
`gtdb.metadata_file`. The pipeline downloads and decompresses the metadata
when the cache is absent. Set `gtdb.run_gtdbtk: true` to run GTDB-Tk instead.
`run_gtdbtk` defaults to `true` for compatibility with existing configurations
that specify only `gtdb.db_path`.

## Main rules and sub-rules

This pipeline is self-contained in its `Snakefile` and does not include separate `.smk` sub-rule files.

| Rule | What it does |
| --- | --- |
| `check_database` | Checks whether the expected database inputs/outputs are present. |
| `get_summary` | Retrieves or prepares metadata summaries used for genome selection. |
| `prepare_download_list` | Creates the list of genomes to download. |
| `download_genomes` | Downloads selected genomes. |
| `extract_genomes` | Extracts downloaded genome archives. |
| `link_run_genomes` | Links genomes into the current run workspace. |
| `add_manual_genome` | Adds manually supplied genomes into the run. |
| `checkm2` | Runs CheckM2 quality assessment. |
| `link_filtered` | Links genomes that pass quality/filter criteria. |
| `skder` | Runs sketching/dereplication support steps for genomes. |
| `gtdbtk` | Runs GTDB-Tk taxonomy classification. |
| `taxa_counts` | Plots taxonomic count summaries. |
| `update_master_checkm2` | Updates the master CheckM2 quality table. |
| `update_master_gtdbtk` | Updates the master GTDB-Tk taxonomy table. |
| `bakta` | Runs Bakta genome annotation. |
| `interproscan` | Runs InterProScan for an annotated genome/proteome. |
| `all_interproscan` | Aggregates completion of InterProScan jobs. |
| `collate_proteins` | Collates protein FASTA/annotation data from run outputs. |
| `make_genome_csv` | Builds the updated genome metadata CSV. |
| `make_protein_csv` | Builds the updated protein metadata CSV. |

## Python scripts called

No standalone Python scripts are invoked directly by the current `Snakefile`. Most processing is performed through shell commands and external tools.

## Non-Python helpers called

The workflow calls external bioinformatics tools such as CheckM2, GTDB-Tk, Bakta, InterProScan, and R script `plot_taxa_counts.r` for taxonomy-count plots.
