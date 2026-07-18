########################################
# Parse IPS
########################################
"""
USED TO RETRIEVE PROTEINS GIVEN A DOMAIN
This smk parses the interproscan collated result parquet and retrieves proteins 
that match the given constraints of the search term
The search term can be in regex or regular search. This searches for the search term in the
sig_desc column of the collated database

"""


rule parse_ips:
    input:
        database     = config['inputs']["database"],      # protein_summary parquet (ipr_desc aggregated)
        raw_database = config['inputs']["raw_database"],  # raw InterPro parquet (all analyses)
    output:
        outfile     = f"{EXPLORATION_DIR}/{PROTEIN}_domain_proteins.tsv",
        protein_ids = f"{EXPLORATION_DIR}/{PROTEIN}.ids",
        itol_dir    = directory(f"{EXPLORATION_DIR}/{PROTEIN}_itol_domains"),
    params:
        search_string = config['modes']['parse_ips']['search']['search_string'], 
        rstring       = config['modes']['parse_ips']['search']['rstring'],
    conda:
        f"{ENV_DIR}/duckdb_handler.yaml"
    message:
        """
        ===============================
        Running parse_ips
        ===============================
        """
    # run:
    #     import sys
    #     print(sys.executable)
    #     print(sys.path)
        
    script:
        f"{CURRENT_DIR}/bin/units/parse_ips4.py"


########################################
# Extract FASTA + CSV
########################################

rule merge_file:
    input:
        protein_file = config['inputs']["protein_file"],
        fasta = config['inputs']["fasta_file"],
        protein_ids = f"{EXPLORATION_DIR}/{PROTEIN}.ids",
        genome_file = config['inputs']['genome_file'],
        
    output:
        outfasta = f"{EXPLORATION_DIR}/{PROTEIN}.unr.fasta",
        protein_csv = f"{EXPLORATION_DIR}/{PROTEIN}.unr.csv"
    params:
        remove_hypotheticals = config.get("modes", {}).get("parse_ips", {}).get("remove_hypotheticals", False),
        protein_name = config['run']['protein_name'],
        taxon_level = config['run'].get('taxon_level', 'phylum'),
        taxon_filter = config['run'].get('taxon_filter', 'Asgardarchaeota')
    conda:
        f"{ENV_DIR}/Reg.yaml"
    message:
        """
        ===============================
        Running merge_file
        ===============================
        """
    script:
        f"{CURRENT_DIR}/bin/units/get_fasta_csv_from_ids.py"
