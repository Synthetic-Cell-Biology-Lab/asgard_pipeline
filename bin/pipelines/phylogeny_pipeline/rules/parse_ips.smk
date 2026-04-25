########################################
# Parse IPS
########################################

rule parse_ips:
    input:
        database     = config["database"],      # protein_summary parquet (ipr_desc aggregated)
        raw_database = config["raw_database"],  # raw InterPro parquet (all analyses)
    output:
        outfile     = f"{EXPLORATION_DIR}/{PROTEIN}_domain_proteins.tsv",
        protein_ids = f"{EXPLORATION_DIR}/{PROTEIN}.ids",
        itol_dir    = directory(f"{EXPLORATION_DIR}/{PROTEIN}_itol_domains"),
    params:
        search_string = config.get("search_string", None),
        rstring       = config.get("rstring", None),
    conda:
        f"{config['env_dir']}/duckdb_handler.yaml"
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
        protein_file = config["protein_file"],
        fasta = config["fasta_file"],
        protein_ids = f"{EXPLORATION_DIR}/{PROTEIN}.ids",
        genome_file = config['genome_file'],
        
    output:
        outfasta = f"{EXPLORATION_DIR}/{PROTEIN}.unr.fasta",
        protein_csv = f"{EXPLORATION_DIR}/{PROTEIN}.unr.csv"
    params:
        remove_hypotheticals = config.get("remove_hypotheticals", False),
        protein_name = config['protein_name']
    conda:
        f"{config['env_dir']}/Reg.yaml"
    message:
        """
        ===============================
        Running merge_file
        ===============================
        """
    script:
        f"{CURRENT_DIR}/bin/units/get_fasta_csv_from_ids.py"
