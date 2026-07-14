rule ips_acc_parse:
    input:
        database     = config["raw_database"],  # raw InterPro parquet (all analyses)
    output:
        outfile     = f"{EXPLORATION_DIR}/{PROTEIN}.ids",
        itol_dir    = directory(f"{EXPLORATION_DIR}/{PROTEIN}_itol_domains"),
        matching_tsv = f"{EXPLORATION_DIR}/{PROTEIN}_domain_proteins.tsv"
    params:
        acc=config['acc_groups']
    conda:
        f"{config['env_dir']}/duckdb_handler.yaml"
    message:
        """
        ===============================
        Running parse_ips
        ===============================
        """
    shell:
        """
        python {CURRENT_DIR}/bin/units/parse_ips_acc.py \
            --parquet {input.database} \
            --acc {params.acc} \
            --protein-ids {output.outfile} \
            --matching-tsv {output.matching_tsv} \
            --itol {output.itol_dir}
        """

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
    shell:
        """        
        python {CURRENT_DIR}/bin/units/get_fasta_csv_from_ids.py \
        --fasta {input.fasta} \
        --csv {input.protein_file} \
        --ids {input.protein_ids} \
        --genome_file {input.genome_file} \
        --outfasta {output.outfasta} \
        --outcsv {output.protein_csv} \
        --protein_name {params.protein_name}\
        --remove_hypotheticals {params.remove_hypotheticals}      
        
        """