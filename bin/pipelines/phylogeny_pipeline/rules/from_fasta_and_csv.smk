rule copy_inputs:
    input:
        csv = config['INPUT_CSV'],
        fasta = config['INPUT_FASTA'],
        domain_files = config['INPUT_DOMAINS']
    output:
        outfasta = f"{EXPLORATION_DIR}/{PROTEIN}.unr.fasta",
        protein_csv = f"{EXPLORATION_DIR}/{PROTEIN}.unr.csv",
        domains_tsv = f"{EXPLORATION_DIR}/{PROTEIN}_domain_proteins.tsv",
        domain_files = directory(f"{EXPLORATION_DIR}/{PROTEIN}_itol_domains")
    conda:
        f"{config['env_dir']}/Reg.yaml"
    params:
        protein = config['protein_name'],
        domain_files = config['INPUT_DOMAINS']
    shell:
        """
        cp {input.fasta} {output.outfasta}
        cp {input.csv} {output.protein_csv}

        # Create domains TSV
        awk -F',' 'NR>1 {{print $1 "\\t" "{params.protein}"}}' {input.csv} \
        | awk 'BEGIN {{print "protein\\tdomains"}} {{print}}' \
        > {output.domains_tsv}

        cp -r {params.domain_files} {output.domain_files}
        """