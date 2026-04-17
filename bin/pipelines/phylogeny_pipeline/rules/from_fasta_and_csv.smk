rule copy_inputs:
    input:
        csv = config['INPUT_CSV'],
        fasta = config['INPUT_FASTA'],
    output:
        outfasta = f"{EXPLORATION_DIR}/{PROTEIN}.unr.fasta",
        protein_csv = f"{EXPLORATION_DIR}/{PROTEIN}.unr.csv"
    conda:
        f"{config['env_dir']}/Reg.yaml"
    shell:
        """
        cp {input.fasta} {output.outfasta}
        cp {input.csv} {output.protein_csv}
        """
