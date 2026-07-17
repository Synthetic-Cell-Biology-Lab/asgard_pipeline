rule from_manual:
    """
    Inputs:
    protein_file: csv file containing locus_tag (protein_id in the fasta), 
    genome_file (the source genome of the sequence), manual_annotation(previously defined annotation of the protein), and others
    genome_file: consists of the genome metadata, particularly the contamination and completetion metrics, and the taxonomy details
    fasta_file: .faa file that collates all the proteins in genome set

    This rule retrieves all the sequences that fall within the constraints that are
    given the params (the level of taxon filtering and which column of manual_annotation in the csv the script should filter on)

    This is the first rule that runs during a phylogeny run (during the exploration phase)
   
    """
    input:
        protein_file = config['inputs']['protein_file'],
        genome_file = config['inputs']['genome_file'],
        fasta_file = config['inputs']['fasta_file']
    params:
        taxonomic_level = config.get("run", []).get('taxon_level', None), # example phylum, class, order etc
        taxon_filter = config.get('run', []).get("taxon_filter", None), # the specific name of the taxa, for eg, if the taxon_level is phylum, this would be asgardarchaeota
        annotation_filter = config['modes']['from_manual'].get("manual_annotation_filter", []) # the name of the column that contains the manual_annotation
    output:
        protein_ids = f"{EXPLORATION_DIR}/{PROTEIN}.ids",
        outfasta = f"{EXPLORATION_DIR}/{PROTEIN}.unr.fasta",
        protein_csv = f"{EXPLORATION_DIR}/{PROTEIN}.unr.csv"

    conda:
        f"{ENV_DIR}/Reg.yaml"
    shell:
        """
        python {CURRENT_DIR}/bin/units/from_manual.py \
            {input.protein_file} \
            {input.genome_file} \
            {input.fasta_file} \
            {output.protein_ids} \
            {output.outfasta} \
            {output.protein_csv} \
            {params.taxonomic_level} \
            {params.taxon_filter} \
            {params.annotation_filter}
        """



rule parse_ips:
    """
    This rule parses the interproscan collated results (parquet) and gives annotation files
    for mapping the domains in the phylogenetic tree visualization
    
    """
    input:
        raw_database = config["inputs"]["raw_database"],  # raw InterPro parquet (all analyses)
        protein_ids = f"{EXPLORATION_DIR}/{PROTEIN}.ids",
    output:
        outfile     = f"{EXPLORATION_DIR}/{PROTEIN}_domain_proteins.tsv",
        itol_dir    = directory(f"{EXPLORATION_DIR}/{PROTEIN}_itol_domains"),
    conda:
        f"{ENV_DIR}/duckdb_handler.yaml"
    message:
        """
        ===============================
        Running parse_ips
        ===============================
        """
    script:
        f"{CURRENT_DIR}/bin/units/parse_ips_manuall_annot_ver.py"
