import os
envvars:
    "ITOL_API_KEY"
########################################
# Sequence Length Histogram
########################################


# Defines if the final tree needs to be rooted using madroot tool
# Yet to be functional, MadRoot requires the sequences to be deduplicated
# which is not done as iqtree handles duplicates
def get_final_tree(wc):
    if config.get("phylogeny", []).get("madroot", False):
        return f"{EXPLORATION_DIR}/{PROTEIN}_unr_fasttree.rooted.treefile"
    else:
        return f"{EXPLORATION_DIR}/{PROTEIN}_unr_fasttree.treefile"


# gives a length analysis of the protein set
rule length_histogram:
    input:
        fasta = f"{EXPLORATION_DIR}/{PROTEIN}.unr.fasta"
    output:
        plot = f"{EXPLORATION_DIR}/{PROTEIN}_length_hist.png",
        split_dir = directory(f"{EXPLORATION_DIR}/{PROTEIN}_length_bins")
    conda:
        f"{ENV_DIR}/Reg.yaml"
    message:
        """
        ==========================================
        📊 Sequence Length Analysis + Binning
        ==========================================
        """
    shell:
        """
        python {CURRENT_DIR}/bin/units/plot_lengths.py \
            {input.fasta} \
            {output.plot} \
            {output.split_dir}
        """



# plots the counts of the genomes in different classes
rule taxa_counts:
    input:
        csv = f"{EXPLORATION_DIR}/{PROTEIN}.unr.csv"
    output:
        plot = f"{EXPLORATION_DIR}/{PROTEIN}_taxa_count.svg",
    conda:
        f"{ENV_DIR}/bio-r.yaml"
    message:
        """
        ==========================================
        📊      Taxonomic Sampling Plot
        ==========================================
        """
    shell:
        """
        Rscript {CURRENT_DIR}/bin/units/plot_taxa_counts.r \
            -i {input.csv} \
            -o {output.plot} \
        """




############################################
# Exploratory FastTree
############################################

rule exploratory_fasttree:
    input:
        fasta = f"{EXPLORATION_DIR}/{PROTEIN}.unr.fasta"
    output:
        tree = f"{EXPLORATION_DIR}/{PROTEIN}_unr_fasttree.treefile",
        msa = f"{EXPLORATION_DIR}/{PROTEIN}_unr_fasttree.aligned.fasta"
    threads: config.get("phylogeny", {}).get("threads", 8)
    conda:
        f"{ENV_DIR}/phylogeny.yaml"
    params:
        prefix = f"{EXPLORATION_DIR}/{PROTEIN}_unr_fasttree"
    message:
        """
        ==========================================
        🌳 Exploratory FastTree
        ==========================================
        """
    shell:
        """
        bash {CURRENT_DIR}/bin/units/fasttree_pipeline.sh \
            {input.fasta} \
            {params.prefix} \
            {threads}
        """


########################################
# Generate iTOL Colorstrip
########################################
# generates color strips for headers of the protein sequences
rule itol_colorstrip:
    input:
        fasta = f"{EXPLORATION_DIR}/{PROTEIN}.unr.fasta"
    output:
        f"{EXPLORATION_DIR}/{PROTEIN}_colorstrip.txt"
    shell:
        """
        python {CURRENT_DIR}/bin/units/generate_colorstrip.py \
            {input.fasta} \
            {output} \
            "header"
        """

# compiles the dataframe that is required to be passed to table2itol to get 
# all the annotation files required for the visualization
rule df_for_annotation:
    input: 
        fasta = f"{EXPLORATION_DIR}/{PROTEIN}.unr.fasta", 
        protein_csv = f"{EXPLORATION_DIR}/{PROTEIN}.unr.csv",
        cluster_csv = f"{SSN_DIR}/{PROTEIN}.clusters.expanded.csv",
        domains = f"{EXPLORATION_DIR}/{PROTEIN}_domain_proteins.tsv"
    output:
        annot_csv = f"{PHYLO_DIR}/{PROTEIN}.annot.csv"
    params:
        seq_id = config.get("run", []).get("LOCUS_TAG", "locus_tag")
    conda:
        f"{ENV_DIR}/Reg.yaml"
    shell:
        """
        python {CURRENT_DIR}/bin/units/get_annotation_csv.py \
                {input.fasta} \
                {input.protein_csv} \
                {input.cluster_csv} \
                {input.domains} \
                {output.annot_csv} \
                {params.seq_id}      
        
        """

rule msa_to_itol:
    input: 
        msa = f"{EXPLORATION_DIR}/{PROTEIN}_unr_fasttree.aligned.fasta"
    output:
        itol_msa = f"{PHYLO_DIR}/itol_msa.txt"
    conda:
        f"{ENV_DIR}/Reg.yaml"
    shell:
        """
        python {CURRENT_DIR}/bin/units/msa_to_itol_dataset.py \
            {input.msa} \
            {output.itol_msa} \
        """

rule table2itol:
    input:
        annot = f"{PHYLO_DIR}/{PROTEIN}.annot.csv"
    output:
        annotation_dir = directory(f"{PHYLO_DIR}/annotation"),    
        done_annot = touch(f"{PHYLO_DIR}/annotation.done.flag")
    conda:
        f"{ENV_DIR}/bio-r.yaml"
    params:
        LOCUS_TAG = config.get('run', []).get("LOCUS_TAG", "locus_tag")
    shell:
        """
        mkdir -p {output.annotation_dir}

        Rscript {CURRENT_DIR}/bin/units/table2itol.R \
            -i {params.LOCUS_TAG} \
            -l {params.LOCUS_TAG} \
            -s "," \
            -D {output.annotation_dir} \
            {input.annot}
        """

# madRoot

rule madroot:
    input:
        tree = f"{EXPLORATION_DIR}/{PROTEIN}_unr_fasttree.treefile"
    output:
        rooted_tree = f"{EXPLORATION_DIR}/{PROTEIN}_unr_fasttree.rooted.treefile"

    threads: config.get("phylogeny", []).get("threads", 8)
    conda:
        f"{ENV_DIR}/phylogeny.yaml"
    params:
        prefix = f"{EXPLORATION_DIR}/{PROTEIN}_unr_fasttree"
    message:
        """
        ==========================================
        🌳 rooting tree using madroot
        ==========================================
        """
    shell:
        """
        madRoot {input.tree} > {output.rooted_tree}
        """


########################################
# Upload to iTOL
########################################

rule upload_to_itol:
    input:
        tree        = get_final_tree,
        colorstrip  = f"{EXPLORATION_DIR}/{PROTEIN}_colorstrip.txt",
        annot_files = lambda wildcards: sorted(
            glob.glob(f"{PHYLO_DIR}/annotation/*.txt")
        ),
        default     = lambda wildcards: config.get("default_annotation", []),
        marker      = f"{PHYLO_DIR}/annotation.done.flag",
        itol_dir    = f"{EXPLORATION_DIR}/{PROTEIN}_itol_domains",
        itol_msa = f"{PHYLO_DIR}/itol_msa.txt"
    output:
        tree_ids = f"{EXPLORATION_DIR}/{PROTEIN}_fast_itol_uploaded.flag"

    params:
        project   = config.get("itol_project", "Asgard"),
        tree_name = f"fast_{config.get('run', []).get('id', 'run')}_{PROTEIN}",
        all_annots = lambda wildcards, input: (
            sorted(glob.glob(f"{PHYLO_DIR}/annotation/*.txt"))
            + ([input.colorstrip])
            + (input.default if isinstance(input.default, list) else [input.default])
        ),
        domain_files = lambda wildcards, input: sorted(
            glob.glob(f"{input.itol_dir}/*.txt")
        ),
    shell:
        """
        bash {CURRENT_DIR}/bin/units/itol_upload.sh \
            {input.tree}            \
            {params.project}        \
            {params.tree_name}      \
            {output.tree_ids}       \
            {params.all_annots:q}   \
            {params.domain_files:q} \
            {input.itol_msa}

        """




