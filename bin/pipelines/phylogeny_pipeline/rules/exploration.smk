import os

########################################
# Sequence Length Histogram
########################################

rule length_histogram:
    input:
        fasta = f"{EXPLORATION_DIR}/{PROTEIN}.unr.fasta"
    output:
        plot = f"{EXPLORATION_DIR}/{PROTEIN}_length_hist.png",
        split_dir = directory(f"{EXPLORATION_DIR}/{PROTEIN}_length_bins")
    conda:
        f"{config['env_dir']}/Reg.yaml"
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






############################################
# Exploratory FastTree
############################################

rule exploratory_fasttree:
    input:
        fasta = f"{EXPLORATION_DIR}/{PROTEIN}.unr.fasta"
    output:
        tree = f"{EXPLORATION_DIR}/{PROTEIN}_unr_fasttree.treefile",
        msa = f"{EXPLORATION_DIR}/{PROTEIN}_unr_fasttree.aligned.fasta"
    threads: config.get("phylogeny_threads", 8)
    conda:
        f"{config['env_dir']}/phylogeny.yaml"
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


rule df_for_annotation:
    input: 
        fasta = f"{EXPLORATION_DIR}/{PROTEIN}.unr.fasta", 
        protein_csv = f"{EXPLORATION_DIR}/{PROTEIN}.unr.csv",
        cluster_csv = f"{SSN_DIR}/{PROTEIN}.clusters.expanded.csv",
        domains = f"{EXPLORATION_DIR}/{PROTEIN}_domain_proteins.tsv"
    output:
        annot_csv = f"{PHYLO_DIR}/{PROTEIN}.annot.csv"
    conda:
        f"{config['env_dir']}/Reg.yaml"
    shell:
        """
        python {CURRENT_DIR}/bin/units/get_annotation_csv.py \
                {input.fasta} \
                {input.protein_csv} \
                {input.cluster_csv} \
                {input.domains} \
                {output.annot_csv}        
        
        """

rule msa_to_itol:
    input: 
        msa = f"{EXPLORATION_DIR}/{PROTEIN}_unr_fasttree.aligned.fasta"
    output:
        itol_msa = f"{PHYLO_DIR}/itol_msa.txt"
    conda:
        f"{config['env_dir']}/Reg.yaml"
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
        f"{config['env_dir']}/bio-r.yaml"
    shell:
        """
        mkdir -p {output.annotation_dir}

        Rscript {CURRENT_DIR}/bin/units/table2itol.R \
            -i locus_tag \
            -l locus_tag \
            -s "," \
            -D {output.annotation_dir} \
            {input.annot}
        """



########################################
# Upload to iTOL
########################################

rule upload_to_itol:
    input:
        tree        = f"{EXPLORATION_DIR}/{PROTEIN}_unr_fasttree.treefile",
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
        tree_name = f"fast_{config.get('run_id', 'run')}_{PROTEIN}",
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




########################################
# Manual Review Gate
########################################

REVIEW_GATE_INPUTS = [
    f"{EXPLORATION_DIR}/{PROTEIN}.unr.fasta",
    f"{EXPLORATION_DIR}/{PROTEIN}_unr_fasttree.treefile"
]

if config.get("run_itol_upload", False):
    REVIEW_GATE_INPUTS.append(f"{EXPLORATION_DIR}/{PROTEIN}_fast_itol_uploaded.flag")



rule review_gate:
    input:
        REVIEW_GATE_INPUTS
    output:
        rev_fasta = f"{RESULT_DIR}/{PROTEIN}.rev.fasta",
        marker = f"{RESULT_DIR}/REVIEW_DONE.flag"
    run:
        import os

        if os.path.exists(output.rev_fasta) and os.path.exists(output.marker):
            print("✅ Manual review already completed.")
            return

        print("\n🔍 MANUAL REVIEW REQUIRED\n")
        print(f"Review FASTA: {input[0]}")
        print(f"Tree file: {input[1]}")

        if len(input) > 2:
            print(f"iTOL upload flag: {input[2]}")
        print(f"Save curated file as: {output.rev_fasta}")
        print(f"Then run: touch {output.marker}\n")

        raise SystemExit(1)

rule make_rev_csv:
    input:
        rev_fasta = f"{RESULT_DIR}/{PROTEIN}.rev.fasta",
        marker = f"{RESULT_DIR}/REVIEW_DONE.flag",
        unr_csv = f"{EXPLORATION_DIR}/{PROTEIN}.unr.csv"
    output:
        rev_csv = f"{EXPLORATION_DIR}/{PROTEIN}.rev.csv"
    
    conda:
        f"{config['env_dir']}/Reg.yaml"

    run:
        import pandas as pd

        # ── Extract locus_tags from FASTA ───────────────────────
        locus_tags = set()

        with open(input.rev_fasta) as f:
            for line in f:
                if line.startswith(">"):
                    locus = line[1:].strip().split()[0]
                    locus_tags.add(locus)

        # ── Load CSV ────────────────────────────────────────────
        df = pd.read_csv(input.unr_csv)

        if "locus_tag" not in df.columns:
            raise ValueError("Column 'locus_tag' not found in CSV")

        # Clean just in case
        df["locus_tag"] = df["locus_tag"].astype(str).str.strip()

        # ── Filter ──────────────────────────────────────────────
        filtered = df[df["locus_tag"].isin(locus_tags)]

        # ── Save ────────────────────────────────────────────────
        filtered.to_csv(output.rev_csv, index=False)

        # ── Debug (optional but useful) ─────────────────────────
        print(f"Extracted {len(locus_tags)} locus tags")
        print(f"Matched {len(filtered)} rows")