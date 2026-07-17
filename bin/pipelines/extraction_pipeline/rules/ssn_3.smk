
"""
This smk file codes for the sequence similarity network part of the 
exploration pipeline. 
"""


# Defines the different thresholds at which the SSNs need to be created.
# Any 2 proteins (nodes) are joined with an edge provided it
# has a bitscore greater than the threshold
BITSCORES = config.get('SSN', []).get("SSN_BITSCORES", [50, 80, 100])

# =============================================================================
# Rule 1: CD-HIT clustering
# =============================================================================
"""
Sequences are first deduplicated to reduce the computational expense as
well as the rendering in cytoscape

CD-hit is used for the deduplication with cut-offs defined in the config file
Defaults are used in the absence of the cut-offs
"""
rule ssn_cdhit:
    input:
        fasta = f"{EXPLORATION_DIR}/{PROTEIN}.unr.fasta",
    output:
        nr = f"{SSN_DIR}/{PROTEIN}.nr.fasta",
        clstr = f"{SSN_DIR}/{PROTEIN}.nr.fasta.clstr",
    params:
        identity = config.get("SSN_IDENTITY", 0.95),
    threads: config.get("SSN", config.get("cores", 16)).get("SSN_CORES", config.get("cores", 16))
    resources:
        mem_mb  = config.get('SSN', []).get("mem_mb", 16000),
        runtime = config.get('SSN', []).get("runtime_min", 120)
    log:
        f"{LOG_DIR}/ssn/{PROTEIN}.cdhit.log"
    benchmark:
        f"{BENCHMARK_DIR}/ssn/{PROTEIN}.cdhit.benchmark.tsv"
    conda:
        f"{ENV_DIR}/ssn.yaml"
    shell:
        """
        bash {CURRENT_DIR}/bin/units/ssn_cdhit.sh \
            {input.fasta}     \
            {output.nr}       \
            {params.identity} \
            {threads}         \
            2>&1 | tee {log}
        """


# =============================================================================
# Rule 2: All-vs-all similarity search
# =============================================================================
"""
Creates a similarities tsv which contains all pairs of nodes
and the BLAST results between the pair. Here the pairs with evalue less
that 1e-5 are removed
"""
rule ssn_search:
    input:
        nr = f"{SSN_DIR}/{PROTEIN}.nr.fasta",
    output:
        similarities = f"{SSN_DIR}/{PROTEIN}.similarities.tsv",
    params:
        evalue = config.get('SSN').get("SSN_EVALUE", 1e-5),
    threads: config.get("SSN", config.get("cores", 16)).get("SSN_CORES", config.get("cores", 16))
    resources:
        mem_mb  = config.get('SSN', []).get("mem_mb", 16000),
        runtime = config.get('SSN', []).get("runtime_min", 120)
    log:
        f"{LOG_DIR}/ssn/{PROTEIN}.search.log"
    benchmark:
        f"{BENCHMARK_DIR}/ssn/{PROTEIN}.search.benchmark.tsv"
    conda:
        f"{ENV_DIR}/ssn.yaml"
    shell:
        """
        bash {CURRENT_DIR}/bin/units/ssn_search.sh \
            {input.nr}            \
            {output.similarities} \
            {threads}             \
            {params.evalue}       \
            2>&1 | tee {log}
        """


# =============================================================================
# Rule 3a: Node metadata — runs once, shared across all bitscore levels
# =============================================================================
"""
Contains node metadata - protein length

"""
rule ssn_nodes:
    input:
        nr = f"{SSN_DIR}/{PROTEIN}.nr.fasta",
    output:
        nodes = f"{SSN_DIR}/{PROTEIN}.nodes.tsv",
    resources:
        mem_mb  = 2000,
        runtime = 10
    log:
        f"{LOG_DIR}/ssn/{PROTEIN}.nodes.log"
    benchmark:
        f"{BENCHMARK_DIR}/ssn/{PROTEIN}.nodes.benchmark.tsv"
    conda:
        f"{ENV_DIR}/ssn.yaml"
    shell:
        """
        bash {CURRENT_DIR}/bin/units/ssn_nodes.sh \
            {input.nr}     \
            {output.nodes} \
            2>&1 | tee {log}
        """


# =============================================================================
# Rule 3b: Edge filtering — one job per bitscore level
# =============================================================================
"""
creates the edges.tsv file based on the similarities file by creating
different filters of bitscores (as defined in the start of the script)

"""
rule ssn_filter:
    input:
        similarities = f"{SSN_DIR}/{PROTEIN}.similarities.tsv",
        nodes = f"{SSN_DIR}/{PROTEIN}.nodes.tsv"
    output:
        edges = f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.edges.tsv",
    params:
        bitscore = lambda wc: wc.bitscore,
        coverage = config.get('SSN', []).get("coverage", 0.6),
        evalue   = config.get("SSN", []).get('evalue',   1e-5),
    resources:
        mem_mb  = 4000,
        runtime = 30
    log:
        f"{LOG_DIR}/ssn/{PROTEIN}.bs{{bitscore}}.filter.log"
    benchmark:
        f"{BENCHMARK_DIR}/ssn/{PROTEIN}.bs{{bitscore}}.filter.benchmark.tsv"
    conda:
        f"{ENV_DIR}/ssn.yaml"
    shell:
        """
        bash {CURRENT_DIR}/bin/units/ssn_filter.sh \
            {input.similarities} \
            {output.edges}       \
            {params.bitscore}    \
            {params.coverage}    \
            {params.evalue}      \
            {input.nodes}        \
            2>&1 | tee {log}
        """


# =============================================================================
# Rule 4: Cytoscape export
# =============================================================================
"""
creates files in formats readable by cytoscape
ea has protein1, protein2, score
sif has protein1, relation(here similarity), protein2

"""
rule ssn_cytoscape:
    input:
        edges = f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.edges.tsv",
        nodes = f"{SSN_DIR}/{PROTEIN}.nodes.tsv",
    output:
        sif = f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.edges.sif",
        ea  = f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.edges.ea",
    resources:
        mem_mb  = 2000,
        runtime = 10
    log:
        f"{LOG_DIR}/ssn/{PROTEIN}.bs{{bitscore}}.cytoscape.log"
    benchmark:
        f"{BENCHMARK_DIR}/ssn/{PROTEIN}.bs{{bitscore}}.cytoscape.benchmark.tsv"
    conda:
        f"{ENV_DIR}/ssn.yaml"
    shell:
        """
        bash {CURRENT_DIR}/bin/units/ssn_cytoscape.sh \
            {input.edges} \
            {output.sif}  \
            {output.ea}   \
            2>&1 | tee {log}
        """


# =============================================================================
# Rule 5: Annotation / taxonomy
# =============================================================================
"""
creates the taxnomy csv for the annotation in cytoscape

"""
rule ssn_annotate:
    input:
        nodes = f"{SSN_DIR}/{PROTEIN}.nodes.tsv",
        csv   = f"{EXPLORATION_DIR}/{PROTEIN}.unr.csv",
    output:
        taxonomy = f"{SSN_DIR}/{PROTEIN}.tax.tsv",
    params:
        locus_tag = config['run'].get("LOCUS_TAG", "locus_tag")
    resources:
        mem_mb  = 2000,
        runtime = 10
    log:
        f"{LOG_DIR}/ssn/{PROTEIN}.annotate.log"
    benchmark:
        f"{BENCHMARK_DIR}/ssn/{PROTEIN}.annotate.benchmark.tsv"
    conda:
        f"{ENV_DIR}/ssn.yaml"
    shell:
        """
        bash {CURRENT_DIR}/bin/units/ssn_annotate.sh \
            {input.csv}       \
            {output.taxonomy} \
            {params.locus_tag} \
            2>&1 | tee {log}
        """


# =============================================================================
# Rule 6: Cluster assignment + per-cluster FASTAs
# =============================================================================
"""
Allots cluster ids to each protein based on their cytoscape cluster
this is used to visualize the clusters during tree visualization in iTOL
"""
rule ssn_cluster:
    input:
        nodes = f"{SSN_DIR}/{PROTEIN}.nodes.tsv",
        nr    = f"{SSN_DIR}/{PROTEIN}.nr.fasta",
        clstr = f"{SSN_DIR}/{PROTEIN}.nr.fasta.clstr",
        edges = expand(
            f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.edges.tsv",
            bitscore=BITSCORES
        )
    output:
        csv       = f"{SSN_DIR}/{PROTEIN}.clusters.csv",
        expanded_csv       = f"{SSN_DIR}/{PROTEIN}.clusters.expanded.csv",
        fasta_dir = directory(f"{SSN_DIR}/cluster_fastas"),
    params:
        bitscores = ",".join(str(b) for b in BITSCORES)
    resources:
        mem_mb  = 8000,
        runtime = 60
    log:
        f"{LOG_DIR}/ssn/{PROTEIN}.cluster.log"
    benchmark:
        f"{BENCHMARK_DIR}/ssn/{PROTEIN}.cluster.benchmark.tsv"
    conda:
        f"{ENV_DIR}/Reg.yaml"
    shell:
        """
        python {CURRENT_DIR}/bin/units/ssn_cluster.py \
            --nodes     {input.nodes}      \
            --fasta     {input.nr}         \
            --edges     {input.edges}      \
            --clstr     {input.clstr}      \
            --bitscores {params.bitscores} \
            --out-csv   {output.csv}       \
            --out-dir   {output.fasta_dir} \
            2>&1 | tee {log}
        """


# =============================================================================
# Aggregator (calls for the outputs for all bitscores in BITSCORES)
# =============================================================================


rule ssn_network:
    input:
        expand(
            [
                f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.edges.tsv",
                f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.edges.sif",
                f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.edges.ea",
            ],
            bitscore=BITSCORES
        ),
        nodes     = f"{SSN_DIR}/{PROTEIN}.nodes.tsv",
        clusters  = f"{SSN_DIR}/{PROTEIN}.clusters.csv",
        expanded_csv       = f"{SSN_DIR}/{PROTEIN}.clusters.expanded.csv",
        fasta_dir = f"{SSN_DIR}/cluster_fastas",
        tax_annotation = f"{SSN_DIR}/{PROTEIN}.tax.tsv",


# =============================================================================
# Optional: copy to Windows
# =============================================================================

"""
Moves it to windows system
The original user was working in WSL therefore this need 
"""
if config.get('SSN', {}).get("copy_to_windows", False):

    rule copy_outs_to_windows:
        input:
            per_bs = expand(
                [
                    f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.edges.ea",
                    f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.edges.sif",
                ],
                bitscore=BITSCORES
            ),
            nodes = f"{SSN_DIR}/{PROTEIN}.nodes.tsv",
            clusters = f"{SSN_DIR}/{PROTEIN}.clusters.expanded.csv",
            fasta_dir = f"{SSN_DIR}/cluster_fastas",
            tax_annotation = f"{SSN_DIR}/{PROTEIN}.tax.tsv",
            full_csv = f"{EXPLORATION_DIR}/{PROTEIN}.unr.csv"

        output:
            flag = f"{SSN_DIR}/copied_to_windows.flag",

        log:
            f"{LOG_DIR}/copy_outs_to_windows.log"

        params:
            outdir = f"{config['SSN']['windows_path']}/{PROTEIN}/{RUN_ID}"

        shell:
            r"""
            exec > >(tee "{log}") 2>&1
            set -euxo pipefail

            mkdir -p "{params.outdir}"

            # Copy all edge files
            rsync -av --progress \
                {input.per_bs} \
                "{params.outdir}/"

            # Copy remaining files
            rsync -av --progress \
                "{input.nodes}" \
                "{input.clusters}" \
                "{input.tax_annotation}" \
                "{input.full_csv}" \
                "{params.outdir}/"

            # Copy FASTA directory
            rsync -av --progress \
                "{input.fasta_dir}/" \
                "{params.outdir}/cluster_fastas/"

            touch "{output.flag}"
            """