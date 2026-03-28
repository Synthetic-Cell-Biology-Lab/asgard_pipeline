# =============================================================================
# Rule 1: CD-HIT clustering (expensive, but only re-runs if input changes)
# =============================================================================
rule ssn_cdhit:
    input:
        fasta = f"{EXPLORATION_DIR}/{PROTEIN}.unr.fasta",
    output:
        nr    = f"{SSN_DIR}/{PROTEIN}.nr.fasta",
    params:
        identity = config.get("SSN_IDENTITY", 0.90),
    threads: config.get("SSN_CORES", config.get("cores", 16))
    resources:
        mem_mb   = config.get("SSN_MEM_MB", 16000),
        runtime  = config.get("SSN_RUNTIME_MIN", 120)
    log:
        f"{LOG_DIR}/ssn/{PROTEIN}.cdhit.log"
    benchmark:
        f"{BENCHMARK_DIR}/ssn/{PROTEIN}.cdhit.benchmark.tsv"
    conda:
        f"{config['env_dir']}/ssn.yaml"
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
# Rule 2: All-vs-all similarity search (most expensive — BLAST or MMseqs2)
# =============================================================================
rule ssn_search:
    input:
        nr = f"{SSN_DIR}/{PROTEIN}.nr.fasta",
    output:
        similarities = f"{SSN_DIR}/{PROTEIN}.similarities.tsv",
    params:
        evalue  = config.get("SSN_EVALUE", 1e-5),
    threads: config.get("SSN_CORES", config.get("cores", 16))
    resources:
        mem_mb  = config.get("SSN_MEM_MB", 16000),
        runtime = config.get("SSN_RUNTIME_MIN", 600)
    log:
        f"{LOG_DIR}/ssn/{PROTEIN}.search.log"
    benchmark:
        f"{BENCHMARK_DIR}/ssn/{PROTEIN}.search.benchmark.tsv"
    conda:
        f"{config['env_dir']}/ssn.yaml"
    shell:
        """
        bash {CURRENT_DIR}/bin/units/ssn_search.sh \
            {input.nr}           \
            {output.similarities}\
            {threads}            \
            {params.evalue}      \
            2>&1 | tee {log}
        """


# =============================================================================
# Rule 3: Edge filtering + node metadata (cheap — tweak bitscore/coverage freely)
# =============================================================================
rule ssn_filter:
    input:
        nr           = f"{SSN_DIR}/{PROTEIN}.nr.fasta",
        similarities = f"{SSN_DIR}/{PROTEIN}.similarities.tsv",
    output:
        edges = f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.edges.tsv",
        nodes = f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.nodes.tsv",
    params:
        bitscore = lambda wc: wc.bitscore,
        coverage = config.get("SSN_COVERAGE",  0.6),
        evalue   = config.get("SSN_EVALUE",    1e-5),
    resources:
        mem_mb  = 4000,
        runtime = 30
    log:
        f"{LOG_DIR}/ssn/{PROTEIN}.bs{{bitscore}}.filter.log"
    benchmark:
        f"{BENCHMARK_DIR}/ssn/{PROTEIN}.bs{{bitscore}}.filter.benchmark.tsv"
    conda:
        f"{config['env_dir']}/ssn.yaml"
    shell:
        """
        bash {CURRENT_DIR}/bin/units/ssn_filter.sh \
            {input.similarities} \
            {input.nr}           \
            {output.edges}       \
            {output.nodes}       \
            {params.bitscore}    \
            {params.coverage}    \
            {params.evalue}      \
            2>&1 | tee {log}
        """


# =============================================================================
# Rule 4: Cytoscape export + network stats (cheap, depends only on edges/nodes)
# =============================================================================
rule ssn_cytoscape:
    input:
        edges = f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.edges.tsv",
        nodes = f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.nodes.tsv",
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
        f"{config['env_dir']}/ssn.yaml"
    shell:
        """
        bash {CURRENT_DIR}/bin/units/ssn_cytoscape.sh \
            {input.edges} \
            {output.sif}  \
            {output.ea}   \
            2>&1 | tee {log}
        """


# =============================================================================
# Rule 5: Annotation / taxonomy (cheap, depends only on nodes + input CSV)
# =============================================================================
rule ssn_annotate:
    input:
        nodes = f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.nodes.tsv",
        csv   = f"{EXPLORATION_DIR}/{PROTEIN}.unr.csv",
    output:
        taxonomy = f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.tax.tsv",
    resources:
        mem_mb  = 2000,
        runtime = 10
    log:
        f"{LOG_DIR}/ssn/{PROTEIN}.bs{{bitscore}}.annotate.log"
    benchmark:
        f"{BENCHMARK_DIR}/ssn/{PROTEIN}.bs{{bitscore}}.annotate.benchmark.tsv"
    conda:
        f"{config['env_dir']}/ssn.yaml"
    shell:
        """
        bash {CURRENT_DIR}/bin/units/ssn_annotate.sh \
            {input.csv}      \
            {output.taxonomy}\
            2>&1 | tee {log}
        """

rule ssn_cluster:
    input:
        nodes = f"{SSN_DIR}/{PROTEIN}.nodes.tsv",
        nr    = f"{SSN_DIR}/{PROTEIN}.nr.fasta",
        edges = expand(
            f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.edges.tsv",
            bitscore=config.get("SSN_BITSCORES", [50, 80, 100])
        )
    output:
        csv       = f"{SSN_DIR}/{PROTEIN}.clusters.csv",
        fasta_dir = directory(f"{SSN_DIR}/cluster_fastas")
    params:
        bitscores = ",".join(str(b) for b in config.get("SSN_BITSCORES", [50, 80, 100]))
    log:
        f"{LOG_DIR}/ssn/{PROTEIN}.cluster.log"
    benchmark:
        f"{BENCHMARK_DIR}/ssn/{PROTEIN}.cluster.benchmark.tsv"
    conda:
        f"{config['env_dir']}/ssn.yaml"
    shell:
        """
        python {CURRENT_DIR}/bin/units/ssn_cluster.py \
            --nodes     {input.nodes}    \
            --fasta     {input.nr}       \
            --edges     {input.edges}    \
            --bitscores {params.bitscores} \
            --out-csv   {output.csv}     \
            --out-dir   {output.fasta_dir} \
            2>&1 | tee {log}
        """

# =============================================================================
# Aggregator rule — request all final outputs
# =============================================================================
rule ssn_network:
    input:
        rules.ssn_filter.output,
        rules.ssn_cytoscape.output,
        rules.ssn_annotate.output,
        expand(
            [
                f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.edges.tsv",
                f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.edges.sif",
                f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.tax.tsv",
            ],
            bitscore=config.get("SSN_BITSCORES", [50, 80, 100])
        )


        
# =============================================================================
# Optional: copy outputs to Windows
# =============================================================================
if config.get("copy_to_windows", False):

    rule copy_outs_to_windows:
        input:
            expand(
                [
                    f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.edges.ea",
                    f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.nodes.tsv",
                    f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.edges.sif",
                    f"{SSN_DIR}/{PROTEIN}.bs{{bitscore}}.tax.tsv",
                ],
                bitscore=config.get("SSN_BITSCORES", [50, 80, 100])
            )
        output:
            flag = f"{SSN_DIR}/copied_to_windows.flag"
        params:
            outdir = config["windows_path"]
        shell:
            """
            mkdir -p {params.outdir}
            cp {input} {params.outdir}/
            touch {output.flag}
            """
