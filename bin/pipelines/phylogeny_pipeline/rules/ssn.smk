rule ssn_network:
    input:
        fasta = f"{EXPLORATION_DIR}/{PROTEIN}.unr.fasta",
        csv = f"{EXPLORATION_DIR}/{PROTEIN}.unr.csv",

    output:
        nr          = f"{SSN_DIR}/{PROTEIN}.nr.fasta",
        similarities = f"{SSN_DIR}/{PROTEIN}.similarities.tsv",
        edges       = f"{SSN_DIR}/{PROTEIN}.edges.tsv",
        nodes       = f"{SSN_DIR}/{PROTEIN}.nodes.tsv",
        sif         = f"{SSN_DIR}/{PROTEIN}.edges.sif",
        ea          = f"{SSN_DIR}/{PROTEIN}.edges.ea",
        taxonomy    = f"{SSN_DIR}/{PROTEIN}.tax.tsv"

    params:
        identity = config.get("SSN_IDENTITY", 0.90),
        bitscore = config.get("SSN_BITSCORES", 50),
        coverage = config.get("SSN_COVERAGE", 0.6),
        evalue   = config.get("SSN_EVALUE",   1e-5),

    threads: config.get("SSN_CORES", config.get("cores", 16))

    resources:
        mem_mb = config.get("SSN_MEM_MB", 16000),
        runtime = config.get("SSN_RUNTIME_MIN", 720)

    log:
        f"{LOG_DIR}/ssn/{PROTEIN}.log"

    benchmark:
        f"{BENCHMARK_DIR}/ssn/{PROTEIN}.benchmark.tsv"

    conda:
        f"{config['env_dir']}/ssn.yaml"

    shell:
        """
        bash {CURRENT_DIR}/bin/units/run_ssn.sh \
            {input.fasta}        \
            {output.nr}          \
            {output.similarities}\
            {output.edges}       \
            {output.nodes}       \
            {threads}            \
            {params.identity}    \
            {params.bitscore}    \
            {params.coverage}    \
            {params.evalue}      \
            {input.csv}          \
            {output.taxonomy}    \
            2>&1 | tee {log}
        """



if config.get("copy_to_windows", False):

    rule copy_outs_to_windows:
        input:
            edges = f"{SSN_DIR}/{PROTEIN}.edges.ea",
            nodes = f"{SSN_DIR}/{PROTEIN}.nodes.tsv",
            sif   = f"{SSN_DIR}/{PROTEIN}.edges.sif",
            taxonomy    = f"{SSN_DIR}/{PROTEIN}.tax.tsv"

        output:
            flag = f"{SSN_DIR}/copied_to_windows.flag"
        params:
            outdir = config["windows_path"]
        shell:
            """
            mkdir -p {params.outdir}

            cp {input.edges} {params.outdir}/
            cp {input.nodes} {params.outdir}/
            cp {input.sif}   {params.outdir}/
            cp {input.taxonomy} {params.outdir}/

            touch {output.flag}
            """