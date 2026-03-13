rule ssn_network:
    input:
        fasta = f"{RESULT_DIR}/{PROTEIN}.rev.fasta"

    output:
        nr = f"{SSN_DIR}/{PROTEIN}.nr.fasta",
        similarities=f"{SSN_DIR}/{PROTEIN}.similarities.tsv",
        edges=f"{SSN_DIR}/{PROTEIN}.edges.tsv",
        nodes=f"{SSN_DIR}/{PROTEIN}.nodes.tsv"

    params:
        identity=config.get("SSN_IDENTITY", 0.98),
        bitscore=config.get("SSN_BITSCORE", 50),
        coverage=config.get("SSN_COVERAGE", 0.6)

    threads: config.get('cores', 16)

    conda:
        f"{config['env_dir']}/ssn.yaml"
    shell:
        """
        bash {CURRENT_DIR}/bin/units/run_ssn.sh \
        {input.fasta} \
        {output.nr} \
        {output.similarities} \
        {output.edges} \
        {output.nodes} \
        {threads} \
        {params.identity} \
        {params.bitscore} \
        {params.coverage}
        """