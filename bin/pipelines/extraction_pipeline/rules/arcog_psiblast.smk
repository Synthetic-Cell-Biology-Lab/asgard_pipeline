
from pathlib import Path
import os


PROTEIN = config['run']['protein_name']
RESULT_DIR = Path(config['run']['base_dir'])/"database"/"protein_sets"/ PROTEIN
RESULT_DIR.mkdir(exist_ok=True)
PSIBLAST = RESULT_DIR / "PSIBLAST"
PSIBLAST.mkdir(exist_ok=True )
ENV_DIR = config['run']['env_dir']
LOG_DIR = RESULT_DIR / "LOGS"

rule all:
    input:
        f"{PSIBLAST}/arcog/candidate_{PROTEIN}_asgard_hits.tsv"


rule create_arcog_id_set:
    input:
        definition_tsv = config['inputs']['arcog']['def_tab']
    output:
        arcog_ids = f"{PSIBLAST}/arcog/arcog_ids_{PROTEIN}.csv"
    params:
        search_terms = config['modes']['arcog']['search_terms'],
        ids = config['modes']['arcog']['ids'],
    conda:
        f"{ENV_DIR}/Reg.yaml"
    run:
        import pandas as pd

        # import the id search terms
        df = pd.read_csv(input.definition_tsv, sep="\t")
        df.columns = df.columns.str.strip()
        
        # get by search terms
        id_list = []
        for search_term in params.search_terms:
            id_list.extend(df.loc[
                df['Description'].str.contains(search_term, case=False, na=False)
            ])
        
        # get by the ids
        id_list.extend([f"arCOG{int(id_):05d}" for id_ in params.ids])

        # filter and export
        out_txt_df = df[df['arCOG'].isin(id_list)]

        out_txt_df.to_csv(output.arcog_ids)


rule retrieve_proteins_from_ids:
    input:
        arcog_ids = f"{PSIBLAST}/arcog/arcog_ids_{PROTEIN}.csv",
        fasta = config['inputs']['arcog']['fasta'],
        arcog_csv = config['inputs']['arcog']['csv']
    output:
        fasta_out = f"{PSIBLAST}/arcog/candidate_{PROTEIN}_arcog_proteins.faa",
        candidate_ids = f"{PSIBLAST}/arcog/candidate_{PROTEIN}_arcog_proteins.ids",
    conda:
        f"{ENV_DIR}/Reg.yaml"
    run:
        from Bio import SeqIO
        import pandas as pd

        arcsv = pd.read_csv(
            input.arcog_csv,
            names=[
                "locus_tag", "genome_id", "protein_id",
                "metric1", "metric2", "metric3",
                "arcog_id", "metric4", "metric5"
            ]
        )

        # Read desired arCOG IDs
        arcog_ids = set(pd.read_csv(input.arcog_ids)["arCOG"])

        # Subset the mapping table
        arcsv_subset = arcsv[arcsv["arcog_id"].isin(arcog_ids)]

        # Get corresponding protein IDs
        protein_ids = set(arcsv_subset["protein_id"])

        with open(output.candidate_ids, "w") as out:
            out.write("\n".join(sorted(protein_ids)))

        # Filter FASTA
        with open(output.fasta_out, "w") as out:
            for record in SeqIO.parse(input.fasta, "fasta"):
                if record.id in protein_ids:
                    SeqIO.write(record, out, "fasta")



# change to our internal database!!
rule make_blastdb:
    input:
        fasta_db = config['inputs']['fasta_file'],
        
    output:
        multiext(
            config['inputs']['fasta_file'].strip(".fasta"),
            ".phr",
            ".pin",
            ".psq"
        )
    params:
        db_path = config['inputs']['fasta_file'].strip(".fasta")
    conda:
        f"{ENV_DIR}/ssn.yaml"

    shell:
        """
        makeblastdb \
            -dbtype prot \
            -in {input} \
            -out {params.db_path}
        """

# database is our internal database and query is the arcog set of proteins
rule psiblast:
    input:
        query=f"{PSIBLAST}/arcog/candidate_{PROTEIN}_arcog_proteins.faa",
        db = ancient(multiext(
            config['inputs']['fasta_file'].strip(".fasta"),
            ".phr",
            ".pin",
            ".psq"
        ))

    output:
        f"{PSIBLAST}/arcog/candidate_{PROTEIN}_asgard_hits.tsv"
    threads: 8

    conda:
        f"{ENV_DIR}/ssn.yaml"


    log:
        f"{LOG_DIR}/blast/{PROTEIN}.psiblast.log"

    params:
        db_path = config['inputs']['fasta_file'].strip(".fasta")

    shell:
        """
        psiblast \
            -query {input.query} \
            -db {params.db_path}\
            -evalue 0.01 \
            -num_iterations 5 \
            -dbsize 20000000 \
            -comp_based_stats F \
            -seg no \
            -outfmt 6 \
            -num_threads {threads} \
            -out {output} \
            2> {log}
        """