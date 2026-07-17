"""
Yet to be finished
under preliminary work
"""
from pathlib import Path
import pandas as pd


OUTDIR = config["output_dir"]
ORGS   = [org["name"].replace(" ", "_").lower() for org in config["organisms"]]
BAKTA_OUTDIR = Path(config.get("cds", {}).get("existing"))



def get_filtered_genomes():
    genomes = []
    for org in ORGS:
        csv_path = f"{OUTDIR}/gtdbtk/{org}_gtdbtk_classification_split.ar53.csv"
        if not os.path.exists(csv_path):
            continue
        df = pd.read_csv(csv_path)
        phyla = config['interproscan']['phylum_filter']
        if isinstance(phyla, str):
            phyla = [phyla]
        genomes.extend(df[df['phylum'].isin(phyla)]['user_genome'].tolist())
    return genomes

FILTERED_GENOMES = get_filtered_genomes()



rule all:
    input:
        f"{OUTDIR}/orthofinder/OrthoFinder.done"

rule link_proteomes:
    input:
        expand(
            f"{BAKTA_OUTDIR}/{{genome}}/{{genome}}.faa",
            genome=FILTERED_GENOMES
        )
    output:
        linked = touch(f"{OUTDIR}/orthofinder/input/.linked")
    params:
        path = f"{OUTDIR}/orthofinder/input"
    run:
        from pathlib import Path
        outdir = Path(params.path)
        outdir.mkdir(parents=True, exist_ok=True)
        for faa in input:
            target = outdir / Path(faa).name
            if target.exists() or target.is_symlink():
                target.unlink()
            os.symlink(os.path.abspath(faa), target)


rule orthofinder:
    input:
        f"{OUTDIR}/orthofinder/input/.linked"
    output:
        touch(f"{OUTDIR}/orthofinder/OrthoFinder.done")
    params:
        indir  = f"{OUTDIR}/orthofinder/input",
        outdir = f"{OUTDIR}/orthofinder/results"
    threads: 8
    conda:
        f"{config.get('env_dir')}/bakta.yaml"
    shell:
        """
        orthofinder \
            -f {params.indir} \
            -o {params.outdir} \
            -t {threads} \
            -a {threads} 
        touch {output}
        """