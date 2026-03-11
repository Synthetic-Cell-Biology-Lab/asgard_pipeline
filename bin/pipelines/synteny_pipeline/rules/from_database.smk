rule get_data_from_database:
    input:
        protein_file = config["protein_file"]
    output:
        csv = f"{SYNTENY_DIR}/{RUN_ID}.synteny.csv"
    conda:
        f"{config['env_dir']}/Reg.yaml"
    params:
        proteins = ",".join(config.get("gene_neighborhood_protein_sets", ['FtsZ1','FtsZ2','CetZ','Tubulin'])),
        protein_column = config.get("protein_column", "Manual_annotation")
    message:
        """
        ==========================================
          Subsetting Data from main Protein file
        ==========================================
        """
    shell:
        """
        python {CURRENT_DIR}/bin/units/subset_dataframe.py \
            {input.protein_file} \
            {params.proteins} \
            {params.protein_column} \
            {output.csv}
        """