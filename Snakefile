import os

configfile: "scripts/setup/config.json"

rule all:
    input:
        expand("data/processed/{dataset}_counts.rds",
               dataset=["heart", "organoid"]),
        expand("data/processed/{dataset}_{suffix}.csv",
               dataset=["heart", "organoid"],
               suffix=["coldata", "gene_names"]),
        "data/processed/contrast_definitions.rds",
        config["qc_report_output"]

rule load_data:
    input:
        heart = "data/raw/" + config["heart_counts_file"],
        organoid = "data/raw/" + config["organoid_counts_file"],
        metadata = "data/raw/" + config["metadata_file"],
        config_file = "scripts/setup/config.json"
    output:
        expand("data/processed/{dataset}_counts.rds",
               dataset=["heart", "organoid"]),
        expand("data/processed/{dataset}_{suffix}.csv",
               dataset=["heart", "organoid"],
               suffix=["coldata", "gene_names"]),
        "data/processed/contrast_definitions.rds"
    resources:
        mem_mb = 4000
    shell:
        "Rscript scripts/intake/00_loadData.R"

rule qc_report:
    input:
        expand("data/processed/{dataset}_counts.rds",
               dataset=["heart", "organoid"]),
        expand("data/processed/{dataset}_coldata.csv",
               dataset=["heart", "organoid"]),
        config_file = "scripts/setup/config.json",
        rmd = "scripts/intake/01_runQc.Rmd"
    output:
        config["qc_report_output"]
    resources:
        mem_mb = 8000
    params:
        output_path = lambda wildcards, output: os.path.abspath(str(output)),
        proj_root = os.getcwd()
    shell:
        "Rscript -e \"rmarkdown::render('{input.rmd}', "
        "params=list(proj_root='{params.proj_root}'), "
        "output_file='{params.output_path}')\""
