import os

configfile: "scripts/setup/config.json"

DECON_REFS = config["decon_references"]

rule all:
    input:
        expand("data/processed/{dataset}_counts.rds",
               dataset=["heart", "organoid"]),
        expand("data/processed/{dataset}_{suffix}.csv",
               dataset=["heart", "organoid"],
               suffix=["coldata", "gene_names"]),
        "data/processed/contrast_definitions.rds",
        config["qc_report_output"],
        expand("data/processed/decon/{ref}_compositions.csv",
               ref=DECON_REFS),
        expand("data/processed/models/{ref}_de_comparison.csv",
               ref=DECON_REFS),
        expand("data/processed/models/{ref}_dirichlet_coefficients.csv",
               ref=DECON_REFS),
        expand("results/figures/{ref}_dirichlet_coeff.png",
               ref=DECON_REFS),
        config["decon_report_output"]

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
        "conda run --live-stream -n orbital_mus "
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
        "conda run --live-stream -n orbital_mus "
        "Rscript -e \"rmarkdown::render('{input.rmd}', "
        "params=list(proj_root='{params.proj_root}'), "
        "output_file='{params.output_path}')\""

rule find_markers:
    input:
        gene_names = "data/processed/heart_gene_names.csv",
        sc_labeled = "data/processed/sc/celltype_labeled.rds",
        sc_igvf = "data/processed/sc/IGVFFI6644FMFS.rds"
    output:
        expand("data/processed/decon/{ref}_{suffix}",
               ref=DECON_REFS,
               suffix=["counts.rds", "metadata.csv", "markers.csv"])
    resources:
        mem_mb = 24000
    shell:
        "conda run --live-stream -n orbital_mus_decon "
        "Rscript scripts/decon/02_find_markers.R"

rule deconvolute:
    input:
        sc_counts = "data/processed/decon/{ref}_counts.rds",
        sc_meta = "data/processed/decon/{ref}_metadata.csv",
        markers = "data/processed/decon/{ref}_markers.csv",
        counts = "data/processed/heart_counts.rds",
        coldata = "data/processed/heart_coldata.csv",
        gene_names = "data/processed/heart_gene_names.csv"
    output:
        compositions = "data/processed/decon/{ref}_compositions.csv",
        weights = "data/processed/decon/{ref}_gene_weights.csv"
    resources:
        mem_mb = 16000
    shell:
        "conda run --live-stream -n orbital_mus_music "
        "Rscript scripts/decon/03_deconvolute.R {wildcards.ref}"

rule differential_expression:
    input:
        counts = "data/processed/heart_counts.rds",
        coldata = "data/processed/heart_coldata.csv",
        gene_names = "data/processed/heart_gene_names.csv",
        contrasts = "data/processed/contrast_definitions.rds",
        compositions = "data/processed/decon/{ref}_compositions.csv"
    output:
        dds_adj = "data/processed/models/{ref}_dds_adjusted.rds",
        dds_unadj = "data/processed/models/{ref}_dds_unadjusted.rds",
        res_adj = "data/processed/models/{ref}_results_adjusted.rds",
        res_unadj = "data/processed/models/{ref}_results_unadjusted.rds",
        comparison = "data/processed/models/{ref}_de_comparison.csv"
    resources:
        mem_mb = 16000
    shell:
        "conda run --live-stream -n orbital_mus "
        "Rscript scripts/decon/04_runDe.R {wildcards.ref}"

rule dirichlet:
    input:
        compositions = "data/processed/decon/{ref}_compositions.csv",
        coldata = "data/processed/heart_coldata.csv",
        config_file = "scripts/setup/config.json"
    output:
        coefficients = "data/processed/models/{ref}_dirichlet_coefficients.csv",
        plot = "results/figures/{ref}_dirichlet_coeff.png"
    resources:
        mem_mb = 4000
    shell:
        "conda run --live-stream -n orbital_mus "
        "Rscript scripts/decon/06_dirichlet.R {wildcards.ref}"

rule decon_summary:
    input:
        compositions = expand("data/processed/decon/{ref}_compositions.csv",
                              ref=DECON_REFS),
        comparisons = expand("data/processed/models/{ref}_de_comparison.csv",
                             ref=DECON_REFS),
        config_file = "scripts/setup/config.json",
        rmd = "scripts/decon/05_summarize.Rmd"
    output:
        config["decon_report_output"]
    resources:
        mem_mb = 8000
    params:
        output_path = lambda wildcards, output: os.path.abspath(str(output)),
        proj_root = os.getcwd()
    shell:
        "conda run --live-stream -n orbital_mus "
        "Rscript -e \"rmarkdown::render('{input.rmd}', "
        "params=list(proj_root='{params.proj_root}'), "
        "output_file='{params.output_path}')\""
