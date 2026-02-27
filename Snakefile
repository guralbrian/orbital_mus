import os
import re

configfile: "scripts/setup/config.json"

DECON_REFS = config["decon_references"]
MODEL_TYPES = ["adjusted", "unadjusted"]

# Constrain wildcards to prevent ambiguous rule matching
wildcard_constraints:
    ref = "|".join([re.escape(r) for r in DECON_REFS]),
    model = "adjusted|unadjusted"

# All 19 contrast slugs (filesystem-safe versions of contrast names).
# Deterministic from 04_runDe.R: pairwise contrasts from metadata +
# hard-coded main effects/marginal means/interactions.
# The _contrast_index.csv maps slugs back to original names.
ALL_CONTRAST_SLUGS = [
    # Pairwise (from Excel metadata)
    "HU_vs_Ctrl",
    "IR_vs_Ctrl",
    "KMP_vs_Ctrl",
    "HU_KMP_vs_Ctrl",
    "IR_KMP_vs_Ctrl",
    "HU_IR_vs_Ctrl",
    "HU_IR_KMP_vs_Ctrl",
    "KMP_effect_in_HU",
    "KMP_effect_in_IR",
    "KMP_effect_in_HU_IR",
    # Main effects
    "HU_vs_NL_MAIN_EFFECT",
    "IR_vs_Sham_MAIN_EFFECT",
    "KMP_vs_Sham_MAIN_EFFECT",
    # Marginal means
    "IR_KMP__HU_NL__vs_Ctrl__HU_NL_",
    "HU_KMP__IR_Sham__vs_Ctrl__IR_Sham_",
    "HU_IR__KMP_Ctrl__vs_Ctrl__KMP_Ctrl_",
    # Two-way interactions
    "Interaction_KMP_x_HU",
    "Interaction_KMP_x_IR",
    "Interaction_HU_x_IR",
]

# Key contrasts: 3 main effects + 3 two-way interactions
# Used for focused plots (UpSet, GSEA dots, adj comparison)
KEY_CONTRAST_SLUGS = [
    "HU_vs_NL_MAIN_EFFECT",
    "IR_vs_Sham_MAIN_EFFECT",
    "KMP_vs_Sham_MAIN_EFFECT",
    "Interaction_KMP_x_HU",
    "Interaction_KMP_x_IR",
    "Interaction_HU_x_IR",
]

# --- GSEA scope (subset comments kept for quick narrowing) ---
GSEA_REFS = DECON_REFS               # subset: ["bulk_decon"]
GSEA_MODELS = MODEL_TYPES            # subset: ["adjusted"]
GSEA_SLUGS = ALL_CONTRAST_SLUGS      # subset: KEY_CONTRAST_SLUGS

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
        expand("data/processed/models/{ref}_clr_covariates.csv",
               ref=DECON_REFS),
        expand("data/processed/models/{ref}_contrast_index.csv",
               ref=DECON_REFS),
        expand("data/processed/models/{ref}_dirichlet_coefficients.csv",
               ref=DECON_REFS),
        expand("results/figures/{ref}_dirichlet_coeff.png",
               ref=DECON_REFS),
        expand("data/processed/gsea/{ref}_{model}_combined.csv",
               ref=GSEA_REFS, model=GSEA_MODELS),
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
        comparison = "data/processed/models/{ref}_de_comparison.csv",
        clr_covs = "data/processed/models/{ref}_clr_covariates.csv",
        contrast_idx = "data/processed/models/{ref}_contrast_index.csv"
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

rule gsea:
    input:
        results = "data/processed/models/{ref}_results_{model}.rds",
        contrast_idx = "data/processed/models/{ref}_contrast_index.csv"
    output:
        go = "data/processed/gsea/{ref}_{model}_{cslug}_go.rds",
        kegg = "data/processed/gsea/{ref}_{model}_{cslug}_kegg.rds",
        summary = "data/processed/gsea/{ref}_{model}_{cslug}_summary.csv"
    resources:
        mem_mb = 8000,
        runtime = 30
    shell:
        "conda run --live-stream -n orbital_mus_gsea "
        "Rscript scripts/decon/07_gsea.R {wildcards.ref} {wildcards.model} {wildcards.cslug}"

rule collect_gsea:
    input:
        summaries = expand("data/processed/gsea/{{ref}}_{{model}}_{cslug}_summary.csv",
                           cslug=GSEA_SLUGS)
    output:
        "data/processed/gsea/{ref}_{model}_combined.csv"
    resources:
        mem_mb = 2000
    shell:
        "conda run --live-stream -n orbital_mus "
        "Rscript scripts/decon/08_collect_gsea.R {wildcards.ref} {wildcards.model}"

rule decon_summary:
    input:
        compositions = expand("data/processed/decon/{ref}_compositions.csv",
                              ref=DECON_REFS),
        comparisons = expand("data/processed/models/{ref}_de_comparison.csv",
                             ref=DECON_REFS),
        gsea_summaries = expand("data/processed/gsea/{ref}_{model}_combined.csv",
                                ref=GSEA_REFS, model=GSEA_MODELS),
        results_adj = expand("data/processed/models/{ref}_results_adjusted.rds",
                             ref=DECON_REFS),
        results_unadj = expand("data/processed/models/{ref}_results_unadjusted.rds",
                               ref=DECON_REFS),
        clr_covs = expand("data/processed/models/{ref}_clr_covariates.csv",
                          ref=DECON_REFS),
        dirichlet = expand("data/processed/models/{ref}_dirichlet_coefficients.csv",
                           ref=DECON_REFS),
        markers = expand("data/processed/decon/{ref}_markers.csv",
                         ref=DECON_REFS),
        contrast_idx = expand("data/processed/models/{ref}_contrast_index.csv",
                              ref=DECON_REFS),
        dds_adj = expand("data/processed/models/{ref}_dds_adjusted.rds",
                         ref=DECON_REFS),
        config_file = "scripts/setup/config.json",
        rmd = "scripts/decon/05_summarize.Rmd"
    output:
        config["decon_report_output"]
    resources:
        mem_mb = 16000
    params:
        output_path = lambda wildcards, output: os.path.abspath(str(output)),
        proj_root = os.getcwd()
    shell:
        "conda run --live-stream -n orbital_mus "
        "Rscript -e \"rmarkdown::render('{input.rmd}', "
        "params=list(proj_root='{params.proj_root}'), "
        "output_file='{params.output_path}')\""
