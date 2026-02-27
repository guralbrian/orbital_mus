# 08_collect_gsea.R
# Aggregate per-contrast GSEA summary CSVs into one file per ref+model.
#
# Usage: Rscript scripts/decon/08_collect_gsea.R <ref_name> <model_type>
#
# Inputs:
#   - data/processed/gsea/{ref}_{model}_*_summary.csv (all matching files)
#
# Outputs:
#   - data/processed/gsea/{ref}_{model}_gsea_summary.csv

library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
ref_name   <- args[1]
model_type <- args[2]

stopifnot(!is.null(ref_name), !is.null(model_type))

gsea_dir <- file.path("data", "processed", "gsea")

# Discover all per-contrast summary CSVs for this ref + model
pattern <- paste0("^", ref_name, "_", model_type, "_.*_summary\\.csv$")
csv_files <- list.files(gsea_dir, pattern = pattern, full.names = TRUE)

# Exclude the aggregated file itself if it already exists
agg_name <- paste0(ref_name, "_", model_type, "_combined.csv")
csv_files <- csv_files[basename(csv_files) != agg_name]

stopifnot("No per-contrast GSEA summary CSVs found" = length(csv_files) > 0)

combined <- lapply(csv_files, read.csv) |> bind_rows()

write.csv(combined,
          file.path(gsea_dir, agg_name),
          row.names = FALSE)
