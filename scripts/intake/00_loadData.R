# 00_loadData.R
# Load and standardize Salmon-quantified gene count matrices (heart tissue
# and cardiac organoid) and sample metadata. Produces DESeq2-ready count
# matrices and colData data frames.
# Data was provided via Google Drive folder "Heart Deconvolution" by
# Dr. Jonathan C. Schisler on December 15th 2025 to Dr. Christoph Rau
# who then shared it with Brian Gural on Febuary 20th, 2026. Lab
# work was conducted in the lab of Dr. Brian Jensen at the University of
# North Carolina at Chapel Hill.
#
# Inputs:
#   - data/raw/heart_salmon.merged.gene_counts.tsv
#   - data/raw/organoid_salmon.merged.gene_counts.tsv
#   - data/raw/Copy of Master miRNA-Seq & RNA-Seq comparison sheet.xlsx
#   - scripts/setup/config.json
#
# Outputs:
#   - data/processed/heart_counts.rds       (integer matrix, genes x 80 samples)
#   - data/processed/heart_coldata.csv       (data.frame with HU/IR/KMP factors)
#   - data/processed/heart_gene_names.csv    (gene_id to gene_name mapping)
#   - data/processed/organoid_counts.rds     (integer matrix, genes x 54 samples)
#   - data/processed/organoid_coldata.csv    (data.frame with beam/treatment/timepoint)
#   - data/processed/organoid_gene_names.csv (gene_id to gene_name mapping)
#   - data/processed/contrast_definitions.rds (tissue + organoid contrast tables)

library(jsonlite)
library(readxl)
library(dplyr)

cfg <- fromJSON("scripts/setup/config.json")

if (!dir.exists(cfg$processed_dir)) {
  dir.create(cfg$processed_dir, recursive = TRUE)
}

# =============================================================================
# Heart tissue count matrix (FR-001, FR-003)
# =============================================================================
heart_raw <- read.delim(
  file.path(cfg$raw_dir, cfg$heart_counts_file),
  check.names = FALSE
)

heart_gene_names <- data.frame(
  gene_id   = heart_raw$gene_id,
  gene_name = heart_raw$gene_name
)

heart_counts <- as.matrix(heart_raw[, -(1:2)])
rownames(heart_counts) <- heart_raw$gene_id

# Salmon can produce fractional estimated counts
heart_counts <- round(heart_counts)
storage.mode(heart_counts) <- "integer"

# =============================================================================
# Heart colData construction (FR-005)
# =============================================================================
heart_sample_ids <- colnames(heart_counts)
heart_sample_numbers <- as.integer(gsub("_L_Heart", "", heart_sample_ids))

group_order <- cfg$heart_group_order
n_per_group <- cfg$heart_samples_per_group

heart_coldata <- data.frame(
  sample_id     = heart_sample_ids,
  sample_number = heart_sample_numbers
)

#! SAMPLE GROUPS WERE INFERRED BASED ON PRIOR PCA PLOTS FOUND
# IN THE GOOGLE DRIVE PROVIDED BY JONATHAN
heart_coldata$group <- NA_character_
for (i in seq_along(group_order)) {
  lo <- (i - 1) * n_per_group + 1
  hi <- i * n_per_group
  idx <- heart_coldata$sample_number >= lo & heart_coldata$sample_number <= hi
  heart_coldata$group[idx] <- group_order[i]
}

# Group format is HU_IR_KMP — split into factorial design columns
parts <- strsplit(heart_coldata$group, "_")
heart_coldata$HU  <- factor(sapply(parts, `[`, 1), levels = c("NL", "HU"))
heart_coldata$IR  <- factor(sapply(parts, `[`, 2), levels = c("Sham", "IR"))
heart_coldata$KMP <- factor(sapply(parts, `[`, 3), levels = c("Ctrl", "KMP"))
heart_coldata$group <- factor(heart_coldata$group, levels = group_order)

rownames(heart_coldata) <- heart_coldata$sample_id

stopifnot("Heart: colData rows must match count matrix columns" =
            all(colnames(heart_counts) == rownames(heart_coldata)))
stopifnot("Heart: expected 80 samples" = ncol(heart_counts) == 80)

# =============================================================================
# Organoid count matrix (FR-002, FR-003)
# =============================================================================
organoid_raw <- read.delim(
  file.path(cfg$raw_dir, cfg$organoid_counts_file),
  check.names = FALSE
)

organoid_gene_names <- data.frame(
  gene_id   = organoid_raw$gene_id,
  gene_name = organoid_raw$gene_name
)

organoid_counts <- as.matrix(organoid_raw[, -(1:2)])
rownames(organoid_counts) <- organoid_raw$gene_id

# Salmon can produce fractional estimated counts
organoid_counts <- round(organoid_counts)
storage.mode(organoid_counts) <- "integer"

# =============================================================================
# Organoid colData construction (FR-006)
# =============================================================================
organoid_sample_ids <- colnames(organoid_counts)

# Naming convention: first char = batch, second = B/N (beam), third = D/K (treatment)
# D9 anywhere in name indicates Day9, otherwise Day2
organoid_coldata <- data.frame(sample_id = organoid_sample_ids) |>
  mutate(
    batch     = as.integer(substr(sample_id, 1, 1)),
    beam      = case_match(substr(sample_id, 2, 2), "B" ~ "Beam", "N" ~ "NoBeam"),
    treatment = case_match(substr(sample_id, 3, 3), "D" ~ "KMP", "K" ~ "Control"),
    timepoint = if_else(grepl("D9", sample_id), "Day9", "Day2")
  )

stopifnot("Unexpected beam codes in organoid sample names" =
            !any(is.na(organoid_coldata$beam)))
stopifnot("Unexpected treatment codes in organoid sample names" =
            !any(is.na(organoid_coldata$treatment)))

organoid_coldata <- organoid_coldata |>
  mutate(
    group     = paste(beam, treatment, timepoint, sep = "_"),
    batch     = factor(batch),
    beam      = factor(beam, levels = c("NoBeam", "Beam")),
    treatment = factor(treatment, levels = c("Control", "KMP")),
    timepoint = factor(timepoint, levels = c("Day2", "Day9")),
    group     = factor(group)
  ) |>
  as.data.frame()

rownames(organoid_coldata) <- organoid_coldata$sample_id

stopifnot("Organoid: colData rows must match count matrix columns" =
            all(colnames(organoid_counts) == rownames(organoid_coldata)))
stopifnot("Organoid: expected 54 samples" = ncol(organoid_counts) == 54)
stopifnot("Organoid: no NA in factor columns" =
            !any(is.na(organoid_coldata$beam)) &&
            !any(is.na(organoid_coldata$treatment)) &&
            !any(is.na(organoid_coldata$timepoint)))

# =============================================================================
# Excel metadata parsing (FR-004)
# =============================================================================
meta_path <- file.path(cfg$raw_dir, cfg$metadata_file)

rna_sheet <- read_excel(meta_path, sheet = cfg$metadata_sheet, col_names = FALSE)
rna_mat <- as.data.frame(rna_sheet)

# Tissue contrast table (rows ~34-58 based on inspection)
tissue_header_idx <- which(rna_mat[, 1] == "Tissue type" & rna_mat[, 2] == "Comparison Name (in filenames)")
if (length(tissue_header_idx) > 0) {
  tissue_header_idx <- tissue_header_idx[1]
  tissue_rows <- (tissue_header_idx + 1):nrow(rna_mat)

  tissue_contrasts <- data.frame(
    tissue_type     = character(),
    comparison_name = character(),
    group1          = character(),
    group2          = character(),
    biological_question = character()
  )

  for (r in tissue_rows) {
    comp_name <- rna_mat[r, 2]
    if (is.na(comp_name) || comp_name == "") next
    # Stop at the organoid section boundary
    if (!is.na(rna_mat[r, 1]) && rna_mat[r, 1] == "Organoid type") break
    if (is.na(rna_mat[r, 3]) || rna_mat[r, 3] == "") next

    tissue_contrasts <- rbind(tissue_contrasts, data.frame(
      tissue_type     = ifelse(is.na(rna_mat[r, 1]), "All", rna_mat[r, 1]),
      comparison_name = comp_name,
      group1          = rna_mat[r, 3],
      group2          = rna_mat[r, 4],
      biological_question = ifelse(is.na(rna_mat[r, 7]), "", rna_mat[r, 7])
    ))
  }
} else {
  tissue_contrasts <- data.frame()
  warning("Could not find tissue contrast header in RNA_update sheet")
}

# Organoid contrast table (rows ~61-127)
organoid_header_idx <- which(rna_mat[, 1] == "Organoid type" & rna_mat[, 2] == "Comparison Name (in filenames)")
if (length(organoid_header_idx) > 0) {
  organoid_header_idx <- organoid_header_idx[1]
  organoid_rows <- (organoid_header_idx + 1):nrow(rna_mat)

  organoid_contrasts <- data.frame(
    comparison_name = character(),
    group1          = character(),
    group2          = character(),
    biological_question = character()
  )

  for (r in organoid_rows) {
    comp_name <- rna_mat[r, 2]
    if (is.na(comp_name) || comp_name == "") next
    if (is.na(rna_mat[r, 3]) || rna_mat[r, 3] == "") next

    organoid_contrasts <- rbind(organoid_contrasts, data.frame(
      comparison_name = comp_name,
      group1          = rna_mat[r, 3],
      group2          = rna_mat[r, 4],
      biological_question = ifelse(is.na(rna_mat[r, 7]), "", rna_mat[r, 7])
    ))
  }
} else {
  organoid_contrasts <- data.frame()
  warning("Could not find organoid contrast header in RNA_update sheet")
}

contrast_definitions <- list(
  tissue   = tissue_contrasts,
  organoid = organoid_contrasts
)

# =============================================================================
# Save outputs (FR-007, FR-008, FR-009)
# =============================================================================
saveRDS(heart_counts, file.path(cfg$processed_dir, "heart_counts.rds"))
write.csv(heart_coldata, file.path(cfg$processed_dir, "heart_coldata.csv"),
          row.names = FALSE)
write.csv(heart_gene_names, file.path(cfg$processed_dir, "heart_gene_names.csv"),
          row.names = FALSE)
saveRDS(organoid_counts, file.path(cfg$processed_dir, "organoid_counts.rds"))
write.csv(organoid_coldata, file.path(cfg$processed_dir, "organoid_coldata.csv"),
          row.names = FALSE)
write.csv(organoid_gene_names, file.path(cfg$processed_dir, "organoid_gene_names.csv"),
          row.names = FALSE)
saveRDS(contrast_definitions, file.path(cfg$processed_dir, "contrast_definitions.rds"))
