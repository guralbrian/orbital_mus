# 00_loadData.R
# Load and standardize Salmon-quantified gene count matrices (heart tissue
# and cardiac organoid) and sample metadata. Produces DESeq2-ready count
# matrices, colData data frames, and gene name mappings.
#
# Inputs:
#   - data/raw/heart_salmon.merged.gene_counts.tsv
#   - data/raw/organoid_salmon.merged.gene_counts.tsv
#   - data/raw/Copy of Master miRNA-Seq & RNA-Seq comparison sheet.xlsx
#   - scripts/setup/config.json
#
# Outputs:
#   - data/processed/heart_counts.rds       (integer matrix, genes x 80 samples)
#   - data/processed/heart_coldata.rds       (data.frame with HU/IR/KMP factors)
#   - data/processed/heart_gene_names.rds    (gene_id to gene_name mapping)
#   - data/processed/organoid_counts.rds     (integer matrix, genes x 54 samples)
#   - data/processed/organoid_coldata.rds    (data.frame with beam/treatment/timepoint)
#   - data/processed/organoid_gene_names.rds (gene_id to gene_name mapping)
#   - data/processed/contrast_definitions.rds (tissue + organoid contrast tables)

# Load Libraries
libs <- c("jsonlite", "readxl", "dplyr")
lapply(libs, require, character.only = TRUE)
rm(libs)

# Load configuration
cfg <- fromJSON("scripts/setup/config.json")

# Ensure output directory exists
if (!dir.exists(cfg$processed_dir)) {
  dir.create(cfg$processed_dir, recursive = TRUE)
}

# =============================================================================
# Heart tissue count matrix (FR-001, FR-003)
# =============================================================================
message("Loading heart count matrix...")
heart_raw <- read.delim(
  file.path(cfg$raw_dir, cfg$heart_counts_file),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Extract gene name mapping
heart_gene_names <- data.frame(
  gene_id   = heart_raw$gene_id,
  gene_name = heart_raw$gene_name,
  stringsAsFactors = FALSE
)

# Build count matrix: genes x samples
heart_counts <- as.matrix(heart_raw[, -(1:2)])
rownames(heart_counts) <- heart_raw$gene_id

# Round to integer if needed (Salmon can produce fractional estimates)
if (any(heart_counts != round(heart_counts))) {
  message("  Rounding non-integer heart counts to nearest integer")
}
heart_counts <- round(heart_counts)
storage.mode(heart_counts) <- "integer"

message(sprintf("  Heart matrix: %d genes x %d samples", nrow(heart_counts), ncol(heart_counts)))

# =============================================================================
# Heart colData construction (FR-005)
# =============================================================================
message("Constructing heart colData...")

heart_sample_ids <- colnames(heart_counts)
heart_sample_numbers <- as.integer(gsub("_L_Heart", "", heart_sample_ids))

# Build group assignment from sequential blocks of 10
group_order <- cfg$heart_group_order
n_per_group <- cfg$heart_samples_per_group

heart_coldata <- data.frame(
  sample_id     = heart_sample_ids,
  sample_number = heart_sample_numbers,
  stringsAsFactors = FALSE
)

# Assign group based on sample number
heart_coldata$group <- NA_character_
for (i in seq_along(group_order)) {
  lo <- (i - 1) * n_per_group + 1
  hi <- i * n_per_group
  idx <- heart_coldata$sample_number >= lo & heart_coldata$sample_number <= hi
  heart_coldata$group[idx] <- group_order[i]
}

# Parse HU, IR, KMP from group name (format: HU_IR_KMP)
parts <- strsplit(heart_coldata$group, "_")
heart_coldata$HU  <- factor(sapply(parts, `[`, 1), levels = c("NL", "HU"))
heart_coldata$IR  <- factor(sapply(parts, `[`, 2), levels = c("Sham", "IR"))
heart_coldata$KMP <- factor(sapply(parts, `[`, 3), levels = c("Ctrl", "KMP"))
heart_coldata$group <- factor(heart_coldata$group, levels = group_order)

rownames(heart_coldata) <- heart_coldata$sample_id

# Validate heart data
stopifnot("Heart: colData rows must match count matrix columns" =
            all(colnames(heart_counts) == rownames(heart_coldata)))
stopifnot("Heart: expected 80 samples" = ncol(heart_counts) == 80)
stopifnot("Heart: no NA in group assignment" = !any(is.na(heart_coldata$group)))
stopifnot("Heart: 10 samples per group" =
            all(table(heart_coldata$group) == n_per_group))

# =============================================================================
# Organoid count matrix (FR-002, FR-003)
# =============================================================================
message("Loading organoid count matrix...")
organoid_raw <- read.delim(
  file.path(cfg$raw_dir, cfg$organoid_counts_file),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Extract gene name mapping
organoid_gene_names <- data.frame(
  gene_id   = organoid_raw$gene_id,
  gene_name = organoid_raw$gene_name,
  stringsAsFactors = FALSE
)

# Build count matrix
organoid_counts <- as.matrix(organoid_raw[, -(1:2)])
rownames(organoid_counts) <- organoid_raw$gene_id

# Round to integer if needed
if (any(organoid_counts != round(organoid_counts))) {
  message("  Rounding non-integer organoid counts to nearest integer")
}
organoid_counts <- round(organoid_counts)
storage.mode(organoid_counts) <- "integer"

message(sprintf("  Organoid matrix: %d genes x %d samples", nrow(organoid_counts), ncol(organoid_counts)))

# =============================================================================
# Organoid colData construction (FR-006)
# =============================================================================
message("Constructing organoid colData...")

organoid_sample_ids <- colnames(organoid_counts)

# Parse sample names: first char = batch, second = B/N (beam), third = D/K (treatment)
# D9 anywhere in name indicates Day9, otherwise Day2
parse_organoid_name <- function(name) {
  batch_num  <- as.integer(substr(name, 1, 1))
  beam_code  <- substr(name, 2, 2)
  treat_code <- substr(name, 3, 3)

  if (!(beam_code %in% c("B", "N"))) {
    warning(sprintf("Unexpected beam code '%s' in sample: %s", beam_code, name))
    return(NULL)
  }
  if (!(treat_code %in% c("D", "K"))) {
    warning(sprintf("Unexpected treatment code '%s' in sample: %s", treat_code, name))
    return(NULL)
  }

  data.frame(
    sample_id = name,
    batch     = batch_num,
    beam      = ifelse(beam_code == "B", "Beam", "NoBeam"),
    treatment = ifelse(treat_code == "D", "KMP", "Control"),
    timepoint = ifelse(grepl("D9", name), "Day9", "Day2"),
    stringsAsFactors = FALSE
  )
}

organoid_coldata <- do.call(rbind, lapply(organoid_sample_ids, parse_organoid_name))

# Check for parsing failures
n_parsed <- nrow(organoid_coldata)
n_expected <- length(organoid_sample_ids)
if (n_parsed != n_expected) {
  warning(sprintf("Organoid: parsed %d of %d samples. Check warnings above.", n_parsed, n_expected))
}

# Create composite group and convert to factors
organoid_coldata$group <- paste(organoid_coldata$beam, organoid_coldata$treatment,
                                organoid_coldata$timepoint, sep = "_")
organoid_coldata$batch     <- factor(organoid_coldata$batch)
organoid_coldata$beam      <- factor(organoid_coldata$beam, levels = c("NoBeam", "Beam"))
organoid_coldata$treatment <- factor(organoid_coldata$treatment, levels = c("Control", "KMP"))
organoid_coldata$timepoint <- factor(organoid_coldata$timepoint, levels = c("Day2", "Day9"))
organoid_coldata$group     <- factor(organoid_coldata$group)

rownames(organoid_coldata) <- organoid_coldata$sample_id

# Validate organoid data
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
message("Parsing Excel metadata...")
meta_path <- file.path(cfg$raw_dir, cfg$metadata_file)

rna_sheet <- read_excel(meta_path, sheet = cfg$metadata_sheet, col_names = FALSE)
rna_mat <- as.data.frame(rna_sheet, stringsAsFactors = FALSE)

# Tissue contrast table (rows ~34-58 based on inspection)
# Find header row for tissue section: "Tissue type | Comparison Name..."
tissue_header_idx <- which(rna_mat[, 1] == "Tissue type" & rna_mat[, 2] == "Comparison Name (in filenames)")
if (length(tissue_header_idx) > 0) {
  tissue_header_idx <- tissue_header_idx[1]
  tissue_rows <- (tissue_header_idx + 1):nrow(rna_mat)

  # Collect rows that have a comparison name in column 2 until we hit organoid section
  tissue_contrasts <- data.frame(
    tissue_type     = character(),
    comparison_name = character(),
    group1          = character(),
    group2          = character(),
    biological_question = character(),
    stringsAsFactors = FALSE
  )

  for (r in tissue_rows) {
    comp_name <- rna_mat[r, 2]
    if (is.na(comp_name) || comp_name == "") next
    # Stop if we hit the organoid section marker
    if (!is.na(rna_mat[r, 1]) && rna_mat[r, 1] == "Organoid type") break
    # Skip section headers (no group1/group2)
    if (is.na(rna_mat[r, 3]) || rna_mat[r, 3] == "") next

    tissue_contrasts <- rbind(tissue_contrasts, data.frame(
      tissue_type     = ifelse(is.na(rna_mat[r, 1]), "All", rna_mat[r, 1]),
      comparison_name = comp_name,
      group1          = rna_mat[r, 3],
      group2          = rna_mat[r, 4],
      biological_question = ifelse(is.na(rna_mat[r, 7]), "", rna_mat[r, 7]),
      stringsAsFactors = FALSE
    ))
  }
  message(sprintf("  Parsed %d tissue contrasts from Excel", nrow(tissue_contrasts)))
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
    biological_question = character(),
    stringsAsFactors = FALSE
  )

  for (r in organoid_rows) {
    comp_name <- rna_mat[r, 2]
    if (is.na(comp_name) || comp_name == "") next
    # Skip if no group info
    if (is.na(rna_mat[r, 3]) || rna_mat[r, 3] == "") next

    organoid_contrasts <- rbind(organoid_contrasts, data.frame(
      comparison_name = comp_name,
      group1          = rna_mat[r, 3],
      group2          = rna_mat[r, 4],
      biological_question = ifelse(is.na(rna_mat[r, 7]), "", rna_mat[r, 7]),
      stringsAsFactors = FALSE
    ))
  }
  message(sprintf("  Parsed %d organoid contrasts from Excel", nrow(organoid_contrasts)))
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
message("Saving outputs to ", cfg$processed_dir, "/...")

saveRDS(heart_counts,      file.path(cfg$processed_dir, "heart_counts.rds"))
saveRDS(heart_coldata,     file.path(cfg$processed_dir, "heart_coldata.rds"))
saveRDS(heart_gene_names,  file.path(cfg$processed_dir, "heart_gene_names.rds"))
saveRDS(organoid_counts,   file.path(cfg$processed_dir, "organoid_counts.rds"))
saveRDS(organoid_coldata,  file.path(cfg$processed_dir, "organoid_coldata.rds"))
saveRDS(organoid_gene_names, file.path(cfg$processed_dir, "organoid_gene_names.rds"))
saveRDS(contrast_definitions, file.path(cfg$processed_dir, "contrast_definitions.rds"))

# Print experimental design summary
message("\n=== Experimental Design Summary ===\n")

message("--- Heart Tissue (2x2x2 factorial) ---")
message(sprintf("  Samples: %d genes x %d samples", nrow(heart_counts), ncol(heart_counts)))
message("  Factors: HU (HU/NL), IR (IR/Sham), KMP (KMP/Ctrl)")
message("  Samples per group:")
print(table(heart_coldata$group))

message("\n--- Cardiac Organoid ---")
message(sprintf("  Samples: %d genes x %d samples", nrow(organoid_counts), ncol(organoid_counts)))
message("  Factors: beam (Beam/NoBeam), treatment (KMP/Control), timepoint (Day2/Day9), batch (1-4)")
message("  Samples per group:")
print(table(organoid_coldata$group))
message("  Samples per batch:")
print(table(organoid_coldata$batch))

message(sprintf("\n  Tissue contrasts defined: %d", nrow(tissue_contrasts)))
message(sprintf("  Organoid contrasts defined: %d", nrow(organoid_contrasts)))

message("\nDone. All outputs saved to ", cfg$processed_dir, "/")
