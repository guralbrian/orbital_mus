# 04_runDe.R
# Differential expression analysis on heart bulk RNAseq, comparing models
# with and without CLR-transformed cell type composition covariates.
# Uses a cell-means parameterization (~ 0 + group) with numeric contrast
# vectors for all 19 pre-defined comparisons. LFC shrinkage via ashr
# (supports contrast vectors, unlike apeglm which requires coef indices).
#
# Usage: Rscript scripts/decon/04_runDe.R <ref_name>
#   ref_name: one of bulk_decon, mortazavi_general, mortazavi_celltype, mortazavi_b6j
#
# Inputs:
#   - data/processed/heart_counts.rds
#   - data/processed/heart_coldata.csv
#   - data/processed/heart_gene_names.csv
#   - data/processed/contrast_definitions.rds
#   - data/processed/decon/{ref_name}_compositions.csv
#   - scripts/setup/config.json
#
# Outputs:
#   - data/processed/models/{ref_name}_dds_adjusted.rds
#   - data/processed/models/{ref_name}_dds_unadjusted.rds
#   - data/processed/models/{ref_name}_results_adjusted.rds
#   - data/processed/models/{ref_name}_results_unadjusted.rds
#   - data/processed/models/{ref_name}_de_comparison.csv
#   - data/processed/models/{ref_name}_clr_covariates.csv
#   - data/processed/models/{ref_name}_contrast_index.csv

library(jsonlite)
library(dplyr)
library(tidyr)
library(compositions)
library(DESeq2)
library(ashr)

set.seed(42)

ref_name <- commandArgs(trailingOnly = TRUE)[1]
cfg <- fromJSON("scripts/setup/config.json")

stopifnot("ref_name CLI arg required" = !is.null(ref_name))
stopifnot("Unknown reference name" = ref_name %in% cfg$decon_references)

out_dir <- file.path(cfg$processed_dir, "models")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

alpha <- cfg$de_alpha

# ---- Section 1: Load inputs ----

bulk_counts <- readRDS(file.path(cfg$processed_dir, "heart_counts.rds"))
coldata     <- read.csv(file.path(cfg$processed_dir, "heart_coldata.csv"),
                        row.names = 1)
gene_map    <- read.csv(file.path(cfg$processed_dir, "heart_gene_names.csv"))
contrasts   <- readRDS(file.path(cfg$processed_dir, "contrast_definitions.rds"))
compositions <- read.csv(
  file.path(cfg$processed_dir, "decon", paste0(ref_name, "_compositions.csv"))
)

tissue_contrasts <- contrasts$tissue

# ---- Section 2: CLR-transform compositions ----

# Pivot to wide: samples x cell types
comp_wide <- compositions |>
  select(sample_id, CellType, Prop) |>
  pivot_wider(names_from = CellType, values_from = Prop) |>
  as.data.frame()
rownames(comp_wide) <- comp_wide$sample_id
comp_wide$sample_id <- NULL

# Impute zeros: replace 0 with min(nonzero) / 2, then renormalize rows
all_props <- as.numeric(as.matrix(comp_wide))
min_nonzero <- min(all_props[all_props > 0])
comp_wide[comp_wide == 0] <- min_nonzero / 2

row_sums <- rowSums(comp_wide)
comp_wide <- comp_wide / row_sums

# CLR transform: maps proportions from the simplex to unconstrained real space,
# making them suitable as continuous covariates in the DESeq2 linear model
comp_clr <- as.data.frame(compositions::clr(as.matrix(comp_wide)))
colnames(comp_clr) <- paste0("clr.", make.names(colnames(comp_wide)))

# Select top 2 most abundant cell types as CLR covariates
mean_props <- colMeans(comp_wide)
top2_types <- names(sort(mean_props, decreasing = TRUE))[1:2]
clr_keep <- paste0("clr.", make.names(top2_types))
comp_clr <- comp_clr[, clr_keep, drop = FALSE]

write.csv(
  data.frame(reference = ref_name, celltype = top2_types, clr_name = clr_keep),
  file.path(out_dir, paste0(ref_name, "_clr_covariates.csv")),
  row.names = FALSE
)

# Merge CLR columns into coldata (match row order)
coldata <- cbind(coldata, comp_clr[rownames(coldata), , drop = FALSE])

# Set reference level for group
coldata$group <- factor(coldata$group, levels = cfg$heart_group_order)
coldata$group <- relevel(coldata$group, ref = "NL_Sham_Ctrl")

# Ensure count matrix columns match coldata rows
bulk_counts <- bulk_counts[, rownames(coldata)]

# ---- Section 3: Build and run DESeq2 models ----

run_deseq <- function(counts, sample_info, design_formula, output_path) {
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = sample_info,
    design    = design_formula
  )
  keep <- rowSums(counts(dds) >= cfg$de_min_count) >= cfg$de_min_samples
  dds <- dds[keep, ]
  dds <- DESeq(dds)
  saveRDS(dds, output_path)
  dds
}

# Build adjusted formula dynamically from CLR column names
adj_formula <- as.formula(
  paste("~ 0 + group +", paste(clr_keep, collapse = " + "))
)
unadj_formula <- ~ 0 + group

dds_adj <- run_deseq(
  bulk_counts, coldata, adj_formula,
  file.path(out_dir, paste0(ref_name, "_dds_adjusted.rds"))
)

dds_unadj <- run_deseq(
  bulk_counts, coldata, unadj_formula,
  file.path(out_dir, paste0(ref_name, "_dds_unadjusted.rds"))
)

# ---- Section 4: Build contrast vectors ----

# Group levels as they appear in the model (prefixed with "group")
group_levels <- cfg$heart_group_order
model_names_adj   <- colnames(model.matrix(design(dds_adj), data = colData(dds_adj)))
model_names_unadj <- colnames(model.matrix(design(dds_unadj), data = colData(dds_unadj)))

n_adj   <- length(model_names_adj)
n_unadj <- length(model_names_unadj)

# Helper: create a numeric contrast vector for the adjusted model.
# group_weights is a named numeric vector (names = group levels, values = weights).
# CLR positions are always 0.
make_contrast_adj <- function(group_weights) {
  vec <- numeric(n_adj)
  names(vec) <- model_names_adj
  for (g in names(group_weights)) {
    coef_name <- paste0("group", g)
    vec[coef_name] <- group_weights[g]
  }
  vec
}

make_contrast_unadj <- function(group_weights) {
  vec <- numeric(n_unadj)
  names(vec) <- model_names_unadj
  for (g in names(group_weights)) {
    coef_name <- paste0("group", g)
    vec[coef_name] <- group_weights[g]
  }
  vec
}

# Build group_weights for each of the 19 contrasts
contrast_weights <- list()

for (i in seq_len(nrow(tissue_contrasts))) {
  cname  <- tissue_contrasts$comparison_name[i]
  group1 <- tissue_contrasts$group1[i]
  group2 <- tissue_contrasts$group2[i]

  # Simple pairwise: both group1 and group2 are single valid group names
  if (group1 %in% group_levels && group2 %in% group_levels) {
    w <- setNames(c(1, -1), c(group1, group2))
    contrast_weights[[cname]] <- w
    next
  }

  # Hard-coded contrasts for main effects, marginal means, and interactions
  w <- switch(cname,
    "HU_vs_NL_MAIN_EFFECT" = setNames(
      c(rep(1/4, 4), rep(-1/4, 4)),
      c("HU_Sham_Ctrl", "HU_Sham_KMP", "HU_IR_Ctrl", "HU_IR_KMP",
        "NL_Sham_Ctrl", "NL_Sham_KMP", "NL_IR_Ctrl", "NL_IR_KMP")
    ),
    "IR_vs_Sham_MAIN_EFFECT" = setNames(
      c(rep(1/4, 4), rep(-1/4, 4)),
      c("NL_IR_Ctrl", "NL_IR_KMP", "HU_IR_Ctrl", "HU_IR_KMP",
        "NL_Sham_Ctrl", "NL_Sham_KMP", "HU_Sham_Ctrl", "HU_Sham_KMP")
    ),
    "KMP_vs_Sham_MAIN_EFFECT" = setNames(
      c(rep(1/4, 4), rep(-1/4, 4)),
      c("NL_Sham_KMP", "NL_IR_KMP", "HU_Sham_KMP", "HU_IR_KMP",
        "NL_Sham_Ctrl", "NL_IR_Ctrl", "HU_Sham_Ctrl", "HU_IR_Ctrl")
    ),
    # Marginal means: IR+KMP averaged over HU/NL vs Ctrl averaged over HU/NL
    "IR+KMP (HU+NL)_vs_Ctrl (HU+NL)" = setNames(
      c(1/2, 1/2, -1/2, -1/2),
      c("NL_IR_KMP", "HU_IR_KMP", "NL_Sham_Ctrl", "HU_Sham_Ctrl")
    ),
    # HU+KMP averaged over IR/Sham vs Ctrl averaged over IR/Sham
    "HU+KMP (IR+Sham)_vs_Ctrl (IR+Sham)" = setNames(
      c(1/2, 1/2, -1/2, -1/2),
      c("HU_Sham_KMP", "HU_IR_KMP", "NL_Sham_Ctrl", "NL_IR_Ctrl")
    ),
    # HU+IR averaged over KMP/Ctrl vs Ctrl averaged over KMP/Ctrl
    "HU+IR (KMP+Ctrl)_vs_Ctrl (KMP+Ctrl)" = setNames(
      c(1/2, 1/2, -1/2, -1/2),
      c("HU_IR_KMP", "HU_IR_Ctrl", "NL_IR_Ctrl", "NL_Sham_Ctrl")
    ),
    # Interaction: KMP x HU = (KMP effect in HU) - (KMP effect in NL)
    #   KMP effect in HU = avg(HU_Sham_KMP, HU_IR_KMP) - avg(HU_Sham_Ctrl, HU_IR_Ctrl)
    #   KMP effect in NL = avg(NL_Sham_KMP, NL_IR_KMP) - avg(NL_Sham_Ctrl, NL_IR_Ctrl)
    "Interaction_KMP_x_HU" = setNames(
      c(1/2, 1/2, -1/2, -1/2, -1/2, -1/2, 1/2, 1/2),
      c("HU_Sham_KMP", "HU_IR_KMP", "HU_Sham_Ctrl", "HU_IR_Ctrl",
        "NL_Sham_KMP", "NL_IR_KMP", "NL_Sham_Ctrl", "NL_IR_Ctrl")
    ),
    # Interaction: KMP x IR = (KMP effect in IR) - (KMP effect in Sham)
    #   KMP effect in IR   = avg(NL_IR_KMP, HU_IR_KMP) - avg(NL_IR_Ctrl, HU_IR_Ctrl)
    #   KMP effect in Sham = avg(NL_Sham_KMP, HU_Sham_KMP) - avg(NL_Sham_Ctrl, HU_Sham_Ctrl)
    "Interaction_KMP_x_IR" = setNames(
      c(1/2, 1/2, -1/2, -1/2, -1/2, -1/2, 1/2, 1/2),
      c("NL_IR_KMP", "HU_IR_KMP", "NL_IR_Ctrl", "HU_IR_Ctrl",
        "NL_Sham_KMP", "HU_Sham_KMP", "NL_Sham_Ctrl", "HU_Sham_Ctrl")
    ),
    # Interaction: HU x IR = (HU effect in IR) - (HU effect in Sham)
    #   HU effect in IR   = avg(HU_IR_Ctrl, HU_IR_KMP) - avg(NL_IR_Ctrl, NL_IR_KMP)
    #   HU effect in Sham = avg(HU_Sham_Ctrl, HU_Sham_KMP) - avg(NL_Sham_Ctrl, NL_Sham_KMP)
    "Interaction_HU_x_IR" = setNames(
      c(1/2, 1/2, -1/2, -1/2, -1/2, -1/2, 1/2, 1/2),
      c("HU_IR_Ctrl", "HU_IR_KMP", "NL_IR_Ctrl", "NL_IR_KMP",
        "HU_Sham_Ctrl", "HU_Sham_KMP", "NL_Sham_Ctrl", "NL_Sham_KMP")
    ),
    stop(paste("Unrecognized contrast:", cname))
  )

  contrast_weights[[cname]] <- w
}

# ---- Section 5: Extract results with shrinkage ----

extract_results <- function(dds, make_fn) {
  res_list <- list()
  for (cname in names(contrast_weights)) {
    cvec <- make_fn(contrast_weights[[cname]])
    res <- results(dds, contrast = cvec)
    res_shrunk <- lfcShrink(dds, contrast = cvec, res = res, type = "ashr")
    res_list[[cname]] <- as.data.frame(res_shrunk)
    res_list[[cname]]$gene_id <- rownames(res_shrunk)
  }
  res_list
}

res_adj   <- extract_results(dds_adj, make_contrast_adj)
res_unadj <- extract_results(dds_unadj, make_contrast_unadj)

saveRDS(res_adj,   file.path(out_dir, paste0(ref_name, "_results_adjusted.rds")))
saveRDS(res_unadj, file.path(out_dir, paste0(ref_name, "_results_unadjusted.rds")))

# ---- Section 6: Merge adjusted vs unadjusted ----

comparison_list <- list()

for (cname in names(contrast_weights)) {
  merged <- merge(
    res_adj[[cname]], res_unadj[[cname]],
    by = "gene_id", suffixes = c(".adj", ".unadj")
  )
  merged$contrast <- cname
  merged$p_diff <- -log10(merged$padj.adj) - (-log10(merged$padj.unadj))
  merged$adj_only_sig   <- merged$padj.adj <= alpha & merged$padj.unadj > alpha
  merged$unadj_only_sig <- merged$padj.adj > alpha  & merged$padj.unadj <= alpha
  merged$both_sig       <- merged$padj.adj <= alpha  & merged$padj.unadj <= alpha

  # Add gene names
  merged <- merge(merged, gene_map, by = "gene_id", all.x = TRUE)

  comparison_list[[cname]] <- merged
}

comparison <- do.call(rbind, comparison_list)
rownames(comparison) <- NULL

write.csv(comparison,
          file.path(out_dir, paste0(ref_name, "_de_comparison.csv")),
          row.names = FALSE)

# ---- Section 7: Write contrast index for GSEA jobs ----

contrast_index <- data.frame(
  index = seq_along(names(contrast_weights)),
  name = names(contrast_weights),
  slug = gsub("[^A-Za-z0-9_]", "_", names(contrast_weights))
)
write.csv(contrast_index,
          file.path(out_dir, paste0(ref_name, "_contrast_index.csv")),
          row.names = FALSE)
