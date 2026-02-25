# 03_deconvolute.R
# Run MuSiC cell type deconvolution on heart bulk RNAseq using portable
# single-cell reference files (dgCMatrix RDS + CSV metadata/markers).
# No Seurat dependency — runs in the minimal orbital_mus_music env.
#
# Usage: Rscript scripts/decon/03_deconvolute.R <ref_name>
#   ref_name: one of bulk_decon, mortazavi_general, mortazavi_celltype, mortazavi_b6j
#
# Inputs:
#   - data/processed/decon/{ref_name}_counts.rds    (dgCMatrix)
#   - data/processed/decon/{ref_name}_metadata.csv  (barcode, celltype, sample)
#   - data/processed/decon/{ref_name}_markers.csv   (marker genes)
#   - data/processed/heart_counts.rds               (bulk count matrix)
#   - data/processed/heart_coldata.csv              (sample metadata)
#   - data/processed/heart_gene_names.csv           (Ensembl -> symbol map)
#
# Outputs:
#   - data/processed/decon/{ref_name}_compositions.csv (cell type proportions)
#   - data/processed/decon/{ref_name}_gene_weights.csv (MuSiC gene weights)

library(jsonlite)
library(dplyr)
library(tibble)
library(Biobase)
library(MuSiC)
library(reshape2)
library(SingleCellExperiment)

set.seed(42)

ref_name <- commandArgs(trailingOnly = TRUE)[1]
cfg <- fromJSON("scripts/setup/config.json")

stopifnot("ref_name CLI arg required" = !is.null(ref_name))
stopifnot("Unknown reference name" = ref_name %in% cfg$decon_references)

out_dir <- file.path(cfg$processed_dir, "decon")

# ---- Load portable single-cell reference ----

sc_counts <- readRDS(file.path(out_dir, paste0(ref_name, "_counts.rds")))
sc_meta   <- read.csv(file.path(out_dir, paste0(ref_name, "_metadata.csv")))
markers   <- read.csv(file.path(out_dir, paste0(ref_name, "_markers.csv")))

# Reconstruct SCE from dgCMatrix + metadata
sc_meta$celltype <- as.factor(sc_meta$celltype)
sc_meta$sample   <- as.factor(sc_meta$sample)
rownames(sc_meta) <- sc_meta$barcode

sce <- SingleCellExperiment(
  assays  = list(counts = sc_counts),
  colData = DataFrame(sc_meta)
)

# ---- Load and convert bulk data ----

bulk_counts <- readRDS(file.path(cfg$processed_dir, "heart_counts.rds"))
gene_map    <- read.csv(file.path(cfg$processed_dir, "heart_gene_names.csv"))
coldata     <- read.csv(file.path(cfg$processed_dir, "heart_coldata.csv"))

# Convert Ensembl IDs to gene symbols — for duplicates, keep highest total count
gene_map <- gene_map |>
  filter(gene_name != "" & !is.na(gene_name))

bulk_df <- as.data.frame(bulk_counts) |>
  tibble::rownames_to_column("gene_id") |>
  inner_join(gene_map, by = "gene_id") |>
  mutate(total = rowSums(across(where(is.numeric)))) |>
  arrange(desc(total)) |>
  filter(!duplicated(gene_name)) |>
  select(-gene_id, -total)

bulk_mat <- as.matrix(bulk_df[, -which(names(bulk_df) == "gene_name")])
rownames(bulk_mat) <- bulk_df$gene_name

# ---- Prepare bulk matrix (marker genes only, integer) ----

marker_genes <- markers$gene[markers$gene %in% rownames(bulk_mat)]
bulk_es_exp <- bulk_mat[marker_genes, ] |>
  apply(2, as.integer)
rownames(bulk_es_exp) <- marker_genes

bulk_es_exp <- bulk_es_exp[!is.na(rowMeans(bulk_es_exp)), ]

# ---- Run MuSiC ----

decon <- music_prop(
  bulk.mtx = bulk_es_exp,
  sc.sce   = sce,
  markers  = rownames(bulk_es_exp),
  clusters = "celltype",
  samples  = "sample"
)

# ---- Format and save results ----

# Gene weights
gene_weights <- as.data.frame(decon$Weight.gene)
gene_weights$gene <- rownames(gene_weights)
write.csv(gene_weights,
          file.path(out_dir, paste0(ref_name, "_gene_weights.csv")),
          row.names = FALSE)

# Cell type proportions (long format) merged with sample metadata
decon_melt <- reshape2::melt(decon$Est.prop.weighted)
colnames(decon_melt) <- c("sample_id", "CellType", "Prop")
decon_melt$sample_id <- as.character(decon_melt$sample_id)

compositions <- decon_melt |>
  left_join(coldata, by = "sample_id")

write.csv(compositions,
          file.path(out_dir, paste0(ref_name, "_compositions.csv")),
          row.names = FALSE)
