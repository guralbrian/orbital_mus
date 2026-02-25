# 02_find_markers.R
# Load all single-cell references, identify cell-type marker genes, and save
# portable outputs (dgCMatrix RDS + CSV) for downstream MuSiC deconvolution
# in a separate conda env.
#
# Usage: Rscript scripts/decon/02_find_markers.R
#   (no CLI args — processes all 4 references in one run)
#
# Inputs:
#   - data/processed/sc/celltype_labeled.rds  (Seurat object)
#   - data/processed/sc/IGVFFI6644FMFS.rds    (Seurat object)
#   - data/processed/heart_gene_names.csv      (gene_id <-> gene_name mapping)
#
# Outputs (per reference):
#   - data/processed/decon/{ref}_counts.rds    (dgCMatrix — portable across R versions)
#   - data/processed/decon/{ref}_metadata.csv  (barcode, celltype, sample)
#   - data/processed/decon/{ref}_markers.csv   (gene, celltype, p.value, FDR, summary.logFC)

library(dplyr)
library(tibble)
library(Seurat)
library(SingleCellExperiment)

set.seed(123)

out_dir <- "data/processed/decon"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Load bulk gene names (shared gene universe) ----

gene_map <- read.csv("data/processed/heart_gene_names.csv")
bulk_genes <- unique(gene_map$gene_name)

# ---- Load single-cell references (each source loaded once) ----

sc_labeled <- readRDS("data/processed/sc/celltype_labeled.rds")
sc_igvf    <- readRDS("data/processed/sc/IGVFFI6644FMFS.rds")

# Handle list-style RDS (from convert_h5seurat.R)
if (is.list(sc_labeled) && !inherits(sc_labeled, "Seurat")) {
  sc_labeled <- CreateSeuratObject(counts = sc_labeled$counts,
                                   meta.data = sc_labeled$meta)
}
if (is.list(sc_igvf) && !inherits(sc_igvf, "Seurat")) {
  sc_igvf <- CreateSeuratObject(counts = sc_igvf$counts,
                                meta.data = sc_igvf$meta)
}

# ---- Reference definitions ----

refs <- list(
  bulk_decon = list(
    sc             = sc_labeled,
    celltype_col   = "active_ident",
    sample_col     = "orig.ident",
    genotype_filter = NULL
  ),
  mortazavi_general = list(
    sc             = sc_igvf,
    celltype_col   = "general_celltype",
    sample_col     = "lab_sample_id",
    genotype_filter = NULL
  ),
  mortazavi_celltype = list(
    sc             = sc_igvf,
    celltype_col   = "celltype",
    sample_col     = "lab_sample_id",
    genotype_filter = NULL
  ),
  mortazavi_b6j = list(
    sc             = sc_igvf,
    celltype_col   = "general_celltype",
    sample_col     = "lab_sample_id",
    genotype_filter = "B6J"
  )
)

# ---- Helper: process one reference ----

process_reference <- function(ref_name, ref) {
  sc <- ref$sc
  ct_col <- ref$celltype_col
  sample_col <- ref$sample_col

  # Remove cells with NA or "low quality" cell type labels
  keep <- !is.na(sc@meta.data[[ct_col]]) &
          sc@meta.data[[ct_col]] != "low quality"
  sc <- sc[, keep]

  # Optional genotype filter (e.g. B6J only)
  if (!is.null(ref$genotype_filter)) {
    geno_keep <- sc$Genotype == ref$genotype_filter
    sc <- sc[, geno_keep]
  }

  Idents(sc) <- ct_col

  # Drop celltypes with <20 cells
  ct_counts <- table(Idents(sc))
  keep_types <- names(ct_counts[ct_counts >= 20])
  sc <- sc[, Idents(sc) %in% keep_types]
  Idents(sc) <- droplevels(Idents(sc))

  # Handle non-integer counts (h5ad conversion may store log1p values)
  counts_mat <- LayerData(sc, layer = "counts")
  if (!all(counts_mat@x == floor(counts_mat@x))) {
    counts_mat@x <- round(expm1(counts_mat@x))
    LayerData(sc, layer = "counts") <- counts_mat
  }

  # Subset to genes measured in bulk
  shared_genes <- intersect(rownames(sc), bulk_genes)
  sc <- subset(sc, features = shared_genes)

  # ---- Find markers with scran ----

  sc_sce <- as.SingleCellExperiment(sc, assay = "RNA")

  scn_markers <- scran::findMarkers(
    sc_sce,
    groups = Idents(sc),
    pval.type = "all",
    assay.type = "counts"
  )

  all_markers <- lapply(levels(Idents(sc)), function(type) {
    scn_markers@listData[[type]] |>
      as.data.frame() |>
      dplyr::select(p.value, FDR, summary.logFC) |>
      mutate(celltype = type) |>
      rownames_to_column(var = "gene")
  }) |>
    bind_rows()

  # Top 15 per celltype, logFC >= 1, remove cross-type duplicates
  top_markers <- all_markers |>
    group_by(celltype) |>
    filter(summary.logFC >= 1) |>
    arrange(p.value) |>
    slice_head(n = 15) |>
    ungroup()

  # Genes appearing as top markers for multiple cell types are ambiguous and
  # would bias MuSiC's weighted deconvolution — remove them entirely
  dup_genes <- top_markers$gene[duplicated(top_markers$gene)]
  top_markers <- top_markers |>
    filter(!(gene %in% dup_genes))

  # ---- Save portable outputs ----

  # dgCMatrix (stable across R 4.3/4.4)
  final_counts <- LayerData(sc, layer = "counts")
  saveRDS(final_counts, file.path(out_dir, paste0(ref_name, "_counts.rds")))

  # Metadata CSV (barcode, celltype, sample)
  meta <- data.frame(
    barcode  = colnames(sc),
    celltype = as.character(Idents(sc)),
    sample   = sc@meta.data[[sample_col]],
    stringsAsFactors = FALSE
  )
  write.csv(meta, file.path(out_dir, paste0(ref_name, "_metadata.csv")),
            row.names = FALSE)

  # Markers CSV
  write.csv(top_markers, file.path(out_dir, paste0(ref_name, "_markers.csv")),
            row.names = FALSE)
}

# ---- Run all references ----

for (ref_name in names(refs)) {
  process_reference(ref_name, refs[[ref_name]])
}
