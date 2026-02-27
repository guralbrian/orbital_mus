# 07_gsea.R
# Per-contrast GSEA using clusterProfiler (GO-BP + KEGG).
# Follows the ranking approach from bulk_decon/scripts/11_GSEA.R:
# signed -log10(pvalue) as the ranking metric.
#
# Usage: Rscript scripts/decon/07_gsea.R <ref_name> <model_type> <contrast_slug>
#   ref_name:      one of bulk_decon, mortazavi_general, mortazavi_celltype, mortazavi_b6j
#   model_type:    "adjusted" or "unadjusted"
#   contrast_slug: filesystem-safe contrast name (from _contrast_index.csv)
#
# Inputs:
#   - data/processed/models/{ref}_results_{model}.rds
#   - data/processed/models/{ref}_contrast_index.csv
#
# Outputs:
#   - data/processed/gsea/{ref}_{model}_{cslug}_go.rds
#   - data/processed/gsea/{ref}_{model}_{cslug}_kegg.rds
#   - data/processed/gsea/{ref}_{model}_{cslug}_summary.csv

library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(enrichplot)

set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
ref_name       <- args[1]
model_type     <- args[2]
contrast_slug  <- args[3]

stopifnot(!is.null(ref_name), !is.null(model_type), !is.null(contrast_slug))

proc_dir <- "data/processed"
out_dir  <- file.path(proc_dir, "gsea")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Load inputs ----

res_all <- readRDS(file.path(proc_dir, "models",
                             paste0(ref_name, "_results_", model_type, ".rds")))

contrast_idx <- read.csv(file.path(proc_dir, "models",
                                   paste0(ref_name, "_contrast_index.csv")))

# Look up original contrast name from slug
contrast_name <- contrast_idx$name[contrast_idx$slug == contrast_slug]
stopifnot("Contrast slug not found in index" = length(contrast_name) == 1)

res <- res_all[[contrast_name]]
stopifnot("Contrast not found in results RDS" = !is.null(res))

# ---- Build ranked gene list ----

# Ranking metric: signed -log10(pvalue), matching bulk_decon approach
pvals <- res$pvalue
pvals[is.na(pvals)] <- 1
pvals[pvals == 0] <- .Machine$double.xmin

gene_list <- -log10(pvals) * sign(res$log2FoldChange)
names(gene_list) <- res$gene_id
gene_list <- sort(gene_list, decreasing = TRUE)

# Remove duplicates and NAs
gene_list <- gene_list[!is.na(names(gene_list))]
gene_list <- gene_list[!duplicated(names(gene_list))]

# ---- GO Biological Process (ENSEMBL IDs) ----

gsea_go <- gseGO(
  geneList     = gene_list,
  OrgDb        = org.Mm.eg.db,
  ont          = "BP",
  keyType      = "ENSEMBL",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  eps          = 0
)

# ---- KEGG (requires ENTREZID) ----

id_map <- bitr(names(gene_list),
               fromType = "ENSEMBL",
               toType   = "ENTREZID",
               OrgDb    = org.Mm.eg.db)

# Build ENTREZ-keyed gene list (deduplicate both directions)
id_map <- id_map[!duplicated(id_map$ENSEMBL), ]
id_map <- id_map[!duplicated(id_map$ENTREZID), ]
entrez_list <- gene_list[id_map$ENSEMBL]
names(entrez_list) <- id_map$ENTREZID
entrez_list <- sort(entrez_list, decreasing = TRUE)

gsea_kegg <- gseKEGG(
  geneList     = entrez_list,
  organism     = "mmu",
  keyType      = "ncbi-geneid",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  eps          = 0
)

# ---- Save RDS objects ----

prefix <- file.path(out_dir, paste0(ref_name, "_", model_type, "_", contrast_slug))

saveRDS(gsea_go,   paste0(prefix, "_go.rds"))
saveRDS(gsea_kegg, paste0(prefix, "_kegg.rds"))

# ---- Write flat summary CSV ----

go_df <- if (nrow(as.data.frame(gsea_go)) > 0) {
  as.data.frame(gsea_go) |>
    dplyr::select(ID, Description, NES, pvalue, p.adjust, setSize) |>
    mutate(database = "GO_BP")
} else {
  data.frame(ID = character(), Description = character(), NES = numeric(),
             pvalue = numeric(), p.adjust = numeric(), setSize = integer(),
             database = character())
}

kegg_df <- if (nrow(as.data.frame(gsea_kegg)) > 0) {
  as.data.frame(gsea_kegg) |>
    dplyr::select(ID, Description, NES, pvalue, p.adjust, setSize) |>
    mutate(database = "KEGG")
} else {
  data.frame(ID = character(), Description = character(), NES = numeric(),
             pvalue = numeric(), p.adjust = numeric(), setSize = integer(),
             database = character())
}

summary_df <- bind_rows(go_df, kegg_df) |>
  mutate(
    reference     = ref_name,
    model         = model_type,
    contrast_slug = contrast_slug,
    contrast_name = contrast_name
  )

write.csv(summary_df, paste0(prefix, "_summary.csv"), row.names = FALSE)
