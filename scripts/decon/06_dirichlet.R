# 06_dirichlet.R
# Dirichlet regression to test whether cell type composition differs across
# experimental groups (HU × IR × KMP, 3-factor design).
# Adapted from ../bulk_decon/scripts/08_dirichlet.R for the 3-factor design.
#
# Usage: Rscript scripts/decon/06_dirichlet.R <ref_name>
#   ref_name: one of bulk_decon, mortazavi_general, mortazavi_celltype, mortazavi_b6j
#
# Inputs:
#   - data/processed/decon/{ref_name}_compositions.csv
#   - data/processed/heart_coldata.csv
#   - scripts/setup/config.json
#
# Outputs:
#   - data/processed/models/{ref_name}_dirichlet_coefficients.csv
#   - results/figures/{ref_name}_dirichlet_coeff.png

library(jsonlite)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DirichletReg)

set.seed(42)

ref_name <- commandArgs(trailingOnly = TRUE)[1]
cfg <- fromJSON("scripts/setup/config.json")

stopifnot("ref_name CLI arg required" = !is.null(ref_name))
stopifnot("Unknown reference name" = ref_name %in% cfg$decon_references)

model_dir <- file.path(cfg$processed_dir, "models")
fig_dir   <- file.path(cfg$results_dir, "figures")
if (!dir.exists(model_dir)) dir.create(model_dir, recursive = TRUE)
if (!dir.exists(fig_dir))   dir.create(fig_dir, recursive = TRUE)

alpha <- cfg$de_alpha

# ---- Section 1: Load and reshape data ----

compositions <- read.csv(
  file.path(cfg$processed_dir, "decon", paste0(ref_name, "_compositions.csv"))
)

# Pivot to wide: one row per sample, columns = cell type proportions
dir_mat <- compositions |>
  pivot_wider(id_cols = c(sample_id, HU, IR, KMP),
              names_from = CellType, values_from = Prop)

# Drop cell types with zero proportions in more than 1 sample
ct_cols <- setdiff(colnames(dir_mat), c("sample_id", "HU", "IR", "KMP"))
zero_counts <- colSums(dir_mat[, ct_cols] == 0)
drop_cts <- names(zero_counts[zero_counts > 1])
if (length(drop_cts) > 0) {
  warning("Dropping cell types with >1 zero-proportion sample: ",
          paste(drop_cts, collapse = ", "))
  ct_cols <- setdiff(ct_cols, drop_cts)
}

# Build DirichletRegData object from remaining proportion columns
dir_mat$CellTypes <- DR_data(dir_mat[, ct_cols])

# Set factor levels with reference categories
dir_mat$HU  <- factor(dir_mat$HU,  levels = c("NL", "HU"))
dir_mat$IR  <- factor(dir_mat$IR,  levels = c("Sham", "IR"))
dir_mat$KMP <- factor(dir_mat$KMP, levels = c("Ctrl", "KMP"))

# ---- Section 2: Fit Dirichlet regression ----

model <- DirichReg(CellTypes ~ HU * IR * KMP, data = dir_mat, model = "common")

# Extract coefficient table
coef_mat <- summary(model)[["coef.mat"]]
dir_results <- as.data.frame(coef_mat, check.names = FALSE)
rownames(dir_results) <- seq_len(nrow(dir_results))

n_coefs <- nrow(coef_mat) / length(ct_cols)
dir_results$CellType <- rep(summary(model)[["varnames"]], each = n_coefs)
dir_results$Feature  <- rownames(coef_mat)
colnames(dir_results)[colnames(dir_results) == "Std. Error"] <- "StdError"

write.csv(
  dir_results,
  file.path(model_dir, paste0(ref_name, "_dirichlet_coefficients.csv")),
  row.names = FALSE
)

# ---- Section 3: Coefficient error-bar plot ----

# Readable labels for the 3-factor terms
dir_results <- dir_results |>
  mutate(Feature_wrap = case_when(
    Feature == "(Intercept)"                  ~ "Intercept",
    Feature == "HUHU"                         ~ "HU",
    Feature == "IRIR"                         ~ "IR",
    Feature == "KMPKMP"                       ~ "KMP",
    Feature == "HUHU:IRIR"                    ~ "HU x IR",
    Feature == "HUHU:KMPKMP"                  ~ "HU x KMP",
    Feature == "IRIR:KMPKMP"                  ~ "IR x KMP",
    Feature == "HUHU:IRIR:KMPKMP"            ~ "HU x IR\n  x KMP",
    TRUE                                      ~ Feature
  ))

# Plot non-intercept coefficients
err_plot <- dir_results |>
  filter(Feature != "(Intercept)") |>
  ggplot(aes(y = Feature_wrap, x = Estimate, color = `Pr(>|z|)`, shape = CellType)) +
  geom_errorbarh(
    aes(xmin = Estimate - StdError, xmax = Estimate + StdError),
    height = 0.5, position = position_dodge(width = 0.7), linewidth = 0.5
  ) +
  geom_point(position = position_dodge(width = 0.7), size = 1.4, color = "#474747") +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.5) +
  binned_scale(
    aesthetics = "color",
    scale_name = "stepsn",
    palette = function(x) c("#FDE725FF", "#73D055FF", "#20A387FF", "#482677FF"),
    breaks = c(0.005, 0.01, 0.05),
    limits = c(0.0005, 0.5),
    show.limits = TRUE,
    guide = "colorsteps"
  ) +
  guides(shape = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(
    axis.text.x     = element_text(color = "black"),
    axis.text.y     = element_text(color = "black"),
    axis.ticks       = element_blank(),
    axis.title.y     = element_blank(),
    plot.title       = element_blank(),
    legend.position  = "right",
    legend.text      = element_text(size = 6),
    legend.key.size  = unit(0.3, "lines"),
    legend.spacing.y = unit(0.5, "mm"),
    panel.grid.major = element_line(colour = "grey", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.background  = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent"),
    text             = element_text(size = 6),
    plot.margin      = unit(c(0, -0.2, 0, 0), "cm")
  ) +
  labs(
    x     = "Coefficients (estimates of effect)",
    shape = "Cell Types",
    color = "p-value"
  )

ggsave(
  filename = file.path(fig_dir, paste0(ref_name, "_dirichlet_coeff.png")),
  plot     = err_plot,
  width    = 6.23, height = 4.70, units = "cm", dpi = 600
)
