# 03_projectils.R
# Classify cells using the ProjecTILs TIL reference atlas.
#
# ProjecTILs projects query cells into a pre-built TIL reference embedding
# and assigns each cell a TIL subtype label (e.g. Tex, Tem, Teff, Tn, Tpex).
# This gives a data-driven functional annotation independent of the CAR-T
# stimulation labels, making it possible to ask whether CR patients are
# enriched for Tpex (progenitor-exhausted) vs. Tex (terminally exhausted).
#
# Prerequisites:
#   remotes::install_github("carmonalab/UCell")
#   remotes::install_github("carmonalab/ProjecTILs")
#
# Input:  data/processed/01_annotated.rds
# Output: data/processed/03_projectils.rds
#         results/plots/03_projectils_umap.pdf
#         results/plots/03_projectils_composition.pdf

library(Seurat)
library(ggplot2)
library(patchwork)

# Fail fast with an informative message if ProjecTILs is not installed.
if (!requireNamespace("ProjecTILs", quietly = TRUE))
  stop("ProjecTILs not installed.\n",
       "Run: remotes::install_github('carmonalab/UCell')\n",
       "     remotes::install_github('carmonalab/ProjecTILs')")
library(ProjecTILs)

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
seu <- readRDS("data/processed/01_annotated.rds")
cat("Loaded object:", ncol(seu), "cells\n")

# ProjecTILs operates on the RNA assay.  Ensure it is the default.
DefaultAssay(seu) <- "RNA"

# ===========================================================================
# 1. Load TIL reference atlas
# ===========================================================================
# CHOICE: use the default human TIL CD8 reference (Carmona et al. 2021).
# This is the most appropriate reference for CD8+ CAR-T cells from ALL patients.
# If a CD4 analysis is also needed, call ProjecTILs again with ref.type = "CD4".
#
# The atlas is downloaded and cached automatically by ProjecTILs on first call.
# Subsequent calls use the local cache (~300 MB, stored in the R temp or home dir).

cat("Fetching TIL reference atlas...\n")
ref_cd8 <- ProjecTILs::load.reference.map()   # default = human CD8 TIL atlas

cat("\n--- Reference cell type composition ---\n")
print(table(ref_cd8$functional.cluster))

# ===========================================================================
# 2. Run ProjecTILs classifier
# ===========================================================================
# CHOICE: classifier = "knn" is the default and is recommended for Seurat v4
# objects without a pre-existing SingleCellExperiment step.
# run.projection = TRUE also returns the projected UMAP coordinates in the
# reference space, which we use to make the atlas-overlay plot below.

cat("\nRunning ProjecTILs classifier...\n")
time_proj <- system.time({
  seu_proj <- ProjecTILs.classifier(
    query     = seu,
    ref       = ref_cd8,
    filter.cells = TRUE,   # remove cells that project poorly (low confidence)
    ncores    = 1L         # increase if running on a multi-core machine
  )
})
cat("ProjecTILs run time:\n"); print(time_proj)

# ===========================================================================
# 3. Inspect classification output
# ===========================================================================
# ProjecTILs adds "functional.cluster" and "functional.cluster.conf" to meta.data.
# conf is a 0–1 confidence score; cells below ~0.3 should be treated with caution.

cat("\n--- ProjecTILs labels assigned ---\n")
print(table(seu_proj$functional.cluster, useNA = "ifany"))

cat("\n--- Confidence score summary ---\n")
if ("functional.cluster.conf" %in% colnames(seu_proj@meta.data)) {
  print(summary(seu_proj$functional.cluster.conf))

  low_conf <- mean(seu_proj$functional.cluster.conf < 0.3, na.rm = TRUE)
  cat(sprintf("%.1f%% of cells have confidence < 0.3 — treat their labels cautiously.\n",
              low_conf * 100))
}

cat("\n--- Label composition by response ---\n")
tbl <- table(seu_proj$functional.cluster, seu_proj$response)
print(tbl)
print(prop.table(tbl, margin = 2))   # column percentages

# ===========================================================================
# 4. Visualise on original UMAP
# ===========================================================================
plot_dir <- "results/plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# TIL type colour scheme following ProjecTILs convention.
til_cols <- c(
  "Tex"   = "#d62728",   # terminally exhausted
  "Tpex"  = "#ff7f0e",   # progenitor-exhausted
  "Tem"   = "#2ca02c",   # effector-memory
  "Teff"  = "#1f77b4",   # effector
  "Tn"    = "#9467bd",   # naive
  "Tcm"   = "#8c564b"    # central-memory
)

p_til <- DimPlot(
  seu_proj,
  reduction = "umap",
  group.by  = "functional.cluster",
  cols      = til_cols,
  pt.size   = 0.3,
  alpha     = 0.6
) +
  ggtitle("ProjecTILs TIL subtype") +
  theme_classic(base_size = 11)

p_til_split <- DimPlot(
  seu_proj,
  reduction = "umap",
  group.by  = "functional.cluster",
  split.by  = "response",
  cols      = til_cols,
  pt.size   = 0.3,
  alpha     = 0.5,
  ncol      = 3
) +
  ggtitle("TIL subtype split by response") +
  theme_classic(base_size = 10)

p_layout <- p_til / p_til_split
ggsave(file.path(plot_dir, "03_projectils_umap.pdf"), p_layout,
       width = 14, height = 12)
message("Saved: results/plots/03_projectils_umap.pdf")

# ===========================================================================
# 5. Composition barplot per patient
# ===========================================================================
# Stacked barplots show whether CR/RL/NR patients differ in TIL subtype mix.
# CHOICE: normalise to proportions (not raw counts) so unequal cell numbers
# per patient do not confound the comparison.

comp_df <- as.data.frame(
  prop.table(table(seu_proj$patient_id, seu_proj$functional.cluster), margin = 1)
)
colnames(comp_df) <- c("patient_id", "TIL_type", "proportion")

# Add response for fill / facet ordering.
resp_df <- unique(seu_proj@meta.data[, c("patient_id", "response")])
comp_df  <- merge(comp_df, resp_df, by = "patient_id")

p_comp <- ggplot(comp_df, aes(x = patient_id, y = proportion, fill = TIL_type)) +
  geom_col() +
  scale_fill_manual(values = til_cols, na.value = "grey70") +
  facet_grid(~ response, scales = "free_x", space = "free_x") +
  labs(
    title = "TIL subtype composition per patient",
    x     = NULL, y = "Proportion"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    strip.background = element_blank()
  )

ggsave(file.path(plot_dir, "03_projectils_composition.pdf"), p_comp,
       width = 10, height = 5)
message("Saved: results/plots/03_projectils_composition.pdf")

# ===========================================================================
# 6. Save
# ===========================================================================
saveRDS(seu_proj, "data/processed/03_projectils.rds")
message("Saved: data/processed/03_projectils.rds")

cat("\n--- Power note ---\n")
cat("n = 12 donors; per-group n = 5 (CR), 5 (RL), 2 (NR).\n")
cat("Treat composition differences as descriptive / hypothesis-generating.\n")
cat("Statistical tests on TIL proportions are severely underpowered at n = 2 (NR).\n")
