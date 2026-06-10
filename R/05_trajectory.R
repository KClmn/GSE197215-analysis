# 05_trajectory.R
# Trajectory inference with Slingshot, rooted at stem-like memory T cells.
#
# Strategy:
#   1. Convert Seurat → SingleCellExperiment (SCE)
#   2. Run Slingshot on existing UMAP embedding, using Seurat clusters as
#      provisional states and the stem-like memory cluster as the root
#   3. Validate trajectory: exhaustion genes should increase, memory genes
#      should decrease along pseudotime
#   4. Add pseudotime to Seurat metadata for downstream use (MOFA, variancePartition)
#
# Input:  data/processed/04_scored.rds
# Output: data/processed/05_trajectory.rds
#         results/plots/05_pseudotime_umap.pdf
#         results/plots/05_trajectory_validation.pdf
#         results/tables/05_pseudotime_summary.csv

library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(ggplot2)
library(patchwork)
library(viridis)

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
seu <- readRDS("data/processed/04_scored.rds")
cat("Loaded:", ncol(seu), "cells\n")

DefaultAssay(seu) <- "RNA"
plot_dir  <- "results/plots";  if (!dir.exists(plot_dir))  dir.create(plot_dir,  recursive = TRUE)
table_dir <- "results/tables"; if (!dir.exists(table_dir)) dir.create(table_dir, recursive = TRUE)

# ===========================================================================
# 1. Identify the stem-like memory cluster (trajectory root)
# ===========================================================================
# Root cells are identified by HIGH expression of:
#   TCF7 (Tcf1 — stem-like TIL transcription factor)
#   CCR7 (homing to lymph nodes, naive/central-memory marker)
#   SELL (CD62L — lymph node homing)
#   IL7R (CD127 — survival signal, naive/memory)
# and LOW expression of exhaustion genes.
#
# CHOICE: we use per-cluster mean expression of these markers to pick the root
# cluster, rather than hard-coding a cluster number, so the script is robust to
# re-clustering with different resolution parameters.

stem_markers <- c("TCF7", "CCR7", "SELL", "IL7R")
exh_markers  <- c("PDCD1", "TOX", "LAG3", "HAVCR2", "TIGIT")

present_stem <- intersect(stem_markers, rownames(seu[["RNA"]]))
present_exh  <- intersect(exh_markers,  rownames(seu[["RNA"]]))

cat("Stem marker genes found:", paste(present_stem, collapse = ", "), "\n")
cat("Exhaustion marker genes found:", paste(present_exh, collapse = ", "), "\n")

# Ensure Idents are seurat_clusters, not orig.ident.
# Idents can revert to orig.ident across save/load cycles if not explicitly set.
if ("seurat_clusters" %in% colnames(seu@meta.data)) {
  Idents(seu) <- "seurat_clusters"
} else {
  stop("seurat_clusters not found in metadata — run 02_visualize.R first.")
}

# Compute per-cluster mean for stem and exhaustion markers.
cluster_ids <- levels(Idents(seu))
mean_stem <- sapply(cluster_ids, function(cl) {
  cells <- WhichCells(seu, idents = cl)
  mean(colMeans(GetAssayData(seu, layer = "data")[present_stem, cells, drop = FALSE]))
})
mean_exh <- sapply(cluster_ids, function(cl) {
  cells <- WhichCells(seu, idents = cl)
  mean(colMeans(GetAssayData(seu, layer = "data")[present_exh, cells, drop = FALSE]))
})

# Stem score = high stem, low exhaustion.
stem_score <- mean_stem - mean_exh
root_cluster <- cluster_ids[which.max(stem_score)]

cat("\n--- Cluster stem scores (high = more stem-like) ---\n")
stem_df <- data.frame(
  cluster    = cluster_ids,
  mean_stem  = round(mean_stem, 3),
  mean_exh   = round(mean_exh,  3),
  stem_score = round(stem_score, 3)
)
print(stem_df[order(-stem_df$stem_score), ])
cat("\nSelected root cluster:", root_cluster, "\n")

# ===========================================================================
# 2. Convert to SCE and run Slingshot
# ===========================================================================
# CHOICE: use the existing UMAP embedding rather than PCA.
# UMAP captures non-linear structure that is important for T cell exhaustion
# trajectories (the Tex branch is often separated from Teff in UMAP but not PCA).
# However, Slingshot's minimum spanning tree is more stable in PCA space for
# highly branched topologies — if the UMAP trajectory looks incorrect, re-run
# with reducedDimName = "PCA" as a sanity check.

umap_key <- if ("umap" %in% Reductions(seu)) "umap" else Reductions(seu)[1]
cat("Using reduction:", umap_key, "\n")

sce <- as.SingleCellExperiment(seu)

# Seurat v5's as.SingleCellExperiment() does not always transfer reductions.
# Inject UMAP embeddings manually so Slingshot can find them.
cat("reducedDims available in SCE:", paste(reducedDimNames(sce), collapse = ", "), "\n")

sce_umap_name <- "UMAP"
if (!sce_umap_name %in% reducedDimNames(sce)) {
  umap_embed <- Embeddings(seu, reduction = umap_key)
  reducedDim(sce, "UMAP") <- umap_embed
  cat("UMAP manually injected into SCE reducedDims\n")
}
cat("Using reducedDim:", sce_umap_name, "\n")

cat("\nRunning Slingshot...\n")
set.seed(42)
time_sling <- system.time({
  sds <- slingshot(
    data          = sce,
    clusterLabels = as.character(colData(sce)$seurat_clusters),
    reducedDim    = sce_umap_name,
    start.clus    = as.character(root_cluster),
    approx_points = 150
  )
})
cat("Slingshot run time:\n"); print(time_sling)

# ===========================================================================
# 3. Extract pseudotime
# ===========================================================================
# slingPseudotime() returns a matrix: cells × lineages.
# Cells not assigned to a lineage have NA.
pt_mat <- slingPseudotime(sds)

cat("\n--- Pseudotime dimensions (cells × lineages) ---\n")
print(dim(pt_mat))

n_lineages <- ncol(pt_mat)
cat("Number of lineages detected:", n_lineages, "\n")

# Summarise NA rate per lineage (cells off that branch).
for (i in seq_len(n_lineages)) {
  cat(sprintf("  Lineage %d: %d cells assigned (%.1f%%), %d NA\n",
              i, sum(!is.na(pt_mat[, i])),
              100 * mean(!is.na(pt_mat[, i])),
              sum(is.na(pt_mat[, i]))))
}

# Add pseudotime columns to Seurat metadata.
# If multiple lineages, also compute the minimum-NA pseudotime across them.
for (i in seq_len(n_lineages)) {
  col_name         <- paste0("pseudotime_L", i)
  seu[[col_name]]  <- pt_mat[colnames(seu), i]
}
if (n_lineages > 1) {
  seu$pseudotime_mean <- rowMeans(pt_mat[colnames(seu), ], na.rm = TRUE)
}
seu$pseudotime <- pt_mat[colnames(seu), 1]   # primary lineage shorthand

# ===========================================================================
# 4. Visualise pseudotime on UMAP
# ===========================================================================
pt_cols <- c("pseudotime", paste0("pseudotime_L", seq_len(min(n_lineages, 3))))
pt_cols <- intersect(pt_cols, colnames(seu@meta.data))

fp_list <- lapply(pt_cols, function(col) {
  FeaturePlot(seu, features = col, reduction = "umap", pt.size = 0.2, order = TRUE) +
    scale_colour_viridis_c(option = "plasma", na.value = "grey90") +
    theme_classic(base_size = 9) +
    ggtitle(col)
})
p_pt <- wrap_plots(fp_list, ncol = min(n_lineages + 1, 3)) +
  plot_annotation(title = "Pseudotime — Slingshot (root: stem-like memory cluster)")

ggsave(file.path(plot_dir, "05_pseudotime_umap.pdf"), p_pt,
       width = 12, height = ceiling(length(fp_list) / 3) * 4)
message("Saved: results/plots/05_pseudotime_umap.pdf")

# ===========================================================================
# 5. Trajectory validation — gene expression along pseudotime
# ===========================================================================
# Validation criterion (specified in the analysis plan):
#   Exhaustion genes (TOX, NR4A1, PDCD1, LAG3, HAVCR2, TIGIT) should INCREASE
#   Memory genes    (TCF7, CCR7, IL7R) should DECREASE along pseudotime

validation_genes <- list(
  exhaustion = c("TOX", "NR4A1", "PDCD1", "LAG3", "HAVCR2", "TIGIT"),
  memory     = c("TCF7", "CCR7", "IL7R")
)

# Bin pseudotime into 30 equal-width bins and compute mean expression per bin.
# This gives a smooth trend without needing tradeSeq.
n_bins <- 30L
cells_with_pt <- !is.na(seu$pseudotime)
pt_vals        <- seu$pseudotime[cells_with_pt]
pt_cuts        <- cut(pt_vals, breaks = n_bins, labels = FALSE)

expr_data <- GetAssayData(seu, layer = "data")[
  intersect(unlist(validation_genes), rownames(seu[["RNA"]])),
  cells_with_pt,
  drop = FALSE
]

bin_means <- sapply(sort(unique(pt_cuts)), function(b) {
  cells_in_bin <- which(pt_cuts == b)
  rowMeans(expr_data[, cells_in_bin, drop = FALSE])
})

bin_df_list <- lapply(rownames(bin_means), function(gene) {
  data.frame(
    gene     = gene,
    bin      = seq_len(n_bins),
    mean_exp = bin_means[gene, ],
    group    = if (gene %in% validation_genes$exhaustion) "exhaustion" else "memory"
  )
})
bin_df <- do.call(rbind, bin_df_list)
# Normalise each gene to [0, 1] so different expression scales are comparable.
bin_df <- do.call(rbind, lapply(split(bin_df, bin_df$gene), function(df) {
  rng <- range(df$mean_exp, na.rm = TRUE)
  if (diff(rng) < 1e-10) df$norm_exp <- 0.5 else
    df$norm_exp <- (df$mean_exp - rng[1]) / diff(rng)
  df
}))

group_cols <- c(exhaustion = "#d73027", memory = "#4575b4")

p_valid <- ggplot(bin_df, aes(x = bin, y = norm_exp, colour = group, group = gene)) +
  geom_smooth(se = FALSE, method = "loess", span = 0.5, linewidth = 0.8, alpha = 0.7) +
  scale_colour_manual(values = group_cols) +
  labs(
    title    = "Gene expression along pseudotime (validation)",
    subtitle = "Exhaustion genes should rise; memory genes should fall",
    x        = "Pseudotime bin (1 = root, 30 = tip)",
    y        = "Normalised mean expression",
    colour   = "Gene group"
  ) +
  facet_wrap(~ gene, ncol = 4) +
  theme_classic(base_size = 10)

ggsave(file.path(plot_dir, "05_trajectory_validation.pdf"), p_valid,
       width = 14, height = ceiling(length(unlist(validation_genes)) / 4) * 3.5)
message("Saved: results/plots/05_trajectory_validation.pdf")

# Pearson correlation of each validation gene with pseudotime.
cat("\n--- Gene–pseudotime correlations ---\n")
cor_df <- data.frame(
  gene  = rownames(bin_means),
  r     = apply(expr_data, 1, function(x) cor(x, pt_vals, use = "complete.obs")),
  group = ifelse(rownames(bin_means) %in% validation_genes$exhaustion, "exhaustion", "memory")
)
cor_df <- cor_df[order(cor_df$r, decreasing = TRUE), ]
print(cor_df)
write.csv(cor_df, file.path(table_dir, "05_pseudotime_gene_correlations.csv"), row.names = FALSE)

# ===========================================================================
# 6. Save
# ===========================================================================
saveRDS(list(seu = seu, sds = sds, pt_mat = pt_mat),
        "data/processed/05_trajectory.rds")
message("Saved: data/processed/05_trajectory.rds")
