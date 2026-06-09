# 02_visualize.R
# UMAP visualisations: response category overview + all 14 ADT marker FeaturePlots.
#
# Input:  data/processed/01_annotated.rds
# Output: results/plots/02_dimplot_response.pdf
#         results/plots/02_featureplot_adt.pdf
#         results/plots/02_violins_exhaustion.pdf
#         results/plots/02_violins_memory.pdf

library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
seu <- readRDS("data/processed/01_annotated.rds")

# Run clustering on the integrated PCA if not already present.
# CHOICE: resolution 0.5 is a moderate starting point for ~28k cells; gives
# roughly 10-20 clusters. Adjust and re-run 02 if clusters look over/under-split.
if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
  cat("seurat_clusters not found — running FindNeighbors + FindClusters on integrated PCA\n")
  DefaultAssay(seu) <- "integrated"
  seu <- FindNeighbors(seu, reduction = "pca", dims = 1:30)
  seu <- FindClusters(seu, resolution = 0.5)
  # Save clusters back so downstream scripts don't repeat this step.
  saveRDS(seu, "data/processed/01_annotated.rds")
  cat("Clusters added and saved back to 01_annotated.rds\n")
}

# Confirm the UMAP exists; stop with a helpful message if missing.
if (!"umap" %in% tolower(Reductions(seu)))
  stop("No UMAP reduction found. Reductions present: ",
       paste(Reductions(seu), collapse = ", "),
       "\nRun Seurat's RunUMAP() or check the object structure in 01_annotate.R.")

# Use the ADT_renorm assay for protein visualisation if present; fall back to ADT.
# ADT_renorm is the DSB- or CLR-normalised version and has better visual contrast.
adt_assay <- if ("ADT_renorm" %in% names(seu@assays)) "ADT_renorm" else "ADT"
cat("Using assay '", adt_assay, "' for ADT FeaturePlots\n", sep = "")

adt_features <- rownames(seu[[adt_assay]])
cat("ADT features:\n"); print(adt_features)

# ---------------------------------------------------------------------------
# Colour palette for response categories
# ---------------------------------------------------------------------------
# CHOICE: green/orange/red maps intuitively to CR/RL/NR clinical severity.
response_cols <- c(CR = "#2ca02c", RL = "#ff7f0e", NR = "#d62728")

plot_dir <- "results/plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ===========================================================================
# 1. DimPlot — response category
# ===========================================================================
p_response <- DimPlot(
  seu,
  reduction = "umap",
  group.by  = "response",
  cols      = response_cols,
  pt.size   = 0.3,
  alpha     = 0.5
) +
  ggtitle("Response category") +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")

p_patient <- DimPlot(
  seu,
  reduction = "umap",
  group.by  = "patient_id",
  pt.size   = 0.3,
  alpha     = 0.4
) +
  ggtitle("Donor (patient_id)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")

p_sort <- DimPlot(
  seu,
  reduction = "umap",
  group.by  = "sort_frac",
  pt.size   = 0.3,
  alpha     = 0.5,
  cols      = c(CD4 = "#4393c3", CD8 = "#d6604d", unknown = "grey70")
) +
  ggtitle("Sort fraction (CD4 / CD8)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")

# CHOICE: split_by response gives the clearest view of cluster composition
# differences — each panel has the same axis limits so spatial differences
# are directly comparable.
p_split <- DimPlot(
  seu,
  reduction = "umap",
  group.by  = "seurat_clusters",
  split.by  = "response",
  pt.size   = 0.3,
  alpha     = 0.5,
  ncol      = 3
) +
  ggtitle("Clusters split by response") +
  theme_classic(base_size = 10)

layout_overview <- (p_response | p_patient | p_sort) / p_split +
  plot_annotation(
    title   = "Bai et al. 2022 - CD19-3T3 stimulation condition",
    caption = "GSE197215"
  )

ggsave(
  file.path(plot_dir, "02_dimplot_response.pdf"),
  layout_overview,
  width = 18, height = 12
)
message("Saved: results/plots/02_dimplot_response.pdf")

# ===========================================================================
# 2. FeaturePlots — all 14 ADT markers
# ===========================================================================
# CHOICE: use the viridis colour scale — it is perceptually uniform and prints
# well in greyscale, both of which matter for publication figures.
# FeaturePlot returns one ggplot per feature; we stitch them with patchwork.

DefaultAssay(seu) <- adt_assay

# ADT feature names confirmed from 01_annotate.R output.
# "ADT-" (empty name) is a ghost barcode — excluded from all plots.
exhaustion_markers <- intersect(
  c("ADT-PD-1", "ADT-LAG-3", "ADT-TIGIT", "ADT-CTLA-4", "ADT-2B4"),
  adt_features
)
memory_markers <- intersect(
  c("ADT-CCR7", "ADT-CD45RO", "ADT-CD45RA", "ADT-CD62L", "ADT-CD27"),
  adt_features
)
# Exclude the empty ghost barcode from all panels
other_markers <- setdiff(adt_features, c(exhaustion_markers, memory_markers, "ADT-"))

cat("\nExhaustion ADT features found: ", paste(exhaustion_markers, collapse = ", "), "\n")
cat("Memory ADT features found: ",      paste(memory_markers,     collapse = ", "), "\n")
cat("Other ADT features: ",             paste(other_markers,      collapse = ", "), "\n")

make_feature_plots <- function(features, assay, title_prefix = "") {
  plots <- lapply(features, function(f) {
    FeaturePlot(
      seu,
      features  = f,
      reduction = "umap",
      pt.size   = 0.2,
      order     = TRUE   # plot highest-expressing cells on top for clarity
    ) +
      scale_colour_viridis_c(option = "magma") +
      theme_classic(base_size = 9) +
      theme(
        legend.key.height = unit(0.4, "cm"),
        plot.title = element_text(size = 9, face = "bold")
      ) +
      ggtitle(paste0(title_prefix, f))
  })
  plots
}

# All 14 ADT markers in one figure (4 columns)
all_plots <- make_feature_plots(adt_features, adt_assay)
n_cols    <- min(4L, length(all_plots))
p_adt     <- wrap_plots(all_plots, ncol = n_cols) +
  plot_annotation(
    title   = paste("ADT markers -", adt_assay, "assay"),
    caption = "GSE197215 | Bai et al. 2022"
  )

ggsave(
  file.path(plot_dir, "02_featureplot_adt.pdf"),
  p_adt,
  width  = n_cols * 3.5,
  height = ceiling(length(all_plots) / n_cols) * 3.2
)
message("Saved: results/plots/02_featureplot_adt.pdf")

# ===========================================================================
# 3. Violin plots — exhaustion and memory markers by response
# ===========================================================================
# VlnPlots split by response give a per-marker statistical overview that
# complements the spatial UMAP view.

# "ADT-" ghost barcode already excluded from exhaustion/memory/other lists above.
DefaultAssay(seu) <- adt_assay

if (length(exhaustion_markers) > 0) {
  p_vln_exh <- VlnPlot(
    seu,
    features  = exhaustion_markers,
    group.by  = "response",
    cols      = response_cols,
    pt.size   = 0,
    ncol      = length(exhaustion_markers)
  ) &
    theme(axis.title.x = element_blank(), plot.title = element_text(size = 9))

  ggsave(
    file.path(plot_dir, "02_violins_exhaustion.pdf"),
    p_vln_exh,
    width  = length(exhaustion_markers) * 2.5,
    height = 4
  )
  message("Saved: results/plots/02_violins_exhaustion.pdf")
}

if (length(memory_markers) > 0) {
  p_vln_mem <- VlnPlot(
    seu,
    features  = memory_markers,
    group.by  = "response",
    cols      = response_cols,
    pt.size   = 0,
    ncol      = length(memory_markers)
  ) &
    theme(axis.title.x = element_blank(), plot.title = element_text(size = 9))

  ggsave(
    file.path(plot_dir, "02_violins_memory.pdf"),
    p_vln_mem,
    width  = length(memory_markers) * 2.5,
    height = 4
  )
  message("Saved: results/plots/02_violins_memory.pdf")
}

# Reset to RNA for downstream scripts
DefaultAssay(seu) <- "RNA"
