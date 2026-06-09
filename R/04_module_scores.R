# 04_module_scores.R
# Compute module scores for exhaustion, TH2 polarisation, and cell cycle.
#
# Scores are added as new metadata columns so all downstream analyses
# (trajectory, MOFA, variancePartition) can use them without re-running.
#
# Input:  data/processed/03_projectils.rds   (or 01_annotated.rds if 03 not run)
# Output: data/processed/04_scored.rds
#         results/plots/04_module_scores.pdf

library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)

# ---------------------------------------------------------------------------
# Load — fall back gracefully if ProjecTILs step was skipped
# ---------------------------------------------------------------------------
input_path <- if (file.exists("data/processed/03_projectils.rds")) {
  "data/processed/03_projectils.rds"
} else {
  warning("03_projectils.rds not found — loading 01_annotated.rds instead.")
  "data/processed/01_annotated.rds"
}
seu <- readRDS(input_path)
DefaultAssay(seu) <- "RNA"
cat("Loaded:", ncol(seu), "cells from", input_path, "\n")

plot_dir <- "results/plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ===========================================================================
# 1. Exhaustion / dysfunction module score
# ===========================================================================
# SOURCE: Good CR, Aznar MA, Kuramitsu S, et al. "An NK-like CAR T cell
#   transition in CAR T cell dysfunction." Cell. 2021;184(24):6081-6100.e26.
#
# The manuscript Supplementary Table S4 lists the top up-regulated genes in
# dysfunctional CAR-T cells relative to functional effectors.
# IMPORTANT: replace the gene list below with the exact genes from Table S4
# of Good et al. 2021 if you have access to the paper's supplement.
# The list here is the canonical CAR-T exhaustion gene set (Fraietta, Lynn,
# Good) used consistently in the Bai et al. 2022 analysis.
#
# CHOICE: nbin = 25 (Seurat default). The score is the average expression of
# genes in the set minus the average of a random background set sampled from
# bins of equal expression — so it is less confounded by library size than
# raw sum expression.

exhaustion_genes <- c(
  # Core transcription factors of exhaustion
  "TOX", "TOX2", "NR4A1", "NR4A2", "NR4A3",
  # Inhibitory receptors
  "PDCD1",    # PD-1
  "LAG3",
  "HAVCR2",   # TIM-3
  "TIGIT",
  "CTLA4",
  "CD244",    # 2B4
  "ENTPD1",   # CD39
  # Transcription factors downstream of exhaustion
  "BATF", "IRF4", "PRDM1",
  # NK-like transition genes (Good et al. 2021 specific)
  "KLRG1", "KLRB1", "KLRC1", "NCAM1",
  # Additional dysfunction markers
  "EOMES", "TBX21"
)

# Keep only genes present in this object's RNA assay.
all_rna_genes   <- rownames(seu[["RNA"]])
exh_genes_found <- intersect(exhaustion_genes, all_rna_genes)
exh_genes_miss  <- setdiff(exhaustion_genes, all_rna_genes)
cat("\n--- Exhaustion genes found:", length(exh_genes_found), "---\n")
cat("Missing from object:", paste(exh_genes_miss, collapse = ", "), "\n")

seu <- AddModuleScore(
  seu,
  features = list(exh_genes_found),
  name     = "exhaustion_score",
  nbin     = 25,
  seed     = 42
)
# Rename the auto-suffixed column (AddModuleScore appends "1")
colnames(seu@meta.data)[colnames(seu@meta.data) == "exhaustion_score1"] <- "exhaustion_score"

cat("Exhaustion score range:", range(seu$exhaustion_score), "\n")

# ===========================================================================
# 2. TH2 module score
# ===========================================================================
# TH2 polarisation in CAR-T infusion products has been associated with
# poor response in ALL (Bai et al. 2022 examined this directly).
# Canonical TH2 genes: cytokines + master regulator.

th2_genes_all   <- c("IL4", "IL5", "IL13", "GATA3", "IL10", "IL33")
th2_genes_found <- intersect(th2_genes_all, all_rna_genes)
th2_genes_miss  <- setdiff(th2_genes_all, all_rna_genes)
cat("\n--- TH2 genes found:", length(th2_genes_found), "---\n")
cat("Missing:", paste(th2_genes_miss, collapse = ", "), "\n")

if (length(th2_genes_found) >= 2) {
  seu <- AddModuleScore(
    seu,
    features = list(th2_genes_found),
    name     = "th2_score",
    nbin     = 25,
    seed     = 42
  )
  colnames(seu@meta.data)[colnames(seu@meta.data) == "th2_score1"] <- "th2_score"
  cat("TH2 score range:", range(seu$th2_score), "\n")
} else {
  warning("Too few TH2 genes found (< 2) — th2_score not computed.")
  seu$th2_score <- NA_real_
}

# ===========================================================================
# 3. Effector / cytotoxicity module score
# ===========================================================================
# Included as a contrast to exhaustion — useful for the MOFA factor correlation.

cytotox_genes_all   <- c("GZMB", "GZMA", "GZMH", "GZMK", "PRF1", "GNLY",
                          "IFNG", "TNF", "FASLG", "NKG7")
cytotox_genes_found <- intersect(cytotox_genes_all, all_rna_genes)
cat("\n--- Cytotoxicity genes found:", length(cytotox_genes_found), "---\n")

if (length(cytotox_genes_found) >= 2) {
  seu <- AddModuleScore(
    seu,
    features = list(cytotox_genes_found),
    name     = "cytotox_score",
    nbin     = 25,
    seed     = 42
  )
  colnames(seu@meta.data)[colnames(seu@meta.data) == "cytotox_score1"] <- "cytotox_score"
}

# ===========================================================================
# 4. Memory module score
# ===========================================================================
memory_genes_all   <- c("TCF7", "CCR7", "IL7R", "SELL",  # CCR7 = CD197, SELL = CD62L
                        "LEF1", "KLF2", "BCL6", "FoxO1")
memory_genes_found <- intersect(memory_genes_all, all_rna_genes)
cat("\n--- Memory genes found:", length(memory_genes_found), "---\n")

if (length(memory_genes_found) >= 2) {
  seu <- AddModuleScore(
    seu,
    features = list(memory_genes_found),
    name     = "memory_score",
    nbin     = 25,
    seed     = 42
  )
  colnames(seu@meta.data)[colnames(seu@meta.data) == "memory_score1"] <- "memory_score"
}

# ===========================================================================
# 5. Cell cycle scoring (Seurat built-in)
# ===========================================================================
# Uses Seurat's built-in cc.genes (Tirosh et al. 2016 Science) to compute
# S-phase and G2M-phase scores, then assigns each cell to G1, S, or G2M.
# CHOICE: cells with high S or G2M scores are cycling and may confound
# trajectory and MOFA analyses — their phase is recorded so we can regress
# it out or subset it away in downstream steps.

data("cc.genes", package = "Seurat")   # loads s.genes and g2m.genes
s_genes   <- intersect(cc.genes$s.genes,   all_rna_genes)
g2m_genes <- intersect(cc.genes$g2m.genes, all_rna_genes)
cat("\n--- Cell cycle genes found: S =", length(s_genes),
    ", G2M =", length(g2m_genes), "---\n")

seu <- CellCycleScoring(
  seu,
  s.features   = s_genes,
  g2m.features = g2m_genes,
  set.ident    = FALSE   # don't overwrite Idents — we want response/cluster Idents
)
cat("Cell cycle phase distribution:\n")
print(table(seu$Phase))

# ===========================================================================
# 6. Visualise all scores
# ===========================================================================
response_cols <- c(CR = "#2ca02c", RL = "#ff7f0e", NR = "#d62728")
score_cols    <- c("exhaustion_score", "th2_score", "cytotox_score",
                   "memory_score", "S.Score", "G2M.Score")
score_cols    <- score_cols[score_cols %in% colnames(seu@meta.data)]

# FeaturePlots on UMAP
fp_list <- lapply(score_cols, function(sc) {
  FeaturePlot(seu, features = sc, reduction = "umap", pt.size = 0.2, order = TRUE) +
    scale_colour_gradientn(colours = c("grey90", "#fee08b", "#d73027")) +
    theme_classic(base_size = 9) +
    ggtitle(sc)
})
p_fp <- wrap_plots(fp_list, ncol = 3) +
  plot_annotation(title = "Module scores on UMAP")
ggsave(file.path(plot_dir, "04_module_scores_umap.pdf"), p_fp,
       width = 15, height = ceiling(length(score_cols) / 3) * 4)

# Violin plots split by response
vln_list <- lapply(score_cols, function(sc) {
  VlnPlot(seu, features = sc, group.by = "response",
          cols = response_cols, pt.size = 0, log = FALSE) +
    theme(axis.title.x = element_blank(), plot.title = element_text(size = 9)) +
    ggtitle(sc)
})
p_vln <- wrap_plots(vln_list, ncol = 3) +
  plot_annotation(title = "Module scores by response category")
ggsave(file.path(plot_dir, "04_module_scores_violin.pdf"), p_vln,
       width = 12, height = ceiling(length(score_cols) / 3) * 4)

message("Saved: results/plots/04_module_scores_umap.pdf")
message("Saved: results/plots/04_module_scores_violin.pdf")

# ===========================================================================
# 7. Save
# ===========================================================================
saveRDS(seu, "data/processed/04_scored.rds")
message("Saved: data/processed/04_scored.rds")

cat("\n--- New metadata columns ---\n")
print(score_cols)
