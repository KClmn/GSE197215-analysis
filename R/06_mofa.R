# 06_mofa.R
# Unsupervised multi-modal factor analysis of RNA + ADT with MOFA2.
#
# Goal: discover latent factors that explain co-variation across the two
# modalities WITHOUT using response labels — then correlate each factor with
# known clinical/biological variables to characterise it, and FLAG factors
# that explain substantial variance but do NOT correlate with any known variable.
# Those unexplained factors are the primary leads.
#
# MOFA2 requires:
#   R: BiocManager::install(c("MOFA2", "basilisk"))
#   Python: handled automatically — use_basilisk = TRUE in run_mofa() creates
#           an isolated conda environment and installs mofapy2 on first run.
#           No manual pip install needed (works out of the box in RStudio).
#
# If MOFA2 itself is not installed, the script falls back to multi-modal PCA.
#
# Input:  data/processed/05_trajectory.rds   (or 04_scored.rds if 05 not run)
# Output: data/processed/06_mofa_model.hdf5
#         data/processed/06_mofa.rds
#         results/plots/06_mofa_variance.pdf
#         results/plots/06_mofa_factor_correlations.pdf
#         results/plots/06_mofa_dimplot.pdf
#         results/tables/06_mofa_factor_correlations.csv

library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
traj_data <- if (file.exists("data/processed/05_trajectory.rds")) {
  readRDS("data/processed/05_trajectory.rds")
} else {
  warning("05_trajectory.rds not found — falling back to 04_scored.rds")
  NULL
}
seu <- if (!is.null(traj_data)) traj_data$seu else readRDS("data/processed/04_scored.rds")
cat("Loaded:", ncol(seu), "cells\n")

plot_dir  <- "results/plots";  if (!dir.exists(plot_dir))  dir.create(plot_dir, recursive = TRUE)
table_dir <- "results/tables"; if (!dir.exists(table_dir)) dir.create(table_dir, recursive = TRUE)

# ===========================================================================
# 1. Prepare input matrices for MOFA2
# ===========================================================================
# MOFA2 expects a list of views: each view is a features × cells matrix.
# CHOICE: use the top 2000 highly variable genes for RNA (too many features
# slow training without improving factors) and all ADT features (only 14).
# Cells are the same in both views — MOFA2 handles any missing entries as NA.

DefaultAssay(seu) <- "RNA"

# Identify highly variable genes if not already computed.
if (!"vst.variance.standardized" %in% colnames(seu[["RNA"]]@meta.features)) {
  cat("Computing HVGs...\n")
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
}
hvg <- VariableFeatures(seu)
cat("HVGs selected:", length(hvg), "\n")

# RNA view: scaled data of HVGs (features × cells)
# CHOICE: use scaled data (zero-centred, unit variance) so MOFA2 factors are
# not dominated by the highest-expressing genes.  MOFA2's Gaussian likelihood
# is appropriate for scaled continuous data.
rna_mat <- GetAssayData(seu, assay = "RNA", layer = "scale.data")
if (nrow(rna_mat) == 0) {
  cat("scale.data slot empty — scaling now...\n")
  seu     <- ScaleData(seu, features = hvg)
  rna_mat <- GetAssayData(seu, assay = "RNA", layer = "scale.data")
}
rna_mat <- rna_mat[intersect(hvg, rownames(rna_mat)), ]
cat("RNA view dimensions:", nrow(rna_mat), "×", ncol(rna_mat), "\n")

# ADT view: use ADT_renorm if available (CLR/DSB normalised), else ADT.
# CHOICE: Gaussian likelihood on CLR-normalised ADT values is standard for
# MOFA2 CITE-seq analyses (Argelaguet et al. 2020 Genome Biology).
adt_assay_name <- if ("ADT_renorm" %in% names(seu@assays)) "ADT_renorm" else "ADT"
adt_mat <- GetAssayData(seu, assay = adt_assay_name, layer = "data")
cat("ADT view dimensions:", nrow(adt_mat), "×", ncol(adt_mat), "\n")

# Guard: cell names must align between views.
shared_cells <- intersect(colnames(rna_mat), colnames(adt_mat))
rna_mat <- rna_mat[, shared_cells]
adt_mat <- adt_mat[, shared_cells]
cat("Shared cells:", length(shared_cells), "\n")

mofa_views <- list(RNA = as.matrix(rna_mat), ADT = as.matrix(adt_mat))

# ===========================================================================
# 2. Run MOFA2 (or fall back to PCA)
# ===========================================================================
# use_basilisk = TRUE means MOFA2 manages its own Python env — no need to
# pre-check for mofapy2; the only hard requirement is the R package itself.
mofa_available <- requireNamespace("MOFA2", quietly = TRUE)

if (mofa_available) {
  # Point reticulate at OSCAR's system Python (which has OpenSSL and mofapy2).
  # basilisk fails on OSCAR because pyenv can't compile Python without OpenSSL headers.
  # Run in terminal first:
  #   module load python && pip install --user mofapy2
  # Then set MOFA2_PYTHON to the output of `which python3`.
  py_bin <- Sys.getenv("MOFA2_PYTHON", unset = "")
  if (nzchar(py_bin)) {
    reticulate::use_python(py_bin, required = TRUE)
    cat("Using Python:", py_bin, "\n")
  }

  library(MOFA2)

  cat("\nRunning MOFA2...\n")

  mofa_obj <- create_mofa(mofa_views)

  # CHOICE: 15 factors with default data and model options.
  # 15 factors is generous for 12 donors; residual factors after ~5–7 meaningful
  # ones provide a natural noise floor to calibrate variance thresholds.
  data_opts  <- get_default_data_options(mofa_obj)
  model_opts <- get_default_model_options(mofa_obj)
  model_opts$num_factors <- 15L

  train_opts <- get_default_training_options(mofa_obj)
  train_opts$seed             <- 42L
  train_opts$maxiter          <- 1000L
  train_opts$convergence_mode <- "medium"

  mofa_obj <- prepare_mofa(
    mofa_obj,
    data_options     = data_opts,
    model_options    = model_opts,
    training_options = train_opts
  )

  hdf5_path <- "data/processed/06_mofa_model.hdf5"
  time_mofa <- system.time({
    use_bas  <- !nzchar(Sys.getenv("MOFA2_PYTHON", unset = ""))
    mofa_obj <- run_mofa(mofa_obj, outfile = hdf5_path, use_basilisk = use_bas)
  })
  cat("MOFA2 run time:\n"); print(time_mofa)

  # Extract factor values (cells × factors)
  factor_mat <- get_factors(mofa_obj, factors = "all")[[1]]
  cat("Factor matrix dimensions:", nrow(factor_mat), "×", ncol(factor_mat), "\n")

  # Add factors to Seurat metadata for downstream use
  stopifnot(all(rownames(factor_mat) %in% colnames(seu)))
  factor_df   <- as.data.frame(factor_mat[colnames(seu), ])
  colnames(factor_df) <- paste0("MOFA_F", seq_len(ncol(factor_df)))
  seu <- AddMetaData(seu, factor_df)

  # Variance explained per factor per view
  var_exp <- calculate_variance_explained(mofa_obj)
  cat("\n--- Variance explained per factor per view ---\n")
  print(var_exp$r2_per_factor)

  # ===========================================================================
  # 3. Plots — variance explained
  # ===========================================================================
  p_var <- plot_variance_explained(mofa_obj, plot_total = TRUE) &
    theme_classic(base_size = 10)
  ggsave(file.path(plot_dir, "06_mofa_variance.pdf"), p_var, width = 10, height = 6)
  message("Saved: results/plots/06_mofa_variance.pdf")

  # ===========================================================================
  # 4. Correlate factors with known variables
  # ===========================================================================
  known_vars_cont <- c(
    "exhaustion_score", "th2_score", "cytotox_score", "memory_score",
    "S.Score", "G2M.Score", "pseudotime", "nCount_RNA", "age_years",
    "time_to_relapse_days"
  )
  known_vars_cont <- intersect(known_vars_cont, colnames(seu@meta.data))

  meta_sub <- seu@meta.data[colnames(seu), ]

  # Pearson r between each factor and each continuous variable.
  cor_list <- lapply(colnames(factor_df), function(fac) {
    lapply(known_vars_cont, function(var) {
      vals <- meta_sub[[var]]
      if (sum(!is.na(vals)) < 10) return(data.frame(factor = fac, variable = var, r = NA_real_, p = NA_real_))
      ct <- cor.test(factor_df[[fac]], vals, use = "complete.obs")
      data.frame(factor = fac, variable = var, r = ct$estimate, p = ct$p.value)
    })
  })
  cor_df <- do.call(rbind, do.call(c, cor_list))

  # Eta-squared for categorical variables (response, patient_id, Phase).
  cat_vars <- c("response", "patient_id", "Phase", "sort_frac")
  cat_vars <- intersect(cat_vars, colnames(seu@meta.data))
  eta_list <- lapply(colnames(factor_df), function(fac) {
    lapply(cat_vars, function(var) {
      vals <- meta_sub[[var]]
      if (length(unique(na.omit(vals))) < 2) return(NULL)
      fit  <- lm(factor_df[[fac]] ~ vals)
      ss   <- anova(fit)
      eta2 <- ss["vals", "Sum Sq"] / sum(ss[, "Sum Sq"])
      data.frame(factor = fac, variable = var, r = sqrt(eta2), p = NA_real_)
    })
  })
  eta_df <- do.call(rbind, Filter(Negate(is.null), do.call(c, eta_list)))

  cor_all <- rbind(cor_df, eta_df)
  write.csv(cor_all, file.path(table_dir, "06_mofa_factor_correlations.csv"), row.names = FALSE)

  # Heatmap of |r| — correlation magnitude
  cor_wide <- reshape(cor_all[, c("factor", "variable", "r")],
                      idvar = "factor", timevar = "variable", direction = "wide")
  rownames(cor_wide) <- cor_wide$factor
  cor_wide$factor    <- NULL
  colnames(cor_wide) <- sub("^r\\.", "", colnames(cor_wide))
  cor_mat <- abs(as.matrix(cor_wide))

  # Identify "mystery factors" — high variance but low correlation with everything.
  # CHOICE: threshold = 0.3 for |r|.  Any factor whose maximum |r| across all
  # known variables is < 0.3 is flagged as unexplained.
  max_r      <- apply(cor_mat, 1, max, na.rm = TRUE)
  total_var  <- rowSums(var_exp$r2_per_factor, na.rm = TRUE)
  mystery_df <- data.frame(
    factor    = rownames(cor_mat),
    max_corr  = round(max_r, 3),
    total_var_pct = round(total_var * 100, 1)
  )
  mystery_df <- mystery_df[order(-mystery_df$total_var_pct), ]
  flagged    <- mystery_df[mystery_df$max_corr < 0.3, ]

  cat("\n=== Mystery factors (total var explained > 0% AND max |r| < 0.3) ===\n")
  print(flagged)
  cat("These factors explain variance NOT captured by known variables — primary leads.\n")

  write.csv(mystery_df, file.path(table_dir, "06_mofa_mystery_factors.csv"), row.names = FALSE)

  # Heatmap plot
  cor_df_long <- cor_all
  cor_df_long$abs_r <- abs(cor_df_long$r)

  p_heat <- ggplot(cor_df_long, aes(x = variable, y = factor, fill = abs_r)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = ifelse(!is.na(r), round(r, 2), "")),
              size = 2.5, colour = "black") +
    scale_fill_viridis_c(option = "magma", limits = c(0, 1), na.value = "grey90",
                         name = "|r|") +
    labs(title = "MOFA2 factor × known variable correlations", x = NULL, y = NULL) +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(plot_dir, "06_mofa_factor_correlations.pdf"), p_heat,
         width = max(8, ncol(cor_mat) * 0.7 + 3),
         height = max(5, nrow(cor_mat) * 0.5 + 2))
  message("Saved: results/plots/06_mofa_factor_correlations.pdf")

  # DimPlot coloured by top factors on the existing UMAP
  top_facs <- colnames(factor_df)[1:min(6, ncol(factor_df))]
  fp_fac   <- lapply(top_facs, function(f) {
    FeaturePlot(seu, features = f, reduction = "umap", pt.size = 0.2, order = TRUE) +
      scale_colour_gradientn(colours = c("#4575b4","grey90","#d73027")) +
      theme_classic(base_size = 9) + ggtitle(f)
  })
  p_fac <- wrap_plots(fp_fac, ncol = 3) +
    plot_annotation(title = "Top MOFA2 factors on UMAP")
  ggsave(file.path(plot_dir, "06_mofa_dimplot.pdf"), p_fac, width = 15, height = 8)
  message("Saved: results/plots/06_mofa_dimplot.pdf")

  saveRDS(list(seu = seu, mofa_obj = mofa_obj, cor_all = cor_all,
               mystery_df = mystery_df),
          "data/processed/06_mofa.rds")

} else {
  # ---------------------------------------------------------------------------
  # Fallback: PCA on combined RNA + ADT using the integrated assay
  # ---------------------------------------------------------------------------
  warning("MOFA2 not installed. Running multi-modal PCA fallback.\n",
          "Install: BiocManager::install(c('MOFA2', 'basilisk'))  then re-run.")
  library(stats)

  cat("\nRunning multi-modal PCA fallback...\n")

  # Scale ADT matrix to unit variance before concatenating.
  adt_scaled <- t(scale(t(as.matrix(adt_mat))))
  rna_sub    <- rna_mat[, shared_cells]

  combined_mat <- rbind(rna_sub, adt_scaled)
  cat("Combined matrix for PCA:", nrow(combined_mat), "×", ncol(combined_mat), "\n")

  set.seed(42)
  pca_res  <- irlba::irlba(t(combined_mat), nv = 20)
  pc_mat   <- pca_res$u
  rownames(pc_mat) <- shared_cells
  colnames(pc_mat) <- paste0("PC", seq_len(ncol(pc_mat)))

  # Variance explained
  pct_var <- (pca_res$d^2 / sum(pca_res$d^2)) * 100
  cat("\nVariance explained per PC:\n")
  print(round(pct_var, 2))

  # Add to Seurat metadata
  pc_df   <- as.data.frame(pc_mat[colnames(seu), ])
  colnames(pc_df) <- paste0("MOFA_F", seq_len(ncol(pc_df)))   # same naming convention
  seu     <- AddMetaData(seu, pc_df)

  saveRDS(list(seu = seu, pca_fallback = pca_res),
          "data/processed/06_mofa.rds")
  cat("Saved fallback PCA results to data/processed/06_mofa.rds\n")
}

message("Saved: data/processed/06_mofa.rds")
