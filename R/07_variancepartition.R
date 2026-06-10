# 07_variancepartition.R
# Variance decomposition: what fraction of pseudotime and exhaustion score
# variance is explained by donor, response, ADT markers, RNA markers,
# module scores, and technical covariates?
#
# variancePartition fits a mixed model per gene (or per cell-level score)
# and partitions R² among predictors.  With only 12 donors, donor is a random
# effect (not a fixed effect) — this is the correct treatment for a small cohort.
#
# IMPORTANT POWER NOTE:
#   n = 12 donors; per-group n = 5 (CR), 5 (RL), 2 (NR).
#   The NR group (n = 2) is severely underpowered for any between-group test.
#   Variance components here are descriptive, not confirmatory.
#   Treat Residuals as your primary output — it bounds unexplained structure.
#
# Input:  data/processed/06_mofa.rds   (or 04_scored.rds if 06 not run)
# Output: data/processed/07_vp.rds
#         results/plots/07_vp_pseudotime.pdf
#         results/plots/07_vp_exhaustion.pdf
#         results/tables/07_vp_pseudotime.csv
#         results/tables/07_vp_exhaustion.csv
#         results/tables/07_vp_top_genes.csv

library(Seurat)
library(variancePartition)
library(BiocParallel)
library(ggplot2)
library(patchwork)

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
mofa_data <- if (file.exists("data/processed/06_mofa.rds")) {
  readRDS("data/processed/06_mofa.rds")
} else {
  warning("06_mofa.rds not found — loading 04_scored.rds")
  NULL
}
seu <- if (!is.null(mofa_data)) mofa_data$seu else readRDS("data/processed/04_scored.rds")
cat("Loaded:", ncol(seu), "cells\n")

# Pseudotime must exist for the primary analysis.
if (!"pseudotime" %in% colnames(seu@meta.data)) {
  if (file.exists("data/processed/05_trajectory.rds")) {
    traj <- readRDS("data/processed/05_trajectory.rds")
    seu$pseudotime <- traj$seu$pseudotime
    cat("Pseudotime loaded from 05_trajectory.rds\n")
  } else {
    warning("Pseudotime not found — pseudotime variance decomposition will be skipped.")
  }
}

plot_dir  <- "results/plots";  if (!dir.exists(plot_dir))  dir.create(plot_dir,  recursive = TRUE)
table_dir <- "results/tables"; if (!dir.exists(table_dir)) dir.create(table_dir, recursive = TRUE)

# ---------------------------------------------------------------------------
# Parallel back-end (conservative: 4 cores; adjust to available CPUs)
# ---------------------------------------------------------------------------
# CHOICE: SnowParam for compatibility across Mac/Linux.  Set workers = 1 for
# debugging or if memory is tight (variancePartition can be memory-hungry with
# many features and predictors).
n_workers <- min(4L, parallel::detectCores() - 1L)
register(SnowParam(workers = n_workers, type = "SOCK"))
cat("Using", n_workers, "parallel workers\n")

# ===========================================================================
# Helper: build the formula and run variancePartition on a response vector
# ===========================================================================
# variancePartition for a single response (pseudotime or exhaustion_score)
# and a full RNA gene expression matrix.
#
# Formula:
#   ~ (1|patient_id)     random effect — accounts for donor differences without
#                        consuming 11 degrees of freedom on 12 donors
#   + response           fixed effect (CR/RL/NR)
#   + sort_frac          CD4/CD8 sort (fixed)
#   + nCount_RNA         technical covariate (sequencing depth)
#   + S.Score + G2M.Score cell cycle scores
#   + ADT markers        continuous protein expression
#   + key RNA markers    as continuous covariates (normalised expression)

build_meta <- function(seu, extra_continuous = character(0)) {
  adt_assay_name <- if ("ADT_renorm" %in% names(seu@assays)) "ADT_renorm" else "ADT"
  adt_df         <- as.data.frame(t(as.matrix(
    GetAssayData(seu, assay = adt_assay_name, layer = "data")
  )))
  colnames(adt_df) <- paste0("ADT_", make.names(colnames(adt_df)))

  key_rna_genes <- c("TCF7", "TOX", "NR4A1", "GZMB", "PRF1", "PDCD1",
                     "LAG3", "HAVCR2", "CCR7", "IL7R")
  key_rna_genes <- intersect(key_rna_genes, rownames(seu[["RNA"]]))
  rna_df <- as.data.frame(t(as.matrix(
    GetAssayData(seu, assay = "RNA", layer = "data")[key_rna_genes, , drop = FALSE]
  )))
  colnames(rna_df) <- paste0("RNA_", make.names(colnames(rna_df)))

  # Core metadata
  # time_to_relapse_days is included here so it is available in the RL-only
  # model (section 3 below); it will be NA for CR/NR cells and is excluded
  # from the full-cohort formula by the >50% NA filter below.
  core_cols <- c("patient_id", "response", "sort_frac",
                 "nCount_RNA", "S.Score", "G2M.Score",
                 "exhaustion_score", "th2_score", "memory_score",
                 "pseudotime", "time_to_relapse_days")
  core_cols <- intersect(core_cols, colnames(seu@meta.data))
  meta <- seu@meta.data[, core_cols, drop = FALSE]

  # Drop columns with >50% NA (insufficient data for model fitting).
  ok_cols <- colnames(meta)[colMeans(is.na(meta)) < 0.5]
  meta    <- meta[, ok_cols, drop = FALSE]

  cbind(meta, adt_df, rna_df)
}

# ===========================================================================
# 1. Variance decomposition across RNA genes
# ===========================================================================
# Run variancePartition on a set of HVGs (faster than all genes, more
# informative than a single score).  This tells us which predictor drives
# gene-expression variance genome-wide.

DefaultAssay(seu) <- "RNA"
hvg <- VariableFeatures(seu)
if (length(hvg) == 0) {
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  hvg <- VariableFeatures(seu)
}
hvg <- hvg[seq_len(min(1000L, length(hvg)))]   # cap at 1000 for speed
cat("Running variancePartition on", length(hvg), "HVGs...\n")

expr_mat <- GetAssayData(seu, assay = "RNA", layer = "data")[hvg, ]

meta_full <- build_meta(seu)

# Drop any predictor column with zero variance (would cause model failure).
numeric_cols  <- colnames(meta_full)[sapply(meta_full, is.numeric)]
zero_var_cols <- numeric_cols[sapply(meta_full[numeric_cols], function(x) var(x, na.rm = TRUE) < 1e-10)]
if (length(zero_var_cols) > 0) {
  cat("Dropping zero-variance predictors:", paste(zero_var_cols, collapse = ", "), "\n")
  meta_full <- meta_full[, !colnames(meta_full) %in% zero_var_cols, drop = FALSE]
}

# Scale all continuous predictors to mean=0, sd=1 so lmer is numerically
# stable across the large range difference between nCount_RNA (~10k) and
# cell cycle scores (~0). Categorical columns are left as-is.
numeric_cols <- colnames(meta_full)[sapply(meta_full, is.numeric)]
meta_full[numeric_cols] <- lapply(meta_full[numeric_cols], function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s < 1e-10) x else as.numeric(scale(x))
})

# Build the formula dynamically from whatever columns are present.
# Random effects: patient_id
# Fixed effects: everything else
rand_effects  <- c("patient_id")
fixed_effects <- setdiff(colnames(meta_full),
                         c(rand_effects, "pseudotime", "exhaustion_score",
                           "memory_score", "th2_score"))

rand_terms  <- paste(sprintf("(1|%s)", rand_effects), collapse = " + ")
fixed_terms <- if (length(fixed_effects) > 0) paste(fixed_effects, collapse = " + ") else ""
form_str    <- if (nzchar(fixed_terms)) paste("~", rand_terms, "+", fixed_terms) else paste("~", rand_terms)
form        <- as.formula(form_str)
cat("\nVariance partition formula:\n", deparse(form), "\n\n")

# Align cells across expr_mat and meta_full.
shared_cells <- intersect(colnames(expr_mat), rownames(meta_full))
expr_mat     <- expr_mat[, shared_cells]
meta_full    <- meta_full[shared_cells, ]

time_vp <- system.time({
  vp_rna <- fitExtractVarPartModel(
    exprObj  = expr_mat,
    formula  = form,
    data     = meta_full,
    BPPARAM  = bpparam()
  )
})
cat("variancePartition (RNA) run time:\n"); print(time_vp)

cat("\n--- Median variance fractions across HVGs ---\n")
print(sort(apply(vp_rna, 2, median), decreasing = TRUE))

# Top genes by residual variance (least explained)
top_resid <- head(vp_rna[order(-vp_rna$Residuals), ], 20)
cat("\n--- Top 20 genes by Residuals (most unexplained structure) ---\n")
print(round(top_resid, 3))
write.csv(as.data.frame(vp_rna), file.path(table_dir, "07_vp_rna_full.csv"))

# ===========================================================================
# 2. Variance decomposition of pseudotime score (single response)
# ===========================================================================
run_vp_single <- function(response_col, label) {
  if (!response_col %in% colnames(meta_full)) {
    cat("Skipping", label, "— column not found in metadata.\n")
    return(invisible(NULL))
  }
  y    <- meta_full[[response_col]]
  ok   <- !is.na(y)
  if (sum(ok) < 20) {
    cat("Skipping", label, "— fewer than 20 non-NA values.\n")
    return(invisible(NULL))
  }

  # Build formula excluding the response column itself as a predictor.
  excl  <- c(rand_effects, response_col, "pseudotime", "exhaustion_score",
             "memory_score", "th2_score")
  fixed <- setdiff(colnames(meta_full), excl)
  rand_t  <- paste(sprintf("(1|%s)", rand_effects), collapse = " + ")
  fixed_t <- if (length(fixed) > 0) paste(fixed, collapse = " + ") else ""
  f_str   <- if (nzchar(fixed_t)) paste("~", rand_t, "+", fixed_t) else paste("~", rand_t)
  f       <- as.formula(f_str)

  # variancePartition on a single numeric response uses fitVarPartModel.
  # Wrap response as a 1-row matrix (features = 1).
  y_mat <- matrix(y[ok], nrow = 1, dimnames = list(response_col, names(y)[ok]))
  meta_ok <- meta_full[ok, ]

  vp_single <- tryCatch(
    fitExtractVarPartModel(y_mat, f, meta_ok, BPPARAM = SerialParam()),
    error = function(e) { cat("Error in", label, ":", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(vp_single)) return(invisible(NULL))

  cat("\n---", label, "variance components ---\n")
  print(round(sort(as.numeric(vp_single), decreasing = TRUE), 3))

  # Bar plot
  vp_df <- data.frame(
    predictor = colnames(vp_single),
    fraction  = as.numeric(vp_single)
  )
  vp_df <- vp_df[order(-vp_df$fraction), ]
  vp_df$predictor <- factor(vp_df$predictor, levels = vp_df$predictor)

  p_bar <- ggplot(vp_df, aes(x = predictor, y = fraction * 100, fill = predictor)) +
    geom_col() +
    scale_fill_manual(values = c(
      Residuals   = "grey60",
      patient_id  = "#8c564b",
      response    = "#2ca02c"
    ), na.value = "#1f77b4") +
    labs(
      title = paste("Variance explained in", label),
      x     = NULL, y = "% variance explained"
    ) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    geom_hline(yintercept = 0, linewidth = 0.3)

  fname_pdf <- file.path(plot_dir, paste0("07_vp_", gsub("[^a-z0-9]", "_", tolower(label)), ".pdf"))
  fname_csv <- file.path(table_dir, paste0("07_vp_", gsub("[^a-z0-9]", "_", tolower(label)), ".csv"))
  ggsave(fname_pdf, p_bar, width = 8, height = 5)
  write.csv(vp_df, fname_csv, row.names = FALSE)
  message("Saved: ", fname_pdf)

  cat("\nPOWER NOTE: n =", sum(ok), "cells from 12 donors.\n",
      "Donor-level effects (patient_id) rely on only 12 levels — treat as rough estimates.\n")

  invisible(vp_single)
}

vp_pt  <- run_vp_single("pseudotime",       "Pseudotime")
vp_exh <- run_vp_single("exhaustion_score", "Exhaustion score")

# ===========================================================================
# 3. time_to_relapse_days — RL patients only
# ===========================================================================
# CHP139 is response = "CR" with an atypical CD19-negative relapse event.
# Its time_to_relapse_days is NA in the metadata (set in 01_annotate.R) so it
# cannot enter this model; all other CR and NR cells have NA by definition.
# Restricting to RL cells means:
#   - only 5 RL donors remain (110, 111, 157, 161, 171)
#   - patient_id random effect has only 5 levels — treat output as descriptive only.
#
# What this model asks: within the RL group, does how quickly a patient relapsed
# (early = 44d, late = 632d) predict cell-level variation in pseudotime,
# exhaustion, ADT markers, or RNA markers, above and beyond donor identity?

vp_rl <- local({
  if (!"time_to_relapse_days" %in% colnames(seu@meta.data)) {
    cat("Skipping RL-only time_to_relapse model — column absent.\n")
    return(invisible(NULL))
  }

  rl_cells <- colnames(seu)[!is.na(seu$response) & seu$response == "RL" &
                              !is.na(seu$time_to_relapse_days)]
  cat("\nRL cells with time_to_relapse_days:", length(rl_cells),
      "from", length(unique(seu$patient_id[rl_cells])), "donors\n")

  if (length(rl_cells) < 50) {
    cat("Skipping RL-only model — fewer than 50 eligible cells.\n")
    return(invisible(NULL))
  }

  seu_rl     <- seu[, rl_cells]
  hvg_rl     <- intersect(hvg, rownames(seu_rl[["RNA"]]))
  expr_rl    <- GetAssayData(seu_rl, assay = "RNA", layer = "data")[hvg_rl, ]
  meta_rl    <- build_meta(seu_rl)

  # time_to_relapse_days is now a meaningful continuous predictor for all cells.
  # Drop zero-variance and >50%-NA columns as before.
  ok_cols_rl <- colnames(meta_rl)[colMeans(is.na(meta_rl)) < 0.5]
  meta_rl    <- meta_rl[, ok_cols_rl, drop = FALSE]
  num_rl     <- colnames(meta_rl)[sapply(meta_rl, is.numeric)]
  zv_rl      <- num_rl[sapply(meta_rl[num_rl], function(x) var(x, na.rm = TRUE) < 1e-10)]
  if (length(zv_rl)) meta_rl <- meta_rl[, !colnames(meta_rl) %in% zv_rl, drop = FALSE]

  # response is constant (all RL) — remove it from the formula.
  excl_rl   <- c("response", "pseudotime", "exhaustion_score", "memory_score", "th2_score")
  fixed_rl  <- setdiff(colnames(meta_rl), c("patient_id", excl_rl))
  rand_t    <- "(1|patient_id)"
  fixed_t   <- if (length(fixed_rl) > 0) paste(fixed_rl, collapse = " + ") else ""
  f_rl      <- as.formula(if (nzchar(fixed_t)) paste("~", rand_t, "+", fixed_t)
                          else paste("~", rand_t))
  cat("RL-only formula:\n", deparse(f_rl), "\n")

  shared_rl <- intersect(colnames(expr_rl), rownames(meta_rl))
  expr_rl   <- expr_rl[, shared_rl]
  meta_rl   <- meta_rl[shared_rl, ]

  vp_rl_fit <- tryCatch(
    fitExtractVarPartModel(expr_rl, f_rl, meta_rl, BPPARAM = bpparam()),
    error = function(e) { cat("Error in RL model:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(vp_rl_fit)) return(invisible(NULL))

  cat("\n--- RL-only median variance fractions ---\n")
  print(sort(apply(vp_rl_fit, 2, median), decreasing = TRUE))

  p_vp_rl <- plotVarPart(vp_rl_fit) +
    labs(title = "Variance explained per gene — RL patients only",
         subtitle = "n = 5 donors; time_to_relapse_days as continuous predictor") +
    theme_classic(base_size = 10)
  ggsave(file.path(plot_dir, "07_vp_rl_time_to_relapse.pdf"), p_vp_rl, width = 10, height = 5)
  write.csv(as.data.frame(vp_rl_fit),
            file.path(table_dir, "07_vp_rl_time_to_relapse.csv"))
  message("Saved: results/plots/07_vp_rl_time_to_relapse.pdf")

  cat("\nPOWER NOTE: n = 5 RL donors only. time_to_relapse_days variance\n",
      "component is descriptive — do not over-interpret effect sizes.\n")

  invisible(vp_rl_fit)
})
vp_rl_time <- vp_rl

# ===========================================================================
# 3. Violin / candle plot of RNA-level variance components
# ===========================================================================
p_vp_rna <- plotVarPart(vp_rna) +
  labs(title = "Variance explained per gene (n = 1000 HVGs)") +
  theme_classic(base_size = 10)

ggsave(file.path(plot_dir, "07_vp_rna_violin.pdf"), p_vp_rna, width = 10, height = 5)
message("Saved: results/plots/07_vp_rna_violin.pdf")

# Genes with the highest residuals — hardest to explain with known predictors.
cat("\n=== SUMMARY: top 10 genes with most residual variance ===\n")
cat("These are candidates for novel regulatory axes.\n")
top_resid_genes <- head(rownames(vp_rna)[order(-vp_rna$Residuals)], 10)
print(round(vp_rna[top_resid_genes, c("patient_id", "response", "Residuals")], 3))

write.csv(as.data.frame(vp_rna[top_resid_genes, ]),
          file.path(table_dir, "07_vp_top_residual_genes.csv"))

# ===========================================================================
# 4. Save
# ===========================================================================
saveRDS(
  list(vp_rna = vp_rna, vp_pseudotime = vp_pt, vp_exhaustion = vp_exh,
       vp_rl_time_to_relapse = vp_rl_time),
  "data/processed/07_vp.rds"
)
message("Saved: data/processed/07_vp.rds")

cat("\n=== POWER NOTE (final reminder) ===\n")
cat("n = 12 donors: 5 CR, 5 RL, 2 NR.\n")
cat("NR group (n = 2) cannot support any statistical comparison.\n")
cat("All variance components should be interpreted as exploratory.\n")
cat("The Residuals component is your primary estimate of unexplained structure.\n")
