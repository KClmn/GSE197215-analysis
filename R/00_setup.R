# 00_setup.R
# Install (if missing) and load all packages needed for the GSE197215 pipeline.
# Run once per machine; the requireNamespace checks make re-runs fast.
#
# Input:  nothing
# Output: nothing (side-effect: packages installed into user library)

# ---------------------------------------------------------------------------
# User library
# ---------------------------------------------------------------------------
user_lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(user_lib) && !dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
if (nzchar(user_lib)) .libPaths(c(user_lib, .libPaths()))

# ---------------------------------------------------------------------------
# BiocManager bootstrap
# ---------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = user_lib)

# ---------------------------------------------------------------------------
# CRAN packages
# ---------------------------------------------------------------------------
cran_pkgs <- c(
  "Seurat",         # single-cell framework; object, DimPlot, FeaturePlot, etc.
  "ggplot2",        # plotting backend
  "patchwork",      # multi-panel layout for ggplot figures
  "dplyr",          # metadata manipulation
  "readxl",         # read Table S1 from the supplement
  "viridis",        # perceptually uniform colour scales
  "ggrepel",        # non-overlapping labels on DimPlots
  "Matrix"          # sparse matrix support (Seurat dependency)
)

# ---------------------------------------------------------------------------
# Bioconductor packages
# ---------------------------------------------------------------------------
bioc_pkgs <- c(
  "SingleCellExperiment",  # SCE class; required by slingshot and MOFA2
  "slingshot",             # trajectory inference (05_trajectory.R)
  "tradeSeq",              # test gene–pseudotime association post-slingshot
  "MOFA2",                 # multi-omics factor analysis (06_mofa.R)
  "basilisk",              # isolated Python env; installs mofapy2 automatically for MOFA2
  "variancePartition",     # variance decomposition (07_variancepartition.R)
  "BiocParallel"           # parallel back-end for variancePartition
)

# ---------------------------------------------------------------------------
# GitHub-only packages
# ---------------------------------------------------------------------------
# ProjecTILs (03_projectils.R) has a chain of carmonalab dependencies.
# Install in this order — each one must succeed before the next:
#
#   BiocManager::install("UCell")                       # Bioconductor
#   remotes::install_github("carmonalab/STACAS")        # required by ProjecTILs
#   remotes::install_github("carmonalab/ProjecTILs")
#
# We do NOT auto-install GitHub packages here — a missing package gives an
# informative error from library() without hiding an opaque install step.

# ---------------------------------------------------------------------------
# Install missing CRAN packages
# ---------------------------------------------------------------------------
missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) {
  message("Installing missing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran, lib = user_lib)
}

# ---------------------------------------------------------------------------
# Install missing Bioconductor packages
# ---------------------------------------------------------------------------
missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  message("Installing missing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, lib = user_lib, update = FALSE, ask = FALSE)
}

# ---------------------------------------------------------------------------
# MOFA2 Python back-end — handled automatically via basilisk.
# basilisk (already in bioc_pkgs above) creates an isolated conda environment
# and installs mofapy2 inside it on first use; no manual pip call needed.
# In 06_mofa.R we pass use_basilisk = TRUE to run_mofa() — that is all.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Load all packages to surface any installation problems now
# ---------------------------------------------------------------------------
invisible(lapply(c(cran_pkgs, bioc_pkgs), library, character.only = TRUE))

# Check ProjecTILs separately — it is optional; warn rather than hard-fail.
if (requireNamespace("ProjecTILs", quietly = TRUE)) {
  library(ProjecTILs)
  message("ProjecTILs loaded successfully.")
} else {
  warning("ProjecTILs not found. Run:\n",
          "  remotes::install_github('carmonalab/UCell')\n",
          "  remotes::install_github('carmonalab/ProjecTILs')\n",
          "Then re-source 00_setup.R before running 03_projectils.R.")
}

# ---------------------------------------------------------------------------
sessionInfo()
