# config.R
# Machine-specific paths. source() this at the top of any script that reads
# raw data. Intermediate outputs (data/processed/) use relative paths and are
# not affected by this file.

if (Sys.info()["sysname"] == "Darwin") {
  # Local Mac
  DATA_PATH <- "~/scratch/GSE197215-analysis/RAW/"
} else {
  # OSCAR (Linux)
  DATA_PATH <- "/users/kclin/GSE197215-analysis/RAW/"
}
