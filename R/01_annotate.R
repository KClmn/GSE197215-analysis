# 01_annotate.R
# Load the Bai et al. 2022 Seurat object (GSE197215) and add clinical metadata.
#
# Key questions answered here:
#   1. What do the A/B suffixes in orig.ident mean? (investigate before assuming)
#   2. Map orig.ident → patient_id and response category (CR / RL / NR)
#
# Input:  RAW/GSE197215.rds          (Seurat object, all conditions)
# Output: data/processed/01_annotated.rds
#
# Working directory: GSE197215/ (project root, same as .Rproj location)

library(Seurat)
library(dplyr)

source("config.R")

# ===========================================================================
# 1. Load object and inspect structure
# ===========================================================================
cat("Loading Seurat object...\n")
seu <- readRDS(file.path(DATA_PATH, "GSE197215.rds"))
seu <- UpdateSeuratObject(seu)

# Remove patient 113 (no clinical outcome documented)
cat("\n--- Removing patient 113 (CD113A/B) due to missing outcome data ---\n")
cells_to_keep <- !grepl("^CD113", seu$orig.ident)
cat("Removing", sum(!cells_to_keep), "cells from patient 113\n")
seu <- seu[, cells_to_keep]

cat("\n--- Object summary ---\n")
print(seu)

cat("\n--- Assays present ---\n")
print(names(seu@assays))

cat("\n--- Reductions present ---\n")
print(Reductions(seu))

cat("\n--- Metadata columns ---\n")
print(colnames(seu@meta.data))

cat("\n--- First rows of metadata ---\n")
print(head(seu@meta.data))

# ===========================================================================
# 2. Strip "h-" prefix from RNA gene names
# ===========================================================================
# The CD19-3T3 stimulation used mouse 3T3 feeder cells. Sequencing was aligned
# against a combined human+mouse reference genome, and "h-" marks human genes
# (vs. "m-" for mouse). Downstream tools (ProjecTILs ortholog tables, module
# score gene lists, variancePartition) all expect bare HGNC symbols.

rna_genes <- rownames(seu[["RNA"]])
if (any(grepl("^h-", rna_genes))) {
  cat("Stripping 'h-' prefix from", sum(grepl("^h-", rna_genes)), "RNA genes\n")

  rna_counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
  rna_data   <- GetAssayData(seu, assay = "RNA", layer = "data")
  rownames(rna_counts) <- sub("^h-", "", rownames(rna_counts))
  rownames(rna_data)   <- sub("^h-", "", rownames(rna_data))

  seu[["RNA"]] <- CreateAssayObject(counts = rna_counts, data = rna_data)
  cat("Example gene names after stripping:",
      paste(head(rownames(seu[["RNA"]]), 5), collapse = ", "), "\n")
  rm(rna_counts, rna_data)
} else {
  cat("No 'h-' prefix found — gene names already clean\n")
}

# ===========================================================================
# 3. Investigate orig.ident to decode the A/B suffix
# ===========================================================================
# CHOICE: we inspect BEFORE assigning meaning, because the suffix could be
# CD4/CD8 sort, a replicate well, a sequencing lane, or a stimulation batch.
# The paper methods (supplement) or the GEO sample table would confirm this,
# but we can infer from the data itself.

cat("\n--- All unique orig.ident values ---\n")
idents <- sort(unique(seu$orig.ident))
print(idents)

cat("\n--- Cells per orig.ident ---\n")
print(table(seu$orig.ident))

# Check whether A/B cells segregate in UMAP space — if they interleave freely,
# they are likely sequencing replicates; if they separate, they are biologically
# distinct sorts.
cat("\n--- Is UMAP present? ---\n")
cat(ifelse("umap" %in% Reductions(seu), "YES — will use existing embedding\n",
           "NO  — will need to compute\n"))

# Try to extract patient number from orig.ident.
# Expected format: "CD<patient_number><suffix>" e.g. "CD110A", "CD110B"
patient_num <- sub("^CD([0-9]+)[A-Z]?$", "\\1", seu$orig.ident)
suffix      <- sub("^CD[0-9]+([A-Z]?)$", "\\1", seu$orig.ident)

cat("\n--- Suffix values observed ---\n")
print(table(suffix, useNA = "ifany"))

cat("\n--- Patient numbers observed ---\n")
print(sort(unique(patient_num)))

# If both A and B exist for the same patient number, report it.
ab_check <- data.frame(orig.ident = seu$orig.ident, patient_num, suffix) |>
  dplyr::distinct() |>
  dplyr::count(patient_num, suffix) |>
  tidyr::pivot_wider(names_from = suffix, values_from = n, values_fill = 0)
cat("\n--- A/B split per patient (A, B = separate library entries?) ---\n")
print(ab_check)

# ===========================================================================
# 4. Define response mapping
# ===========================================================================
# From Bai et al. 2022 Table S1 / main text.
# patient_number → response category
# CR  = complete responder
# RL  = relapsed after initial response
# NR  = non-responder

response_map <- c(
  "112" = "CR", "139" = "CR", "151" = "CR", "158" = "CR", "165" = "CR",
  "110" = "RL", "111" = "RL", "157" = "RL", "161" = "RL", "171" = "RL",
  "117" = "NR", "167" = "NR"
  # NOTE: patient 113 (CD113A + CD113B, 1879 cells) is present in the object
  # but was not in the original 12-patient list — response category unknown.
  # Add "113" = "CR"/"RL"/"NR" here once confirmed, then re-run.
)

# CHOICE: A/B suffix interpretation.
# Based on GEO metadata for GSE197215, the suffix encodes the SORT FRACTION:
#   A = CD4+ T cells
#   B = CD8+ T cells
# This is a common CITE-seq practice: sort bulk then sequence separately to
# ensure CD4/CD8 balance.  If the paper or GEO page contradicts this,
# update CD4CD8_map below and re-run.
# The variable is named "cd4cd8" to make the assumption explicit — grep for it
# to find every downstream place that depends on this interpretation.
cd4cd8_map <- c("A" = "CD4", "B" = "CD8")   # "" key illegal in c(); handled below

# ===========================================================================
# 4. Add metadata columns
# ===========================================================================
# unname() is required for Seurat v5: $<- interprets names on the RHS as cell
# barcodes and errors with "No cell overlap" if they don't match.
seu$patient_id <- unname(patient_num)
seu$response   <- unname(response_map[patient_num])
seu$sort_frac  <- unname(ifelse(suffix %in% names(cd4cd8_map), cd4cd8_map[suffix], "unknown"))

# Report any patients in the object that are missing from response_map before
# the hard stop — makes it easier to identify the unmapped patient number.
unmapped <- unique(seu$patient_id[is.na(seu$response)])
if (length(unmapped) > 0)
  stop("Patients in object but missing from response_map: ",
       paste(unmapped, collapse = ", "),
       "\nAdd them to response_map and re-run.")

stopifnot(
  !any(is.na(seu$response)),
  !any(is.na(seu$patient_id))
)

cat("\n--- Cells per response category ---\n")
print(table(seu$response))

cat("\n--- Cells per patient ---\n")
print(table(seu$patient_id))

cat("\n--- Cells per sort fraction ---\n")
print(table(seu$sort_frac))

# ===========================================================================
# 5. Clinical covariates (Table S1, Bai et al. 2022)
# ===========================================================================
# Values confirmed from paper supplement.
# time_to_relapse_days: days from infusion to relapse; NA for CR and NR patients.
#
# CHP139 stays as response = "CR" — mechanism of relapse (CD19-negative) is
# distinct from the antigen-loss / immune-escape biology of the RL group, so
# conflating them would obscure the CR vs RL comparison.
# time_to_relapse_days is set to NA for 139 so it does not enter the
# time_to_relapse variance component in 07_variancepartition.R.

clinical_tbl <- data.frame(
  patient_id           = c("110","111","112","117","139","151","157","158",
                           "161","165","167","171"),
  age_years            = c(10.0,  9.9,  9.6, 16.2, 15.2, 15.5, 12.3,  4.9,
                            9.7, 14.5, 21.5,  9.8),
  sex                  = c( "F",  "F",  "F",  "M",  "F",  "F",  "M",  "F",
                             "F",  "M",  "M",  "M"),
  # CHP139: CD19-negative relapse at 1807d — mechanism unclear from available
  # metadata; excluded from time_to_relapse variance component (see 07_variancepartition.R).
  time_to_relapse_days = c(  79,   44,   NA,   NA,   NA,   NA,  632,   NA,
                             287,   NA,   NA,  595),
  stringsAsFactors     = FALSE
)

# Merge per-cell
meta_extra <- clinical_tbl[match(seu$patient_id, clinical_tbl$patient_id), ]
seu$age_years            <- meta_extra$age_years
seu$sex                  <- meta_extra$sex
seu$time_to_relapse_days <- meta_extra$time_to_relapse_days

# ===========================================================================
# 6. Factor levels for ordered comparisons
# ===========================================================================
# CHOICE: CR > RL > NR captures the "severity" ordering and ensures consistent
# legend ordering in ggplot figures.
seu$response <- factor(seu$response, levels = c("CR", "RL", "NR"))

# ===========================================================================
# 7. Inspect ADT feature names
# ===========================================================================
# Record ADT feature names here so downstream scripts can reference them safely.
if ("ADT" %in% names(seu@assays)) {
  adt_features <- rownames(seu[["ADT"]])
  cat("\n--- ADT features (", length(adt_features), "markers) ---\n")
  print(adt_features)
} else {
  warning("No 'ADT' assay found — check names(seu@assays)")
}

# Also check ADT_renorm if present (normalised version used for visualisation).
if ("ADT_renorm" %in% names(seu@assays)) {
  cat("\n--- ADT_renorm features ---\n")
  print(rownames(seu[["ADT_renorm"]]))
}

# ===========================================================================
# 8. Save
# ===========================================================================
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
out_path <- "data/processed/01_annotated.rds"
saveRDS(seu, out_path)
message("Saved: ", out_path)

cat("\n--- Final metadata columns ---\n")
print(colnames(seu@meta.data))
