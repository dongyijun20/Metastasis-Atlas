#!/usr/bin/env Rscript
#
# Installs the R package stack for the Metastasis-Atlas single-cell scripts.
# Idempotent: packages already present are skipped. CRAN packages come from
# Posit Package Manager as Ubuntu binaries; Bioconductor and a few GitHub-only
# packages are installed on top.
#
# CRAN and Bioconductor repos are combined into a single `repos` set so
# cross-repo dependencies resolve automatically (e.g. Signac -> Rsamtools,
# NMF -> Biobase, which live on Bioconductor even though the packages are on
# CRAN).
#
# The CORE set must install for the environment to be considered healthy; the
# script errors if any core package is not loadable at the end. Heavier
# optional packages (ArchR, CellChat, DoubletFinder, ...) are best-effort so a
# single upstream breakage does not block the whole environment, and any misses
# are reported clearly.

options(
  repos = c(P3M = "https://packagemanager.posit.co/cran/__linux__/noble/latest"),
  # bioconductor.org is intermittently unavailable at bleeding-edge releases;
  # Posit's mirror is a reliable CDN that also hosts annotation data packages.
  BioC_mirror = "https://packagemanager.posit.co/bioconductor",
  Ncpus = max(1L, parallel::detectCores()),
  timeout = 3600,
  warn = 1
)

installed <- function(pkg) requireNamespace(pkg, quietly = TRUE)

## Install any missing packages from the combined CRAN + Bioconductor repos.
ensure <- function(pkgs, label) {
  need <- pkgs[!vapply(pkgs, installed, logical(1))]
  if (!length(need)) {
    message(sprintf("==> [%s] all present", label))
    return(invisible())
  }
  message(sprintf("\n==> [%s] installing: %s", label, paste(need, collapse = ", ")))
  for (p in need) {
    tryCatch(install.packages(p),
             error = function(e) message(sprintf("  !! %s failed: %s", p, conditionMessage(e))))
    if (!installed(p)) message(sprintf("  !! %s still not loadable after install", p))
  }
}

# ---------------------------------------------------------------------------
# Bootstrap: BiocManager + remotes, then fold Bioconductor repos into `repos`
# so install.packages() resolves CRAN and Bioconductor dependencies together.
# ---------------------------------------------------------------------------
ensure(c("BiocManager", "remotes"), "bootstrap")
options(repos = BiocManager::repositories())
message("Active repositories:")
print(getOption("repos"))

# ---------------------------------------------------------------------------
# CORE packages (CRAN + Bioconductor) — backbone of the analysis scripts.
# ---------------------------------------------------------------------------
core_pkgs <- c(
  # CRAN: Seurat/Signac ecosystem, plotting, utilities
  "Seurat", "Signac", "tidyverse", "harmony", "patchwork", "cowplot",
  "ggrepel", "pheatmap", "ggalluvial", "ggsci", "gridExtra", "data.table",
  "Matrix", "hdf5r", "RcppML", "NMF", "msigdbr", "reshape2", "viridis",
  "scales", "uwot", "R.utils", "future", "GeneNMF", "forcats",
  # Bioconductor: genomics, annotation, enrichment, motif analysis
  "GenomicRanges", "GenomeInfoDb", "Biobase", "org.Hs.eg.db",
  "ComplexHeatmap", "fgsea", "UCell", "clusterProfiler",
  "EnsDb.Hsapiens.v86", "BSgenome.Hsapiens.UCSC.hg38", "JASPAR2020",
  "TFBSTools", "motifmatchr", "chromVAR"
)
ensure(core_pkgs, "core")

# ---------------------------------------------------------------------------
# Optional GitHub-only packages (best effort).
# ---------------------------------------------------------------------------
gh_specs <- c(
  DoubletFinder = "chris-mcginnis-ucsf/DoubletFinder",
  presto        = "immunogenomics/presto",
  CellChat      = "jinworks/CellChat",
  ArchR         = "GreenleafLab/ArchR"
)
for (nm in names(gh_specs)) {
  if (installed(nm)) {
    message(sprintf("==> [github] %s present", nm))
    next
  }
  message(sprintf("\n==> [github] installing %s (%s)", nm, gh_specs[[nm]]))
  tryCatch(remotes::install_github(gh_specs[[nm]], upgrade = "never", dependencies = TRUE),
           error = function(e) message(sprintf("  !! %s failed: %s", nm, conditionMessage(e))))
}

# ---------------------------------------------------------------------------
# Final health check (recomputed after every tier is installed).
# ---------------------------------------------------------------------------
miss_core <- core_pkgs[!vapply(core_pkgs, installed, logical(1))]
miss_gh   <- names(gh_specs)[!vapply(names(gh_specs), installed, logical(1))]

message("\n================ INSTALL SUMMARY ================")
message("core missing:              ", if (length(miss_core)) paste(miss_core, collapse = ", ") else "none")
message("optional (github) missing: ", if (length(miss_gh)) paste(miss_gh, collapse = ", ") else "none")
message("=================================================")

if (length(miss_core)) {
  stop("Core packages failed to install: ", paste(miss_core, collapse = ", "))
}
message("Core package set is complete.")
