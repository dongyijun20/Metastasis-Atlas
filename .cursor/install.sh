#!/usr/bin/env bash
#
# Per-agent install step for the Metastasis-Atlas single-cell R codebase.
#
# Refreshes the R package stack used by the analysis scripts on the `dev`
# (scATAC-seq: ArchR / Signac) and `Feng` (scRNA-seq: Seurat / CellChat /
# GeneNMF) branches. Idempotent: already-installed packages are skipped, so it
# runs fast against a snapshot/build that already contains them.
#
# R itself, the system libraries, and the binary-repo configuration are
# provided by the base environment (see .cursor/setup_system.sh, baked into the
# snapshot). This script only manages R package state.

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "R is not installed. Provision the base environment first:" >&2
  echo "  bash .cursor/setup_system.sh" >&2
  exit 1
fi

# Safety net: if the base image did not hand the site library to this user,
# claim it so package installation can proceed without root.
SITE_LIB="$(Rscript -e 'cat(.Library.site[1])')"
if [ ! -w "$SITE_LIB" ]; then
  sudo chown -R "$(id -u):$(id -g)" "$SITE_LIB" || true
fi

echo "==> R version:"
R --version | head -2

Rscript "$HERE/install_packages.R"

echo "==> install.sh complete"
