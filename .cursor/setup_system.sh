#!/usr/bin/env bash
#
# Base-environment bootstrap for the Metastasis-Atlas single-cell R codebase.
#
# Installs R and the system libraries required to build the Seurat / Signac /
# ArchR / CellChat package stack, makes the R site-library writable by the
# agent user, and configures fast binary package repositories. This provisions
# the *base image*; it is slow and stable, so it belongs in the environment
# snapshot/build rather than the per-boot install step.
#
# Idempotent and safe to re-run.

set -euxo pipefail
export DEBIAN_FRONTEND=noninteractive

# --- CRAN apt repository for a current R on Ubuntu 24.04 (noble) ------------
sudo apt-get update -y
sudo apt-get install -y --no-install-recommends \
  ca-certificates gnupg curl wget software-properties-common dirmngr

sudo install -d -m 0755 /etc/apt/keyrings
if [ ! -s /etc/apt/keyrings/cran.gpg ]; then
  wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
    | sudo gpg --dearmor -o /etc/apt/keyrings/cran.gpg
fi
echo "deb [signed-by=/etc/apt/keyrings/cran.gpg] https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/" \
  | sudo tee /etc/apt/sources.list.d/cran.list >/dev/null

sudo apt-get update -y
sudo apt-get install -y --no-install-recommends r-base r-base-dev

# --- System libraries needed to build the single-cell R package stack -------
sudo apt-get install -y --no-install-recommends \
  build-essential gfortran git pkg-config \
  libcurl4-openssl-dev libssl-dev libxml2-dev \
  libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
  libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
  libhdf5-dev libglpk-dev libgit2-dev \
  libgsl-dev libmagick++-dev \
  libgeos-dev libgdal-dev libudunits2-dev libproj-dev \
  libbz2-dev liblzma-dev zlib1g-dev libncurses-dev \
  libcairo2-dev libxt-dev cmake poppler-utils

# --- Let the agent user install R packages into the site library -----------
sudo chown -R "$(id -u):$(id -g)" "$(Rscript -e 'cat(.Library.site[1])')"

# --- Configure fast binary package repositories for all R sessions ----------
RPROFILE_SITE="$(R RHOME)/etc/Rprofile.site"
MARKER="# >>> metastasis-atlas package repos >>>"
if ! grep -qF "$MARKER" "$RPROFILE_SITE" 2>/dev/null; then
  sudo tee -a "$RPROFILE_SITE" >/dev/null <<EOF

$MARKER
local({
  options(
    # Precompiled binary CRAN packages for Ubuntu 24.04 (noble).
    repos = c(P3M = "https://packagemanager.posit.co/cran/__linux__/noble/latest"),
    # Reliable Bioconductor mirror (bioconductor.org can be flaky at new releases).
    BioC_mirror = "https://packagemanager.posit.co/bioconductor",
    HTTPUserAgent = sprintf(
      "R/%s R (%s)",
      getRversion(),
      paste(getRversion(), R.version[["platform"]], R.version[["arch"]], R.version[["os"]])
    ),
    Ncpus = max(1L, parallel::detectCores())
  )
})
# <<< metastasis-atlas package repos <<<
EOF
fi

R --version | head -2
echo "=== setup_system.sh complete ==="
