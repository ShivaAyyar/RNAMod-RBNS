#!/bin/bash
#SBATCH --job-name=download_rbns
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=/mnt/vstor/SOM_CCCC_JGS25/sayyar/RNAMod-RBNS/logs/download_rbns_%j.out
#SBATCH --error=/mnt/vstor/SOM_CCCC_JGS25/sayyar/RNAMod-RBNS/logs/download_rbns_%j.err

# ============================================================================
# Automated RBNS 5-mer Enrichment Download Script for RNAMod-RBNS Pipeline
#
# Downloads pre-computed 5-mer R-value enrichment files from ENCODE for
# the 15 RBPs that have both RBNS and eCLIP data available.
#
# Data source: ENCODE Portal (RNA Bind-n-Seq experiments)
# Output: TSV files with R-values across protein concentrations
#
# Usage:
#   bash scripts/download_rbns.sh      # Run interactively
#   sbatch scripts/download_rbns.sh    # Submit as SLURM job
# ============================================================================

set -euo pipefail

# ==========================
# Configuration
# ==========================

BASE_DIR="/mnt/vstor/SOM_CCCC_JGS25/sayyar/RNAMod-RBNS"
DATA_DIR="${BASE_DIR}/data"
RBNS_DIR="${DATA_DIR}/rbns"
ENCODE_BASE="https://www.encodeproject.org"

# ==========================
# Helper Functions
# ==========================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

download_rbns_file() {
    local rbp=$1
    local accession=$2
    local output="${RBNS_DIR}/${rbp}_enrichment.tsv"

    if [[ -f "${output}" && -s "${output}" ]]; then
        log "EXISTS: ${rbp} ($(wc -l < "${output}") k-mers)"
        return 0
    fi

    local url="${ENCODE_BASE}/files/${accession}/@@download/${accession}.tsv"

    log "Downloading: ${rbp} (${accession})..."

    if curl -sL -o "${output}" "${url}" 2>/dev/null; then
        if [[ -s "${output}" ]]; then
            local lines=$(wc -l < "${output}")
            log "SUCCESS: ${rbp} (${lines} lines)"
            return 0
        else
            log "ERROR: Downloaded file is empty for ${rbp}"
            rm -f "${output}"
            return 1
        fi
    else
        log "ERROR: Failed to download ${rbp}"
        rm -f "${output}"
        return 1
    fi
}

# ==========================
# Create Directory
# ==========================

mkdir -p "${RBNS_DIR}"
mkdir -p "${BASE_DIR}/logs"

log "=============================================="
log "Automated RBNS 5-mer Enrichment Download"
log "=============================================="
log "Output directory: ${RBNS_DIR}"
log ""

# ==========================
# RBNS 5-mer Enrichment File Accessions
# ==========================
# Pre-verified from ENCODE: output_type=enrichment, 5-mer files
# Source: experiment_report_2026_2_16_23h_10m.tsv + ENCODE API verification
#
# These are the 15 RBPs with BOTH eCLIP and RBNS data available.
# 8 RBPs from original list lack RBNS data: PTBP1, QKI, SRSF1, MATR3,
# HNRNPA1, HNRNPM, NONO, U2AF2

declare -A RBNS_FILES

# Verified 5-mer enrichment file accessions from ENCODE API
# Each file contains R-values across protein concentrations for ~970 5-mers

# Positive controls (m6A readers)
RBNS_FILES["IGF2BP1"]="ENCFF157MJN"
RBNS_FILES["IGF2BP2"]="ENCFF446BRC"

# Other RBPs with both RBNS and eCLIP data
RBNS_FILES["HNRNPC"]="ENCFF592FMJ"
RBNS_FILES["TIA1"]="ENCFF936GUV"
RBNS_FILES["HNRNPK"]="ENCFF715DIL"
RBNS_FILES["PCBP2"]="ENCFF761GIK"
RBNS_FILES["RBFOX2"]="ENCFF002DFE"
RBNS_FILES["TARDBP"]="ENCFF022PDN"
RBNS_FILES["SRSF9"]="ENCFF518QQO"
RBNS_FILES["RBM22"]="ENCFF688EDA"
RBNS_FILES["TRA2A"]="ENCFF191NLM"
RBNS_FILES["HNRNPL"]="ENCFF422OAC"
RBNS_FILES["LIN28B"]="ENCFF189NSX"
RBNS_FILES["FUS"]="ENCFF835DNX"
RBNS_FILES["EWSR1"]="ENCFF100OEJ"

# ==========================
# Download All Files
# ==========================

success_count=0
fail_count=0

log "Starting downloads..."
log ""

for rbp in "${!RBNS_FILES[@]}"; do
    accession="${RBNS_FILES[$rbp]}"
    if download_rbns_file "${rbp}" "${accession}"; then
        ((++success_count))
    else
        ((++fail_count))
    fi
done

# ==========================
# Summary
# ==========================

log ""
log "=============================================="
log "Download Complete"
log "=============================================="
log "Successful: ${success_count}"
log "Failed: ${fail_count}"
log ""

# List downloaded files
log "Downloaded files:"
for f in "${RBNS_DIR}"/*.tsv; do
    if [[ -f "$f" ]]; then
        lines=$(wc -l < "$f")
        log "  $(basename "$f"): ${lines} lines"
    fi
done

log ""
log "Next step: Run scripts/process_rbns_enrichment.py to convert R-values to Z-scores"
log ""

# Note about RBPs without RBNS data
log "=============================================="
log "NOTE: RBPs without RBNS data (cannot analyze):"
log "=============================================="
log "  - PTBP1 (only PTBP3 in RBNS)"
log "  - QKI"
log "  - SRSF1"
log "  - MATR3"
log "  - HNRNPA1 (only HNRNPA1L2 in RBNS)"
log "  - HNRNPM"
log "  - NONO"
log "  - U2AF2"
log ""
log "Total analyzable RBPs: 15"
