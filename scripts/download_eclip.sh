#!/bin/bash
#SBATCH --job-name=download_eclip
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=/mnt/vstor/SOM_CCCC_JGS25/sayyar/RNAMod-RBNS/logs/download_eclip_%j.out
#SBATCH --error=/mnt/vstor/SOM_CCCC_JGS25/sayyar/RNAMod-RBNS/logs/download_eclip_%j.err

# ============================================================================
# Automated eCLIP Download Script for RNAMod-RBNS Pipeline
#
# Downloads IDR-thresholded eCLIP peaks from ENCODE for all 23 RBPs.
# Uses pre-filtered ENCODE file accessions (GRCh38, merged replicates, released).
#
# Usage:
#   bash scripts/download_eclip.sh      # Run interactively
#   sbatch scripts/download_eclip.sh    # Submit as SLURM job
# ============================================================================

set -euo pipefail

# ==========================
# Configuration
# ==========================

BASE_DIR="/mnt/vstor/SOM_CCCC_JGS25/sayyar/RNAMod-RBNS"
DATA_DIR="${BASE_DIR}/data"
ECLIP_DIR="${DATA_DIR}/eclip"
ENCODE_BASE="https://www.encodeproject.org"

# ==========================
# Helper Functions
# ==========================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

download_eclip() {
    local rbp=$1
    local cell_line=$2
    local accession=$3
    local output="${ECLIP_DIR}/${rbp}_${cell_line}.bed"

    if [[ -f "${output}" && -s "${output}" ]]; then
        log "EXISTS: ${rbp} ${cell_line} ($(wc -l < "${output}") peaks)"
        return 0
    fi

    local url="${ENCODE_BASE}/files/${accession}/@@download/${accession}.bed.gz"
    local temp_file="${output}.gz"

    log "Downloading: ${rbp} ${cell_line} (${accession})..."

    if wget -q --show-progress -O "${temp_file}" "${url}" 2>/dev/null; then
        gunzip -f "${temp_file}"
        if [[ -s "${output}" ]]; then
            log "SUCCESS: ${rbp} ${cell_line} ($(wc -l < "${output}") peaks)"
            return 0
        else
            log "ERROR: Downloaded file is empty for ${rbp} ${cell_line}"
            rm -f "${output}"
            return 1
        fi
    else
        log "ERROR: Failed to download ${rbp} ${cell_line}"
        rm -f "${temp_file}"
        return 1
    fi
}

# ==========================
# Create Directory
# ==========================

mkdir -p "${ECLIP_DIR}"
mkdir -p "${BASE_DIR}/logs"

log "=============================================="
log "Automated eCLIP Download from ENCODE"
log "=============================================="
log "Output directory: ${ECLIP_DIR}"
log ""

# ==========================
# eCLIP File Accessions
# ==========================
# Pre-filtered from ENCODE: GRCh38, IDR peaks (merged bio reps 1,2), released
# Source: file_report_2026_2_16_4h_30m.tsv
#
# Format: RBP CELL_LINE ACCESSION

declare -A ECLIP_FILES

# K562 files (22 RBPs - PCBP2 only has HepG2)
ECLIP_FILES["IGF2BP1_K562"]="ENCFF650LMV"
ECLIP_FILES["IGF2BP2_K562"]="ENCFF524ZZB"
ECLIP_FILES["HNRNPC_K562"]="ENCFF167CDB"
ECLIP_FILES["TIA1_K562"]="ENCFF918KMT"
ECLIP_FILES["HNRNPK_K562"]="ENCFF318RSO"
# PCBP2_K562 not available (HepG2 only)
ECLIP_FILES["RBFOX2_K562"]="ENCFF206RIM"
ECLIP_FILES["PTBP1_K562"]="ENCFF907HNN"
ECLIP_FILES["TARDBP_K562"]="ENCFF037TVC"
ECLIP_FILES["QKI_K562"]="ENCFF786UOW"
ECLIP_FILES["SRSF1_K562"]="ENCFF223KVR"
ECLIP_FILES["SRSF9_K562"]="ENCFF781BNS"
ECLIP_FILES["RBM22_K562"]="ENCFF972ZMJ"
ECLIP_FILES["TRA2A_K562"]="ENCFF726PFJ"
ECLIP_FILES["HNRNPL_K562"]="ENCFF917CBK"
ECLIP_FILES["LIN28B_K562"]="ENCFF061XNA"
ECLIP_FILES["FUS_K562"]="ENCFF861KMV"
ECLIP_FILES["MATR3_K562"]="ENCFF246EPM"
ECLIP_FILES["HNRNPA1_K562"]="ENCFF392AEV"
ECLIP_FILES["HNRNPM_K562"]="ENCFF445ENC"
ECLIP_FILES["NONO_K562"]="ENCFF730QRI"
ECLIP_FILES["U2AF2_K562"]="ENCFF290DFO"
ECLIP_FILES["EWSR1_K562"]="ENCFF607ZRF"

# HepG2 files (19 RBPs - IGF2BP2, NONO, EWSR1 only have K562)
ECLIP_FILES["IGF2BP1_HepG2"]="ENCFF442USD"
# IGF2BP2_HepG2 not available (K562 only)
ECLIP_FILES["HNRNPC_HepG2"]="ENCFF440ROZ"
ECLIP_FILES["TIA1_HepG2"]="ENCFF759KCD"
ECLIP_FILES["HNRNPK_HepG2"]="ENCFF855CPQ"
ECLIP_FILES["PCBP2_HepG2"]="ENCFF642GNE"
ECLIP_FILES["RBFOX2_HepG2"]="ENCFF871NYM"
ECLIP_FILES["PTBP1_HepG2"]="ENCFF726SQU"
ECLIP_FILES["TARDBP_HepG2"]="ENCFF673QBV"
ECLIP_FILES["QKI_HepG2"]="ENCFF704OCI"
ECLIP_FILES["SRSF1_HepG2"]="ENCFF934ANS"
ECLIP_FILES["SRSF9_HepG2"]="ENCFF765PIF"
ECLIP_FILES["RBM22_HepG2"]="ENCFF293IZG"
ECLIP_FILES["TRA2A_HepG2"]="ENCFF766OCH"
ECLIP_FILES["HNRNPL_HepG2"]="ENCFF266TKW"
ECLIP_FILES["LIN28B_HepG2"]="ENCFF341XMP"
ECLIP_FILES["FUS_HepG2"]="ENCFF972DFZ"
ECLIP_FILES["MATR3_HepG2"]="ENCFF587KKM"
ECLIP_FILES["HNRNPA1_HepG2"]="ENCFF797GSK"
ECLIP_FILES["HNRNPM_HepG2"]="ENCFF752JNY"
# NONO_HepG2 not available (K562 only)
ECLIP_FILES["U2AF2_HepG2"]="ENCFF721PWF"
# EWSR1_HepG2 not available (K562 only)

# ==========================
# Download All Files
# ==========================

success_count=0
fail_count=0
skip_count=0

log "Starting downloads..."
log ""

for key in "${!ECLIP_FILES[@]}"; do
    rbp="${key%_*}"
    cell_line="${key#*_}"
    accession="${ECLIP_FILES[$key]}"

    if download_eclip "${rbp}" "${cell_line}" "${accession}"; then
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
ls -la "${ECLIP_DIR}"/*.bed 2>/dev/null | while read line; do
    log "  ${line}"
done

# Note cell-line specific availability
log ""
log "Note: Cell-line availability varies by RBP:"
log ""
log "  K562 only (no HepG2):"
log "    - IGF2BP2"
log "    - NONO"
log "    - EWSR1"
log ""
log "  HepG2 only (no K562):"
log "    - PCBP2"
log ""
log "Total: 41 files (22 K562 + 19 HepG2)"
