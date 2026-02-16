#!/bin/bash
#SBATCH --job-name=download_data
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/mnt/vstor/SOM_CCCC_JGS25/sayyar/RNAMod-RBNS/logs/download_data_%j.out
#SBATCH --error=/mnt/vstor/SOM_CCCC_JGS25/sayyar/RNAMod-RBNS/logs/download_data_%j.err

# ============================================================================
# Data Download Script for RNAMod-RBNS Pipeline
#
# This script handles ONLY automated downloads:
# - Reference genome (hg38) from UCSC
# - Chromosome sizes
# - FASTA index creation
#
# For eCLIP, RBNS, and modification data, follow the manual download
# instructions in: docs/DATA_DOWNLOAD.md
#
# Usage:
#   bash scripts/download_data.sh      # Run interactively
#   sbatch scripts/download_data.sh    # Submit as SLURM job
#
# Note: This script requires wget, gunzip, and ~3GB disk space for genome
# ============================================================================

set -euo pipefail

# ==========================
# Configuration
# ==========================

BASE_DIR="/mnt/vstor/SOM_CCCC_JGS25/sayyar/RNAMod-RBNS"
DATA_DIR="${BASE_DIR}/data"
CONDA_ENV="rbp_mod_env"

# Cell lines for directory structure
CELL_LINES=("K562" "HepG2")

# ==========================
# Helper Functions
# ==========================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

download_file() {
    local url=$1
    local output=$2

    if [[ -f "${output}" ]]; then
        log "File exists, skipping: ${output}"
        return 0
    fi

    log "Downloading: ${url}"
    wget -q --show-progress -O "${output}" "${url}" || {
        log "ERROR: Failed to download ${url}"
        rm -f "${output}"
        return 1
    }
}

# ==========================
# Load Environment for Tools
# ==========================

log "Setting up environment..."

# Try to load samtools module (HPC systems)
if command -v module &> /dev/null; then
    module load samtools 2>/dev/null && log "Loaded samtools module" || true
fi

# Try to activate conda environment if it exists (for samtools from bioconda)
if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)" 2>/dev/null || true
    if conda env list 2>/dev/null | grep -q "^${CONDA_ENV} "; then
        conda activate "${CONDA_ENV}" 2>/dev/null && log "Activated conda environment: ${CONDA_ENV}" || true
    fi
elif [[ -f "${HOME}/miniconda3/bin/conda" ]]; then
    eval "$(${HOME}/miniconda3/bin/conda shell.bash hook)" 2>/dev/null || true
    if ${HOME}/miniconda3/bin/conda env list 2>/dev/null | grep -q "^${CONDA_ENV} "; then
        conda activate "${CONDA_ENV}" 2>/dev/null && log "Activated conda environment: ${CONDA_ENV}" || true
    fi
elif [[ -f "${HOME}/anaconda3/bin/conda" ]]; then
    eval "$(${HOME}/anaconda3/bin/conda shell.bash hook)" 2>/dev/null || true
    if ${HOME}/anaconda3/bin/conda env list 2>/dev/null | grep -q "^${CONDA_ENV} "; then
        conda activate "${CONDA_ENV}" 2>/dev/null && log "Activated conda environment: ${CONDA_ENV}" || true
    fi
fi

# Check what tools are available
log "Checking available tools..."
command -v wget &> /dev/null && log "  wget: available" || log "  wget: NOT FOUND (required)"
command -v gunzip &> /dev/null && log "  gunzip: available" || log "  gunzip: NOT FOUND (required)"
command -v samtools &> /dev/null && log "  samtools: available ($(samtools --version | head -1))" || log "  samtools: not found (optional, for FASTA indexing)"

# ==========================
# Create Directory Structure
# ==========================

log "=============================================="
log "Creating Directory Structure"
log "=============================================="

mkdir -p "${DATA_DIR}/genome"
mkdir -p "${DATA_DIR}/rbns"
mkdir -p "${DATA_DIR}/eclip"
for cell in "${CELL_LINES[@]}"; do
    mkdir -p "${DATA_DIR}/mods/${cell}"
done
mkdir -p "${BASE_DIR}/results"
mkdir -p "${BASE_DIR}/logs"

log "Directory structure created at ${DATA_DIR}"
log ""
log "  ${DATA_DIR}/"
log "  ├── genome/     (auto-downloaded)"
log "  ├── rbns/       (manual download required)"
log "  ├── eclip/      (manual download required)"
log "  └── mods/"
log "      ├── K562/   (manual download required)"
log "      └── HepG2/  (manual download required)"

# ==========================
# Download Reference Genome
# ==========================

log ""
log "=============================================="
log "Downloading Reference Genome (hg38)"
log "=============================================="

GENOME_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
CHROM_SIZES_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"

if [[ ! -f "${DATA_DIR}/genome/hg38.fa" ]]; then
    download_file "${GENOME_URL}" "${DATA_DIR}/genome/hg38.fa.gz"
    log "Decompressing genome (this may take a few minutes)..."
    gunzip "${DATA_DIR}/genome/hg38.fa.gz"
    log "Genome decompressed"
else
    log "Genome already exists, skipping"
fi

download_file "${CHROM_SIZES_URL}" "${DATA_DIR}/genome/hg38.chrom.sizes"

# Create FASTA index if samtools is available
if command -v samtools &> /dev/null; then
    if [[ ! -f "${DATA_DIR}/genome/hg38.fa.fai" ]]; then
        log "Creating FASTA index..."
        samtools faidx "${DATA_DIR}/genome/hg38.fa"
        log "FASTA index created"
    else
        log "FASTA index already exists"
    fi
else
    log "WARNING: samtools not found. FASTA index not created."
    log "Run manually after installing samtools:"
    log "  samtools faidx ${DATA_DIR}/genome/hg38.fa"
fi

# ==========================
# Create Download Summary
# ==========================

log ""
log "=============================================="
log "Creating Download Summary"
log "=============================================="

SUMMARY_FILE="${DATA_DIR}/download_summary.txt"

cat > "${SUMMARY_FILE}" << EOF
RNAMod-RBNS Data Download Summary
Generated: $(date)
Base Directory: ${DATA_DIR}

============================================
Automated Downloads (Complete)
============================================

Genome files:
  hg38.fa            $([ -f "${DATA_DIR}/genome/hg38.fa" ] && echo "[DOWNLOADED - $(du -h "${DATA_DIR}/genome/hg38.fa" | cut -f1)]" || echo "[MISSING]")
  hg38.fa.fai        $([ -f "${DATA_DIR}/genome/hg38.fa.fai" ] && echo "[CREATED]" || echo "[MISSING - run samtools faidx]")
  hg38.chrom.sizes   $([ -f "${DATA_DIR}/genome/hg38.chrom.sizes" ] && echo "[DOWNLOADED]" || echo "[MISSING]")

============================================
Manual Downloads Required
============================================

The following data must be downloaded manually. See docs/DATA_DOWNLOAD.md
for detailed instructions.

1. eCLIP Peak Data (ENCODE)
   Location: ${DATA_DIR}/eclip/
   Format: {RBP}_{CellLine}.bed (e.g., IGF2BP1_K562.bed)
   Source: https://www.encodeproject.org/

2. RBNS Z-scores (Dominguez et al. 2018)
   Location: ${DATA_DIR}/rbns/
   Format: {RBP}_zscores.csv
   Source: Molecular Cell Supplementary Table S2

3. RNA Modification Sites
   Location: ${DATA_DIR}/mods/{CellLine}/
   Files needed:
   - m6A.bed      (Source: REPIC)
   - pseudoU.bed  (Source: RMBase v2.0)
   - m5C.bed      (Source: RMBase v2.0)
   - ac4C.bed     (Source: RMBase v2.0)

============================================
Next Steps
============================================

1. Read the detailed download guide:
   cat docs/DATA_DOWNLOAD.md

2. Download required data files manually

3. Verify all files are in place:
   ls -la ${DATA_DIR}/eclip/
   ls -la ${DATA_DIR}/rbns/
   ls -la ${DATA_DIR}/mods/K562/

4. Create/update conda environment:
   conda env create -f environment.yml
   # OR
   conda env update -f environment.yml

5. Run the analysis:
   sbatch scripts/submit_job.sh
EOF

log "Summary written to: ${SUMMARY_FILE}"

# ==========================
# Final Status Report
# ==========================

log ""
log "=============================================="
log "Download Script Complete"
log "=============================================="
log ""
log "Automated downloads completed:"
log "  [$([ -f "${DATA_DIR}/genome/hg38.fa" ] && echo "OK" || echo "FAIL")] Reference genome (hg38.fa)"
log "  [$([ -f "${DATA_DIR}/genome/hg38.chrom.sizes" ] && echo "OK" || echo "FAIL")] Chromosome sizes"
log "  [$([ -f "${DATA_DIR}/genome/hg38.fa.fai" ] && echo "OK" || echo "SKIP")] FASTA index"
log ""
log "Manual downloads required:"
log "  [ ] eCLIP data (24 RBPs × 2 cell lines)"
log "  [ ] RBNS Z-scores (24 RBPs)"
log "  [ ] Modification data (4 types × 2 cell lines)"
log ""
log "For detailed manual download instructions, see:"
log "  ${BASE_DIR}/docs/DATA_DOWNLOAD.md"
log ""
log "Summary file: ${SUMMARY_FILE}"
log ""
