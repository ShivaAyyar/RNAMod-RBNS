#!/bin/bash
#SBATCH --job-name=download_data
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/mnt/vstor/SOM_CCCC_JGS25/sayyar/RNAMod-RBNS/logs/download_data_%j.out
#SBATCH --error=/mnt/vstor/SOM_CCCC_JGS25/sayyar/RNAMod-RBNS/logs/download_data_%j.err

# ============================================================================
# Data Download Script for RNAMod-RBNS Pipeline
#
# Downloads all prerequisite data:
# 1. Reference genome (hg38)
# 2. eCLIP peaks from ENCODE (24 Gold Standard RBPs)
# 3. RBNS Z-scores from ENCODE
# 4. RNA modification sites from REPIC and RMBase
#
# Usage:
#   bash scripts/download_data.sh      # Run interactively
#   sbatch scripts/download_data.sh    # Submit as SLURM job
#
# Note: This script requires wget, gunzip, and ~50GB disk space
# ============================================================================

set -euo pipefail

# ==========================
# Configuration
# ==========================

BASE_DIR="/mnt/vstor/SOM_CCCC_JGS25/sayyar/RNAMod-RBNS"
DATA_DIR="${BASE_DIR}/data"
CONDA_ENV="rbp_mod_env"

# Cell lines to download
CELL_LINES=("K562" "HepG2")

# 24 Gold Standard RBPs
RBPS=(
    "IGF2BP1" "IGF2BP2" "HNRNPC" "TIA1" "HNRNPK" "PCBP2"
    "RBFOX2" "PTBP3" "TARDBP" "QKI" "SRSF1" "SRSF9"
    "RBM22" "TRA2A" "HNRNPL" "LIN28B" "ZNF326" "FUS"
    "MATR3" "HNRNPA1" "HNRNPM" "NONO" "U2AF2" "EWSR1"
)

# ==========================
# Helper Functions
# ==========================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
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
# Create Directory Structure
# ==========================

log "Creating directory structure..."

mkdir -p "${DATA_DIR}/genome"
mkdir -p "${DATA_DIR}/rbns"
mkdir -p "${DATA_DIR}/eclip"
for cell in "${CELL_LINES[@]}"; do
    mkdir -p "${DATA_DIR}/mods/${cell}"
done
mkdir -p "${BASE_DIR}/results"
mkdir -p "${BASE_DIR}/logs"

log "Directory structure created at ${DATA_DIR}"

# ==========================
# Download Reference Genome
# ==========================

log "=============================================="
log "Downloading Reference Genome (hg38)"
log "=============================================="

GENOME_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
CHROM_SIZES_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"

if [[ ! -f "${DATA_DIR}/genome/hg38.fa" ]]; then
    download_file "${GENOME_URL}" "${DATA_DIR}/genome/hg38.fa.gz"
    log "Decompressing genome..."
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
    fi
else
    log "WARNING: samtools not found. FASTA index not created."
    log "Run 'samtools faidx ${DATA_DIR}/genome/hg38.fa' manually."
fi

# ==========================
# Download eCLIP Data from ENCODE
# ==========================

log "=============================================="
log "Downloading eCLIP Data from ENCODE"
log "=============================================="

# eCLIP data download using ENCODE API
# The script will search for IDR-filtered narrowPeak files for each RBP

download_eclip_via_api() {
    local rbp=$1
    local cell_line=$2
    local output="${DATA_DIR}/eclip/${rbp}_${cell_line}.bed"

    if [[ -f "${output}" && -s "${output}" ]]; then
        log "eCLIP exists: ${rbp} ${cell_line}, skipping"
        return 0
    fi

    log "Searching ENCODE for eCLIP: ${rbp} ${cell_line}..."

    # Use ENCODE API to find the correct file
    # Looking for: eCLIP, specific RBP target, specific cell line, IDR peaks, bed narrowPeak
    local api_url="https://www.encodeproject.org/search/?type=File&assay_title=eCLIP&target.label=${rbp}&biosample_ontology.term_name=${cell_line}&file_type=bed+narrowPeak&output_type=IDR+thresholded+peaks&status=released&format=json&limit=1"

    local json_response
    json_response=$(wget -q -O - "${api_url}" 2>/dev/null) || {
        log "WARNING: Could not query ENCODE API for ${rbp} ${cell_line}"
        create_eclip_readme "${rbp}" "${cell_line}" "${output}"
        return 1
    }

    # Extract the file accession from JSON response
    # Using grep/sed since jq may not be available
    local file_path
    file_path=$(echo "${json_response}" | grep -o '"href": "/files/[^"]*"' | head -1 | sed 's/"href": "\/files\/\([^/]*\)\/.*"/\1/')

    if [[ -z "${file_path}" || "${file_path}" == "null" ]]; then
        log "WARNING: No eCLIP data found for ${rbp} ${cell_line} - may need manual search"
        create_eclip_readme "${rbp}" "${cell_line}" "${output}"
        return 1
    fi

    local download_url="https://www.encodeproject.org/files/${file_path}/@@download/${file_path}.bed.gz"
    local temp_file="${output}.gz"

    log "Downloading eCLIP: ${rbp} ${cell_line} (${file_path})"

    if wget -q --show-progress -O "${temp_file}" "${download_url}" 2>/dev/null; then
        gunzip -f "${temp_file}"
        if [[ -s "${output}" ]]; then
            log "Downloaded: ${rbp} ${cell_line} ($(wc -l < "${output}") peaks)"
        else
            log "WARNING: Downloaded file is empty for ${rbp} ${cell_line}"
            rm -f "${output}"
            create_eclip_readme "${rbp}" "${cell_line}" "${output}"
        fi
    else
        log "WARNING: Could not download ${rbp} ${cell_line}"
        rm -f "${temp_file}"
        create_eclip_readme "${rbp}" "${cell_line}" "${output}"
    fi
}

create_eclip_readme() {
    local rbp=$1
    local cell_line=$2
    local output=$3

    cat > "${output}.README" << EOF
# Manual download required for ${rbp} eCLIP data (${cell_line})
#
# Steps:
# 1. Go to: https://www.encodeproject.org/search/?type=Experiment&assay_title=eCLIP&target.label=${rbp}&biosample_ontology.term_name=${cell_line}
# 2. Click on the experiment
# 3. Under "Files", find "IDR thresholded peaks" with file type "bed narrowPeak"
# 4. Download the .bed.gz file
# 5. Extract and save as: ${output}
#
# Alternative: Use ENCODE batch download
# https://www.encodeproject.org/batch_download/?type=File&assay_title=eCLIP&target.label=${rbp}&biosample_ontology.term_name=${cell_line}&file_type=bed+narrowPeak&output_type=IDR+thresholded+peaks
EOF
}

for rbp in "${RBPS[@]}"; do
    for cell_line in "${CELL_LINES[@]}"; do
        download_eclip_via_api "${rbp}" "${cell_line}"
    done
done

# ==========================
# Download RBNS Z-scores
# ==========================

log "=============================================="
log "Downloading RBNS Z-scores"
log "=============================================="

# RBNS Z-scores must be downloaded from Dominguez et al. (2018) supplementary data
# The ENCODE RBNS portal has raw data but not pre-computed Z-scores
#
# Primary source: Dominguez et al. (2018) Molecular Cell
# Supplementary Table S2 contains 5-mer Z-scores for all RBPs
# DOI: 10.1016/j.molcel.2018.06.012
# Direct link: https://www.cell.com/molecular-cell/supplemental/S1097-2765(18)30437-3

log "Creating template RBNS files..."
log "NOTE: RBNS Z-scores require manual download from:"
log "  1. ENCODE RBNS experiments: https://www.encodeproject.org/rbns/"
log "  2. Dominguez et al. (2018) Supplementary Table S2"
log "     DOI: 10.1016/j.molcel.2018.06.012"

# Create template files with instructions
for rbp in "${RBPS[@]}"; do
    output="${DATA_DIR}/rbns/${rbp}_zscores.csv"
    if [[ ! -f "${output}" ]]; then
        cat > "${output}" << EOF
kmer,z_score,r_value
# Template file for ${rbp} RBNS Z-scores
# Download from ENCODE or Dominguez et al. (2018) supplementary data
# Required columns: kmer, z_score (and optionally r_value)
# Example:
# AAAAA,0.23,1.05
# AAAAC,0.45,1.10
EOF
        log "Created template: ${output}"
    fi
done

# ==========================
# Download RNA Modification Data
# ==========================

log "=============================================="
log "Downloading RNA Modification Data"
log "=============================================="

# REPIC m6A data URLs (cell-line specific)
# Note: REPIC requires accepting terms - may need manual download
REPIC_BASE="https://repicmod.uchicago.edu/download"

log "NOTE: RNA modification data requires manual download due to database terms."
log ""
log "m6A data from REPIC (https://repicmod.uchicago.edu/):"
log "  1. Go to https://repicmod.uchicago.edu/download"
log "  2. Select 'Human' -> 'K562' or 'HepG2'"
log "  3. Download m6A consensus peaks"
log "  4. Save to: ${DATA_DIR}/mods/{cell_line}/m6A.bed"
log ""
log "Pseudouridine, m5C, ac4C from RMBase (http://rna.sysu.edu.cn/rmbase/):"
log "  1. Go to http://rna.sysu.edu.cn/rmbase/download.html"
log "  2. Select species: Human (hg38)"
log "  3. Select cell line: K562 or HepG2"
log "  4. Download each modification type"
log "  5. Save to: ${DATA_DIR}/mods/{cell_line}/{modification}.bed"

# Create placeholder files with download instructions
for cell in "${CELL_LINES[@]}"; do
    for mod in "m6A" "pseudoU" "m5C" "ac4C"; do
        output="${DATA_DIR}/mods/${cell}/${mod}.bed"
        if [[ ! -f "${output}" ]]; then
            cat > "${output}" << EOF
# Placeholder for ${mod} modification sites in ${cell}
#
# Download instructions:
#
# For m6A:
#   Source: REPIC (https://repicmod.uchicago.edu/)
#   1. Navigate to Download section
#   2. Select Human -> ${cell}
#   3. Download m6A consensus peaks (BED format)
#
# For pseudoU, m5C, ac4C:
#   Source: RMBase v2.0 (http://rna.sysu.edu.cn/rmbase/)
#   1. Go to download page
#   2. Select: Human, hg38, ${cell}
#   3. Select modification type: ${mod}
#   4. Download BED file
#
# Expected BED format (tab-separated):
# chr1	100	101	${mod}_site_1	.	+
# chr1	200	201	${mod}_site_2	.	-
EOF
            log "Created placeholder: ${output}"
        fi
    done
done

# ==========================
# Alternative: Download from Supplementary Data
# ==========================

log ""
log "=============================================="
log "Alternative Data Sources"
log "=============================================="
log ""
log "If ENCODE downloads fail, try these alternatives:"
log ""
log "1. Dominguez et al. (2018) Supplementary Data:"
log "   https://www.cell.com/molecular-cell/supplemental/S1097-2765(18)30437-3"
log "   - Table S2: RBNS Z-scores for all RBPs"
log "   - Table S1: RBP metadata"
log ""
log "2. ENCODE Bulk Download:"
log "   https://www.encodeproject.org/batch_download/"
log "   - Use filters: Assay=eCLIP, Biosample=K562/HepG2"
log ""
log "3. RMBase v3.0 (updated):"
log "   https://rmbase.org/"
log ""

# ==========================
# Create Download Summary
# ==========================

SUMMARY_FILE="${DATA_DIR}/download_summary.txt"

cat > "${SUMMARY_FILE}" << EOF
RNAMod-RBNS Data Download Summary
Generated: $(date)
Base Directory: ${DATA_DIR}

============================================
Directory Structure
============================================
${DATA_DIR}/
├── genome/
│   ├── hg38.fa              $([ -f "${DATA_DIR}/genome/hg38.fa" ] && echo "[DOWNLOADED]" || echo "[MISSING]")
│   ├── hg38.fa.fai          $([ -f "${DATA_DIR}/genome/hg38.fa.fai" ] && echo "[DOWNLOADED]" || echo "[MISSING]")
│   └── hg38.chrom.sizes     $([ -f "${DATA_DIR}/genome/hg38.chrom.sizes" ] && echo "[DOWNLOADED]" || echo "[MISSING]")
├── rbns/
│   └── {RBP}_zscores.csv    [TEMPLATES CREATED - MANUAL DOWNLOAD REQUIRED]
├── eclip/
│   └── {RBP}_{cell}.bed     [CHECK INDIVIDUAL FILES]
└── mods/
    ├── K562/
    │   ├── m6A.bed          [MANUAL DOWNLOAD REQUIRED]
    │   ├── pseudoU.bed      [MANUAL DOWNLOAD REQUIRED]
    │   ├── m5C.bed          [MANUAL DOWNLOAD REQUIRED]
    │   └── ac4C.bed         [MANUAL DOWNLOAD REQUIRED]
    └── HepG2/
        └── ...              [MANUAL DOWNLOAD REQUIRED]

============================================
eCLIP Download Status
============================================
EOF

for rbp in "${RBPS[@]}"; do
    for cell in "${CELL_LINES[@]}"; do
        file="${DATA_DIR}/eclip/${rbp}_${cell}.bed"
        if [[ -f "${file}" && ! -s "${file}" ]]; then
            status="[EMPTY]"
        elif [[ -f "${file}" ]]; then
            status="[OK]"
        else
            status="[MISSING]"
        fi
        echo "${rbp} ${cell}: ${status}" >> "${SUMMARY_FILE}"
    done
done

cat >> "${SUMMARY_FILE}" << EOF

============================================
Manual Download Required
============================================
1. RBNS Z-scores: See template files in ${DATA_DIR}/rbns/
2. RNA modifications: See placeholder files in ${DATA_DIR}/mods/

============================================
Next Steps
============================================
1. Check download_summary.txt for missing files
2. Download RBNS data from Dominguez et al. (2018) supplementary
3. Download modification data from REPIC and RMBase
4. Run: bash scripts/setup_conda.sh
5. Submit: sbatch scripts/submit_job.sh
EOF

log ""
log "=============================================="
log "Download Complete"
log "=============================================="
log "Summary written to: ${SUMMARY_FILE}"
log ""
log "IMPORTANT: Several data files require manual download."
log "Please review: ${SUMMARY_FILE}"
log ""
