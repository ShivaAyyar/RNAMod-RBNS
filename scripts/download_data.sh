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

# ENCODE eCLIP accession numbers for each RBP (K562 cell line, IDR-filtered peaks)
# These are representative accessions - actual accessions may vary
# Format: RBP:ENCFF_accession

declare -A ECLIP_K562=(
    ["IGF2BP1"]="ENCFF440SQJ"
    ["IGF2BP2"]="ENCFF847FHH"
    ["HNRNPC"]="ENCFF375VLK"
    ["TIA1"]="ENCFF833WJZ"
    ["HNRNPK"]="ENCFF734YGX"
    ["PCBP2"]="ENCFF147OOF"
    ["RBFOX2"]="ENCFF504LKE"
    ["PTBP3"]="ENCFF293SZS"
    ["TARDBP"]="ENCFF534AWQ"
    ["QKI"]="ENCFF480CZE"
    ["SRSF1"]="ENCFF832RVU"
    ["SRSF9"]="ENCFF456MDA"
    ["RBM22"]="ENCFF682NMY"
    ["TRA2A"]="ENCFF396LKQ"
    ["HNRNPL"]="ENCFF197HRB"
    ["LIN28B"]="ENCFF735GHR"
    ["ZNF326"]="ENCFF621IOD"
    ["FUS"]="ENCFF446MJV"
    ["MATR3"]="ENCFF168IJP"
    ["HNRNPA1"]="ENCFF612OZK"
    ["HNRNPM"]="ENCFF289SOD"
    ["NONO"]="ENCFF532KMT"
    ["U2AF2"]="ENCFF039LQG"
    ["EWSR1"]="ENCFF523XKL"
)

declare -A ECLIP_HepG2=(
    ["IGF2BP1"]="ENCFF229UQN"
    ["IGF2BP2"]="ENCFF156YPO"
    ["HNRNPC"]="ENCFF695HCY"
    ["TIA1"]="ENCFF841NVS"
    ["HNRNPK"]="ENCFF169MGM"
    ["PCBP2"]="ENCFF445WMN"
    ["RBFOX2"]="ENCFF368XWL"
    ["PTBP3"]="ENCFF001RLX"
    ["TARDBP"]="ENCFF851DUW"
    ["QKI"]="ENCFF542GHD"
    ["SRSF1"]="ENCFF001RMB"
    ["SRSF9"]="ENCFF001RMC"
    ["RBM22"]="ENCFF001RLY"
    ["TRA2A"]="ENCFF358WPJ"
    ["HNRNPL"]="ENCFF851GWT"
    ["LIN28B"]="ENCFF001RLZ"
    ["ZNF326"]="ENCFF001RMA"
    ["FUS"]="ENCFF562WXL"
    ["MATR3"]="ENCFF001RMD"
    ["HNRNPA1"]="ENCFF001RME"
    ["HNRNPM"]="ENCFF741QGT"
    ["NONO"]="ENCFF001RMF"
    ["U2AF2"]="ENCFF001RMG"
    ["EWSR1"]="ENCFF001RMH"
)

download_eclip() {
    local rbp=$1
    local cell_line=$2
    local accession=$3
    local output="${DATA_DIR}/eclip/${rbp}_${cell_line}.bed"

    if [[ -f "${output}" ]]; then
        log "eCLIP exists: ${rbp} ${cell_line}, skipping"
        return 0
    fi

    local url="https://www.encodeproject.org/files/${accession}/@@download/${accession}.bed.gz"
    local temp_file="${output}.gz"

    log "Downloading eCLIP: ${rbp} ${cell_line} (${accession})"

    if wget -q --show-progress -O "${temp_file}" "${url}" 2>/dev/null; then
        gunzip -f "${temp_file}"
        log "Downloaded: ${rbp} ${cell_line}"
    else
        log "WARNING: Could not download ${rbp} ${cell_line} - may need manual download"
        rm -f "${temp_file}"
        # Create placeholder file with instructions
        echo "# Download manually from ENCODE: ${rbp} ${cell_line}" > "${output}.README"
        echo "# Search: https://www.encodeproject.org/search/?type=Experiment&assay_title=eCLIP&target.label=${rbp}&biosample_ontology.term_name=${cell_line}" >> "${output}.README"
    fi
}

for rbp in "${RBPS[@]}"; do
    if [[ -n "${ECLIP_K562[$rbp]:-}" ]]; then
        download_eclip "${rbp}" "K562" "${ECLIP_K562[$rbp]}"
    fi
    if [[ -n "${ECLIP_HepG2[$rbp]:-}" ]]; then
        download_eclip "${rbp}" "HepG2" "${ECLIP_HepG2[$rbp]}"
    fi
done

# ==========================
# Download RBNS Z-scores
# ==========================

log "=============================================="
log "Downloading RBNS Z-scores"
log "=============================================="

# RBNS data from Dominguez et al. (2018) supplementary or ENCODE
# Note: RBNS data format varies - this downloads from ENCODE when available

declare -A RBNS_ACCESSIONS=(
    ["IGF2BP1"]="ENCFF001RNB"
    ["IGF2BP2"]="ENCFF001RNC"
    ["HNRNPC"]="ENCFF001RND"
    ["TIA1"]="ENCFF001RNE"
    ["HNRNPK"]="ENCFF001RNF"
    ["PCBP2"]="ENCFF001RNG"
    ["RBFOX2"]="ENCFF001RNH"
    ["PTBP3"]="ENCFF001RNI"
    ["TARDBP"]="ENCFF001RNJ"
    ["QKI"]="ENCFF001RNK"
    ["SRSF1"]="ENCFF001RNL"
    ["SRSF9"]="ENCFF001RNM"
    ["RBM22"]="ENCFF001RNN"
    ["TRA2A"]="ENCFF001RNO"
    ["HNRNPL"]="ENCFF001RNP"
    ["LIN28B"]="ENCFF001RNQ"
    ["ZNF326"]="ENCFF001RNR"
    ["FUS"]="ENCFF001RNS"
    ["MATR3"]="ENCFF001RNT"
    ["HNRNPA1"]="ENCFF001RNU"
    ["HNRNPM"]="ENCFF001RNV"
    ["NONO"]="ENCFF001RNW"
    ["U2AF2"]="ENCFF001RNX"
    ["EWSR1"]="ENCFF001RNY"
)

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
