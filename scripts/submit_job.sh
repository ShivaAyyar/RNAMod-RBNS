#!/bin/bash
#SBATCH --job-name=epitrans_rbp
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/mnt/vstor/SOM_CCCC_JGS25/sayyar/RNAMod-RBNS/logs/%x_%A_%a.out
#SBATCH --error=/mnt/vstor/SOM_CCCC_JGS25/sayyar/RNAMod-RBNS/logs/%x_%A_%a.err
#SBATCH --array=0-14

# ============================================================================
# Epitranscriptomic RBP Analysis Pipeline - SLURM Batch Submission Script
#
# Processes 15 RBPs that have BOTH RBNS and ENCODE eCLIP data available.
#
# Note: 8 RBPs from the original 23 lack RBNS data and cannot be analyzed:
#       PTBP1, QKI, SRSF1, MATR3, HNRNPA1, HNRNPM, NONO, U2AF2
#
# Directory Structure:
#   /mnt/vstor/SOM_CCCC_JGS25/sayyar/RNAMod-RBNS/
#   ├── data/
#   │   ├── genome/
#   │   │   ├── hg38.fa
#   │   │   ├── hg38.fa.fai
#   │   │   └── hg38.chrom.sizes
#   │   ├── rbns/
#   │   │   └── {RBP}_zscores.csv
#   │   ├── eclip/
#   │   │   └── {RBP}_{cell_line}.bed
#   │   └── mods/
#   │       ├── K562/
#   │       │   ├── m6A.bed
#   │       │   ├── pseudoU.bed
#   │       │   ├── m5C.bed
#   │       │   └── ac4C.bed
#   │       └── HepG2/
#   │           └── ...
#   ├── results/
#   │   └── {RBP}/
#   ├── logs/
#   └── src/
#
# Usage:
#   sbatch scripts/submit_job.sh                    # All 15 RBPs
#   sbatch --array=0 scripts/submit_job.sh          # IGF2BP1 only
#   sbatch --array=0-1 scripts/submit_job.sh        # Positive controls only
#   CELL_LINE=HepG2 sbatch scripts/submit_job.sh    # HepG2 cell line
# ============================================================================

set -euo pipefail

# ==========================
# Configuration
# ==========================

# Base directory for the project
BASE_DIR="/mnt/vstor/SOM_CCCC_JGS25/sayyar/RNAMod-RBNS"

# Data directory
DATA_DIR="${BASE_DIR}/data"

# Output directory (results will be in ${OUTPUT_BASE}/${RBP}/)
OUTPUT_BASE="${BASE_DIR}/results"

# Cell line to analyze (K562 or HepG2) - can be overridden via environment variable
CELL_LINE="${CELL_LINE:-K562}"

# Conda environment name
CONDA_ENV="rbp_mod_env"

# Path to source code
SRC_DIR="${BASE_DIR}/src"

# ==========================
# 15 RBPs with BOTH eCLIP AND RBNS Data
# ==========================
# Based on cross-referencing ENCODE eCLIP and RBNS experiments
# IGF2BP1 and IGF2BP2 are known m6A readers (positive controls)
#
# RBPs WITHOUT RBNS data (excluded from analysis):
#   PTBP1, QKI, SRSF1, MATR3, HNRNPA1, HNRNPM, NONO, U2AF2

RBPS=(
    "IGF2BP1"   # Index 0  - m6A reader (POSITIVE CONTROL)
    "IGF2BP2"   # Index 1  - m6A reader (POSITIVE CONTROL)
    "HNRNPC"    # Index 2
    "TIA1"      # Index 3
    "HNRNPK"    # Index 4
    "PCBP2"     # Index 5
    "RBFOX2"    # Index 6
    "TARDBP"    # Index 7
    "SRSF9"     # Index 8
    "RBM22"     # Index 9
    "TRA2A"     # Index 10
    "HNRNPL"    # Index 11
    "LIN28B"    # Index 12
    "FUS"       # Index 13
    "EWSR1"     # Index 14
)

# Get current RBP from array index
RBP="${RBPS[$SLURM_ARRAY_TASK_ID]}"

# ==========================
# Setup Environment
# ==========================

echo "=============================================="
echo "Epitranscriptomic RBP Analysis Pipeline"
echo "=============================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "RBP: ${RBP}"
echo "Cell Line: ${CELL_LINE}"
echo "Start Time: $(date)"
echo "Hostname: $(hostname)"
echo "Base Directory: ${BASE_DIR}"
echo "=============================================="

# Create log directory if needed
mkdir -p "${BASE_DIR}/logs"

# Load required modules (adjust for your HPC system)
# Uncomment and modify as needed for your cluster
# module purge
# module load anaconda3
# module load bedtools/2.30.0

# Activate conda environment
if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate ${CONDA_ENV}
elif [[ -f "${HOME}/miniconda3/bin/conda" ]]; then
    eval "$(${HOME}/miniconda3/bin/conda shell.bash hook)"
    conda activate ${CONDA_ENV}
elif [[ -f "${HOME}/anaconda3/bin/conda" ]]; then
    eval "$(${HOME}/anaconda3/bin/conda shell.bash hook)"
    conda activate ${CONDA_ENV}
else
    echo "ERROR: conda not found. Please ensure conda is in PATH."
    echo "Try: module load anaconda3"
    exit 1
fi

# Verify Python environment
echo "Python: $(which python)"
echo "Python version: $(python --version)"

# ==========================
# Define Input/Output Paths
# ==========================

RBNS_FILE="${DATA_DIR}/rbns/${RBP}_zscores.csv"
ECLIP_FILE="${DATA_DIR}/eclip/${RBP}_${CELL_LINE}.bed"
GENOME_FILE="${DATA_DIR}/genome/hg38.fa"
CHROM_SIZES="${DATA_DIR}/genome/hg38.chrom.sizes"
OUTPUT_DIR="${OUTPUT_BASE}/${RBP}"

# Modification files (cell-line specific)
MOD_DIR="${DATA_DIR}/mods/${CELL_LINE}"
M6A_FILE="${MOD_DIR}/m6A.bed"
PSEUDOU_FILE="${MOD_DIR}/pseudoU.bed"
M5C_FILE="${MOD_DIR}/m5C.bed"
AC4C_FILE="${MOD_DIR}/ac4C.bed"

# ==========================
# Validate Inputs
# ==========================

echo ""
echo "Validating input files..."

missing_files=0

check_file() {
    if [[ ! -f "$1" ]]; then
        echo "  ERROR: Missing file: $1"
        missing_files=$((missing_files + 1))
    else
        echo "  OK: $1"
    fi
}

check_file "${RBNS_FILE}"
check_file "${ECLIP_FILE}"
check_file "${GENOME_FILE}"
check_file "${CHROM_SIZES}"
check_file "${M6A_FILE}"
check_file "${PSEUDOU_FILE}"
check_file "${M5C_FILE}"
check_file "${AC4C_FILE}"

if [[ ${missing_files} -gt 0 ]]; then
    echo ""
    echo "ERROR: ${missing_files} required file(s) missing. Exiting."
    echo "Run 'bash scripts/download_data.sh' to download required data."
    exit 1
fi

echo ""
echo "All input files validated."

# ==========================
# Create Output Directory
# ==========================

mkdir -p "${OUTPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"

# ==========================
# Run Analysis
# ==========================

echo ""
echo "Starting analysis..."
echo ""

cd "${SRC_DIR}"

python main.py \
    --rbp "${RBP}" \
    --rbns "${RBNS_FILE}" \
    --eclip "${ECLIP_FILE}" \
    --genome "${GENOME_FILE}" \
    --chrom-sizes "${CHROM_SIZES}" \
    --mods "${M6A_FILE}" "${PSEUDOU_FILE}" "${M5C_FILE}" "${AC4C_FILE}" \
    --mod-names m6A pseudoU m5C ac4C \
    --output "${OUTPUT_DIR}" \
    --cell-line "${CELL_LINE}"

# ==========================
# Summary
# ==========================

echo ""
echo "=============================================="
echo "Analysis Complete"
echo "=============================================="
echo "RBP: ${RBP}"
echo "Output: ${OUTPUT_DIR}"
echo "End Time: $(date)"
echo "=============================================="

# Check for significant results
if [[ -f "${OUTPUT_DIR}/enrichment_results.csv" ]]; then
    echo ""
    echo "Enrichment Results Summary:"
    head -20 "${OUTPUT_DIR}/enrichment_results.csv"
fi

# Special message for positive controls
if [[ "${RBP}" == "IGF2BP1" ]] || [[ "${RBP}" == "IGF2BP2" ]]; then
    echo ""
    echo "*** NOTE: ${RBP} is a known m6A reader (positive control) ***"
    echo "*** Check that m6A shows significant enrichment in discrepant peaks ***"
fi

exit 0
