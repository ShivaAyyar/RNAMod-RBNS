#!/bin/bash
# ============================================================================
# Conda Environment Setup Script for RNAMod-RBNS Pipeline
#
# Creates and configures the conda environment for the analysis pipeline.
#
# Usage:
#   bash scripts/setup_conda.sh
#
# For SLURM submission:
#   sbatch scripts/setup_conda.sh
#
# ============================================================================

#SBATCH --job-name=setup_conda
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/mnt/vstor/SOM_CCCC_JGS25/sayyar/RNAMod-RBNS/logs/setup_conda_%j.out
#SBATCH --error=/mnt/vstor/SOM_CCCC_JGS25/sayyar/RNAMod-RBNS/logs/setup_conda_%j.err

set -euo pipefail

# ==========================
# Configuration
# ==========================

BASE_DIR="/mnt/vstor/SOM_CCCC_JGS25/sayyar/RNAMod-RBNS"
ENV_NAME="rbp_mod_env"
ENV_FILE="${BASE_DIR}/environment.yml"

# ==========================
# Helper Functions
# ==========================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# ==========================
# Setup
# ==========================

log "=============================================="
log "RNAMod-RBNS Conda Environment Setup"
log "=============================================="
log "Base Directory: ${BASE_DIR}"
log "Environment Name: ${ENV_NAME}"
log "Environment File: ${ENV_FILE}"
log ""

# Create logs directory
mkdir -p "${BASE_DIR}/logs"

# Change to base directory
cd "${BASE_DIR}"

# ==========================
# Check for Conda
# ==========================

log "Checking for conda installation..."

# Try to find conda
if command -v conda &> /dev/null; then
    log "Found conda: $(which conda)"
elif [[ -f "${HOME}/miniconda3/bin/conda" ]]; then
    log "Found conda at ~/miniconda3"
    export PATH="${HOME}/miniconda3/bin:${PATH}"
elif [[ -f "${HOME}/anaconda3/bin/conda" ]]; then
    log "Found conda at ~/anaconda3"
    export PATH="${HOME}/anaconda3/bin:${PATH}"
elif [[ -f "/opt/conda/bin/conda" ]]; then
    log "Found conda at /opt/conda"
    export PATH="/opt/conda/bin:${PATH}"
else
    # Try loading module on HPC systems
    if command -v module &> /dev/null; then
        log "Attempting to load conda module..."
        module load anaconda3 2>/dev/null || \
        module load miniconda3 2>/dev/null || \
        module load conda 2>/dev/null || {
            log "ERROR: Could not find or load conda."
            log "Please install Miniconda or load the appropriate module."
            log "Install Miniconda: https://docs.conda.io/en/latest/miniconda.html"
            exit 1
        }
    else
        log "ERROR: conda not found in PATH"
        log "Please install Miniconda: https://docs.conda.io/en/latest/miniconda.html"
        exit 1
    fi
fi

# Initialize conda for this shell
log "Initializing conda..."
eval "$(conda shell.bash hook)"

log "Conda version: $(conda --version)"

# ==========================
# Check Environment File
# ==========================

if [[ ! -f "${ENV_FILE}" ]]; then
    log "ERROR: Environment file not found: ${ENV_FILE}"
    log "Please ensure environment.yml exists in the repository root."
    exit 1
fi

log "Environment file found: ${ENV_FILE}"

# ==========================
# Remove Existing Environment (if requested)
# ==========================

if [[ "${1:-}" == "--force" ]] || [[ "${1:-}" == "-f" ]]; then
    log "Force flag detected, removing existing environment..."
    conda env remove -n "${ENV_NAME}" -y 2>/dev/null || true
fi

# ==========================
# Create/Update Environment
# ==========================

if conda env list | grep -q "^${ENV_NAME} "; then
    log "Environment '${ENV_NAME}' exists. Updating..."
    conda env update -n "${ENV_NAME}" -f "${ENV_FILE}" --prune
else
    log "Creating new environment '${ENV_NAME}'..."
    conda env create -f "${ENV_FILE}"
fi

# ==========================
# Verify Installation
# ==========================

log ""
log "Verifying installation..."

# Activate environment
conda activate "${ENV_NAME}"

log "Python: $(which python)"
log "Python version: $(python --version)"

# Test imports
log ""
log "Testing package imports..."

python << 'EOF'
import sys

packages = [
    ('pandas', 'pd'),
    ('numpy', 'np'),
    ('pybedtools', 'pybedtools'),
    ('pyfaidx', 'pyfaidx'),
    ('scipy', 'scipy'),
    ('scipy.stats', 'stats'),
    ('statsmodels', 'statsmodels'),
    ('statsmodels.stats.multitest', 'multitest'),
    ('matplotlib', 'plt'),
    ('seaborn', 'sns'),
    ('tqdm', 'tqdm'),
]

failed = []
for pkg, alias in packages:
    try:
        exec(f"import {pkg}")
        print(f"  [OK] {pkg}")
    except ImportError as e:
        print(f"  [FAIL] {pkg}: {e}")
        failed.append(pkg)

if failed:
    print(f"\nERROR: Failed to import: {', '.join(failed)}")
    sys.exit(1)
else:
    print("\nAll packages imported successfully!")
EOF

if [[ $? -ne 0 ]]; then
    log "ERROR: Package verification failed!"
    exit 1
fi

# ==========================
# Check bedtools
# ==========================

log ""
log "Checking bedtools..."

if command -v bedtools &> /dev/null; then
    log "bedtools version: $(bedtools --version)"
else
    log "WARNING: bedtools not found in PATH"
    log "bedtools should be installed via conda, checking..."

    # Check if it's in the conda environment
    if [[ -f "${CONDA_PREFIX}/bin/bedtools" ]]; then
        log "bedtools found in conda environment: ${CONDA_PREFIX}/bin/bedtools"
        log "bedtools version: $(${CONDA_PREFIX}/bin/bedtools --version)"
    else
        log "ERROR: bedtools not found. Please install manually:"
        log "  conda install -c bioconda bedtools"
    fi
fi

# ==========================
# Create Activation Script
# ==========================

ACTIVATE_SCRIPT="${BASE_DIR}/activate_env.sh"

cat > "${ACTIVATE_SCRIPT}" << EOF
#!/bin/bash
# Source this file to activate the RNAMod-RBNS environment
# Usage: source activate_env.sh

# Initialize conda
eval "\$(conda shell.bash hook)"

# Activate environment
conda activate ${ENV_NAME}

echo "Activated environment: ${ENV_NAME}"
echo "Python: \$(which python)"
EOF

chmod +x "${ACTIVATE_SCRIPT}"
log "Created activation script: ${ACTIVATE_SCRIPT}"

# ==========================
# Summary
# ==========================

log ""
log "=============================================="
log "Setup Complete!"
log "=============================================="
log ""
log "Environment: ${ENV_NAME}"
log "Location: $(conda env list | grep "^${ENV_NAME} " | awk '{print $2}')"
log ""
log "To activate the environment:"
log "  source ${ACTIVATE_SCRIPT}"
log ""
log "Or manually:"
log "  conda activate ${ENV_NAME}"
log ""
log "Next steps:"
log "  1. Download data: bash scripts/download_data.sh"
log "  2. Submit analysis: sbatch scripts/submit_job.sh"
log ""
