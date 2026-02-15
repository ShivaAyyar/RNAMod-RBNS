# RNAMod-RBNS: Epitranscriptomic Determinants of RBP Specificity

A computational pipeline for investigating the "Specificity Paradox" - where RNA Binding Proteins (RBPs) bind sequences in vivo (eCLIP) that they show low affinity for in vitro (RBNS).

**Hypothesis:** RNA modifications (m6A, Ψ, m5C, ac4C) present in cellular RNA but absent in synthetic RBNS libraries explain binding discrepancies.

## Overview

This pipeline implements the analysis framework from [Dominguez et al. (2018)](https://doi.org/10.1016/j.molcel.2018.06.012) to:

1. Score eCLIP peaks using pre-computed RBNS k-mer Z-scores
2. Classify peaks as **Canonical** (high RBNS affinity) or **Discrepant** (low RBNS affinity)
3. Test for enrichment of RNA modifications in discrepant peaks
4. Generate publication-quality visualizations

## Quick Start (HPC)

For SLURM-based HPC environments (e.g., CWRU HPC):

```bash
# 1. Clone the repository to your HPC storage
cd /mnt/vstor/SOM_CCCC_JGS25/sayyar
git clone https://github.com/ShivaAyyar/RNAMod-RBNS.git
cd RNAMod-RBNS

# 2. Setup conda environment (submit as SLURM job)
sbatch scripts/setup_conda.sh

# 3. Download all prerequisite data (submit as SLURM job)
sbatch scripts/download_data.sh

# 4. After data download, manually download RBNS Z-scores and modification data
#    (see download_summary.txt for instructions)

# 5. Run analysis for all 24 RBPs
sbatch scripts/submit_job.sh

# Or run just the positive controls (IGF2BP1, IGF2BP2)
sbatch --array=0-1 scripts/submit_job.sh
```

## Quick Start (Local)

For local development or testing:

```bash
# 1. Clone the repository
git clone https://github.com/ShivaAyyar/RNAMod-RBNS.git
cd RNAMod-RBNS

# 2. Create conda environment
conda env create -f environment.yml
conda activate rbp_mod_env

# 3. Download data (run interactively)
bash scripts/download_data.sh

# 4. Run analysis for a single RBP
python src/main.py \
    --rbp IGF2BP1 \
    --rbns data/rbns/IGF2BP1_zscores.csv \
    --eclip data/eclip/IGF2BP1_K562.bed \
    --genome data/genome/hg38.fa \
    --chrom-sizes data/genome/hg38.chrom.sizes \
    --mods data/mods/K562/m6A.bed data/mods/K562/pseudoU.bed \
    --mod-names m6A pseudoU \
    --output results/IGF2BP1
```

## Installation

### Prerequisites

- Python 3.8+
- Conda (recommended) or pip
- bedtools 2.30.0+
- wget (for data download)
- ~50GB disk space for all data

### Setup via SLURM (Recommended for HPC)

```bash
# Submit conda setup as SLURM job
sbatch scripts/setup_conda.sh

# Or run interactively
bash scripts/setup_conda.sh
```

### Manual Setup

```bash
# Create and activate conda environment
conda env create -f environment.yml
conda activate rbp_mod_env

# Verify installation
python -c "import pybedtools; import scipy; import statsmodels; print('OK')"
```

## Data Download

The pipeline requires several data files. Use the automated download script:

```bash
# Submit as SLURM job (recommended for large downloads)
sbatch scripts/download_data.sh

# Or run interactively
bash scripts/download_data.sh
```

The script will:
1. Download the hg38 reference genome (~3GB)
2. Attempt to download eCLIP peaks from ENCODE
3. Create template files for RBNS Z-scores (manual download required)
4. Create placeholder files for modification data (manual download required)

After running, check `data/download_summary.txt` for status and manual download instructions.

### Manual Downloads Required

Some data requires manual download due to database terms:

1. **RBNS Z-scores**: Download from [Dominguez et al. (2018) Supplementary Table S2](https://www.cell.com/molecular-cell/supplemental/S1097-2765(18)30437-3)

2. **m6A modification sites**: Download from [REPIC](https://repicmod.uchicago.edu/)
   - Select Human → K562 or HepG2
   - Download consensus m6A peaks

3. **Ψ, m5C, ac4C sites**: Download from [RMBase v2.0](http://rna.sysu.edu.cn/rmbase/)
   - Select species: Human (hg38)
   - Select cell line and modification type

## Data Requirements

### Directory Structure

```
data/
├── rbns/
│   └── {RBP}_zscores.csv           # Pre-computed from RBNS pipeline
├── eclip/
│   ├── {RBP}_K562.bed              # ENCODE narrowPeak (IDR-filtered)
│   └── {RBP}_HepG2.bed
├── genome/
│   ├── hg38.fa                     # Reference genome
│   ├── hg38.fa.fai                 # FASTA index
│   └── hg38.chrom.sizes            # Chromosome sizes
└── mods/
    ├── K562/
    │   ├── m6A.bed                 # From REPIC
    │   ├── pseudoU.bed             # From RMBase
    │   ├── m5C.bed
    │   └── ac4C.bed
    └── HepG2/
        └── ...
```

### Data Sources

| Data Type | Source | URL |
|-----------|--------|-----|
| RBNS Z-scores | ENCODE / Dominguez et al. | [ENCODE RBNS](https://www.encodeproject.org/rbns/) |
| eCLIP peaks | ENCODE | [ENCODE eCLIP](https://www.encodeproject.org/eclip/) |
| m6A sites | REPIC | [repicmod.uchicago.edu](https://repicmod.uchicago.edu/) |
| Ψ, m5C, ac4C | RMBase v2.0 | [rna.sysu.edu.cn/rmbase](http://rna.sysu.edu.cn/rmbase/) |
| Reference genome | UCSC | [hg38.fa](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/) |

### RBNS Z-scores Format

Pre-computed k-mer enrichment Z-scores from the [RBNS pipeline](https://github.com/cburgelab/RBNS_pipeline):

```csv
kmer,z_score,r_value
AAAAA,0.23,1.05
AAAAC,0.45,1.10
...
GCAUG,8.45,15.2
```

## Usage

### Single RBP Analysis

```bash
python src/main.py \
    --rbp IGF2BP1 \
    --rbns data/rbns/IGF2BP1_zscores.csv \
    --eclip data/eclip/IGF2BP1_K562.bed \
    --genome data/genome/hg38.fa \
    --chrom-sizes data/genome/hg38.chrom.sizes \
    --mods data/mods/K562/m6A.bed data/mods/K562/pseudoU.bed data/mods/K562/m5C.bed data/mods/K562/ac4C.bed \
    --mod-names m6A pseudoU m5C ac4C \
    --output results/IGF2BP1 \
    --cell-line K562
```

### Batch Analysis (SLURM)

For HPC environments, use the provided SLURM script to analyze all 24 Gold Standard RBPs:

```bash
# Submit array job (all 24 RBPs)
sbatch scripts/submit_job.sh

# Submit subset (positive controls only)
sbatch --array=0-1 scripts/submit_job.sh

# Use HepG2 cell line instead of K562
CELL_LINE=HepG2 sbatch scripts/submit_job.sh
```

The script is pre-configured for:
- Base directory: `/mnt/vstor/SOM_CCCC_JGS25/sayyar/RNAMod-RBNS`
- 4 CPUs, 16GB memory per job
- 4 hour time limit

### Command-Line Options

| Option | Required | Description |
|--------|----------|-------------|
| `--rbp` | Yes | RBP name (e.g., IGF2BP1) |
| `--rbns` | Yes | Path to RBNS Z-scores CSV |
| `--eclip` | Yes | Path to eCLIP narrowPeak BED |
| `--genome` | Yes | Path to reference genome FASTA |
| `--chrom-sizes` | Yes | Path to chromosome sizes file |
| `--mods` | Yes | Paths to modification BED files |
| `--mod-names` | Yes | Names for modifications (same order) |
| `--output` | Yes | Output directory |
| `--cell-line` | No | Cell line (K562 or HepG2, default: K562) |
| `--canonical-threshold` | No | Z-score for canonical (default: 3.0) |
| `--discrepant-threshold` | No | Z-score for discrepant (default: 1.5) |
| `--min-enrichment` | No | eCLIP enrichment threshold (default: 2.0) |

## Output Files

```
results/{RBP}/
├── peaks_filtered.bed          # After extension + enrichment filter
├── scored_peaks.csv            # All peaks with Score_max and Score_sum
├── canonical_peaks.bed         # Score_max >= 3.0
├── discrepant_peaks.bed        # Score_max < 1.5
├── intermediate_peaks.bed      # 1.5 <= Score_max < 3.0
├── enrichment_results.csv      # Statistical test results
├── summary.csv                 # Analysis summary
├── analysis.log                # Detailed run log
└── figures/
    ├── zscore_distribution.png
    ├── classification_summary.png
    ├── enrichment_barplot.png
    └── score_comparison.png
```

## The 24 Gold Standard RBPs

RBPs with both RBNS and ENCODE eCLIP data:

| RBP | Domain | Notes |
|-----|--------|-------|
| **IGF2BP1*** | KH | m6A reader (positive control) |
| **IGF2BP2*** | KH | m6A reader (positive control) |
| HNRNPC | RRM | Poly-U binding |
| TIA1 | RRM | Stress granule component |
| HNRNPK | KH | Poly-C binding |
| PCBP2 | KH | Poly-C binding |
| RBFOX2 | RRM | GCAUG motif |
| PTBP3 | RRM | Splicing regulator |
| TARDBP | RRM | TDP-43, ALS-associated |
| QKI | KH | Bipartite motif |
| SRSF1 | RRM | SR protein |
| SRSF9 | RRM | SR protein |
| RBM22 | RRM | Stem-loop binding |
| TRA2A | RRM | Splicing regulator |
| HNRNPL | RRM | AC repeats |
| LIN28B | CSD/ZnF | Let-7 regulation |
| ZNF326 | ZnF | Structure-dependent |
| FUS | RRM/ZnF | ALS-associated |
| MATR3 | RRM | Nuclear matrix |
| HNRNPA1 | RRM | Splicing repressor |
| HNRNPM | RRM | G/U-rich binding |
| NONO | RRM | Paraspeckle component |
| U2AF2 | RRM | 3' splice site |
| EWSR1 | RRM/ZnF | FET family |

*Known m6A readers - used as positive controls

## Methodology

### Peak Classification

Based on [Dominguez et al. (2018)](https://doi.org/10.1016/j.molcel.2018.06.012):

1. **Score_max**: Maximum 5-mer Z-score in peak (single high-affinity sites)
2. **Score_sum**: Sum of all 5-mer Z-scores (avidity/multivalent binding)

| Category | Definition | Interpretation |
|----------|------------|----------------|
| Canonical | Score_max >= 3.0 | Standard binding (explained by sequence) |
| Discrepant | Score_max < 1.5 | Modification-enabled binding |
| Intermediate | 1.5 <= Score_max < 3.0 | Ambiguous |

### Statistical Testing

- **Fisher's exact test**: Compares modification frequency in discrepant vs. canonical peaks
- **Multiple testing correction**: Benjamini-Hochberg FDR (adjusted p < 0.05)
- **Effect size**: Odds ratio > 1 indicates enrichment in discrepant peaks

### Validation

The pipeline includes positive control validation:
- IGF2BP1 and IGF2BP2 are known m6A readers
- Analysis should show significant m6A enrichment in their discrepant peaks
- Warning logged if positive controls fail validation

## References

- Dominguez D, et al. (2018). Sequence, Structure, and Context Preferences of Human RNA Binding Proteins. *Molecular Cell*, 70(5):854-867. [DOI: 10.1016/j.molcel.2018.06.012](https://doi.org/10.1016/j.molcel.2018.06.012)

- RBNS Pipeline: [github.com/cburgelab/RBNS_pipeline](https://github.com/cburgelab/RBNS_pipeline)

- ENCODE eCLIP: [encodeproject.org/eclip](https://www.encodeproject.org/eclip/)

- REPIC Database: [repicmod.uchicago.edu](https://repicmod.uchicago.edu/)

- RMBase v2.0: [rna.sysu.edu.cn/rmbase](http://rna.sysu.edu.cn/rmbase/)

## License

MIT License

## Citation

If you use this pipeline, please cite:

```
RNAMod-RBNS: Epitranscriptomic Determinants of RBP Specificity
https://github.com/ShivaAyyar/RNAMod-RBNS
```
