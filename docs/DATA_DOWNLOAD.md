# Data Download Instructions

This document provides step-by-step instructions for downloading all required data for the RNAMod-RBNS pipeline.

## Overview

| Data Type | Source | Auto-Download | Manual Required |
|-----------|--------|---------------|-----------------|
| Reference Genome (hg38) | UCSC | Yes | No |
| Chromosome Sizes | UCSC | Yes | No |
| eCLIP Peaks | ENCODE | **Yes** | No |
| RBNS R-values/Z-scores | ENCODE | **Yes** | No |
| m6A Sites | REPIC | No | Yes |
| Ψ, m5C, ac4C Sites | RMBase | No | Yes |

## The 15 Analyzable RBPs

Only **15 RBPs** have BOTH eCLIP and RBNS data available in ENCODE:

| RBP | K562 | HepG2 | Notes |
|-----|------|-------|-------|
| **IGF2BP1** | ✓ | ✓ | m6A reader (positive control) |
| **IGF2BP2** | ✓ | - | m6A reader (positive control), K562 only |
| HNRNPC | ✓ | ✓ | |
| TIA1 | ✓ | ✓ | |
| HNRNPK | ✓ | ✓ | |
| PCBP2 | - | ✓ | HepG2 only |
| RBFOX2 | ✓ | ✓ | |
| TARDBP | ✓ | ✓ | |
| SRSF9 | ✓ | ✓ | |
| RBM22 | ✓ | ✓ | |
| TRA2A | ✓ | ✓ | |
| HNRNPL | ✓ | ✓ | |
| LIN28B | ✓ | ✓ | |
| FUS | ✓ | ✓ | |
| EWSR1 | ✓ | - | K562 only |

### RBPs WITHOUT RBNS Data (Cannot Analyze)

The following 8 RBPs have eCLIP data but **lack RBNS data** in ENCODE:
- PTBP1 (only PTBP3 has RBNS)
- QKI
- SRSF1
- MATR3
- HNRNPA1 (only HNRNPA1L2 has RBNS)
- HNRNPM
- NONO
- U2AF2

## Step 1: Download Reference Genome (Automated)

Run the genome download script:

```bash
# Submit as SLURM job
sbatch scripts/download_data.sh

# Or run interactively
bash scripts/download_data.sh
```

This downloads:
- hg38 reference genome (~3GB)
- Chromosome sizes file
- Creates FASTA index

## Step 2: Download eCLIP Data (Automated)

**Use the automated download script:**

```bash
# Submit as SLURM job (recommended)
sbatch scripts/download_eclip.sh

# Or run interactively
bash scripts/download_eclip.sh
```

This script contains pre-verified ENCODE file accessions for all 23 RBPs and will download:
- 22 K562 eCLIP files
- 19 HepG2 eCLIP files
- Total: 41 BED files with IDR-thresholded peaks

### eCLIP File Accessions Reference

For reference, here are the ENCODE accessions used:

**K562:**
| RBP | Accession |
|-----|-----------|
| IGF2BP1 | ENCFF650LMV |
| IGF2BP2 | ENCFF524ZZB |
| HNRNPC | ENCFF167CDB |
| TIA1 | ENCFF918KMT |
| HNRNPK | ENCFF318RSO |
| RBFOX2 | ENCFF206RIM |
| PTBP1 | ENCFF907HNN |
| TARDBP | ENCFF037TVC |
| QKI | ENCFF786UOW |
| SRSF1 | ENCFF223KVR |
| SRSF9 | ENCFF781BNS |
| RBM22 | ENCFF972ZMJ |
| TRA2A | ENCFF726PFJ |
| HNRNPL | ENCFF917CBK |
| LIN28B | ENCFF061XNA |
| FUS | ENCFF861KMV |
| MATR3 | ENCFF246EPM |
| HNRNPA1 | ENCFF392AEV |
| HNRNPM | ENCFF445ENC |
| NONO | ENCFF730QRI |
| U2AF2 | ENCFF290DFO |
| EWSR1 | ENCFF607ZRF |

**HepG2:**
| RBP | Accession |
|-----|-----------|
| IGF2BP1 | ENCFF442USD |
| HNRNPC | ENCFF440ROZ |
| TIA1 | ENCFF759KCD |
| HNRNPK | ENCFF855CPQ |
| PCBP2 | ENCFF642GNE |
| RBFOX2 | ENCFF871NYM |
| PTBP1 | ENCFF726SQU |
| TARDBP | ENCFF673QBV |
| QKI | ENCFF704OCI |
| SRSF1 | ENCFF934ANS |
| SRSF9 | ENCFF765PIF |
| RBM22 | ENCFF293IZG |
| TRA2A | ENCFF766OCH |
| HNRNPL | ENCFF266TKW |
| LIN28B | ENCFF341XMP |
| FUS | ENCFF972DFZ |
| MATR3 | ENCFF587KKM |
| HNRNPA1 | ENCFF797GSK |
| HNRNPM | ENCFF752JNY |
| U2AF2 | ENCFF721PWF |

## Step 3: Download RBNS Data (Automated)

**NEW: Use the automated ENCODE RBNS download script:**

```bash
# Submit as SLURM job (recommended)
sbatch scripts/download_rbns.sh

# Or run interactively
bash scripts/download_rbns.sh
```

This script downloads pre-computed 5-mer R-value enrichment files from ENCODE for all 15 RBPs with RBNS data.

### Process RBNS Files

After downloading, convert R-values to Z-scores:

```bash
# Convert all downloaded enrichment files to Z-scores
python scripts/process_rbns_enrichment.py

# Or process a single file
python scripts/process_rbns_enrichment.py --single data/rbns/IGF2BP1_enrichment.tsv data/rbns/IGF2BP1_zscores.csv

# Validate output files
python scripts/process_rbns_enrichment.py --validate data/rbns
```

### RBNS File Accessions Reference

| RBP | Experiment | 5-mer File |
|-----|------------|------------|
| IGF2BP1 | ENCSR928XOW | ENCFF157MJN |
| IGF2BP2 | ENCSR588GYZ | ENCFF446BRC |
| HNRNPC | ENCSR569UIU | ENCFF592FMJ |
| TIA1 | ENCSR064NOY | ENCFF936GUV |
| HNRNPK | ENCSR368NMO | ENCFF715DIL |
| PCBP2 | ENCSR673FLQ | ENCFF761GIK |
| RBFOX2 | ENCSR441HLP | ENCFF002DFE |
| TARDBP | ENCSR466JPT | ENCFF022PDN |
| SRSF9 | ENCSR724HZI | ENCFF518QQO |
| RBM22 | ENCSR006TPX | ENCFF688EDA |
| TRA2A | ENCSR741VUK | ENCFF191NLM |
| HNRNPL | ENCSR954TYO | ENCFF422OAC |
| LIN28B | ENCSR369RLA | ENCFF189NSX |
| FUS | ENCSR936LOF | ENCFF835DNX |
| EWSR1 | ENCSR063HQO | ENCFF100OEJ |

### ENCODE RBNS File Format

The downloaded enrichment TSV files have this format:
```
[RBP_NAME]    0    5    20    80    320    1300
AAAAA    0.98    1.00    1.02    1.05    1.08    1.10
AAAAC    0.97    0.99    1.01    1.04    1.07    1.09
...
```

- First row: header with RBP name and protein concentrations (nM)
- First column: 5-mer sequence
- Values: R-values (enrichment ratios relative to input)
- ~970 5-mers per file (1024 minus adapter-containing k-mers)

### Z-score Calculation

The processing script converts R-values to Z-scores using:

```
Z = (R - mean_R) / std_R
```

Where R is the R-value at the highest protein concentration (e.g., 1300 nM).

Output CSV format:
```csv
kmer,z_score,r_value
GCAUG,8.45,3.25
GCACG,7.23,2.89
...
```

## Step 4: Download and Process RNA Modification Sites

All modification data comes from **RMBase v3.0** (http://rna.sysu.edu.cn/rmbase/).

### 4a. Download RMBase Files

1. Go to: http://rna.sysu.edu.cn/rmbase/download.html

2. Download the following files for **Human (hg38)**:
   - `hg38.m6A.tar.gz` (~63 MB)
   - `hg38.Pseudo.tar.gz` (~360 KB) - Pseudouridine (Ψ)
   - `hg38.m5C.tar.gz` (~2.3 MB)
   - `hg38.ac4C.tar.gz` (~195 KB)

3. Save all files to: `data/mods/`

### 4b. Process RMBase Files

Run the processing script to extract cell-line-specific BED files:

```bash
# Process all modifications for K562 and HepG2
python scripts/process_rmbase_mods.py --setup-all
```

This script will:
- Extract the tar.gz archives
- Filter m6A for HepG2-specific sites (178,800 sites)
- For K562: use all m6A sites with support ≥2 (558,360 sites) since K562 m6A is not profiled in RMBase
- Use all Ψ, m5C, ac4C sites for both cell lines (not cell-line-specific in RMBase)

### Data Availability Summary

| Modification | K562 Sites | HepG2 Sites | Notes |
|--------------|------------|-------------|-------|
| **m6A** | 558,360 | 178,800 | HepG2 cell-specific; K562 uses all sites |
| **Ψ (Pseudo)** | 5,705 | 5,705 | Same sites (not cell-specific) |
| **m5C** | 46,025 | 46,025 | Same sites (not cell-specific) |
| **ac4C** | 1,861 | 1,861 | Same sites (not cell-specific) |

### Output BED6 Format

```
chr1    14414    14415    m6A_site_1    674    -
chr1    14626    14627    m6A_site_4    728    -
```

Columns: chrom, start, end, modID, score (motif_score * 200), strand

## Step 5: Verify Downloads

After completing all downloads, verify the directory structure:

```bash
cd /mnt/vstor/SOM_CCCC_JGS25/sayyar/RNAMod-RBNS

# Check eCLIP files
ls -la data/eclip/*.bed | wc -l  # Should be 41

# Check RBNS enrichment files (downloaded)
ls -la data/rbns/*_enrichment.tsv | wc -l  # Should be 15

# Check RBNS Z-score files (after processing)
ls -la data/rbns/*_zscores.csv | wc -l     # Should be 15

# Check modification files
ls -la data/mods/K562/*.bed      # Should be 4
ls -la data/mods/HepG2/*.bed     # Should be 4
```

Expected structure:
```
data/
├── genome/
│   ├── hg38.fa
│   ├── hg38.fa.fai
│   └── hg38.chrom.sizes
├── eclip/
│   ├── IGF2BP1_K562.bed
│   ├── IGF2BP1_HepG2.bed
│   ├── IGF2BP2_K562.bed      # K562 only
│   └── ... (41 files total)
├── rbns/
│   ├── IGF2BP1_enrichment.tsv  # Downloaded from ENCODE
│   ├── IGF2BP1_zscores.csv     # Processed Z-scores
│   ├── IGF2BP2_enrichment.tsv
│   ├── IGF2BP2_zscores.csv
│   └── ... (15 RBPs with both files)
└── mods/
    ├── K562/
    │   ├── m6A.bed
    │   ├── pseudoU.bed
    │   ├── m5C.bed
    │   └── ac4C.bed
    └── HepG2/
        ├── m6A.bed
        ├── pseudoU.bed
        ├── m5C.bed
        └── ac4C.bed
```

## Step 6: Run Analysis

Once all data is in place:

```bash
# Run for all 15 RBPs (K562)
sbatch scripts/submit_job.sh

# Run just the positive controls first
sbatch --array=0-1 scripts/submit_job.sh

# Run with HepG2 (note: some RBPs K562-only)
CELL_LINE=HepG2 sbatch scripts/submit_job.sh
```

## Troubleshooting

### "File not found" errors
- Check that filenames match exactly (case-sensitive)
- Verify BED files are uncompressed (not .gz)

### "No peaks" or empty results
- Ensure eCLIP files contain IDR-thresholded peaks, not raw reads
- Check that modification BED files use hg38 coordinates

### RBNS Z-score parsing errors
- Ensure CSV has header row with `kmer`, `z_score`, and `r_value` columns
- Run the validation: `python scripts/process_rbns_enrichment.py --validate data/rbns`
- Check for extra whitespace or special characters

### Cell-line mismatch
- PCBP2 only has HepG2 data
- IGF2BP2, EWSR1 only have K562 data
- Analysis will skip missing combinations

### RBNS data not available
- 8 RBPs (PTBP1, QKI, SRSF1, MATR3, HNRNPA1, HNRNPM, NONO, U2AF2) lack RBNS data
- These cannot be analyzed with this pipeline

## References

- Dominguez D, et al. (2018). Sequence, Structure, and Context Preferences of Human RNA Binding Proteins. *Molecular Cell*, 70(5):854-867. [DOI: 10.1016/j.molcel.2018.06.012](https://doi.org/10.1016/j.molcel.2018.06.012)
- ENCODE Project: https://www.encodeproject.org/
- REPIC Database: https://repicmod.uchicago.edu/
- RMBase v2.0: http://rna.sysu.edu.cn/rmbase/
