# Data Download Instructions

This document provides step-by-step instructions for downloading all required data for the RNAMod-RBNS pipeline.

## Overview

| Data Type | Source | Auto-Download | Manual Required |
|-----------|--------|---------------|-----------------|
| Reference Genome (hg38) | UCSC | Yes | No |
| Chromosome Sizes | UCSC | Yes | No |
| eCLIP Peaks | ENCODE | **Yes** | No |
| RBNS Z-scores | Dominguez et al. (2018) | No | Yes |
| m6A Sites | REPIC | No | Yes |
| Ψ, m5C, ac4C Sites | RMBase | No | Yes |

## The 23 RBPs

Based on Dominguez et al. (2018) with ENCODE availability verification:

| RBP | K562 | HepG2 | Notes |
|-----|------|-------|-------|
| **IGF2BP1** | ✓ | ✓ | m6A reader (positive control) |
| **IGF2BP2** | ✓ | - | m6A reader (positive control), K562 only |
| HNRNPC | ✓ | ✓ | |
| TIA1 | ✓ | ✓ | |
| HNRNPK | ✓ | ✓ | |
| PCBP2 | - | ✓ | HepG2 only |
| RBFOX2 | ✓ | ✓ | |
| PTBP1 | ✓ | ✓ | Replaces PTBP3 (not in ENCODE) |
| TARDBP | ✓ | ✓ | |
| QKI | ✓ | ✓ | |
| SRSF1 | ✓ | ✓ | |
| SRSF9 | ✓ | ✓ | |
| RBM22 | ✓ | ✓ | |
| TRA2A | ✓ | ✓ | |
| HNRNPL | ✓ | ✓ | |
| LIN28B | ✓ | ✓ | |
| FUS | ✓ | ✓ | |
| MATR3 | ✓ | ✓ | |
| HNRNPA1 | ✓ | ✓ | |
| HNRNPM | ✓ | ✓ | |
| NONO | ✓ | - | K562 only |
| U2AF2 | ✓ | ✓ | |
| EWSR1 | ✓ | - | K562 only |

**Note:** ZNF326 and PTBP3 are not available in ENCODE and have been excluded.

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

**NEW: Use the automated download script:**

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

## Step 3: Download RBNS Z-scores (Manual)

### Source: Dominguez et al. (2018) Supplementary Data

1. Go to: https://www.cell.com/molecular-cell/supplemental/S1097-2765(18)30437-3

2. Download **Table S2** (Excel file containing 5-mer Z-scores)

3. For each RBP, extract the k-mer Z-scores and save as CSV:

**Required CSV format:**
```csv
kmer,z_score
AAAAA,0.23
AAAAC,0.45
AAAAG,0.12
...
```

4. Save each RBP's data to: `data/rbns/{RBP}_zscores.csv`

### Notes on RBNS Data

- The supplementary table contains Z-scores for 5-mers (1024 possible k-mers)
- Z-scores represent enrichment relative to input library
- Higher Z-scores indicate stronger binding affinity
- Use DNA notation (T, not U) - the pipeline converts automatically

## Step 4: Download RNA Modification Sites (Manual)

### 4a. m6A Sites from REPIC

1. Go to: https://repicmod.uchicago.edu/download

2. Select:
   - Species: **Human**
   - Cell Line: **K562** (and repeat for **HepG2**)
   - Modification: **m6A**

3. Download the consensus peaks (BED format)

4. Save to: `data/mods/K562/m6A.bed` and `data/mods/HepG2/m6A.bed`

### 4b. Pseudouridine (Ψ), m5C, ac4C from RMBase

1. Go to: http://rna.sysu.edu.cn/rmbase/download.html

2. For each modification type, select:
   - Species: **Homo sapiens**
   - Assembly: **hg38**
   - Cell Line: **K562** or **HepG2**
   - Modification: **pseudoU**, **m5C**, or **ac4C**

3. Download the BED file

4. Save to appropriate location:
   - `data/mods/K562/pseudoU.bed`
   - `data/mods/K562/m5C.bed`
   - `data/mods/K562/ac4C.bed`
   - (Repeat for HepG2)

### Expected BED Format

Modification BED files should have at least 6 columns:
```
chr1    1000    1001    mod_site_1    .    +
chr1    2000    2001    mod_site_2    .    -
```

## Step 5: Verify Downloads

After completing all downloads, verify the directory structure:

```bash
cd /mnt/vstor/SOM_CCCC_JGS25/sayyar/RNAMod-RBNS

# Check eCLIP files
ls -la data/eclip/*.bed | wc -l  # Should be 41

# Check RBNS files
ls -la data/rbns/*.csv | wc -l   # Should be 23

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
│   ├── IGF2BP1_zscores.csv
│   ├── IGF2BP2_zscores.csv
│   └── ... (23 files total)
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
# Run for all 23 RBPs (K562)
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
- Ensure CSV has header row with `kmer` and `z_score` columns
- Check for extra whitespace or special characters

### Cell-line mismatch
- PCBP2 only has HepG2 data
- IGF2BP2, NONO, EWSR1 only have K562 data
- Analysis will skip missing combinations

## References

- Dominguez D, et al. (2018). Sequence, Structure, and Context Preferences of Human RNA Binding Proteins. *Molecular Cell*, 70(5):854-867. [DOI: 10.1016/j.molcel.2018.06.012](https://doi.org/10.1016/j.molcel.2018.06.012)
- ENCODE Project: https://www.encodeproject.org/
- REPIC Database: https://repicmod.uchicago.edu/
- RMBase v2.0: http://rna.sysu.edu.cn/rmbase/
