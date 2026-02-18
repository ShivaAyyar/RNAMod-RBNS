# RNAMod-RBNS Analysis Results Summary

**Analysis Date:** 2026-02-16
**Cell Line:** K562
**RBPs Analyzed:** 14 of 15 (PCBP2 pending)

---

## Overview

This analysis investigates the "Specificity Paradox" - where RNA Binding Proteins (RBPs) bind sequences in vivo (eCLIP) that show low affinity in vitro (RBNS). The hypothesis is that RNA modifications (m6A, pseudoU, m5C, ac4C) present in cellular RNA but absent in synthetic RBNS libraries explain this discrepancy.

### Classification Criteria

| Category | Score_max Threshold | Interpretation |
|----------|---------------------|----------------|
| **Canonical** | >= 3.0 | High eCLIP + High RBNS (binding explained by sequence) |
| **Intermediate** | 1.5 - 3.0 | Moderate RBNS affinity |
| **Discrepant** | < 1.5 | High eCLIP + Low RBNS (potential modification-dependent) |

---

## Summary Statistics

| RBP | Peaks | Canonical | Intermediate | Discrepant | % Discrepant |
|-----|-------|-----------|--------------|------------|--------------|
| **IGF2BP1*** | 5,706 | 0 | 5,597 | 109 | 1.9% |
| **IGF2BP2*** | 3,834 | 1,494 | 2,286 | 54 | 1.4% |
| EWSR1 | 8,803 | 7,643 | 1,152 | 8 | 0.1% |
| FUS | 4,473 | 0 | 4,461 | 12 | 0.3% |
| HNRNPC | 516 | 510 | 5 | 1 | 0.2% |
| HNRNPK | 4,919 | 4,440 | 471 | 8 | 0.2% |
| HNRNPL | 5,929 | 5,167 | 672 | 90 | 1.5% |
| LIN28B | 6,916 | 3,733 | 2,236 | 947 | 13.7% |
| RBFOX2 | 3,525 | 2,476 | 914 | 135 | 3.8% |
| RBM22 | 991 | 508 | 159 | 324 | 32.7% |
| SRSF9 | 270 | 143 | 120 | 7 | 2.6% |
| TARDBP | 7,010 | 6,915 | 94 | 1 | 0.0% |
| TIA1 | 5,885 | 5,211 | 294 | 380 | 6.5% |
| TRA2A | 1,741 | 1,680 | 46 | 15 | 0.9% |

\* Known m6A readers (positive controls)

### Key Observations

- **RBM22** has the highest discrepancy rate (32.7%), suggesting strong modification-dependent binding
- **LIN28B** shows 13.7% discrepant peaks - notable for potential epitranscriptomic regulation
- **IGF2BP1/2** (positive controls) show relatively low discrepancy rates (~1-2%), suggesting their binding may be well-explained by RBNS even without m6A
- **TARDBP** and **HNRNPC** show extremely high canonical binding (>98%), indicating sequence-specific binding dominates

---

## Figure Descriptions

Each RBP analysis generates the following visualizations:

### 1. Z-score Distribution (`zscore_distribution.png`)
Histogram showing the distribution of maximum 5-mer RBNS Z-scores across all eCLIP peaks.
- **Green dashed line:** Canonical threshold (Z >= 3.0)
- **Red dashed line:** Discrepant threshold (Z < 1.5)
- Peaks between lines are classified as intermediate

### 2. Classification Summary (`classification_summary.png`)
Pie chart showing the proportion of peaks in each category (canonical, intermediate, discrepant).

### 3. Score Comparison (`score_comparison.png`)
Scatter plot of Score_max vs Score_sum, colored by peak category.
- **Score_max:** Maximum Z-score among all 5-mers (identifies single high-affinity sites)
- **Score_sum:** Sum of all Z-scores (captures multivalent/avidity effects)

### 4. Enrichment Barplot (`enrichment_barplot.png`)
*(Only generated when enrichment analysis completes)*
Bar chart showing odds ratios for each modification type (m6A, pseudoU, m5C, ac4C).
- **Green bars:** Statistically significant enrichment (FDR < 0.05)
- **Gray bars:** Not significant

---

## Individual RBP Results

### IGF2BP1 (Positive Control - m6A Reader)

![Z-score Distribution](IGF2BP1/figures/zscore_distribution.png)

**Statistics:** 5,706 peaks | 0 canonical | 5,597 intermediate | 109 discrepant

Notable: Zero canonical peaks indicates IGF2BP1's binding may be dominated by non-sequence factors or the RBNS data doesn't capture its true binding preferences well.

![Classification](IGF2BP1/figures/classification_summary.png)
![Score Comparison](IGF2BP1/figures/score_comparison.png)

---

### IGF2BP2 (Positive Control - m6A Reader)

![Z-score Distribution](IGF2BP2/figures/zscore_distribution.png)

**Statistics:** 3,834 peaks | 1,494 canonical | 2,286 intermediate | 54 discrepant

![Classification](IGF2BP2/figures/classification_summary.png)
![Score Comparison](IGF2BP2/figures/score_comparison.png)
![Enrichment](IGF2BP2/figures/enrichment_barplot.png)

---

### EWSR1

![Z-score Distribution](EWSR1/figures/zscore_distribution.png)

**Statistics:** 8,803 peaks | 7,643 canonical | 1,152 intermediate | 8 discrepant

Highly sequence-specific binding with 86.8% canonical peaks.

![Classification](EWSR1/figures/classification_summary.png)
![Score Comparison](EWSR1/figures/score_comparison.png)
![Enrichment](EWSR1/figures/enrichment_barplot.png)

---

### FUS

![Z-score Distribution](FUS/figures/zscore_distribution.png)

**Statistics:** 4,473 peaks | 0 canonical | 4,461 intermediate | 12 discrepant

Similar to IGF2BP1, no canonical peaks - binding poorly explained by RBNS Z-scores.

![Classification](FUS/figures/classification_summary.png)
![Score Comparison](FUS/figures/score_comparison.png)

---

### HNRNPC

![Z-score Distribution](HNRNPC/figures/zscore_distribution.png)

**Statistics:** 516 peaks | 510 canonical | 5 intermediate | 1 discrepant

Extremely sequence-specific (98.8% canonical).

![Classification](HNRNPC/figures/classification_summary.png)
![Score Comparison](HNRNPC/figures/score_comparison.png)
![Enrichment](HNRNPC/figures/enrichment_barplot.png)

---

### HNRNPK

![Z-score Distribution](HNRNPK/figures/zscore_distribution.png)

**Statistics:** 4,919 peaks | 4,440 canonical | 471 intermediate | 8 discrepant

![Classification](HNRNPK/figures/classification_summary.png)
![Score Comparison](HNRNPK/figures/score_comparison.png)
![Enrichment](HNRNPK/figures/enrichment_barplot.png)

---

### HNRNPL

![Z-score Distribution](HNRNPL/figures/zscore_distribution.png)

**Statistics:** 5,929 peaks | 5,167 canonical | 672 intermediate | 90 discrepant

![Classification](HNRNPL/figures/classification_summary.png)
![Score Comparison](HNRNPL/figures/score_comparison.png)
![Enrichment](HNRNPL/figures/enrichment_barplot.png)

---

### LIN28B

![Z-score Distribution](LIN28B/figures/zscore_distribution.png)

**Statistics:** 6,916 peaks | 3,733 canonical | 2,236 intermediate | 947 discrepant

**High discrepancy rate (13.7%)** - significant candidate for modification-dependent binding.

![Classification](LIN28B/figures/classification_summary.png)
![Score Comparison](LIN28B/figures/score_comparison.png)
![Enrichment](LIN28B/figures/enrichment_barplot.png)

---

### RBFOX2

![Z-score Distribution](RBFOX2/figures/zscore_distribution.png)

**Statistics:** 3,525 peaks | 2,476 canonical | 914 intermediate | 135 discrepant

![Classification](RBFOX2/figures/classification_summary.png)
![Score Comparison](RBFOX2/figures/score_comparison.png)
![Enrichment](RBFOX2/figures/enrichment_barplot.png)

---

### RBM22

![Z-score Distribution](RBM22/figures/zscore_distribution.png)

**Statistics:** 991 peaks | 508 canonical | 159 intermediate | 324 discrepant

**Highest discrepancy rate (32.7%)** - strong candidate for modification-dependent binding mechanism.

![Classification](RBM22/figures/classification_summary.png)
![Score Comparison](RBM22/figures/score_comparison.png)
![Enrichment](RBM22/figures/enrichment_barplot.png)

---

### SRSF9

![Z-score Distribution](SRSF9/figures/zscore_distribution.png)

**Statistics:** 270 peaks | 143 canonical | 120 intermediate | 7 discrepant

![Classification](SRSF9/figures/classification_summary.png)
![Score Comparison](SRSF9/figures/score_comparison.png)
![Enrichment](SRSF9/figures/enrichment_barplot.png)

---

### TARDBP

![Z-score Distribution](TARDBP/figures/zscore_distribution.png)

**Statistics:** 7,010 peaks | 6,915 canonical | 94 intermediate | 1 discrepant

Extremely sequence-specific (98.6% canonical) - TDP-43 binding well-explained by UG-rich motifs.

![Classification](TARDBP/figures/classification_summary.png)
![Score Comparison](TARDBP/figures/score_comparison.png)
![Enrichment](TARDBP/figures/enrichment_barplot.png)

---

### TIA1

![Z-score Distribution](TIA1/figures/zscore_distribution.png)

**Statistics:** 5,885 peaks | 5,211 canonical | 294 intermediate | 380 discrepant

![Classification](TIA1/figures/classification_summary.png)
![Score Comparison](TIA1/figures/score_comparison.png)
![Enrichment](TIA1/figures/enrichment_barplot.png)

---

### TRA2A

![Z-score Distribution](TRA2A/figures/zscore_distribution.png)

**Statistics:** 1,741 peaks | 1,680 canonical | 46 intermediate | 15 discrepant

![Classification](TRA2A/figures/classification_summary.png)
![Score Comparison](TRA2A/figures/score_comparison.png)
![Enrichment](TRA2A/figures/enrichment_barplot.png)

---

## Modification Enrichment Analysis

**Note:** Current enrichment analysis shows no significant overlaps between peaks and modification sites. This may indicate:
1. Coordinate system mismatch between peak BED files and modification BED files
2. Need for window-based overlap rather than exact intersection
3. Sparse modification coverage in K562 cell line data

Further investigation is needed to validate the enrichment analysis pipeline.

---

## Methods Summary

### Pipeline Steps
1. **Load RBNS Z-scores** - Pre-computed 5-mer enrichment Z-scores from ENCODE RBNS
2. **Extend eCLIP peaks** - 50nt 5' extension (strand-aware) to capture binding context
3. **Filter peaks** - Enrichment threshold >= 2.0 (signalValue)
4. **Score peaks** - Calculate Score_max and Score_sum from RBNS Z-scores
5. **Classify peaks** - Apply thresholds (canonical >= 3.0, discrepant < 1.5)
6. **Enrichment analysis** - Fisher's exact test comparing modification overlap in discrepant vs canonical peaks

### Data Sources
- **eCLIP:** ENCODE Project (K562 cell line)
- **RBNS:** ENCODE Project (5-mer enrichment data)
- **Modifications:** RMBase v3.0 (m6A, pseudoU, m5C, ac4C)

---

## File Structure

```
results/
├── RESULTS_SUMMARY.md          # This file
├── {RBP}/
│   ├── summary.csv             # Single-row summary statistics
│   ├── scored_peaks.csv        # All peaks with scores and categories
│   ├── canonical_peaks.bed     # High RBNS affinity peaks
│   ├── discrepant_peaks.bed    # Low RBNS affinity peaks (of interest)
│   ├── intermediate_peaks.bed  # Moderate RBNS affinity peaks
│   ├── peaks_filtered.bed      # All filtered eCLIP peaks
│   ├── enrichment_results.csv  # Modification enrichment statistics
│   └── figures/
│       ├── zscore_distribution.png
│       ├── classification_summary.png
│       ├── score_comparison.png
│       └── enrichment_barplot.png
```

---

*Generated by RNAMod-RBNS Pipeline*
