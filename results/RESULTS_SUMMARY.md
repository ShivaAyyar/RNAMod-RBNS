# RNAMod-RBNS Analysis Results Summary

**Analysis Date:** 2026-02-25 (re-run after pipeline bug fix)
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

| RBP | Peaks | Canonical | Intermediate | Discrepant | % Discrepant | Significant Enrichments |
|-----|-------|-----------|--------------|------------|--------------|------------------------|
| **IGF2BP1*** | 5,706 | 0 | 5,597 | 109 | 1.9% | — (no canonical peaks) |
| **IGF2BP2*** | 3,834 | 1,494 | 2,286 | 54 | 1.4% | none |
| EWSR1 | 8,803 | 7,643 | 1,152 | 8 | 0.1% | none |
| FUS | 4,473 | 0 | 4,461 | 12 | 0.3% | — (no canonical peaks) |
| HNRNPC | 516 | 510 | 5 | 1 | 0.2% | none |
| HNRNPK | 4,919 | 4,440 | 471 | 8 | 0.2% | none |
| HNRNPL | 5,929 | 5,167 | 672 | 90 | 1.5% | none |
| LIN28B | 6,916 | 3,733 | 2,236 | 947 | 13.7% | none |
| **RBFOX2** | 3,525 | 2,476 | 914 | 135 | 3.8% | **m6A (OR=2.6, p_adj=6.7e-5)** |
| **RBM22** | 991 | 508 | 159 | 324 | 32.7% | **pseudoU (OR=16.6, p_adj=1.7e-6)** |
| SRSF9 | 270 | 143 | 120 | 7 | 2.6% | none |
| TARDBP | 7,010 | 6,915 | 94 | 1 | 0.0% | none |
| TIA1 | 5,885 | 5,211 | 294 | 380 | 6.5% | none |
| TRA2A | 1,741 | 1,680 | 46 | 15 | 0.9% | none |

\* Known m6A readers (positive controls)

### Key Observations

- **RBFOX2** shows significant m6A enrichment in discrepant peaks (24.4% disc vs 11.0% canon, OR=2.6, FDR=6.7e-5), supporting modification-dependent binding
- **RBM22** shows the strongest signal: highly significant pseudoU enrichment in discrepant peaks (6.2% disc vs 0.4% canon, OR=16.6, FDR=1.7e-6)
- **RBM22** also has the highest discrepancy rate (32.7%), consistent with broad modification-dependent binding
- **IGF2BP1** (positive control) has 0 canonical peaks — all peaks fall below the canonical threshold (Score_max < 3.0), suggesting IGF2BP1's RBNS 5-mer affinities are systematically lower than expected. Enrichment analysis requires at least one canonical peak and was not run
- **IGF2BP2** (positive control) shows m6A strongly enriched in *canonical* peaks (64.7% vs 18.5%), not discrepant — suggesting IGF2BP2's consensus sequence motif overlaps with m6A sites (the DRACH motif contains the CAU/UAC context IGF2BP2 prefers), so its canonical peaks already capture m6A-proximal binding
- **FUS** has 0 canonical peaks (same issue as IGF2BP1)
- **TARDBP** and **HNRNPC** show >98% canonical binding — sequence affinity fully explains their eCLIP patterns

---

## Modification Enrichment Results

### Significant Findings (FDR < 0.05)

| RBP | Modification | Canon% | Disc% | Odds Ratio | p_adj |
|-----|-------------|--------|-------|------------|-------|
| **RBFOX2** | m6A | 11.0% | 24.4% | 2.61 | 6.7e-5 |
| **RBM22** | pseudoU | 0.4% | 6.2% | 16.6 | 1.7e-6 |

### All Enrichment Results

| RBP | m6A OR | pseudoU OR | m5C OR | ac4C OR | Notes |
|-----|--------|------------|--------|---------|-------|
| IGF2BP1 | — | — | — | — | No canonical peaks; not run |
| IGF2BP2 | 0.12 | 1.11 | 2.21 | 0 | m6A enriched in canonical (inverted) |
| EWSR1 | 2.21 | 0 | 0 | 0 | Only 8 discrepant peaks (low power) |
| FUS | — | — | — | — | No canonical peaks; not run |
| HNRNPC | 0 | 0 | 0 | — | Only 1 discrepant peak |
| HNRNPK | 0 | 0 | 4.74 | 0 | Only 8 discrepant peaks |
| HNRNPL | 0.22 | 0 | 0 | 0 | Low discrepant overlap |
| LIN28B | 0.89 | 0.71 | 0.82 | 0.32 | High absolute overlap; no enrichment |
| **RBFOX2** | **2.61 ✓** | 0 | 3.75 | 0 | **m6A significant** |
| **RBM22** | 1.01 | **16.6 ✓** | 0.62 | 0.52 | **pseudoU significant** |
| SRSF9 | 0.43 | 0 | 0 | 0 | Only 7 discrepant peaks |
| TARDBP | 0 | 0 | 0 | 0 | Only 1 discrepant peak |
| TIA1 | 0.27 | 0 | 0.68 | 0 | m6A depleted in discrepant |
| TRA2A | 0.055 | 0 | 0 | 0 | Only 15 discrepant peaks |

✓ = FDR < 0.05

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
*(Not generated for IGF2BP1 or FUS — no canonical peaks)*
Bar chart showing odds ratios for each modification type (m6A, pseudoU, m5C, ac4C).
- **Green bars:** Statistically significant enrichment (FDR < 0.05)
- **Gray bars:** Not significant

---

## Individual RBP Results

### IGF2BP1 (Positive Control - m6A Reader)

![Z-score Distribution](IGF2BP1/figures/zscore_distribution.png)

**Statistics:** 5,706 peaks | 0 canonical | 5,597 intermediate | 109 discrepant

**Note:** Zero canonical peaks — all 5,706 eCLIP peaks score below the canonical threshold (Score_max < 3.0). Enrichment analysis requires canonical peaks and was not run. IGF2BP1's RBNS Z-scores may be systematically lower than expected, or the Z ≥ 3.0 threshold may need recalibration for this protein.

![Classification](IGF2BP1/figures/classification_summary.png)
![Score Comparison](IGF2BP1/figures/score_comparison.png)

---

### IGF2BP2 (Positive Control - m6A Reader)

![Z-score Distribution](IGF2BP2/figures/zscore_distribution.png)

**Statistics:** 3,834 peaks | 1,494 canonical | 2,286 intermediate | 54 discrepant

**Enrichment:** m6A is strongly enriched in canonical peaks (64.7%) relative to discrepant (18.5%), OR=0.12. This is the inverse of the hypothesis, suggesting IGF2BP2's sequence motif overlaps with m6A modification sites — its preferred binding context (CAU/UAC-containing sequences) frequently co-occurs with DRACH m6A sites.

![Classification](IGF2BP2/figures/classification_summary.png)
![Score Comparison](IGF2BP2/figures/score_comparison.png)
![Enrichment](IGF2BP2/figures/enrichment_barplot.png)

---

### EWSR1

![Z-score Distribution](EWSR1/figures/zscore_distribution.png)

**Statistics:** 8,803 peaks | 7,643 canonical | 1,152 intermediate | 8 discrepant

Highly sequence-specific binding (86.8% canonical). Only 8 discrepant peaks — insufficient power for enrichment analysis.

![Classification](EWSR1/figures/classification_summary.png)
![Score Comparison](EWSR1/figures/score_comparison.png)
![Enrichment](EWSR1/figures/enrichment_barplot.png)

---

### FUS

![Z-score Distribution](FUS/figures/zscore_distribution.png)

**Statistics:** 4,473 peaks | 0 canonical | 4,461 intermediate | 12 discrepant

**Note:** Zero canonical peaks (same situation as IGF2BP1). Enrichment analysis not run. FUS Z-scores cluster below the canonical threshold; the RBNS data may underrepresent FUS's sequence preferences, or FUS binding may be broadly context-dependent rather than driven by a single high-affinity motif.

![Classification](FUS/figures/classification_summary.png)
![Score Comparison](FUS/figures/score_comparison.png)

---

### HNRNPC

![Z-score Distribution](HNRNPC/figures/zscore_distribution.png)

**Statistics:** 516 peaks | 510 canonical | 5 intermediate | 1 discrepant

Extremely sequence-specific (98.8% canonical). Only 1 discrepant peak.

![Classification](HNRNPC/figures/classification_summary.png)
![Score Comparison](HNRNPC/figures/score_comparison.png)
![Enrichment](HNRNPC/figures/enrichment_barplot.png)

---

### HNRNPK

![Z-score Distribution](HNRNPK/figures/zscore_distribution.png)

**Statistics:** 4,919 peaks | 4,440 canonical | 471 intermediate | 8 discrepant

Highly sequence-specific (90.3% canonical). Only 8 discrepant peaks.

![Classification](HNRNPK/figures/classification_summary.png)
![Score Comparison](HNRNPK/figures/score_comparison.png)
![Enrichment](HNRNPK/figures/enrichment_barplot.png)

---

### HNRNPL

![Z-score Distribution](HNRNPL/figures/zscore_distribution.png)

**Statistics:** 5,929 peaks | 5,167 canonical | 672 intermediate | 90 discrepant

Predominantly sequence-specific (87.1% canonical). Low modification overlap in discrepant peaks.

![Classification](HNRNPL/figures/classification_summary.png)
![Score Comparison](HNRNPL/figures/score_comparison.png)
![Enrichment](HNRNPL/figures/enrichment_barplot.png)

---

### LIN28B

![Z-score Distribution](LIN28B/figures/zscore_distribution.png)

**Statistics:** 6,916 peaks | 3,733 canonical | 2,236 intermediate | 947 discrepant

High discrepancy rate (13.7%) with 947 discrepant peaks. m6A overlap is high in both canonical (45.7%) and discrepant (42.9%) peaks, suggesting LIN28B's eCLIP peaks broadly co-occur with m6A-modified transcripts but without enrichment in the discrepant set specifically.

![Classification](LIN28B/figures/classification_summary.png)
![Score Comparison](LIN28B/figures/score_comparison.png)
![Enrichment](LIN28B/figures/enrichment_barplot.png)

---

### RBFOX2 ✓

![Z-score Distribution](RBFOX2/figures/zscore_distribution.png)

**Statistics:** 3,525 peaks | 2,476 canonical | 914 intermediate | 135 discrepant

**Significant m6A enrichment in discrepant peaks** (24.4% vs 11.0%, OR=2.61, FDR=6.7e-5). Discrepant RBFOX2 peaks are approximately twice as likely to overlap an m6A site as canonical peaks, consistent with m6A-dependent binding at a subset of sites.

![Classification](RBFOX2/figures/classification_summary.png)
![Score Comparison](RBFOX2/figures/score_comparison.png)
![Enrichment](RBFOX2/figures/enrichment_barplot.png)

---

### RBM22 ✓

![Z-score Distribution](RBM22/figures/zscore_distribution.png)

**Statistics:** 991 peaks | 508 canonical | 159 intermediate | 324 discrepant

**Strongest signal in the dataset: highly significant pseudoU enrichment in discrepant peaks** (6.2% vs 0.4%, OR=16.6, FDR=1.7e-6). Discrepant RBM22 peaks are ~17× more likely to overlap pseudouridine sites than canonical peaks. RBM22 also has the highest discrepancy rate (32.7%) among all RBPs. m6A overlap is similar in canonical and discrepant peaks (OR≈1), suggesting pseudoU rather than m6A is the relevant modification.

![Classification](RBM22/figures/classification_summary.png)
![Score Comparison](RBM22/figures/score_comparison.png)
![Enrichment](RBM22/figures/enrichment_barplot.png)

---

### SRSF9

![Z-score Distribution](SRSF9/figures/zscore_distribution.png)

**Statistics:** 270 peaks | 143 canonical | 120 intermediate | 7 discrepant

Small dataset. Only 7 discrepant peaks; insufficient power for enrichment analysis.

![Classification](SRSF9/figures/classification_summary.png)
![Score Comparison](SRSF9/figures/score_comparison.png)
![Enrichment](SRSF9/figures/enrichment_barplot.png)

---

### TARDBP

![Z-score Distribution](TARDBP/figures/zscore_distribution.png)

**Statistics:** 7,010 peaks | 6,915 canonical | 94 intermediate | 1 discrepant

Extremely sequence-specific (98.6% canonical). TDP-43 binding is almost entirely explained by its UG-rich sequence motif. Only 1 discrepant peak.

![Classification](TARDBP/figures/classification_summary.png)
![Score Comparison](TARDBP/figures/score_comparison.png)
![Enrichment](TARDBP/figures/enrichment_barplot.png)

---

### TIA1

![Z-score Distribution](TIA1/figures/zscore_distribution.png)

**Statistics:** 5,885 peaks | 5,211 canonical | 294 intermediate | 380 discrepant

Predominantly sequence-specific (88.5% canonical). m6A is notably depleted in discrepant relative to canonical peaks (9.7% vs 28.9%, OR=0.27), suggesting TIA1's discrepant peaks occur in regions with low m6A density.

![Classification](TIA1/figures/classification_summary.png)
![Score Comparison](TIA1/figures/score_comparison.png)
![Enrichment](TIA1/figures/enrichment_barplot.png)

---

### TRA2A

![Z-score Distribution](TRA2A/figures/zscore_distribution.png)

**Statistics:** 1,741 peaks | 1,680 canonical | 46 intermediate | 15 discrepant

Highly sequence-specific (96.5% canonical). m6A strongly enriched in canonical peaks (73.8%), consistent with TRA2A's GAA-repeat binding preference co-occurring with m6A. Only 15 discrepant peaks.

![Classification](TRA2A/figures/classification_summary.png)
![Score Comparison](TRA2A/figures/score_comparison.png)
![Enrichment](TRA2A/figures/enrichment_barplot.png)

---

## Methods Summary

### Pipeline Steps
1. **Load RBNS Z-scores** - Pre-computed 5-mer enrichment Z-scores from ENCODE RBNS
2. **Extend eCLIP peaks** - 50nt 5' extension (strand-aware) to capture binding context
3. **Filter peaks** - Enrichment threshold >= 2.0 (signalValue)
4. **Score peaks** - Calculate Score_max and Score_sum from RBNS Z-scores
5. **Classify peaks** - Apply thresholds (canonical >= 3.0, discrepant < 1.5)
6. **Enrichment analysis** - Fisher's exact test (one-sided, greater) comparing modification overlap in discrepant vs canonical peaks; Benjamini-Hochberg FDR correction

### Data Sources
- **eCLIP:** ENCODE Project (K562 cell line)
- **RBNS:** ENCODE Project (5-mer enrichment data)
- **Modifications:** RMBase v3.0 — all four modification types use cumulative multi-study sites (hg38). m6A for HepG2 is cell-line-specific (178,800 sites); m6A for K562 is cumulative multi-study with support ≥2 (558,360 sites) as K562-specific m6A data is absent from RMBase v3.0. Ψ, m5C, and ac4C use the same cumulative multi-study sites for both cell lines.

### Known Limitations
- IGF2BP1 and FUS have no canonical peaks under the current Z ≥ 3.0 threshold; enrichment analysis requires a canonical reference set and was not performed for these two RBPs
- Several RBPs (TARDBP, HNRNPC, EWSR1, SRSF9) have very few discrepant peaks (1–8), giving low statistical power for enrichment tests
- K562 m6A data is cumulative multi-study rather than K562-specific, which may reduce specificity of the m6A enrichment signal

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

*Generated by RNAMod-RBNS Pipeline — updated 2026-02-25 after chromosome name bug fix*
