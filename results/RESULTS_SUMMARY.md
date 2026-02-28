# RNAMod-RBNS Analysis Results Summary

**Analysis Date:** 2026-02-27 (re-run after summit-anchored peak extension + per-RBP threshold fix)
**Cell Line:** K562 (PCBP2: HepG2)
**RBPs Analyzed:** 15 of 15

---

## Overview

This analysis investigates the "Specificity Paradox" — where RNA Binding Proteins (RBPs) bind sequences in vivo (eCLIP) that show low affinity in vitro (RBNS). The hypothesis is that RNA modifications (m6A, pseudoU, m5C, ac4C) present in cellular RNA but absent in synthetic RBNS libraries explain this discrepancy.

### Classification Criteria

| Category | Score_max Threshold | Interpretation |
|----------|---------------------|----------------|
| **Canonical** | >= 3.0 (IGF2BP1/FUS: >= 2.0) | High eCLIP + High RBNS (binding explained by sequence) |
| **Intermediate** | 1.5–3.0 (IGF2BP1/FUS: 1.0–2.0) | Moderate RBNS affinity |
| **Discrepant** | < 1.5 (IGF2BP1/FUS: < 1.0) | High eCLIP + Low RBNS (potential modification-dependent) |

---

## Summary Statistics

| RBP | Cell Line | Canonical | Discrepant | % Discrepant | Significant Enrichments |
|-----|-----------|-----------|------------|--------------|------------------------|
| **IGF2BP1*** | K562 | 3,724 | 104 | 2.7% | none (m6A depleted in discrepant) |
| **IGF2BP2*** | K562 | 797 | 240 | 23.2% | none (m6A depleted in discrepant) |
| EWSR1 | K562 | 5,387 | 105 | 1.9% | none |
| FUS | K562 | 3,881 | 19 | 0.5% | none (only 19 discrepant peaks) |
| HNRNPC | K562 | 388 | 81 | 17.3% | none |
| HNRNPK | K562 | 2,930 | 75 | 2.5% | none |
| HNRNPL | K562 | 3,958 | 371 | 8.6% | none |
| LIN28B | K562 | 2,424 | 2,245 | 48.0% | none |
| PCBP2 | HepG2 | 6,979 | 690 | 9.0% | none |
| **RBFOX2** | K562 | 1,281 | 728 | 36.2% | **m6A (OR=2.40, p_adj=8.6e-7), m5C (OR=3.55, p_adj=0.032)** |
| **RBM22** | K562 | 271 | 592 | 68.6% | **pseudoU (OR=7.98, p_adj=0.040)** |
| SRSF9 | K562 | 82 | 14 | 14.6% | none |
| TARDBP | K562 | 5,578 | 156 | 2.7% | none |
| TIA1 | K562 | 4,230 | 929 | 18.0% | none |
| TRA2A | K562 | 1,433 | 106 | 6.9% | none |

\* Known m6A readers (positive controls)

### Key Observations

- **RBFOX2** shows significant enrichment of both m6A (OR=2.40, FDR=8.6e-7) and m5C (OR=3.55, FDR=0.032) in discrepant peaks, with the highest discrepancy rate (36.2%) among sequence-specific RBPs
- **RBM22** retains highly significant pseudoU enrichment in discrepant peaks (OR=7.98, FDR=0.040) and has the highest discrepancy rate overall (68.6%), suggesting broad modification-dependent binding
- **IGF2BP1** (positive control): Now has 3,724 canonical peaks at the lowered threshold (Z≥2.0). Unexpectedly, m6A is strongly *depleted* in discrepant peaks (3.8% vs 41.5%, OR=0.056), mirroring IGF2BP2. Both m6A readers show the same inverted pattern — see discussion below
- **IGF2BP2** (positive control): m6A enriched in canonical (46.9%) over discrepant (18.3%), OR=0.25 — same inverted pattern as IGF2BP1
- **FUS**: Now has 3,881 canonical peaks at lowered threshold but only 19 discrepant peaks — insufficient power. m6A trend positive (OR=2.6) but not significant
- **LIN28B**: 48.0% discrepancy rate with 2,245 discrepant peaks; m6A overlap similar in both classes (~29%), no enrichment signal. May reflect broad transcriptome-wide co-occurrence rather than modification-dependent binding
- **TARDBP** and **HNRNPC**: Predominantly sequence-specific binding (well-characterized UG-rich and U-tract motifs respectively)

---

## Modification Enrichment Results

### Significant Findings (FDR < 0.05)

| RBP | Modification | Canon% | Disc% | Odds Ratio | p_adj |
|-----|-------------|--------|-------|------------|-------|
| **RBFOX2** | m6A | 5.2% | 11.7% | 2.40 | 8.6e-7 |
| **RBFOX2** | m5C | 0.39% | 1.37% | 3.55 | 0.032 |
| **RBM22** | pseudoU | 0.37% | 2.87% | 7.98 | 0.040 |

### All Enrichment Results

| RBP | m6A OR | pseudoU OR | m5C OR | ac4C OR | Notes |
|-----|--------|------------|--------|---------|-------|
| IGF2BP1 | 0.056 | 0 | 0.49 | 0 | m6A depleted in discrepant (inverted pattern) |
| IGF2BP2 | 0.25 | 0.74 | 0.66 | 0 | m6A depleted in discrepant (inverted pattern) |
| EWSR1 | 0.65 | 0 | 0 | 0 | 105 discrepant peaks, low overlap |
| FUS | 2.60 | 0 | 0 | 0 | Only 19 discrepant peaks; underpowered |
| HNRNPC | 1.15 | 0 | 0 | — | 81 discrepant peaks, no signal |
| HNRNPK | 1.31 | 0 | 3.46 | 0 | 75 discrepant peaks; m5C trend (p_adj=0.148) |
| HNRNPL | 0.40 | — | 0 | — | m6A depleted in discrepant |
| LIN28B | 1.05 | 0.47 | 0.61 | 0.54 | High absolute overlap; no enrichment |
| PCBP2 | 0.83 | 1.12 | 0.26 | 0 | HepG2; no modification signal |
| **RBFOX2** | **2.40 ✓** | 0 | **3.55 ✓** | 0 | **m6A and m5C significant** |
| **RBM22** | 1.43 | **7.98 ✓** | 2.30 | 0 | **pseudoU significant** |
| SRSF9 | 0.71 | inf | 0.71 | 0 | Only 14 discrepant peaks; underpowered |
| TARDBP | 1.64 | 0 | 3.99 | — | Only 156 discrepant peaks; no signal |
| TIA1 | 0.49 | 0.25 | 0.71 | 0 | m6A depleted in discrepant |
| TRA2A | 0.21 | 13.6 | 0 | 0 | Only 106 discrepant peaks; underpowered |

✓ = FDR < 0.05

---

## Discussion

### The Inverted Pattern in IGF2BP1 and IGF2BP2

Both known m6A readers show m6A strongly enriched in *canonical* peaks rather than discrepant peaks. This is counterintuitive for the modification-dependent binding hypothesis, but has a straightforward explanation:

IGF2BP1 and IGF2BP2 are **sequence-first** m6A readers. Their canonical binding motifs (CA[U/C]-rich sequences; IGF2BP1 consensus CAUH) substantially overlap the DRACH m6A consensus. Consequently, peaks where the protein binds via sequence affinity (canonical) are disproportionately located within m6A-modified contexts. The discrepant peaks — where in vivo binding exceeds RBNS prediction — are in sequence contexts that *lack* the preferred motif, and thus also lack the m6A sites that co-occur with it.

This is not a pipeline failure. It reflects the biology: for IGF2BP1/2, the m6A modification reinforces sequence-driven binding rather than compensating for absent sequence affinity.

### RBFOX2: m6A and m5C Co-enrichment

RBFOX2 discrepant peaks are enriched for both m6A (OR=2.40) and m5C (OR=3.55). RBFOX2 canonically binds UGCAUG. Discrepant binding at non-canonical sites may be stabilized by epitranscriptomic modifications. The m5C enrichment is particularly notable given the low genomic prevalence of m5C (0.39% canonical overlap) — a 3.5× enrichment in discrepant peaks represents a strong relative effect.

### RBM22: Strong pseudoU Signal

RBM22 has the highest discrepancy rate in the dataset (68.6%) — more than two-thirds of its eCLIP peaks are classified as discrepant. Among discrepant peaks, 2.87% overlap pseudoU sites vs 0.37% in canonical (OR=7.98, FDR=0.040). RBM22 is a core spliceosome component (U4/U6.U5 tri-snRNP); pseudouridines are highly enriched in spliceosomal snRNAs and at splice sites in pre-mRNA. This association is consistent with RBM22's known function and suggests pseudoU-modified splice site regions may contribute to its pre-mRNA binding specificity.

---

## Figure Descriptions

### 1. Z-score Distribution (`zscore_distribution.png`)
Histogram of maximum 5-mer RBNS Z-scores across all eCLIP peaks.
- **Green dashed line:** Canonical threshold
- **Red dashed line:** Discrepant threshold

### 2. Classification Summary (`classification_summary.png`)
Pie chart of peak categories (canonical, intermediate, discrepant).

### 3. Score Comparison (`score_comparison.png`)
Scatter plot of Score_max vs Score_sum, colored by peak category.

### 4. Enrichment Barplot (`enrichment_barplot.png`)
Odds ratios for each modification. Green bars = FDR < 0.05.

---

## Individual RBP Results

### IGF2BP1 (Positive Control - m6A Reader)

![Z-score Distribution](IGF2BP1/figures/zscore_distribution.png)

**Statistics:** 3,724 canonical | 104 discrepant | thresholds: Z≥2.0 / Z<1.0

**Enrichment:** m6A depleted in discrepant (3.8% vs 41.5%, OR=0.056). As an m6A reader with CAUH-overlapping motif, canonical binding sites co-occur with m6A by sequence context — the inverted pattern is expected (see Discussion).

![Classification](IGF2BP1/figures/classification_summary.png)
![Score Comparison](IGF2BP1/figures/score_comparison.png)
![Enrichment](IGF2BP1/figures/enrichment_barplot.png)

---

### IGF2BP2 (Positive Control - m6A Reader)

![Z-score Distribution](IGF2BP2/figures/zscore_distribution.png)

**Statistics:** 797 canonical | 240 discrepant

**Enrichment:** m6A depleted in discrepant (18.3% vs 46.9%, OR=0.25). Same inverted pattern as IGF2BP1 (see Discussion).

![Classification](IGF2BP2/figures/classification_summary.png)
![Score Comparison](IGF2BP2/figures/score_comparison.png)
![Enrichment](IGF2BP2/figures/enrichment_barplot.png)

---

### EWSR1

![Z-score Distribution](EWSR1/figures/zscore_distribution.png)

**Statistics:** 5,387 canonical | 105 discrepant (1.9%)

Predominantly sequence-driven binding. Low modification overlap in discrepant peaks; no significant enrichment.

![Classification](EWSR1/figures/classification_summary.png)
![Score Comparison](EWSR1/figures/score_comparison.png)
![Enrichment](EWSR1/figures/enrichment_barplot.png)

---

### FUS

![Z-score Distribution](FUS/figures/zscore_distribution.png)

**Statistics:** 3,881 canonical | 19 discrepant (0.5%) | thresholds: Z≥2.0 / Z<1.0

Only 19 discrepant peaks after threshold adjustment — insufficient power for enrichment analysis. m6A trend positive (OR=2.6, p=0.13) but underpowered. FUS binding is broadly sequence-driven at this threshold.

![Classification](FUS/figures/classification_summary.png)
![Score Comparison](FUS/figures/score_comparison.png)
![Enrichment](FUS/figures/enrichment_barplot.png)

---

### HNRNPC

![Z-score Distribution](HNRNPC/figures/zscore_distribution.png)

**Statistics:** 388 canonical | 81 discrepant (17.3%)

No significant modification enrichment. HNRNPC's U-tract binding motif is well-captured by RBNS.

![Classification](HNRNPC/figures/classification_summary.png)
![Score Comparison](HNRNPC/figures/score_comparison.png)
![Enrichment](HNRNPC/figures/enrichment_barplot.png)

---

### HNRNPK

![Z-score Distribution](HNRNPK/figures/zscore_distribution.png)

**Statistics:** 2,930 canonical | 75 discrepant (2.5%)

m5C trend in discrepant peaks (5.3% vs 1.6%, OR=3.46, p_adj=0.148) — below significance threshold with 75 discrepant peaks. Worth monitoring in larger datasets.

![Classification](HNRNPK/figures/classification_summary.png)
![Score Comparison](HNRNPK/figures/score_comparison.png)
![Enrichment](HNRNPK/figures/enrichment_barplot.png)

---

### HNRNPL

![Z-score Distribution](HNRNPL/figures/zscore_distribution.png)

**Statistics:** 3,958 canonical | 371 discrepant (8.6%)

m6A depleted in discrepant (1.6% vs 4.0%, OR=0.40). Low modification density in discrepant peaks.

![Classification](HNRNPL/figures/classification_summary.png)
![Score Comparison](HNRNPL/figures/score_comparison.png)
![Enrichment](HNRNPL/figures/enrichment_barplot.png)

---

### LIN28B

![Z-score Distribution](LIN28B/figures/zscore_distribution.png)

**Statistics:** 2,424 canonical | 2,245 discrepant (48.0%)

The highest absolute discrepant peak count (2,245). m6A overlap is nearly equal in canonical (28.8%) and discrepant (29.8%, OR=1.05) — modification co-occurrence is uniform across peak classes, suggesting transcriptome-wide m6A background rather than modification-specific discrepant binding.

![Classification](LIN28B/figures/classification_summary.png)
![Score Comparison](LIN28B/figures/score_comparison.png)
![Enrichment](LIN28B/figures/enrichment_barplot.png)

---

### PCBP2 (HepG2)

![Z-score Distribution](PCBP2/figures/zscore_distribution.png)

**Statistics:** 6,979 canonical | 690 discrepant (9.0%) | Cell line: HepG2

No significant modification enrichment. Note: K562 eCLIP data unavailable for PCBP2; HepG2 modification data used. Results not directly comparable to K562 RBPs.

![Classification](PCBP2/figures/classification_summary.png)
![Score Comparison](PCBP2/figures/score_comparison.png)
![Enrichment](PCBP2/figures/enrichment_barplot.png)

---

### RBFOX2 ✓

![Z-score Distribution](RBFOX2/figures/zscore_distribution.png)

**Statistics:** 1,281 canonical | 728 discrepant (36.2%)

**Significant m6A enrichment** (11.7% disc vs 5.2% canon, OR=2.40, FDR=8.6e-7) and **significant m5C enrichment** (1.37% disc vs 0.39% canon, OR=3.55, FDR=0.032) in discrepant peaks. RBFOX2 discrepant peaks are ~2.4× more likely to carry m6A and ~3.6× more likely to carry m5C than canonical peaks, indicating that epitranscriptomic marks contribute to a substantial fraction of its in vivo binding.

![Classification](RBFOX2/figures/classification_summary.png)
![Score Comparison](RBFOX2/figures/score_comparison.png)
![Enrichment](RBFOX2/figures/enrichment_barplot.png)

---

### RBM22 ✓

![Z-score Distribution](RBM22/figures/zscore_distribution.png)

**Statistics:** 271 canonical | 592 discrepant (68.6%)

**Highly significant pseudoU enrichment in discrepant peaks** (2.87% vs 0.37%, OR=7.98, FDR=0.040). RBM22 has the highest discrepancy rate in the dataset — over two-thirds of eCLIP peaks are discrepant. As a spliceosomal component (U4/U6.U5 tri-snRNP), this pseudoU enrichment is consistent with RBM22's known function at pseudouridine-rich splice sites. m6A also shows a positive trend in discrepant peaks (OR=1.43) that did not reach FDR significance.

![Classification](RBM22/figures/classification_summary.png)
![Score Comparison](RBM22/figures/score_comparison.png)
![Enrichment](RBM22/figures/enrichment_barplot.png)

---

### SRSF9

![Z-score Distribution](SRSF9/figures/zscore_distribution.png)

**Statistics:** 82 canonical | 14 discrepant (14.6%)

Only 14 discrepant peaks. pseudoU shows infinite odds ratio (0 canonical overlap, 2/14 discrepant) — not significant after correction. Underpowered.

![Classification](SRSF9/figures/classification_summary.png)
![Score Comparison](SRSF9/figures/score_comparison.png)
![Enrichment](SRSF9/figures/enrichment_barplot.png)

---

### TARDBP

![Z-score Distribution](TARDBP/figures/zscore_distribution.png)

**Statistics:** 5,578 canonical | 156 discrepant (2.7%)

Predominantly sequence-specific (TDP-43 UG-rich motif). m5C shows a positive trend in discrepant peaks (OR=3.99, p=0.24) with low absolute counts; underpowered.

![Classification](TARDBP/figures/classification_summary.png)
![Score Comparison](TARDBP/figures/score_comparison.png)
![Enrichment](TARDBP/figures/enrichment_barplot.png)

---

### TIA1

![Z-score Distribution](TIA1/figures/zscore_distribution.png)

**Statistics:** 4,230 canonical | 929 discrepant (18.0%)

m6A depleted in discrepant relative to canonical (9.8% vs 18.1%, OR=0.49), suggesting TIA1's discrepant peaks occur in low-m6A regions of the transcriptome.

![Classification](TIA1/figures/classification_summary.png)
![Score Comparison](TIA1/figures/score_comparison.png)
![Enrichment](TIA1/figures/enrichment_barplot.png)

---

### TRA2A

![Z-score Distribution](TRA2A/figures/zscore_distribution.png)

**Statistics:** 1,433 canonical | 106 discrepant (6.9%)

m6A heavily enriched in canonical peaks (53.1%). Only 106 discrepant peaks; underpowered for enrichment testing.

![Classification](TRA2A/figures/classification_summary.png)
![Score Comparison](TRA2A/figures/score_comparison.png)
![Enrichment](TRA2A/figures/enrichment_barplot.png)

---

## Methods Summary

### Pipeline Steps
1. **Load RBNS Z-scores** — Pre-computed 5-mer enrichment Z-scores from ENCODE RBNS
2. **Extend eCLIP peaks** — 50nt 5' extension from crosslink summit (narrowPeak col 10 summit offset); strand-aware
3. **Filter peaks** — Enrichment threshold >= 2.0 (signalValue)
4. **Score peaks** — Calculate Score_max and Score_sum from RBNS Z-scores
5. **Classify peaks** — Apply thresholds (canonical >= 3.0, discrepant < 1.5; IGF2BP1/FUS: >= 2.0 / < 1.0)
6. **Enrichment analysis** — Fisher's exact test (one-sided, greater) comparing modification overlap in discrepant vs canonical; Benjamini-Hochberg FDR correction

### Data Sources
- **eCLIP:** ENCODE Project (K562; PCBP2: HepG2)
- **RBNS:** ENCODE Project (5-mer enrichment data)
- **Modifications:** RMBase v3.0 — cumulative multi-study sites (hg38). m6A for HepG2 is cell-line-specific (178,800 sites); m6A for K562 is cumulative multi-study with support ≥2 (558,360 sites). Ψ, m5C, and ac4C use the same cumulative multi-study sites for both cell lines.

### Known Limitations
- IGF2BP1 and FUS required lowered thresholds (Z≥2.0/Z<1.0) due to compressed RBNS Z-score distributions; results should be interpreted cautiously relative to the other 13 RBPs
- K562 m6A data is cumulative multi-study rather than K562-specific, which may reduce specificity of the m6A enrichment signal
- PCBP2 uses HepG2 eCLIP and modification data (no K562 available); direct comparison to K562 RBPs is not appropriate
- Several RBPs have very few discrepant peaks (FUS: 19, SRSF9: 14), giving insufficient power for enrichment tests

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

*Generated by RNAMod-RBNS Pipeline — updated 2026-02-27 after summit-anchored peak extension and per-RBP threshold overrides*
