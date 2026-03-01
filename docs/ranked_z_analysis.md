# Ranked Z-score Enrichment Analysis

## Overview

Replace the current categorical Fisher's exact test with a ranked continuous approach
using the Mann-Whitney U test. This uses all peaks (not just canonical/discrepant
subsets), eliminates arbitrary thresholds, and detects enrichment in either direction.

---

## Motivation

The Fisher's exact test compares modification overlap rates between two peak subsets
defined by arbitrary Z-score thresholds. This has two limitations:

1. **Threshold sensitivity** — the canonical/discrepant boundary is arbitrary. A peak
   at Z=1.4 vs Z=1.6 is treated categorically differently.
2. **Directionality** — the one-tailed Fisher test only detects enrichment in
   discrepant peaks. It cannot detect the inverted pattern seen in IGF2BP1/2 (m6A
   enriched in *canonical* peaks) without additional tests.

The Mann-Whitney approach treats RBNS Z-score as a continuous variable and asks:
*Do peaks overlapping modification sites have systematically different Z-scores than
peaks that don't?* A significant negative median difference means modification-bearing
peaks have lower Z-scores (enrichment at discrepant/modification-dependent sites). A
positive difference means the opposite (e.g., IGF2BP1's sequence-co-occurrence pattern).

---

## Design Decisions

1. **Use ALL peaks** (`peaks_filtered.bed`) — maximizes power, no threshold sensitivity
2. **Two-tailed test** — catches enrichment in either direction; direction reported via
   median difference sign
3. **Effect size primary** — median Z-score difference is interpretable (units = Z-score);
   Spearman ρ shows linear trend strength; p-value used only for significance threshold
4. **Binned analysis** — visualizes U-shaped or non-linear enrichment patterns that
   the median-based test would miss
5. **Backward compatibility** — Fisher/GAT analyses unchanged; ranked runs in parallel
   via `--enrichment-method both` (default)

---

## Phase 1: Enrichment Analysis Function

**File:** `src/enrichment_analysis.py`

### New Function: `run_ranked_enrichment_analysis()`

**Inputs:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `peaks_bed` | `str` | Path to ALL filtered peaks (`peaks_filtered.bed`) |
| `scored_df` | `pd.DataFrame` | DataFrame with `peak_id` and `score_max` columns |
| `mod_beds` | `list[str]` | List of modification BED file paths |
| `mod_names` | `list[str]` | List of modification names (same order) |
| `chrom_sizes` | `str` | Path to chromosome sizes file |
| `window` | `int` | Window size for overlap (0 = exact, default) |

**Processing steps per modification:**

1. Use `bedtools intersect` to identify which peaks overlap modification sites
2. Annotate `scored_df` with binary column `has_{mod_name}` (True/False)
3. Split Z-scores into two groups:
   - `z_mod_positive`: Z-scores for peaks **with** modification
   - `z_mod_negative`: Z-scores for peaks **without** modification
4. Run **Mann-Whitney U test** (two-tailed):
   `scipy.stats.mannwhitneyu(z_mod_positive, z_mod_negative, alternative='two-sided')`
5. Calculate effect sizes:
   - Median Z for mod+ peaks, median Z for mod- peaks
   - Median difference (mod+ minus mod-)
   - Direction: `enriched_at_low_Z` if diff < 0, `enriched_at_high_Z` if diff > 0
6. Run **Spearman correlation**:
   `scipy.stats.spearmanr(scored_df['score_max'], has_mod_column)`
   (negative ρ = low-Z enrichment, positive ρ = high-Z enrichment)
7. **Binned analysis**: divide peaks into 5 quantiles by Z-score; calculate
   modification frequency per bin

**Output DataFrame columns:**

| Column | Description |
|--------|-------------|
| `modification` | Modification name |
| `n_total_peaks` | Total peaks analyzed |
| `n_mod_positive` | Count of peaks overlapping modification |
| `n_mod_negative` | Count of peaks NOT overlapping modification |
| `pct_mod_positive` | Percentage with modification |
| `median_z_mod_pos` | Median Z-score for mod+ peaks |
| `median_z_mod_neg` | Median Z-score for mod- peaks |
| `median_z_diff` | Difference (mod+ minus mod-) |
| `direction` | `enriched_at_low_Z`, `enriched_at_high_Z`, or `no_difference` |
| `mann_whitney_u` | U-statistic |
| `mann_whitney_p` | Raw p-value |
| `mann_whitney_p_adj` | BH-corrected p-value |
| `significant` | Boolean (adj p < 0.05) |
| `spearman_rho` | Spearman correlation coefficient |
| `spearman_p` | P-value for correlation |
| `z_bin_frequencies` | JSON dict of modification frequency per Z-score quantile |

---

## Phase 2: Visualization Functions

**File:** `src/visualization.py`

### `plot_ranked_enrichment_analysis()`

Three-panel figure per modification:

**Panel 1 — Cumulative Distribution Function (CDF)**
- X: RBNS Z-score; Y: Cumulative fraction (0–1)
- Two lines: mod+ peaks vs mod- peaks
- Medians marked with vertical dashed lines
- Legend showing median values

**Panel 2 — Violin Plot**
- X: Two categories (mod+, mod-)
- Y: RBNS Z-score
- Medians marked; annotation shows Δmedian, p-value, direction

**Panel 3 — Binned Frequency Plot**
- X: Z-score quantile bins (Q1–Q5)
- Y: Modification frequency (%)
- Dashed horizontal line at overall modification frequency (baseline)
- Visualizes non-linear enrichment patterns (e.g., U-shaped)

Saved as: `{output_dir}/figures/{mod_name}_ranked_enrichment.png`

### `plot_all_mods_ranked_summary()`

Two-panel summary figure across all modifications:

**Panel 1 — Forest Plot of Median Differences**
- Y: Modification names
- X: Median Z-score difference (mod+ minus mod-)
- Vertical line at 0 (no difference)
- Color by significance (green = adj p < 0.05, gray = not significant)
- Significance stars (*, **, ***)

**Panel 2 — Spearman ρ Heatmap**
- Single-column heatmap of Spearman ρ values
- Color scale: red (negative ρ = low-Z enrichment) → white (0) → blue (positive ρ)
- ρ value annotated in each cell

Saved as: `{output_dir}/figures/ranked_enrichment_summary.png`

---

## Phase 3: Main Pipeline Integration

**File:** `src/main.py`

Add after Step 4 (Classify peaks), new **Step 4b**:

```python
if args.enrichment_method in ['ranked', 'both']:
    ranked_results = ea.run_ranked_enrichment_analysis(
        peaks_bed=filtered_bed,   # ALL peaks
        scored_df=scored_df,
        mod_beds=args.mods,
        mod_names=args.mod_names,
        chrom_sizes=args.chrom_sizes,
        window=args.overlap_window
    )
    ranked_results['rbp'] = args.rbp
    ranked_results['cell_line'] = args.cell_line
    ranked_results.to_csv(outdir / 'ranked_enrichment_results.csv', index=False)
```

Add to Step 6 (Visualizations):

```python
for mod_name in args.mod_names:
    plot_ranked_enrichment_analysis(scored_df, mod_name, ...)

plot_all_mods_ranked_summary(ranked_results, ...)
```

---

## Phase 4: CLI Arguments

**File:** `src/main.py` argument parser

```
--enrichment-method   {categorical,ranked,both}   default: both
--overlap-window      int                          default: 0 (exact overlap)
```

Execution logic:

```python
if args.enrichment_method in ['categorical', 'both']:
    enrichment_results = run_enrichment_analysis(...)   # Fisher + GAT

if args.enrichment_method in ['ranked', 'both']:
    ranked_results = run_ranked_enrichment_analysis(...)
```

---

## Phase 5: Output Files

```
results/{RBP}/
├── ranked_enrichment_results.csv       NEW: Mann-Whitney results
├── enrichment_results.csv              EXISTING: Fisher's test results
├── gat_results.csv                     EXISTING: GAT permutation results
├── scored_peaks.csv                    EXISTING: All peaks with scores
└── figures/
    ├── m6A_ranked_enrichment.png       NEW: 3-panel CDF/violin/binned
    ├── pseudoU_ranked_enrichment.png   NEW
    ├── m5C_ranked_enrichment.png       NEW
    ├── ac4C_ranked_enrichment.png      NEW
    ├── ranked_enrichment_summary.png   NEW: Forest plot + rho heatmap
    ├── zscore_distribution.png         EXISTING
    ├── classification_summary.png      EXISTING
    └── enrichment_barplot.png          EXISTING
```

---

## Phase 6: Summary Statistics Updates

Add to `summary.csv`:

| New Column | Description |
|------------|-------------|
| `n_significant_ranked` | Count of significant modifications (Mann-Whitney) |
| `top_mod_ranked` | Modification with largest `|median_z_diff|` |
| `top_mod_direction` | Direction of top modification effect |

---

## Expected Outcomes

**RBFOX2** (modification-enabled binding at low-affinity sites):
```
m6A:   median_z_diff ≈ -0.7, direction = enriched_at_low_Z, p_adj < 0.001
m5C:   median_z_diff < 0,    direction = enriched_at_low_Z, p_adj < 0.05
```

**IGF2BP1** (sequence motif co-occurs with m6A context):
```
m6A:   median_z_diff ≈ +1.4, direction = enriched_at_high_Z, p_adj < 0.001
```

**Most RBPs:**
```
All modifications: |median_z_diff| < 0.5, p_adj > 0.05
```

---

## Success Criteria

- [ ] All 15 RBPs run without errors
- [ ] RBFOX2 shows significant low-Z m6A enrichment
- [ ] IGF2BP1 shows significant high-Z m6A enrichment (inverted pattern)
- [ ] Figures clearly show distributional differences
- [ ] Both ranked and categorical results saved for comparison
- [ ] No information lost from current pipeline (all existing files still generated)