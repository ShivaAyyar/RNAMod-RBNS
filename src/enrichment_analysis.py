#!/usr/bin/env python3
"""
Epitranscriptomic Enrichment Analysis Module
Implements Module B from the research document Section 3.2

This module handles:
1. Intersection of peaks with modification sites (m6A, Ψ, m5C, ac4C)
2. Fisher's exact test for enrichment analysis
3. FDR correction for multiple testing
4. Results aggregation across RBPs and modifications

References:
- Fisher's exact test for contingency tables
- Benjamini-Hochberg FDR correction
- REPIC database: https://repicmod.uchicago.edu/
- RMBase v2.0: http://rna.sysu.edu.cn/rmbase/
"""

import os
import subprocess
import tempfile

import pandas as pd
import pybedtools
import logging
from pathlib import Path
from typing import List, Tuple, Optional, Dict
import json
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests

# Configure logging
logger = logging.getLogger(__name__)


def count_modification_overlaps(
    peaks_bed: str,
    mod_bed: str,
    min_overlap: int = 1,
    window: int = 0,
    chrom_sizes: Optional[str] = None
) -> Tuple[int, int]:
    """
    Count peaks overlapping modification sites.

    Uses bedtools intersect with -u flag for unique peaks only.

    Parameters
    ----------
    peaks_bed : str
        Path to peaks BED file
    mod_bed : str
        Path to modification sites BED file
    min_overlap : int
        Minimum overlap in bp (default: 1)
    window : int
        Extend peaks by this many bp on each side before intersection.
        0 = exact intersection (default). Useful for single-nucleotide
        modification sites (e.g. RMBase) where exact overlap may miss
        sites at peak edges.
    chrom_sizes : str, optional
        Path to chromosome sizes file (required when window > 0)

    Returns
    -------
    tuple
        (overlapping_count, total_count)
    """
    peaks = pybedtools.BedTool(peaks_bed)
    mods = pybedtools.BedTool(mod_bed)

    total = len(peaks)

    if window > 0 and chrom_sizes:
        peaks = peaks.slop(b=window, g=chrom_sizes)

    overlapping = len(peaks.intersect(mods, u=True))

    logger.debug(f"Overlap: {overlapping}/{total} peaks with {mod_bed}")
    return overlapping, total


def fisher_enrichment_test(
    disc_overlap: int,
    disc_total: int,
    canon_overlap: int,
    canon_total: int,
    alternative: str = 'greater'
) -> Tuple[float, float]:
    """
    Fisher's exact test for modification enrichment.

    Tests whether discrepant peaks are enriched for modifications
    compared to canonical peaks.

    Contingency table:
                        Mod+                    Mod-
    Discrepant         disc_overlap    disc_total - disc_overlap
    Canonical          canon_overlap   canon_total - canon_overlap

    Parameters
    ----------
    disc_overlap : int
        Number of discrepant peaks overlapping modification
    disc_total : int
        Total number of discrepant peaks
    canon_overlap : int
        Number of canonical peaks overlapping modification
    canon_total : int
        Total number of canonical peaks
    alternative : str
        Alternative hypothesis: 'two-sided', 'greater', or 'less'
        Default 'greater' tests if discrepant has MORE modifications

    Returns
    -------
    tuple
        (odds_ratio, p_value)

    Notes
    -----
    - Odds ratio > 1 indicates enrichment in discrepant peaks
    - Uses scipy.stats.fisher_exact
    - Handles edge cases (zero counts) gracefully
    """
    # Handle edge cases
    disc_no_mod = disc_total - disc_overlap
    canon_no_mod = canon_total - canon_overlap

    # Ensure non-negative values
    disc_no_mod = max(0, disc_no_mod)
    canon_no_mod = max(0, canon_no_mod)

    table = [
        [disc_overlap, disc_no_mod],
        [canon_overlap, canon_no_mod]
    ]

    try:
        odds_ratio, pvalue = fisher_exact(table, alternative=alternative)
    except ValueError as e:
        logger.warning(f"Fisher's test failed: {e}. Returning NaN.")
        return float('nan'), float('nan')

    return odds_ratio, pvalue


def run_enrichment_analysis(
    canonical_bed: str,
    discrepant_bed: str,
    mod_beds: List[str],
    mod_names: List[str],
    apply_fdr: bool = True,
    fdr_method: str = 'fdr_bh',
    window: int = 0,
    chrom_sizes: Optional[str] = None,
    run_gat: bool = False,
    workspace_bed: Optional[str] = None,
    gat_n_samples: int = 10000,
    output_dir: Optional[str] = None
) -> pd.DataFrame:
    """
    Run enrichment analysis for all modification types.

    Parameters
    ----------
    canonical_bed : str
        Path to canonical peaks BED file
    discrepant_bed : str
        Path to discrepant peaks BED file
    mod_beds : list
        List of paths to modification BED files
    mod_names : list
        List of modification names (same order as mod_beds)
    apply_fdr : bool
        Whether to apply FDR correction (default: True)
    fdr_method : str
        Method for multipletests: 'fdr_bh' (Benjamini-Hochberg),
        'bonferroni', 'holm', etc. (default: 'fdr_bh')

    Returns
    -------
    pd.DataFrame
        Results with columns:
        - modification: name of modification
        - canonical_overlap: count in canonical peaks
        - canonical_total: total canonical peaks
        - discrepant_overlap: count in discrepant peaks
        - discrepant_total: total discrepant peaks
        - odds_ratio: enrichment odds ratio
        - pvalue: raw p-value
        - pvalue_adj: FDR-adjusted p-value (if apply_fdr=True)
        - significant: boolean (if apply_fdr=True)
    """
    logger.info(f"Running enrichment analysis for {len(mod_names)} modifications")

    # Validate inputs
    if len(mod_beds) != len(mod_names):
        raise ValueError("mod_beds and mod_names must have same length")

    # Check files exist
    for bed_file in [canonical_bed, discrepant_bed] + mod_beds:
        if not Path(bed_file).exists():
            raise FileNotFoundError(f"File not found: {bed_file}")

    # Get total counts
    canon_total = len(pybedtools.BedTool(canonical_bed))
    disc_total = len(pybedtools.BedTool(discrepant_bed))

    logger.info(f"Canonical peaks: {canon_total}, Discrepant peaks: {disc_total}")

    if canon_total == 0 or disc_total == 0:
        logger.warning("One or both peak sets are empty!")

    results = []
    pvalues = []

    for mod_bed, mod_name in zip(mod_beds, mod_names):
        logger.info(f"Processing {mod_name}...")

        canon_overlap, _ = count_modification_overlaps(
            canonical_bed, mod_bed, window=window, chrom_sizes=chrom_sizes)
        disc_overlap, _ = count_modification_overlaps(
            discrepant_bed, mod_bed, window=window, chrom_sizes=chrom_sizes)

        odds_ratio, pvalue = fisher_enrichment_test(
            disc_overlap, disc_total,
            canon_overlap, canon_total
        )

        # Calculate percentages
        canon_pct = (canon_overlap / canon_total * 100) if canon_total > 0 else 0
        disc_pct = (disc_overlap / disc_total * 100) if disc_total > 0 else 0

        results.append({
            'modification': mod_name,
            'canonical_overlap': canon_overlap,
            'canonical_total': canon_total,
            'canonical_pct': canon_pct,
            'discrepant_overlap': disc_overlap,
            'discrepant_total': disc_total,
            'discrepant_pct': disc_pct,
            'odds_ratio': odds_ratio,
            'pvalue': pvalue
        })
        pvalues.append(pvalue)

        logger.info(
            f"  {mod_name}: OR={odds_ratio:.2f}, p={pvalue:.2e}, "
            f"disc={disc_overlap}/{disc_total} ({disc_pct:.1f}%), "
            f"canon={canon_overlap}/{canon_total} ({canon_pct:.1f}%)"
        )

    # Apply FDR correction
    if apply_fdr and len(pvalues) > 0:
        # Handle NaN p-values
        valid_pvalues = [p if not pd.isna(p) else 1.0 for p in pvalues]

        reject, pvals_adj, _, _ = multipletests(
            valid_pvalues,
            method=fdr_method,
            alpha=0.05
        )

        for i, result in enumerate(results):
            result['pvalue_adj'] = pvals_adj[i]
            result['significant'] = reject[i]

    fisher_df = pd.DataFrame(results)

    # Optionally run GAT on discrepant peaks as a background-controlled
    # complement to the Fisher test. Results are saved to gat_results.csv
    # but do not affect the returned Fisher DataFrame.
    if run_gat and workspace_bed:
        try:
            gat_df = run_gat_analysis(
                discrepant_bed, mod_beds, mod_names,
                workspace_bed, chrom_sizes,
                n_samples=gat_n_samples,
                output_dir=output_dir
            )
            if output_dir:
                gat_df.to_csv(
                    Path(output_dir) / 'gat_results.csv', index=False)
                logger.info("GAT results saved to gat_results.csv")
        except Exception as e:
            logger.error(f"GAT analysis failed (Fisher results unaffected): {e}")

    return fisher_df


def run_gat_analysis(
    segments_bed: str,
    mod_beds: List[str],
    mod_names: List[str],
    workspace_bed: str,
    chrom_sizes: str,
    n_samples: int = 10000,
    output_dir: Optional[str] = None
) -> pd.DataFrame:
    """
    Run GAT (Genomic Association Tester) for modification enrichment.

    Tests whether peaks overlap modification sites more than expected
    by chance, controlling for genomic region biases by randomly
    relocating peaks within the workspace (filtered eCLIP peaks).

    GAT is invoked as a subprocess (gat-run.py must be in PATH).

    Parameters
    ----------
    segments_bed : str
        Path to peaks BED file to test (discrepant peaks)
    mod_beds : list
        Paths to modification BED files (used as annotations)
    mod_names : list
        Names of modifications (same order as mod_beds)
    workspace_bed : str
        BED file defining regions within which peaks are shuffled.
        Should be peaks_filtered.bed (all eCLIP peaks for the RBP)
        to constrain shuffling to the transcriptomic regions where
        this protein actually binds.
    chrom_sizes : str
        Path to chromosome sizes file
    n_samples : int
        Number of random shuffles for null distribution (default: 10000)
    output_dir : str, optional
        Directory for GAT output files (tsv and log)

    Returns
    -------
    pd.DataFrame
        GAT results with columns: modification, observed, expected,
        fold, pvalue, qvalue
    """
    work_dir = output_dir or tempfile.gettempdir()
    out_tsv = os.path.join(work_dir, 'gat_results.tsv')
    gat_log = os.path.join(work_dir, 'gat.log')

    # GAT's BED parser is strict: it expects exactly 3 or 4 columns and
    # chokes on float scores, long names, or a trailing strand '-' character
    # (which it misreads as an incomplete negative number). Write clean
    # 3-column BEDs for segments and workspace, and 4-column for annotations.

    def _write_3col_bed(src_bed: str, fh) -> None:
        """Write chrom/start/end only — strip all extra columns."""
        with open(src_bed) as f:
            for line in f:
                line = line.rstrip()
                if not line or line.startswith('#'):
                    continue
                fields = line.split('\t')
                fh.write(f"{fields[0]}\t{fields[1]}\t{fields[2]}\n")

    # Segments: 3-column clean BED
    seg_fh = tempfile.NamedTemporaryFile(
        mode='w', suffix='.bed', delete=False)
    seg_path = seg_fh.name
    _write_3col_bed(segments_bed, seg_fh)
    seg_fh.close()

    # Workspace: 3-column clean BED
    ws_fh = tempfile.NamedTemporaryFile(
        mode='w', suffix='.bed', delete=False)
    ws_path = ws_fh.name
    _write_3col_bed(workspace_bed, ws_fh)
    ws_fh.close()

    # Annotations: 4-column BED (chrom, start, end, mod_name).
    # Concatenate all mod beds, replacing col 4 with the mod name.
    ann_fh = tempfile.NamedTemporaryFile(
        mode='w', suffix='.bed', delete=False)
    try:
        ann_path = ann_fh.name
        for mod_bed, mod_name in zip(mod_beds, mod_names):
            with open(mod_bed) as f:
                for line in f:
                    line = line.rstrip()
                    if not line or line.startswith('#'):
                        continue
                    fields = line.split('\t')
                    ann_fh.write(
                        f"{fields[0]}\t{fields[1]}\t{fields[2]}\t{mod_name}\n")
        ann_fh.close()

        # GAT writes results to stdout (single counter = nucleotide-overlap).
        # Capture stdout directly; stderr goes to the log file.
        cmd = [
            'gat-run.py',
            '--segments', seg_path,
            '--annotations', ann_path,
            '--workspace', ws_path,
            '--num-samples', str(n_samples),
        ]

        logger.info(f"Running GAT with {n_samples} permutations...")
        logger.debug(f"GAT command: {' '.join(cmd)}")

        with open(gat_log, 'w') as log_fh:
            result = subprocess.run(
                cmd, stdout=subprocess.PIPE, stderr=log_fh,
                check=True, text=True)

        # GAT stdout mixes progress log lines with the TSV results table.
        # The results table starts at the line beginning with "track\t".
        # Extract only that portion before writing/parsing.
        lines = result.stdout.splitlines()
        tsv_start = next(
            (i for i, l in enumerate(lines) if l.startswith('track\t')),
            None
        )
        if tsv_start is None:
            raise ValueError(
                "Could not find TSV header in GAT output. "
                f"See {gat_log} for details.")
        tsv_content = '\n'.join(lines[tsv_start:])

        with open(out_tsv, 'w') as f:
            f.write(tsv_content)

        # Parse GAT output table (24 columns confirmed from source)
        df = pd.read_csv(out_tsv, sep='\t')

        # GAT output columns include: track, annotation, observed,
        # expected, CI95low, CI95high, stddev, fold, pvalue, qvalue, ...
        # 'annotation' is the mod name we embedded in col 4.
        keep = ['annotation', 'observed', 'expected', 'fold',
                'pvalue', 'qvalue']
        df = df[[c for c in keep if c in df.columns]].copy()
        df = df.rename(columns={'annotation': 'modification'})

        logger.info(f"GAT complete. Results for {len(df)} modifications.")
        return df

    except subprocess.CalledProcessError as e:
        logger.error(f"GAT failed (exit {e.returncode}). See {gat_log} for details.")
        raise
    finally:
        for tmp in (ann_path, seg_path, ws_path):
            try:
                os.unlink(tmp)
            except OSError:
                pass


def run_ranked_enrichment_analysis(
    peaks_bed: str,
    scored_df: pd.DataFrame,
    mod_beds: List[str],
    mod_names: List[str],
    chrom_sizes: Optional[str] = None,
    window: int = 0,
    n_bins: int = 5,
    apply_fdr: bool = True,
    fdr_method: str = 'fdr_bh'
) -> pd.DataFrame:
    """
    Ranked continuous enrichment analysis using Mann-Whitney U test.

    Tests whether peaks overlapping modification sites have systematically
    different RBNS Z-scores than peaks that don't, using all filtered peaks
    (not just canonical/discrepant subsets). This eliminates threshold
    sensitivity and detects enrichment in either direction.

    Parameters
    ----------
    peaks_bed : str
        Path to ALL filtered peaks BED file (peaks_filtered.bed)
    scored_df : pd.DataFrame
        DataFrame with 'peak_id' and 'score_max' columns (all filtered peaks)
    mod_beds : list
        List of paths to modification BED files
    mod_names : list
        List of modification names (same order as mod_beds)
    chrom_sizes : str, optional
        Path to chromosome sizes file (required if window > 0)
    window : int
        Extend peaks by this many bp on each side before intersecting.
        0 = exact intersection (default)
    n_bins : int
        Number of Z-score quantile bins for binned frequency analysis (default: 5)
    apply_fdr : bool
        Whether to apply BH FDR correction (default: True)
    fdr_method : str
        Method for multipletests (default: 'fdr_bh')

    Returns
    -------
    pd.DataFrame
        Results with columns: modification, n_total_peaks, n_mod_positive,
        n_mod_negative, pct_mod_positive, median_z_mod_pos, median_z_mod_neg,
        median_z_diff, direction, mann_whitney_u, mann_whitney_p,
        mann_whitney_p_adj, significant, spearman_rho, spearman_p,
        z_bin_frequencies
    """
    logger.info(
        f"Running ranked enrichment analysis on {len(scored_df)} peaks "
        f"for {len(mod_names)} modifications"
    )

    peaks = pybedtools.BedTool(peaks_bed)

    if window > 0 and chrom_sizes:
        peaks_for_intersect = peaks.slop(b=window, g=chrom_sizes)
    else:
        peaks_for_intersect = peaks

    # Build a set of peak IDs present in scored_df for fast lookup
    # peak_id format matches pa.scored_df_to_bed: "{chrom}:{start}-{end}"
    scored_index = scored_df.set_index('peak_id')['score_max']

    results = []
    pvalues = []

    for mod_bed, mod_name in zip(mod_beds, mod_names):
        logger.info(f"  Ranked analysis: {mod_name}...")

        # Find which peaks overlap this modification
        mods = pybedtools.BedTool(mod_bed)
        overlapping = peaks_for_intersect.intersect(mods, u=True)

        # Build set of overlapping peak IDs (chrom:start-end)
        # Use original (non-slopped) coordinates from scored_df
        overlapping_ids = set()
        for feat in overlapping:
            if window > 0:
                # The slopped peak coords differ from original; re-intersect
                # using -wa to get original peak coords instead
                pass
            else:
                overlapping_ids.add(f"{feat.chrom}:{feat.start}-{feat.end}")

        # For windowed mode, use -wa (write original A record) to recover IDs
        if window > 0 and chrom_sizes:
            overlapping_wa = peaks.intersect(mods, u=True, wa=True)
            overlapping_ids = set()
            for feat in overlapping_wa:
                overlapping_ids.add(f"{feat.chrom}:{feat.start}-{feat.end}")

        # Annotate scored_df
        has_mod = scored_df['peak_id'].isin(overlapping_ids)
        z_mod_pos = scored_df.loc[has_mod, 'score_max'].values
        z_mod_neg = scored_df.loc[~has_mod, 'score_max'].values

        n_total = len(scored_df)
        n_pos = int(has_mod.sum())
        n_neg = n_total - n_pos
        pct_pos = 100.0 * n_pos / n_total if n_total > 0 else 0.0

        # Mann-Whitney U test (two-tailed)
        if n_pos >= 3 and n_neg >= 3:
            u_stat, mw_p = mannwhitneyu(
                z_mod_pos, z_mod_neg, alternative='two-sided')
        else:
            logger.warning(
                f"  {mod_name}: insufficient peaks for Mann-Whitney "
                f"(mod+={n_pos}, mod-={n_neg}). Skipping."
            )
            u_stat, mw_p = float('nan'), float('nan')

        # Median effect sizes
        med_pos = float(pd.Series(z_mod_pos).median()) if n_pos > 0 else float('nan')
        med_neg = float(pd.Series(z_mod_neg).median()) if n_neg > 0 else float('nan')
        med_diff = (med_pos - med_neg) if not (
            pd.isna(med_pos) or pd.isna(med_neg)) else float('nan')

        if pd.isna(med_diff):
            direction = 'insufficient_data'
        elif abs(med_diff) < 1e-9:
            direction = 'no_difference'
        elif med_diff < 0:
            direction = 'enriched_at_low_Z'
        else:
            direction = 'enriched_at_high_Z'

        # Spearman correlation: score_max vs binary has_mod (0/1)
        if n_pos >= 3 and n_neg >= 3:
            rho, sp_p = spearmanr(
                scored_df['score_max'].values,
                has_mod.astype(int).values
            )
        else:
            rho, sp_p = float('nan'), float('nan')

        # Binned analysis: modification frequency across Z-score quantiles
        try:
            bin_labels = [f'Q{i+1}' for i in range(n_bins)]
            scored_df_copy = scored_df.copy()
            scored_df_copy['_has_mod'] = has_mod
            scored_df_copy['_zbin'] = pd.qcut(
                scored_df_copy['score_max'],
                q=n_bins,
                labels=bin_labels,
                duplicates='drop'
            )
            bin_freqs = (
                scored_df_copy.groupby('_zbin', observed=True)['_has_mod']
                .mean()
                .multiply(100)
                .round(2)
                .to_dict()
            )
        except Exception:
            bin_freqs = {}

        results.append({
            'modification': mod_name,
            'n_total_peaks': n_total,
            'n_mod_positive': n_pos,
            'n_mod_negative': n_neg,
            'pct_mod_positive': round(pct_pos, 3),
            'median_z_mod_pos': round(med_pos, 4) if not pd.isna(med_pos) else None,
            'median_z_mod_neg': round(med_neg, 4) if not pd.isna(med_neg) else None,
            'median_z_diff': round(med_diff, 4) if not pd.isna(med_diff) else None,
            'direction': direction,
            'mann_whitney_u': u_stat,
            'mann_whitney_p': mw_p,
            'spearman_rho': round(rho, 4) if not pd.isna(rho) else None,
            'spearman_p': sp_p,
            'z_bin_frequencies': json.dumps(bin_freqs)
        })
        pvalues.append(mw_p)

        logger.info(
            f"    {mod_name}: Δmedian={med_diff:+.3f} ({direction}), "
            f"MW p={mw_p:.2e}, ρ={rho:.3f}"
            if not pd.isna(mw_p) else
            f"    {mod_name}: insufficient data"
        )

    # FDR correction across modifications
    df = pd.DataFrame(results)
    if apply_fdr and len(pvalues) > 0:
        valid_p = [p if not pd.isna(p) else 1.0 for p in pvalues]
        reject, p_adj, _, _ = multipletests(valid_p, method=fdr_method, alpha=0.05)
        df['mann_whitney_p_adj'] = p_adj
        df['significant'] = reject
    else:
        df['mann_whitney_p_adj'] = df['mann_whitney_p']
        df['significant'] = False

    return df


def aggregate_rbp_results(
    results_dir: str,
    rbp_list: List[str],
    output_path: Optional[str] = None
) -> pd.DataFrame:
    """
    Aggregate enrichment results across multiple RBPs.

    Parameters
    ----------
    results_dir : str
        Directory containing per-RBP result directories
    rbp_list : list
        List of RBP names to aggregate
    output_path : str, optional
        Path to save aggregated results

    Returns
    -------
    pd.DataFrame
        Aggregated results with RBP column
    """
    results_dir = Path(results_dir)
    all_results = []

    for rbp in rbp_list:
        result_file = results_dir / rbp / 'enrichment_results.csv'
        if result_file.exists():
            df = pd.read_csv(result_file)
            df['rbp'] = rbp
            all_results.append(df)
        else:
            logger.warning(f"Results not found for {rbp}: {result_file}")

    if not all_results:
        logger.error("No results found to aggregate")
        return pd.DataFrame()

    aggregated = pd.concat(all_results, ignore_index=True)

    if output_path:
        aggregated.to_csv(output_path, index=False)
        logger.info(f"Saved aggregated results to {output_path}")

    return aggregated


def identify_top_candidates(
    results_df: pd.DataFrame,
    min_odds_ratio: float = 2.0,
    max_pvalue: float = 0.05,
    use_adjusted: bool = True
) -> pd.DataFrame:
    """
    Identify top modification-dependent binding candidates.

    Parameters
    ----------
    results_df : pd.DataFrame
        Enrichment results DataFrame
    min_odds_ratio : float
        Minimum odds ratio threshold (default: 2.0)
    max_pvalue : float
        Maximum p-value threshold (default: 0.05)
    use_adjusted : bool
        Use adjusted p-values if available (default: True)

    Returns
    -------
    pd.DataFrame
        Filtered results sorted by odds ratio
    """
    pval_col = 'pvalue_adj' if use_adjusted and 'pvalue_adj' in results_df.columns else 'pvalue'

    candidates = results_df[
        (results_df['odds_ratio'] >= min_odds_ratio) &
        (results_df[pval_col] <= max_pvalue)
    ].copy()

    candidates = candidates.sort_values('odds_ratio', ascending=False)

    logger.info(f"Identified {len(candidates)} significant candidates")
    return candidates


def create_contingency_table(
    results_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Create a formatted contingency table for display.

    Parameters
    ----------
    results_df : pd.DataFrame
        Single row or filtered enrichment results

    Returns
    -------
    pd.DataFrame
        Formatted 2x2 contingency table
    """
    tables = []
    for _, row in results_df.iterrows():
        disc_mod = row['discrepant_overlap']
        disc_no_mod = row['discrepant_total'] - row['discrepant_overlap']
        canon_mod = row['canonical_overlap']
        canon_no_mod = row['canonical_total'] - row['canonical_overlap']

        table = pd.DataFrame({
            'Mod+': [disc_mod, canon_mod],
            'Mod-': [disc_no_mod, canon_no_mod]
        }, index=['Discrepant', 'Canonical'])

        tables.append({
            'modification': row['modification'],
            'table': table
        })

    return tables


if __name__ == '__main__':
    # Example usage
    import sys
    logging.basicConfig(level=logging.INFO)

    if len(sys.argv) < 5:
        print("Usage: python enrichment_analysis.py <canonical.bed> <discrepant.bed> <mod.bed> <mod_name>")
        sys.exit(1)

    canonical, discrepant, mod_bed, mod_name = sys.argv[1:5]

    results = run_enrichment_analysis(
        canonical, discrepant,
        [mod_bed], [mod_name]
    )

    print("\nEnrichment Results:")
    print(results.to_string())
