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

import pandas as pd
import pybedtools
import logging
from pathlib import Path
from typing import List, Tuple, Optional, Dict
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# Configure logging
logger = logging.getLogger(__name__)


def count_modification_overlaps(
    peaks_bed: str,
    mod_bed: str,
    min_overlap: int = 1
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

    Returns
    -------
    tuple
        (overlapping_count, total_count)
    """
    peaks = pybedtools.BedTool(peaks_bed)
    mods = pybedtools.BedTool(mod_bed)

    total = len(peaks)
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
    fdr_method: str = 'fdr_bh'
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

        canon_overlap, _ = count_modification_overlaps(canonical_bed, mod_bed)
        disc_overlap, _ = count_modification_overlaps(discrepant_bed, mod_bed)

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

    return pd.DataFrame(results)


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
