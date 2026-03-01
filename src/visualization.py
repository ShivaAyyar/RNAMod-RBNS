#!/usr/bin/env python3
"""
Visualization Module for Epitranscriptomic RBP Analysis

This module provides plotting functions for:
1. Z-score distribution histograms
2. Enrichment analysis barplots
3. Odds ratio heatmaps across RBPs
4. Peak classification summary plots

Uses matplotlib and seaborn for publication-quality figures.
"""

import json

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from pathlib import Path
from typing import Optional, List, Tuple

# Configure logging
logger = logging.getLogger(__name__)

# Set default style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette('colorblind')

# Default figure parameters
FIGURE_DPI = 150
FIGURE_FORMAT = 'png'


def plot_zscore_distribution(
    scored_df: pd.DataFrame,
    output_path: str,
    canonical_threshold: float = 3.0,
    discrepant_threshold: float = 1.5,
    title: Optional[str] = None,
    figsize: Tuple[int, int] = (10, 6)
) -> None:
    """
    Plot Z-score distribution with classification thresholds.

    Creates a histogram of Score_max values with vertical lines
    indicating canonical and discrepant thresholds.

    Parameters
    ----------
    scored_df : pd.DataFrame
        DataFrame with 'score_max' column
    output_path : str
        Path to save the figure
    canonical_threshold : float
        Threshold for canonical classification (default: 3.0)
    discrepant_threshold : float
        Threshold for discrepant classification (default: 1.5)
    title : str, optional
        Plot title
    figsize : tuple
        Figure size (width, height)
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Filter out infinite values
    plot_data = scored_df[scored_df['score_max'] > float('-inf')]['score_max']

    # Create histogram
    ax.hist(plot_data, bins=50, edgecolor='black', alpha=0.7, color='steelblue')

    # Add threshold lines
    ax.axvline(
        x=canonical_threshold,
        color='green',
        linestyle='--',
        linewidth=2,
        label=f'Canonical (≥{canonical_threshold})'
    )
    ax.axvline(
        x=discrepant_threshold,
        color='red',
        linestyle='--',
        linewidth=2,
        label=f'Discrepant (<{discrepant_threshold})'
    )

    # Add shading for regions
    ymax = ax.get_ylim()[1]
    ax.axvspan(canonical_threshold, plot_data.max() + 1, alpha=0.1, color='green')
    ax.axvspan(plot_data.min() - 1, discrepant_threshold, alpha=0.1, color='red')

    # Labels
    ax.set_xlabel('Maximum 5-mer RBNS Z-score', fontsize=12)
    ax.set_ylabel('Number of Peaks', fontsize=12)

    if title:
        ax.set_title(title, fontsize=14)
    else:
        ax.set_title('Distribution of RBNS Z-scores in eCLIP Peaks', fontsize=14)

    ax.legend(loc='upper right')

    # Add counts annotation
    n_canonical = len(scored_df[scored_df['score_max'] >= canonical_threshold])
    n_discrepant = len(scored_df[scored_df['score_max'] < discrepant_threshold])
    n_intermediate = len(scored_df) - n_canonical - n_discrepant

    text = f'Canonical: {n_canonical}\nIntermediate: {n_intermediate}\nDiscrepant: {n_discrepant}'
    ax.text(
        0.95, 0.95, text,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment='top',
        horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved Z-score distribution plot to {output_path}")


def plot_enrichment_barplot(
    results_df: pd.DataFrame,
    output_path: str,
    title: Optional[str] = None,
    figsize: Tuple[int, int] = (10, 6),
    show_significance: bool = True
) -> None:
    """
    Plot enrichment odds ratios as a barplot.

    Parameters
    ----------
    results_df : pd.DataFrame
        Enrichment results with 'modification', 'odds_ratio', 'pvalue_adj'
    output_path : str
        Path to save the figure
    title : str, optional
        Plot title
    figsize : tuple
        Figure size (width, height)
    show_significance : bool
        Show significance stars (default: True)
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Sort by odds ratio
    plot_df = results_df.sort_values('odds_ratio', ascending=True)

    # Create color based on significance
    if 'pvalue_adj' in plot_df.columns:
        colors = ['forestgreen' if sig else 'gray'
                  for sig in plot_df['significant']]
    else:
        colors = 'steelblue'

    # Horizontal bar plot
    bars = ax.barh(
        plot_df['modification'],
        plot_df['odds_ratio'],
        color=colors,
        edgecolor='black'
    )

    # Add reference line at OR=1
    ax.axvline(x=1, color='black', linestyle='-', linewidth=1, label='No enrichment')

    # Add significance stars
    if show_significance and 'pvalue_adj' in plot_df.columns:
        for i, (_, row) in enumerate(plot_df.iterrows()):
            if row['pvalue_adj'] < 0.001:
                star = '***'
            elif row['pvalue_adj'] < 0.01:
                star = '**'
            elif row['pvalue_adj'] < 0.05:
                star = '*'
            else:
                star = ''

            if star:
                ax.text(
                    row['odds_ratio'] + 0.1,
                    i,
                    star,
                    va='center',
                    fontsize=12,
                    fontweight='bold'
                )

    # Labels
    ax.set_xlabel('Odds Ratio (Discrepant / Canonical)', fontsize=12)
    ax.set_ylabel('Modification Type', fontsize=12)

    if title:
        ax.set_title(title, fontsize=14)
    else:
        ax.set_title('Modification Enrichment in Discrepant Peaks', fontsize=14)

    # Legend
    if 'pvalue_adj' in results_df.columns:
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='forestgreen', edgecolor='black', label='Significant (adj. p < 0.05)'),
            Patch(facecolor='gray', edgecolor='black', label='Not significant')
        ]
        ax.legend(handles=legend_elements, loc='lower right')

    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved enrichment barplot to {output_path}")


def plot_classification_summary(
    scored_df: pd.DataFrame,
    output_path: str,
    title: Optional[str] = None,
    figsize: Tuple[int, int] = (8, 8)
) -> None:
    """
    Plot pie chart of peak classifications.

    Parameters
    ----------
    scored_df : pd.DataFrame
        DataFrame with 'category' column
    output_path : str
        Path to save the figure
    title : str, optional
        Plot title
    figsize : tuple
        Figure size (width, height)
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Count categories
    if 'category' in scored_df.columns:
        counts = scored_df['category'].value_counts()
    else:
        # Calculate from score_max
        canonical = len(scored_df[scored_df['score_max'] >= 3.0])
        discrepant = len(scored_df[scored_df['score_max'] < 1.5])
        intermediate = len(scored_df) - canonical - discrepant
        counts = pd.Series({
            'canonical': canonical,
            'intermediate': intermediate,
            'discrepant': discrepant
        })

    # Colors
    colors = {
        'canonical': 'forestgreen',
        'intermediate': 'gray',
        'discrepant': 'indianred'
    }

    # Plot
    wedges, texts, autotexts = ax.pie(
        counts.values,
        labels=counts.index,
        colors=[colors.get(c, 'steelblue') for c in counts.index],
        autopct='%1.1f%%',
        startangle=90,
        explode=[0.05 if c == 'discrepant' else 0 for c in counts.index]
    )

    # Style autotexts
    for autotext in autotexts:
        autotext.set_fontsize(11)
        autotext.set_fontweight('bold')

    if title:
        ax.set_title(title, fontsize=14)
    else:
        ax.set_title('Peak Classification Summary', fontsize=14)

    # Add count annotation
    total = counts.sum()
    text = f'Total: {total:,} peaks'
    ax.text(
        0.5, -0.1, text,
        transform=ax.transAxes,
        fontsize=11,
        ha='center'
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved classification summary to {output_path}")


def plot_heatmap_across_rbps(
    aggregated_df: pd.DataFrame,
    output_path: str,
    metric: str = 'odds_ratio',
    title: Optional[str] = None,
    figsize: Tuple[int, int] = (12, 8),
    cmap: str = 'RdYlGn'
) -> None:
    """
    Plot heatmap of enrichment metrics across RBPs and modifications.

    Parameters
    ----------
    aggregated_df : pd.DataFrame
        Aggregated results with 'rbp', 'modification', and metric columns
    output_path : str
        Path to save the figure
    metric : str
        Column to plot: 'odds_ratio', 'pvalue_adj', 'discrepant_pct'
    title : str, optional
        Plot title
    figsize : tuple
        Figure size (width, height)
    cmap : str
        Colormap name
    """
    # Pivot to matrix form
    matrix = aggregated_df.pivot(
        index='rbp',
        columns='modification',
        values=metric
    )

    fig, ax = plt.subplots(figsize=figsize)

    # For p-values, use log scale
    if 'pvalue' in metric:
        matrix = -np.log10(matrix + 1e-10)
        cbar_label = f'-log10({metric})'
    else:
        cbar_label = metric.replace('_', ' ').title()

    # Create heatmap
    sns.heatmap(
        matrix,
        ax=ax,
        cmap=cmap,
        center=1 if metric == 'odds_ratio' else None,
        annot=True,
        fmt='.2f',
        linewidths=0.5,
        cbar_kws={'label': cbar_label}
    )

    if title:
        ax.set_title(title, fontsize=14)
    else:
        ax.set_title(f'{metric.replace("_", " ").title()} Across RBPs and Modifications', fontsize=14)

    ax.set_xlabel('Modification', fontsize=12)
    ax.set_ylabel('RBP', fontsize=12)

    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved heatmap to {output_path}")


def plot_score_vs_enrichment(
    scored_df: pd.DataFrame,
    output_path: str,
    x_col: str = 'score_max',
    title: Optional[str] = None,
    figsize: Tuple[int, int] = (10, 6)
) -> None:
    """
    Plot relationship between RBNS score and eCLIP enrichment.

    Parameters
    ----------
    scored_df : pd.DataFrame
        DataFrame with score and enrichment columns
    output_path : str
        Path to save the figure
    x_col : str
        Column for x-axis (default: 'score_max')
    title : str, optional
        Plot title
    figsize : tuple
        Figure size (width, height)
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Filter valid data
    plot_df = scored_df[scored_df[x_col] > float('-inf')].copy()

    # Color by category if available
    if 'category' in plot_df.columns:
        colors = {
            'canonical': 'forestgreen',
            'intermediate': 'gray',
            'discrepant': 'indianred'
        }
        for cat in ['discrepant', 'intermediate', 'canonical']:
            subset = plot_df[plot_df['category'] == cat]
            ax.scatter(
                subset[x_col],
                subset['score_sum'],
                c=colors[cat],
                alpha=0.5,
                label=cat.capitalize(),
                s=20
            )
        ax.legend()
    else:
        ax.scatter(plot_df[x_col], plot_df['score_sum'], alpha=0.5, s=20)

    ax.set_xlabel('Maximum 5-mer Z-score (Score_max)', fontsize=12)
    ax.set_ylabel('Sum of Z-scores (Score_sum)', fontsize=12)

    if title:
        ax.set_title(title, fontsize=14)
    else:
        ax.set_title('RBNS Score Metrics Comparison', fontsize=14)

    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved score comparison plot to {output_path}")


def generate_all_plots(
    scored_df: pd.DataFrame,
    enrichment_df: pd.DataFrame,
    output_dir: str,
    rbp_name: Optional[str] = None
) -> None:
    """
    Generate all standard plots for a single RBP analysis.

    Parameters
    ----------
    scored_df : pd.DataFrame
        Scored peaks DataFrame
    enrichment_df : pd.DataFrame
        Enrichment results DataFrame
    output_dir : str
        Output directory for plots
    rbp_name : str, optional
        RBP name for plot titles
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    title_prefix = f"{rbp_name}: " if rbp_name else ""

    # Z-score distribution
    plot_zscore_distribution(
        scored_df,
        str(output_dir / 'zscore_distribution.png'),
        title=f"{title_prefix}RBNS Z-score Distribution"
    )

    # Classification pie chart
    plot_classification_summary(
        scored_df,
        str(output_dir / 'classification_summary.png'),
        title=f"{title_prefix}Peak Classification"
    )

    # Enrichment barplot
    if len(enrichment_df) > 0:
        plot_enrichment_barplot(
            enrichment_df,
            str(output_dir / 'enrichment_barplot.png'),
            title=f"{title_prefix}Modification Enrichment"
        )

    # Score comparison
    if 'score_sum' in scored_df.columns:
        plot_score_vs_enrichment(
            scored_df,
            str(output_dir / 'score_comparison.png'),
            title=f"{title_prefix}Score Metrics"
        )

    logger.info(f"Generated all plots in {output_dir}")


def plot_ranked_enrichment_analysis(
    scored_df: pd.DataFrame,
    mod_name: str,
    output_path: str,
    stats: dict,
    title_prefix: str = '',
    figsize: Tuple[int, int] = (15, 5)
) -> None:
    """
    Three-panel ranked enrichment figure for one modification.

    Panels: (1) CDF of Z-scores for mod+ vs mod- peaks,
            (2) Violin plot of Z-score distributions,
            (3) Modification frequency across Z-score quantile bins.

    Parameters
    ----------
    scored_df : pd.DataFrame
        Scored peaks with 'score_max' and 'has_{mod_name}' boolean columns
    mod_name : str
        Modification name (used to find 'has_{mod_name}' column)
    output_path : str
        Path to save the figure
    stats : dict
        Row from ranked_enrichment_results as a dict (for annotation)
    title_prefix : str
        Optional prefix for figure title (e.g. RBP name)
    figsize : tuple
        Figure size (width, height)
    """
    col = f'has_{mod_name}'
    if col not in scored_df.columns:
        logger.warning(f"Column {col} not in scored_df; skipping ranked plot.")
        return

    fig, axes = plt.subplots(1, 3, figsize=figsize)
    prefix = f"{title_prefix}: " if title_prefix else ""

    pos_z = scored_df.loc[scored_df[col], 'score_max'].dropna().sort_values()
    neg_z = scored_df.loc[~scored_df[col], 'score_max'].dropna().sort_values()

    med_pos = stats.get('median_z_mod_pos')
    med_neg = stats.get('median_z_mod_neg')
    mw_p_adj = stats.get('mann_whitney_p_adj', float('nan'))
    direction = stats.get('direction', '')
    med_diff = stats.get('median_z_diff')

    # -- Panel 1: CDF --
    ax = axes[0]
    for z_vals, label, color in [
        (pos_z, f'mod+ (n={len(pos_z)})', '#e74c3c'),
        (neg_z, f'mod- (n={len(neg_z)})', '#3498db'),
    ]:
        if len(z_vals) > 0:
            y = np.arange(1, len(z_vals) + 1) / len(z_vals)
            ax.plot(z_vals, y, label=label, color=color, linewidth=1.5)
            med = z_vals.median()
            ax.axvline(med, color=color, linestyle='--', linewidth=1, alpha=0.7)

    ax.set_xlabel('RBNS Z-score', fontsize=10)
    ax.set_ylabel('Cumulative fraction', fontsize=10)
    ax.set_title(f'{prefix}{mod_name}\nCDF', fontsize=11)
    ax.legend(fontsize=8)

    # -- Panel 2: Violin --
    ax = axes[1]
    plot_data = []
    group_labels = []
    for z_vals, label in [(pos_z, f'mod+\n(n={len(pos_z)})'),
                          (neg_z, f'mod-\n(n={len(neg_z)})')]:
        plot_data.append(z_vals.values)
        group_labels.append(label)

    parts = ax.violinplot(
        [d for d in plot_data if len(d) > 0],
        positions=range(len([d for d in plot_data if len(d) > 0])),
        showmedians=True
    )
    for pc in parts.get('bodies', []):
        pc.set_alpha(0.6)
    ax.set_xticks(range(len(group_labels)))
    ax.set_xticklabels(group_labels, fontsize=9)
    ax.set_ylabel('RBNS Z-score', fontsize=10)
    ax.set_title(f'Z-score by Modification', fontsize=11)

    # Annotation box
    if med_diff is not None and not pd.isna(mw_p_adj):
        annot = f'Δmedian={med_diff:+.3f}\n{direction}\np_adj={mw_p_adj:.2e}'
        ax.text(0.97, 0.97, annot,
                transform=ax.transAxes, fontsize=8,
                va='top', ha='right',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    # -- Panel 3: Binned frequency --
    ax = axes[2]
    bin_freqs_raw = stats.get('z_bin_frequencies', '{}')
    try:
        bin_freqs = json.loads(bin_freqs_raw) if isinstance(bin_freqs_raw, str) else bin_freqs_raw
    except (json.JSONDecodeError, TypeError):
        bin_freqs = {}

    overall_pct = stats.get('pct_mod_positive', None)

    if bin_freqs:
        bins = list(bin_freqs.keys())
        freqs = [bin_freqs[b] for b in bins]
        ax.bar(bins, freqs, color='steelblue', edgecolor='black', alpha=0.7)
        if overall_pct is not None:
            ax.axhline(overall_pct, color='red', linestyle='--',
                       linewidth=1.5, label=f'Overall {overall_pct:.1f}%')
            ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, 'No bin data', transform=ax.transAxes,
                ha='center', va='center', fontsize=10)

    ax.set_xlabel('Z-score Quantile', fontsize=10)
    ax.set_ylabel('Modification frequency (%)', fontsize=10)
    ax.set_title('Freq. by Z-score Quantile', fontsize=11)

    fig.suptitle(f'{prefix}{mod_name} Ranked Enrichment', fontsize=12, y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved ranked enrichment plot to {output_path}")


def plot_all_mods_ranked_summary(
    results_df: pd.DataFrame,
    output_path: str,
    title_prefix: str = '',
    figsize: Tuple[int, int] = (10, max(4, 1))
) -> None:
    """
    Two-panel summary figure: forest plot of median differences + Spearman rho.

    Parameters
    ----------
    results_df : pd.DataFrame
        Ranked enrichment results DataFrame (all modifications for one RBP)
    output_path : str
        Path to save the figure
    title_prefix : str
        Optional prefix for figure title (e.g. RBP name)
    figsize : tuple
        Figure size (width, height) — height auto-scales with n modifications
    """
    n_mods = len(results_df)
    if n_mods == 0:
        logger.warning("Empty results_df; skipping ranked summary plot.")
        return

    fig_height = max(4, n_mods * 0.8 + 1)
    fig, axes = plt.subplots(1, 2, figsize=(10, fig_height),
                             gridspec_kw={'width_ratios': [3, 1]})

    prefix = f"{title_prefix}: " if title_prefix else ""
    mods = results_df['modification'].tolist()
    y_pos = np.arange(n_mods)

    # -- Panel 1: Forest plot of median differences --
    ax = axes[0]
    colors = [
        'forestgreen' if sig else 'gray'
        for sig in results_df.get('significant', [False] * n_mods)
    ]
    diffs = results_df['median_z_diff'].fillna(0).tolist()
    ax.barh(y_pos, diffs, color=colors, edgecolor='black', height=0.5)
    ax.axvline(0, color='black', linewidth=1.2)

    # Significance stars
    for i, (_, row) in enumerate(results_df.iterrows()):
        p_adj = row.get('mann_whitney_p_adj', 1.0)
        diff = row.get('median_z_diff', 0) or 0
        if pd.isna(p_adj):
            continue
        if p_adj < 0.001:
            star = '***'
        elif p_adj < 0.01:
            star = '**'
        elif p_adj < 0.05:
            star = '*'
        else:
            star = ''
        if star:
            x_offset = 0.05 if diff >= 0 else -0.05
            ha = 'left' if diff >= 0 else 'right'
            ax.text(diff + x_offset, i, star, va='center',
                    ha=ha, fontsize=11, fontweight='bold')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(mods, fontsize=10)
    ax.set_xlabel('Median Z-score difference\n(mod+ minus mod-)', fontsize=10)
    ax.set_title(f'{prefix}Effect Size (Δ median Z)', fontsize=11)

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='forestgreen', edgecolor='black', label='Significant (p_adj<0.05)'),
        Patch(facecolor='gray', edgecolor='black', label='Not significant'),
    ]
    ax.legend(handles=legend_elements, fontsize=8, loc='lower right')

    # -- Panel 2: Spearman rho heatmap --
    ax = axes[1]
    rho_vals = results_df['spearman_rho'].fillna(0).values.reshape(-1, 1)
    vmax = max(abs(rho_vals).max(), 0.01)

    im = ax.imshow(rho_vals, cmap='RdBu_r', aspect='auto',
                   vmin=-vmax, vmax=vmax)
    ax.set_xticks([0])
    ax.set_xticklabels(['Spearman ρ'], fontsize=9)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([])

    for i, rho in enumerate(rho_vals.flatten()):
        ax.text(0, i, f'{rho:.3f}', ha='center', va='center',
                fontsize=9, color='black' if abs(rho) < vmax * 0.6 else 'white')

    plt.colorbar(im, ax=ax, fraction=0.15, pad=0.05)
    ax.set_title('Spearman ρ\n(neg = low-Z enrich)', fontsize=9)

    fig.suptitle(f'{prefix}Ranked Enrichment Summary', fontsize=12)
    plt.tight_layout()
    plt.savefig(output_path, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved ranked enrichment summary to {output_path}")


if __name__ == '__main__':
    # Example/test usage
    import sys
    logging.basicConfig(level=logging.INFO)

    # Create test data
    np.random.seed(42)
    n = 1000

    test_scored = pd.DataFrame({
        'peak_id': [f'peak_{i}' for i in range(n)],
        'score_max': np.random.normal(2, 1.5, n),
        'score_sum': np.random.normal(10, 5, n)
    })

    test_enrichment = pd.DataFrame({
        'modification': ['m6A', 'pseudoU', 'm5C', 'ac4C'],
        'odds_ratio': [2.5, 1.8, 0.9, 1.2],
        'pvalue': [0.001, 0.02, 0.5, 0.1],
        'pvalue_adj': [0.004, 0.04, 0.5, 0.15],
        'significant': [True, True, False, False],
        'discrepant_overlap': [50, 30, 10, 15],
        'canonical_overlap': [20, 15, 12, 13],
        'discrepant_total': [100, 100, 100, 100],
        'canonical_total': [100, 100, 100, 100]
    })

    # Generate test plots
    generate_all_plots(
        test_scored,
        test_enrichment,
        './test_plots',
        rbp_name='TEST_RBP'
    )

    print("Test plots generated in ./test_plots/")
