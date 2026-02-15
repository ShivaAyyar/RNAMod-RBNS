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
