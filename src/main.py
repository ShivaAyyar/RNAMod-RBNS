#!/usr/bin/env python3
"""
Epitranscriptomic RBP Analysis Pipeline - Main Orchestration Script

This script orchestrates the complete analysis pipeline for investigating
the "Specificity Paradox" - where RNA Binding Proteins bind sequences
in vivo (eCLIP) that they show low affinity for in vitro (RBNS).

Hypothesis: RNA modifications (m6A, Ψ, m5C, ac4C) present in cellular RNA
but absent in synthetic RBNS libraries explain this discrepancy.

Usage:
    python main.py --rbp IGF2BP1 --rbns data/rbns/IGF2BP1_zscores.csv \
                   --eclip data/eclip/IGF2BP1_K562.bed \
                   --genome data/genome/hg38.fa \
                   --chrom-sizes data/genome/hg38.chrom.sizes \
                   --mods data/mods/K562/m6A.bed data/mods/K562/pseudoU.bed \
                   --mod-names m6A pseudoU \
                   --output results/IGF2BP1

References:
- Dominguez et al. (2018) Molecular Cell - RBNS methodology
- ENCODE eCLIP data standards
- REPIC/RMBase for modification data
"""

import argparse
import logging
import sys
from pathlib import Path
from datetime import datetime

import pandas as pd

# Import analysis modules
import peak_analysis as pa
import enrichment_analysis as ea
import visualization as viz

# The 24 "Gold Standard" RBPs with both RBNS and eCLIP data
GOLD_STANDARD_RBPS = [
    'IGF2BP1', 'IGF2BP2', 'HNRNPC', 'TIA1', 'HNRNPK', 'PCBP2',
    'RBFOX2', 'PTBP3', 'TARDBP', 'QKI', 'SRSF1', 'SRSF9',
    'RBM22', 'TRA2A', 'HNRNPL', 'LIN28B', 'ZNF326', 'FUS',
    'MATR3', 'HNRNPA1', 'HNRNPM', 'NONO', 'U2AF2', 'EWSR1'
]

# Known m6A readers (positive controls)
M6A_READERS = ['IGF2BP1', 'IGF2BP2']

# Standard modification types
MODIFICATION_TYPES = ['m6A', 'pseudoU', 'm5C', 'ac4C']


def setup_logging(output_dir: Path, rbp_name: str) -> logging.Logger:
    """Configure logging to both file and console."""
    log_file = output_dir / 'analysis.log'

    # Create logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Clear existing handlers
    logger.handlers = []

    # Console handler
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.INFO)
    console_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console.setFormatter(console_format)
    logger.addHandler(console)

    # File handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_format)
    logger.addHandler(file_handler)

    logger.info(f"Analysis started for {rbp_name}")
    logger.info(f"Log file: {log_file}")

    return logger


def validate_inputs(args) -> bool:
    """Validate that all required input files exist."""
    logger = logging.getLogger(__name__)
    valid = True

    # Check RBNS file
    if not Path(args.rbns).exists():
        logger.error(f"RBNS file not found: {args.rbns}")
        valid = False

    # Check eCLIP file
    if not Path(args.eclip).exists():
        logger.error(f"eCLIP file not found: {args.eclip}")
        valid = False

    # Check genome files
    if not Path(args.genome).exists():
        logger.error(f"Genome file not found: {args.genome}")
        valid = False

    if not Path(args.chrom_sizes).exists():
        logger.error(f"Chromosome sizes file not found: {args.chrom_sizes}")
        valid = False

    # Check modification files
    for mod_file in args.mods:
        if not Path(mod_file).exists():
            logger.error(f"Modification file not found: {mod_file}")
            valid = False

    return valid


def run_analysis(args) -> dict:
    """
    Run the complete analysis pipeline for a single RBP.

    Returns
    -------
    dict
        Summary statistics and results paths
    """
    logger = logging.getLogger(__name__)

    # Create output directory
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    # Setup logging
    setup_logging(outdir, args.rbp)

    logger.info("=" * 60)
    logger.info(f"Epitranscriptomic Analysis Pipeline")
    logger.info(f"RBP: {args.rbp}")
    logger.info(f"Cell Line: {args.cell_line}")
    logger.info(f"Timestamp: {datetime.now().isoformat()}")
    logger.info("=" * 60)

    # Check if this is a positive control
    if args.rbp in M6A_READERS:
        logger.info(f"*** {args.rbp} is a known m6A reader - POSITIVE CONTROL ***")

    # Validate inputs
    if not validate_inputs(args):
        logger.error("Input validation failed. Exiting.")
        sys.exit(1)

    results = {
        'rbp': args.rbp,
        'cell_line': args.cell_line,
        'timestamp': datetime.now().isoformat()
    }

    # =========================================================================
    # Step 1: Load RBNS Z-scores (pre-computed from RBNS pipeline)
    # =========================================================================
    logger.info("\n--- Step 1: Loading RBNS Z-scores ---")
    rbns_dict = pa.load_rbns_zscores(args.rbns)
    results['n_kmers'] = len(rbns_dict)
    logger.info(f"Loaded {len(rbns_dict)} k-mer Z-scores")

    # =========================================================================
    # Step 2: Pre-process eCLIP peaks
    # =========================================================================
    logger.info("\n--- Step 2: Pre-processing eCLIP peaks ---")

    # Extend peaks 50nt in 5' direction
    peaks = pa.extend_peaks_5prime(args.eclip, args.chrom_sizes)
    results['n_peaks_initial'] = len(peaks)

    # Filter by enrichment
    peaks = pa.filter_by_enrichment(peaks, args.min_enrichment)
    results['n_peaks_filtered'] = len(peaks)

    # Save filtered peaks
    filtered_bed = str(outdir / 'peaks_filtered.bed')
    peaks.saveas(filtered_bed)
    logger.info(f"Saved {len(peaks)} filtered peaks to {filtered_bed}")

    # =========================================================================
    # Step 3: Extract sequences and score
    # =========================================================================
    logger.info("\n--- Step 3: Extracting sequences and scoring ---")
    sequences = pa.extract_sequences(filtered_bed, args.genome)

    scored_data = []
    for peak_id, seq in sequences.items():
        score_max, score_sum = pa.calculate_peak_scores(seq, rbns_dict)
        scored_data.append({
            'peak_id': peak_id,
            'sequence': seq,
            'length': len(seq),
            'score_max': score_max,
            'score_sum': score_sum
        })

    scored_df = pd.DataFrame(scored_data)

    # =========================================================================
    # Step 4: Classify peaks
    # =========================================================================
    logger.info("\n--- Step 4: Classifying peaks ---")
    canonical, discrepant, intermediate = pa.classify_peaks(
        scored_df,
        canonical_threshold=args.canonical_threshold,
        discrepant_threshold=args.discrepant_threshold
    )

    # Add category column
    scored_df['category'] = 'intermediate'
    scored_df.loc[scored_df['score_max'] >= args.canonical_threshold, 'category'] = 'canonical'
    scored_df.loc[scored_df['score_max'] < args.discrepant_threshold, 'category'] = 'discrepant'

    # Save scored peaks
    scored_df.to_csv(outdir / 'scored_peaks.csv', index=False)

    # Save classified peaks as BED files
    canonical_bed = str(outdir / 'canonical_peaks.bed')
    discrepant_bed = str(outdir / 'discrepant_peaks.bed')
    intermediate_bed = str(outdir / 'intermediate_peaks.bed')

    if len(canonical) > 0:
        pa.scored_df_to_bed(canonical, canonical_bed)
    else:
        # Create empty file
        Path(canonical_bed).touch()

    if len(discrepant) > 0:
        pa.scored_df_to_bed(discrepant, discrepant_bed)
    else:
        Path(discrepant_bed).touch()

    if len(intermediate) > 0:
        pa.scored_df_to_bed(intermediate, intermediate_bed)
    else:
        Path(intermediate_bed).touch()

    results['n_canonical'] = len(canonical)
    results['n_discrepant'] = len(discrepant)
    results['n_intermediate'] = len(intermediate)

    logger.info(f"Canonical peaks: {len(canonical)}")
    logger.info(f"Discrepant peaks: {len(discrepant)}")
    logger.info(f"Intermediate peaks: {len(intermediate)}")

    # =========================================================================
    # Step 5: Enrichment analysis
    # =========================================================================
    logger.info("\n--- Step 5: Running enrichment analysis ---")

    if len(canonical) > 0 and len(discrepant) > 0:
        enrichment_results = ea.run_enrichment_analysis(
            canonical_bed,
            discrepant_bed,
            args.mods,
            args.mod_names,
            chrom_sizes=args.chrom_sizes
        )

        # Add metadata
        enrichment_results['rbp'] = args.rbp
        enrichment_results['cell_line'] = args.cell_line

        # Save results
        enrichment_results.to_csv(outdir / 'enrichment_results.csv', index=False)

        # Check for significant findings
        if 'significant' in enrichment_results.columns:
            sig_results = enrichment_results[enrichment_results['significant']]
            results['n_significant'] = len(sig_results)

            if len(sig_results) > 0:
                logger.info("\n*** SIGNIFICANT ENRICHMENTS ***")
                for _, row in sig_results.iterrows():
                    logger.info(
                        f"  {row['modification']}: OR={row['odds_ratio']:.2f}, "
                        f"p_adj={row['pvalue_adj']:.2e}"
                    )

                # Special check for positive controls
                if args.rbp in M6A_READERS:
                    m6a_sig = sig_results[sig_results['modification'] == 'm6A']
                    if len(m6a_sig) > 0:
                        logger.info(f"\n*** POSITIVE CONTROL VALIDATED: {args.rbp} shows m6A enrichment ***")
                    else:
                        logger.warning(f"\n*** WARNING: Positive control {args.rbp} did NOT show m6A enrichment ***")
        else:
            results['n_significant'] = 0
    else:
        logger.warning("Insufficient peaks for enrichment analysis")
        enrichment_results = pd.DataFrame()
        results['n_significant'] = 0

    # =========================================================================
    # Step 6: Generate visualizations
    # =========================================================================
    logger.info("\n--- Step 6: Generating visualizations ---")

    figures_dir = outdir / 'figures'
    figures_dir.mkdir(exist_ok=True)

    viz.generate_all_plots(
        scored_df,
        enrichment_results if len(enrichment_results) > 0 else pd.DataFrame(),
        str(figures_dir),
        rbp_name=args.rbp
    )

    # =========================================================================
    # Summary
    # =========================================================================
    logger.info("\n" + "=" * 60)
    logger.info("ANALYSIS COMPLETE")
    logger.info("=" * 60)
    logger.info(f"RBP: {args.rbp}")
    logger.info(f"Total peaks analyzed: {len(scored_df)}")
    logger.info(f"Canonical: {results['n_canonical']} ({100*results['n_canonical']/len(scored_df):.1f}%)")
    logger.info(f"Discrepant: {results['n_discrepant']} ({100*results['n_discrepant']/len(scored_df):.1f}%)")
    logger.info(f"Significant enrichments: {results['n_significant']}")
    logger.info(f"Output directory: {outdir}")

    # Save summary
    summary_df = pd.DataFrame([results])
    summary_df.to_csv(outdir / 'summary.csv', index=False)

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Epitranscriptomic analysis of RBP binding discrepancies',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single RBP analysis
  python main.py --rbp IGF2BP1 --rbns data/rbns/IGF2BP1.csv \\
                 --eclip data/eclip/IGF2BP1_K562.bed \\
                 --genome data/genome/hg38.fa \\
                 --chrom-sizes data/genome/hg38.chrom.sizes \\
                 --mods data/mods/K562/m6A.bed data/mods/K562/pseudoU.bed \\
                 --mod-names m6A pseudoU \\
                 --output results/IGF2BP1

Gold Standard RBPs (24 total):
  IGF2BP1*, IGF2BP2*, HNRNPC, TIA1, HNRNPK, PCBP2, RBFOX2, PTBP3,
  TARDBP, QKI, SRSF1, SRSF9, RBM22, TRA2A, HNRNPL, LIN28B,
  ZNF326, FUS, MATR3, HNRNPA1, HNRNPM, NONO, U2AF2, EWSR1

  * = Known m6A readers (positive controls)
        """
    )

    # Required arguments
    parser.add_argument(
        '--rbp', required=True,
        help='RBP name (e.g., IGF2BP1)'
    )
    parser.add_argument(
        '--rbns', required=True,
        help='Path to RBNS Z-scores CSV (pre-computed from RBNS pipeline)'
    )
    parser.add_argument(
        '--eclip', required=True,
        help='Path to eCLIP narrowPeak BED file'
    )
    parser.add_argument(
        '--genome', required=True,
        help='Path to reference genome FASTA file'
    )
    parser.add_argument(
        '--chrom-sizes', required=True,
        help='Path to chromosome sizes file'
    )
    parser.add_argument(
        '--mods', nargs='+', required=True,
        help='Paths to modification BED files'
    )
    parser.add_argument(
        '--mod-names', nargs='+', required=True,
        help='Names of modification types (same order as --mods)'
    )
    parser.add_argument(
        '--output', required=True,
        help='Output directory'
    )

    # Optional arguments
    parser.add_argument(
        '--cell-line', default='K562',
        choices=['K562', 'HepG2'],
        help='Cell line (default: K562)'
    )
    parser.add_argument(
        '--canonical-threshold', type=float, default=3.0,
        help='Z-score threshold for canonical classification (default: 3.0)'
    )
    parser.add_argument(
        '--discrepant-threshold', type=float, default=1.5,
        help='Z-score threshold for discrepant classification (default: 1.5)'
    )
    parser.add_argument(
        '--min-enrichment', type=float, default=2.0,
        help='Minimum eCLIP enrichment (signalValue) threshold (default: 2.0)'
    )
    parser.add_argument(
        '--extension', type=int, default=50,
        help="5' extension for peaks in nucleotides (default: 50)"
    )

    args = parser.parse_args()

    # Validate mod arguments
    if len(args.mods) != len(args.mod_names):
        parser.error("--mods and --mod-names must have the same number of arguments")

    # Run analysis
    try:
        results = run_analysis(args)
        sys.exit(0)
    except Exception as e:
        logging.error(f"Analysis failed: {e}", exc_info=True)
        sys.exit(1)


if __name__ == '__main__':
    main()
