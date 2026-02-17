#!/usr/bin/env python3
"""
Process RBNS Enrichment Files: Convert R-values to Z-scores

This script converts ENCODE RBNS 5-mer enrichment TSV files (R-values) to
Z-score CSV format compatible with the RNAMod-RBNS pipeline.

Input format (ENCODE enrichment TSV):
    [RBP_NAME]    0    5    20    80    320    1300
    AAAAA    0.98    1.00    1.02    1.05    1.08    1.10
    AAAAC    0.97    0.99    1.01    1.04    1.07    1.09
    ...

Output format (Z-scores CSV):
    kmer,z_score,r_value
    GCAUG,8.45,3.25
    GCACG,7.23,2.89
    ...

Z-score calculation per Dominguez et al.:
    Z = (R - mean_R) / std_R

Usage:
    python scripts/process_rbns_enrichment.py [--input-dir data/rbns] [--output-dir data/rbns]

    # Or process a single file:
    python scripts/process_rbns_enrichment.py --single IGF2BP1_enrichment.tsv IGF2BP1_zscores.csv
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import sys
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def parse_encode_enrichment_tsv(tsv_path):
    """
    Parse ENCODE RBNS enrichment TSV file.

    The file format is:
    - First row: header with [RBP_NAME] followed by protein concentrations
    - Subsequent rows: k-mer followed by R-values at each concentration

    Args:
        tsv_path: Path to ENCODE enrichment TSV file

    Returns:
        tuple: (kmer_col_name, concentration_columns, dataframe)
    """
    # Read the TSV file
    df = pd.read_csv(tsv_path, sep='\t', header=0)

    # First column contains k-mers, rest are concentrations
    kmer_col = df.columns[0]
    conc_cols = list(df.columns[1:])

    # Clean up k-mer column - remove any whitespace
    df[kmer_col] = df[kmer_col].str.strip().str.upper()

    # Convert concentration columns to numeric
    for col in conc_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    return kmer_col, conc_cols, df


def convert_rvalues_to_zscores(df, kmer_col, conc_cols, use_max_conc=True):
    """
    Convert R-values to Z-scores.

    Per Dominguez et al. and the RBNS_pipeline:
        Z = (R - mean_R) / std_R

    Args:
        df: DataFrame with k-mers and R-values
        kmer_col: Name of the k-mer column
        conc_cols: List of concentration column names
        use_max_conc: If True, use highest concentration; else use max across all

    Returns:
        DataFrame with columns: kmer, z_score, r_value
    """
    if use_max_conc:
        # Use the highest protein concentration (last column)
        # This represents maximum binding/enrichment conditions
        r_values = df[conc_cols[-1]].copy()
        conc_used = conc_cols[-1]
    else:
        # Use maximum R-value across all concentrations per k-mer
        r_values = df[conc_cols].max(axis=1)
        conc_used = "max_across_all"

    # Remove NaN values for statistics
    valid_mask = ~r_values.isna()
    valid_r = r_values[valid_mask]

    # Calculate Z-scores: Z = (R - mean) / std
    mean_r = valid_r.mean()
    std_r = valid_r.std()

    if std_r == 0 or np.isnan(std_r):
        logger.warning("Standard deviation is zero or NaN - cannot calculate Z-scores")
        z_scores = pd.Series([0.0] * len(r_values))
    else:
        z_scores = (r_values - mean_r) / std_r

    # Create output DataFrame
    result = pd.DataFrame({
        'kmer': df[kmer_col],
        'z_score': z_scores,
        'r_value': r_values
    })

    # Remove rows with NaN values
    result = result.dropna()

    # Sort by Z-score descending (highest affinity first)
    result = result.sort_values('z_score', ascending=False).reset_index(drop=True)

    return result, mean_r, std_r, conc_used


def process_single_file(input_path, output_path, use_max_conc=True):
    """
    Process a single ENCODE enrichment file and save Z-scores.

    Args:
        input_path: Path to input TSV file
        output_path: Path to output CSV file
        use_max_conc: Use highest concentration (True) or max across all (False)

    Returns:
        dict: Statistics about the conversion
    """
    input_path = Path(input_path)
    output_path = Path(output_path)

    if not input_path.exists():
        logger.error(f"Input file not found: {input_path}")
        return None

    logger.info(f"Processing: {input_path.name}")

    # Parse the enrichment file
    kmer_col, conc_cols, df = parse_encode_enrichment_tsv(input_path)

    logger.info(f"  Found {len(df)} k-mers, {len(conc_cols)} concentrations")
    logger.info(f"  Concentrations: {conc_cols}")

    # Convert to Z-scores
    result, mean_r, std_r, conc_used = convert_rvalues_to_zscores(
        df, kmer_col, conc_cols, use_max_conc
    )

    # Save output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, index=False)

    # Calculate statistics
    stats = {
        'input_file': input_path.name,
        'output_file': output_path.name,
        'n_kmers': len(result),
        'concentration_used': conc_used,
        'mean_r': mean_r,
        'std_r': std_r,
        'z_min': result['z_score'].min(),
        'z_max': result['z_score'].max(),
        'top_kmer': result.iloc[0]['kmer'] if len(result) > 0 else 'N/A',
        'top_zscore': result.iloc[0]['z_score'] if len(result) > 0 else 0
    }

    logger.info(f"  Z-score range: {stats['z_min']:.2f} to {stats['z_max']:.2f}")
    logger.info(f"  Top k-mer: {stats['top_kmer']} (Z={stats['top_zscore']:.2f})")
    logger.info(f"  Output: {output_path}")

    return stats


def process_all_files(input_dir, output_dir=None, use_max_conc=True):
    """
    Process all ENCODE enrichment files in a directory.

    Args:
        input_dir: Directory containing *_enrichment.tsv files
        output_dir: Directory for output files (default: same as input)
        use_max_conc: Use highest concentration (True) or max across all (False)

    Returns:
        list: Statistics for each processed file
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_dir) if output_dir else input_dir

    # Find all enrichment TSV files
    enrichment_files = list(input_dir.glob("*_enrichment.tsv"))

    if not enrichment_files:
        logger.warning(f"No *_enrichment.tsv files found in {input_dir}")
        return []

    logger.info(f"Found {len(enrichment_files)} enrichment files")
    logger.info("")

    all_stats = []
    success_count = 0
    fail_count = 0

    for tsv_path in sorted(enrichment_files):
        # Extract RBP name from filename (e.g., IGF2BP1_enrichment.tsv -> IGF2BP1)
        rbp_name = tsv_path.stem.replace('_enrichment', '')
        output_path = output_dir / f"{rbp_name}_zscores.csv"

        stats = process_single_file(tsv_path, output_path, use_max_conc)

        if stats:
            all_stats.append(stats)
            success_count += 1
        else:
            fail_count += 1

        logger.info("")

    # Print summary
    logger.info("=" * 60)
    logger.info("PROCESSING SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Successful: {success_count}")
    logger.info(f"Failed: {fail_count}")
    logger.info("")

    if all_stats:
        logger.info("Top motifs per RBP:")
        for stats in sorted(all_stats, key=lambda x: x['input_file']):
            rbp = stats['input_file'].replace('_enrichment.tsv', '')
            logger.info(f"  {rbp}: {stats['top_kmer']} (Z={stats['top_zscore']:.2f})")

    return all_stats


def validate_zscore_file(csv_path):
    """
    Validate a Z-score CSV file has the expected format.

    Args:
        csv_path: Path to Z-score CSV file

    Returns:
        bool: True if valid, False otherwise
    """
    try:
        df = pd.read_csv(csv_path)

        required_cols = {'kmer', 'z_score', 'r_value'}
        if not required_cols.issubset(df.columns):
            logger.error(f"Missing columns. Expected: {required_cols}, Found: {set(df.columns)}")
            return False

        if len(df) == 0:
            logger.error("File is empty")
            return False

        # Check k-mer length consistency
        kmer_lengths = df['kmer'].str.len().unique()
        if len(kmer_lengths) != 1:
            logger.warning(f"Inconsistent k-mer lengths: {kmer_lengths}")

        # Check for NaN values
        nan_count = df['z_score'].isna().sum()
        if nan_count > 0:
            logger.warning(f"Found {nan_count} NaN Z-scores")

        logger.info(f"Validated: {csv_path}")
        logger.info(f"  K-mers: {len(df)}, Z-score range: [{df['z_score'].min():.2f}, {df['z_score'].max():.2f}]")

        return True

    except Exception as e:
        logger.error(f"Validation failed: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description='Convert ENCODE RBNS R-value enrichment files to Z-scores',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Process all enrichment files in data/rbns/
    python scripts/process_rbns_enrichment.py

    # Process with explicit directories
    python scripts/process_rbns_enrichment.py --input-dir data/rbns --output-dir data/rbns

    # Process a single file
    python scripts/process_rbns_enrichment.py --single data/rbns/IGF2BP1_enrichment.tsv data/rbns/IGF2BP1_zscores.csv

    # Validate output files
    python scripts/process_rbns_enrichment.py --validate data/rbns
        """
    )

    parser.add_argument(
        '--input-dir',
        default='data/rbns',
        help='Directory containing *_enrichment.tsv files (default: data/rbns)'
    )
    parser.add_argument(
        '--output-dir',
        default=None,
        help='Directory for output *_zscores.csv files (default: same as input)'
    )
    parser.add_argument(
        '--single',
        nargs=2,
        metavar=('INPUT', 'OUTPUT'),
        help='Process a single file: --single input.tsv output.csv'
    )
    parser.add_argument(
        '--validate',
        metavar='DIR',
        help='Validate all *_zscores.csv files in directory'
    )
    parser.add_argument(
        '--use-max-across-concentrations',
        action='store_true',
        help='Use max R-value across all concentrations instead of highest concentration'
    )

    args = parser.parse_args()

    # Determine whether to use max concentration or max across all
    use_max_conc = not args.use_max_across_concentrations

    if args.validate:
        # Validation mode
        validate_dir = Path(args.validate)
        zscore_files = list(validate_dir.glob("*_zscores.csv"))

        if not zscore_files:
            logger.error(f"No *_zscores.csv files found in {validate_dir}")
            sys.exit(1)

        logger.info(f"Validating {len(zscore_files)} Z-score files")
        logger.info("")

        all_valid = True
        for csv_path in sorted(zscore_files):
            if not validate_zscore_file(csv_path):
                all_valid = False

        sys.exit(0 if all_valid else 1)

    elif args.single:
        # Single file mode
        input_path, output_path = args.single
        stats = process_single_file(input_path, output_path, use_max_conc)
        sys.exit(0 if stats else 1)

    else:
        # Batch processing mode
        logger.info("=" * 60)
        logger.info("RBNS R-value to Z-score Conversion")
        logger.info("=" * 60)
        logger.info(f"Input directory: {args.input_dir}")
        logger.info(f"Output directory: {args.output_dir or args.input_dir}")
        logger.info(f"Using: {'highest concentration' if use_max_conc else 'max across concentrations'}")
        logger.info("")

        stats = process_all_files(args.input_dir, args.output_dir, use_max_conc)

        if not stats:
            logger.error("No files were processed successfully")
            sys.exit(1)

        logger.info("")
        logger.info("Next step: Run the main analysis pipeline")
        logger.info("  python src/main.py --rbp IGF2BP1 --rbns data/rbns/IGF2BP1_zscores.csv ...")


if __name__ == '__main__':
    main()
