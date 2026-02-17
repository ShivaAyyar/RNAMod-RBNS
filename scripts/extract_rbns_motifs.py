#!/usr/bin/env python3
"""
Extract RBNS motif data from Dominguez et al. (2018) supplementary files.

This script extracts the TOP ENRICHED MOTIFS from the supplementary data.
NOTE: This is NOT comprehensive k-mer Z-score data (all 1024 5-mers).
      The supplementary only provides top ~5-10 motifs per RBP.

For comprehensive Z-scores, you would need to:
1. Download raw RBNS data from ENCODE using the accession IDs
2. Process with the RBNS_pipeline (github.com/cburgelab/RBNS_pipeline)

Usage:
    python scripts/extract_rbns_motifs.py
"""

import pandas as pd
import os
import re

# Configuration
DATA_DIR = "data/rbns"
OUTPUT_DIR = "data/rbns"

# Our 23 target RBPs
TARGET_RBPS = [
    'IGF2BP1', 'IGF2BP2', 'HNRNPC', 'TIA1', 'HNRNPK', 'PCBP2', 'RBFOX2',
    'PTBP1', 'TARDBP', 'QKI', 'SRSF1', 'SRSF9', 'RBM22', 'TRA2A',
    'HNRNPL', 'LIN28B', 'FUS', 'MATR3', 'HNRNPA1', 'HNRNPM',
    'NONO', 'U2AF2', 'EWSR1'
]


def parse_motif_string(motif_str):
    """
    Parse motif string like 'UUUUU_1_15.38' into (kmer, logo_num, zscore).
    """
    if pd.isna(motif_str):
        return None, None, None

    parts = str(motif_str).split('_')
    if len(parts) >= 3:
        kmer = parts[0]
        logo_num = int(parts[1])
        zscore = float(parts[2])
        return kmer, logo_num, zscore
    return None, None, None


def extract_top_motifs_from_table_s3():
    """
    Extract top motifs from Table S3 (mmc4.xlsx).
    Returns dict: {RBP: [(kmer, zscore), ...]}
    """
    xlsx_path = os.path.join(DATA_DIR, 'mmc4.xlsx')
    if not os.path.exists(xlsx_path):
        print(f"ERROR: {xlsx_path} not found")
        return {}

    df = pd.read_excel(xlsx_path, sheet_name=0, header=0)

    # Motif columns are: 'Motif 5mer_logonum_stepwiseRminus1' and Unnamed: 6-20
    motif_cols = [col for col in df.columns if 'Motif' in str(col) or 'Unnamed' in str(col)]

    results = {}

    for _, row in df.iterrows():
        rbp = row['RBP']
        if rbp not in TARGET_RBPS:
            continue

        motifs = []
        for col in motif_cols:
            kmer, logo_num, zscore = parse_motif_string(row.get(col))
            if kmer and zscore:
                # Convert RNA (U) to DNA (T) for consistency
                kmer_dna = kmer.replace('U', 'T')
                motifs.append((kmer_dna, zscore))

        if motifs:
            results[rbp] = motifs

    return results


def extract_logo_proportions():
    """
    Extract 5-mer logo proportions from mmc4.xlsx.
    This gives us which k-mers appear in the consensus motif logo.
    Returns dict: {RBP: [(kmer, proportion), ...]}
    """
    xlsx_path = os.path.join(DATA_DIR, 'mmc4.xlsx')
    if not os.path.exists(xlsx_path):
        return {}

    df = pd.read_excel(xlsx_path, sheet_name='logo_5mers.prop_in_logo', header=None)

    # Structure: alternating columns of RBP name and proportions
    # Row 0: NAME, Stepwise_R-1, blank, blank, ...
    # Row 1: A1CF, NaN, BOLL, NaN, CELF1, ...
    # Row 2+: kmer, proportion, kmer, proportion, ...

    results = {}

    # Get RBP names from row 1
    rbp_row = df.iloc[1]
    rbp_cols = {}  # col_idx -> RBP name

    for i, val in enumerate(rbp_row):
        if pd.notna(val) and val in TARGET_RBPS:
            rbp_cols[i] = val

    # Extract k-mers and proportions for each RBP
    for col_idx, rbp in rbp_cols.items():
        motifs = []
        for row_idx in range(2, len(df)):
            kmer = df.iloc[row_idx, col_idx]
            prop = df.iloc[row_idx, col_idx + 1]

            if pd.notna(kmer) and pd.notna(prop):
                # Convert T to T (already DNA in this sheet)
                kmer_str = str(kmer).upper()
                try:
                    prop_float = float(prop)
                    motifs.append((kmer_str, prop_float))
                except:
                    pass

        if motifs:
            results[rbp] = motifs

    return results


def write_csv_for_rbp(rbp, motifs, output_dir):
    """
    Write motif data to CSV in expected format.

    Note: This creates a PARTIAL z-score file with only the top motifs.
    The pipeline expects all 1024 5-mers, so this may need adaptation.
    """
    output_path = os.path.join(output_dir, f'{rbp}_zscores.csv')

    with open(output_path, 'w') as f:
        f.write('kmer,z_score\n')
        for kmer, zscore in motifs:
            f.write(f'{kmer},{zscore:.4f}\n')

    print(f"  Wrote {output_path} ({len(motifs)} motifs)")


def main():
    print("=" * 60)
    print("Extracting RBNS motif data from Dominguez et al. (2018)")
    print("=" * 60)
    print()
    print("WARNING: The supplementary data only contains TOP ENRICHED MOTIFS,")
    print("         NOT comprehensive k-mer Z-scores for all 1024 5-mers.")
    print()

    # Extract top motifs from Table S3
    print("Extracting top motifs from Table S3 (mmc4.xlsx)...")
    top_motifs = extract_top_motifs_from_table_s3()

    print(f"\nFound data for {len(top_motifs)}/{len(TARGET_RBPS)} target RBPs:")
    for rbp in TARGET_RBPS:
        if rbp in top_motifs:
            n_motifs = len(top_motifs[rbp])
            top_kmer = top_motifs[rbp][0][0] if top_motifs[rbp] else 'N/A'
            top_z = top_motifs[rbp][0][1] if top_motifs[rbp] else 0
            print(f"  {rbp:12} - {n_motifs:2} motifs (top: {top_kmer} Z={top_z:.2f})")
        else:
            print(f"  {rbp:12} - NOT FOUND in supplementary data")

    # Also extract logo proportions as alternative data
    print("\n\nExtracting logo proportions from 5-mer sheet...")
    logo_props = extract_logo_proportions()
    print(f"Found logo data for {len(logo_props)} RBPs")

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    missing_rbps = [rbp for rbp in TARGET_RBPS if rbp not in top_motifs]
    if missing_rbps:
        print(f"\nMissing RBPs ({len(missing_rbps)}):")
        for rbp in missing_rbps:
            print(f"  - {rbp}")

    print("\n" + "-" * 60)
    print("OPTIONS TO OBTAIN COMPREHENSIVE K-MER DATA:")
    print("-" * 60)
    print("""
1. ENCODE RBNS Raw Data:
   Download raw RBNS FASTQ files from ENCODE using accession IDs
   and process with RBNS_pipeline to compute all k-mer Z-scores.

   Example accessions (from Table S3):
   - IGF2BP1: ENCSR928XOW
   - IGF2BP2: ENCSR588GYZ
   - HNRNPC:  ENCSR569UIU

2. Use Top Motifs Only:
   Modify the pipeline to work with the available top ~5 motifs
   per RBP. Less comprehensive but immediately available.

3. Contact Authors:
   The full k-mer enrichment matrices may be available upon request
   from the Burge lab (burgelab@mit.edu).

4. RBNS Pipeline:
   https://github.com/cburgelab/RBNS_pipeline
   This tool processes raw RBNS data to compute k-mer Z-scores.
""")

    # Ask user what to do
    print("\nWould you like to create partial Z-score files using the top motifs?")
    print("(These will only contain ~5-10 k-mers instead of all 1024)")


if __name__ == '__main__':
    main()
