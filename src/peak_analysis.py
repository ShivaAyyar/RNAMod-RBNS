#!/usr/bin/env python3
"""
Peak Analysis Module - Uses pre-computed RBNS Z-scores
Implements Discrepancy Hunter algorithm from Dominguez et al. (2018)

This module handles:
1. Loading pre-computed RBNS k-mer Z-scores
2. Extending eCLIP peaks 50nt in 5' direction (crosslink site offset)
3. Filtering peaks by enrichment threshold
4. Extracting sequences and scoring with RBNS Z-scores
5. Classifying peaks as Canonical, Discrepant, or Intermediate

References:
- Dominguez et al. (2018) Molecular Cell
- RBNS_pipeline: https://github.com/cburgelab/RBNS_pipeline
"""

import pandas as pd
import pybedtools
import logging
from pathlib import Path
from typing import Dict, Tuple, Optional

# Configure logging
logger = logging.getLogger(__name__)

# Thresholds from Dominguez et al. and research document
CANONICAL_THRESHOLD = 3.0      # Z-score for significant RBNS motif
DISCREPANT_THRESHOLD = 1.5     # Below baseline enrichment
ECLIP_ENRICHMENT_MIN = 2.0     # narrowPeak signalValue column
PEAK_EXTENSION = 50            # 5' extension in nucleotides


def load_rbns_zscores(filepath: str) -> Dict[str, float]:
    """
    Load pre-computed Z-scores from RBNS pipeline output.

    The RBNS pipeline (Freese et al., Burge Lab) computes k-mer enrichment
    values and Z-scores. This function loads those pre-computed values.

    Parameters
    ----------
    filepath : str
        Path to CSV file with columns: kmer, z_score (and optionally r_value)

    Returns
    -------
    dict
        Dictionary mapping k-mers to Z-scores {kmer: z_score}

    Notes
    -----
    - Converts DNA (T) to RNA (U) automatically
    - Expected format: kmer,z_score,r_value (r_value optional)
    """
    logger.info(f"Loading RBNS Z-scores from {filepath}")

    df = pd.read_csv(filepath)

    # Handle different column naming conventions
    kmer_col = None
    zscore_col = None

    for col in df.columns:
        col_lower = col.lower()
        if 'kmer' in col_lower or 'sequence' in col_lower:
            kmer_col = col
        if 'z_score' in col_lower or 'zscore' in col_lower or col_lower == 'z':
            zscore_col = col

    if kmer_col is None or zscore_col is None:
        raise ValueError(
            f"Could not identify kmer and z_score columns. "
            f"Found columns: {list(df.columns)}"
        )

    # Convert DNA to RNA (T -> U)
    df['kmer_rna'] = df[kmer_col].str.upper().str.replace('T', 'U')

    result = dict(zip(df['kmer_rna'], df[zscore_col]))
    logger.info(f"Loaded {len(result)} k-mer Z-scores")

    return result


def extend_peaks_5prime(
    peaks_bed: str,
    chrom_sizes: str,
    extension: int = PEAK_EXTENSION
) -> pybedtools.BedTool:
    """
    Extend peaks in 5' direction per Dominguez methodology.

    The 5' start of eCLIP peaks corresponds to the UV crosslink site.
    Extending 50nt in the 5' direction captures the relevant binding context.

    Parameters
    ----------
    peaks_bed : str
        Path to eCLIP narrowPeak BED file
    chrom_sizes : str
        Path to chromosome sizes file (for boundary handling)
    extension : int
        Number of nucleotides to extend in 5' direction (default: 50)

    Returns
    -------
    pybedtools.BedTool
        Extended peaks (handles chromosome boundaries automatically)

    Notes
    -----
    Uses bedtools slop with:
    - l=extension (left/5' extension)
    - r=0 (no 3' extension)
    - s=True (strand-aware: 5' depends on strand)
    - g=chrom_sizes (genome file for boundary handling)
    """
    logger.info(f"Extending peaks {extension}nt in 5' direction")

    peaks = pybedtools.BedTool(peaks_bed)
    extended = peaks.slop(l=extension, r=0, s=True, g=chrom_sizes)

    logger.info(f"Extended {len(peaks)} peaks")
    return extended


def filter_by_enrichment(
    peaks: pybedtools.BedTool,
    min_enrichment: float = ECLIP_ENRICHMENT_MIN
) -> pybedtools.BedTool:
    """
    Filter narrowPeak by signalValue (fold enrichment over input).

    ENCODE narrowPeak format column 7 (0-indexed as 6) contains the
    signalValue, which represents fold enrichment over input control.

    Parameters
    ----------
    peaks : pybedtools.BedTool
        eCLIP peaks in narrowPeak format
    min_enrichment : float
        Minimum signalValue threshold (default: 2.0)

    Returns
    -------
    pybedtools.BedTool
        Filtered peaks
    """
    logger.info(f"Filtering peaks by enrichment >= {min_enrichment}")

    initial_count = len(peaks)
    filtered = peaks.filter(lambda x: float(x[6]) >= min_enrichment).saveas()
    final_count = len(filtered)

    logger.info(f"Retained {final_count}/{initial_count} peaks after filtering")
    return filtered


def extract_sequences(
    peaks_bed: str,
    genome_fasta: str
) -> Dict[str, str]:
    """
    Extract sequences for peaks with strand awareness.

    Uses bedtools getfasta via pybedtools for efficient extraction.
    Automatically handles reverse complementation for minus-strand peaks.

    Parameters
    ----------
    peaks_bed : str
        Path to BED file with peaks
    genome_fasta : str
        Path to reference genome FASTA file

    Returns
    -------
    dict
        Dictionary mapping peak IDs to RNA sequences {peak_id: sequence}

    Notes
    -----
    - Uses s=True for strand-aware extraction
    - Converts DNA (T) to RNA (U)
    - Peak IDs are in format: chrom:start-end(strand)
    """
    logger.info(f"Extracting sequences from {genome_fasta}")

    peaks = pybedtools.BedTool(peaks_bed)

    # Use bedtools getfasta with strand awareness
    seqs = peaks.sequence(fi=genome_fasta, s=True, tab=True, name=True)

    sequences = {}
    with open(seqs.seqfn) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                name = parts[0]
                seq = parts[1].upper().replace('T', 'U')
                sequences[name] = seq

    logger.info(f"Extracted {len(sequences)} sequences")
    return sequences


def calculate_peak_scores(
    sequence: str,
    rbns_dict: Dict[str, float],
    k: int = 5
) -> Tuple[float, float]:
    """
    Calculate Score_max and Score_sum per Dominguez methodology.

    Two scoring metrics are computed:
    - Score_max: Maximum Z-score among all k-mers (for single high-affinity sites)
    - Score_sum: Sum of all Z-scores (for avidity/multivalent binding)

    Parameters
    ----------
    sequence : str
        RNA sequence (should contain U, not T)
    rbns_dict : dict
        Dictionary mapping k-mers to Z-scores
    k : int
        K-mer size (default: 5)

    Returns
    -------
    tuple
        (score_max, score_sum)

    Notes
    -----
    - Returns (float('-inf'), 0.0) for sequences shorter than k
    - Unknown k-mers are scored as 0.0
    """
    if len(sequence) < k:
        return float('-inf'), 0.0

    scores = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        scores.append(rbns_dict.get(kmer, 0.0))

    if not scores:
        return float('-inf'), 0.0

    return max(scores), sum(scores)


def classify_peaks(
    scored_df: pd.DataFrame,
    canonical_threshold: float = CANONICAL_THRESHOLD,
    discrepant_threshold: float = DISCREPANT_THRESHOLD
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Classify peaks into Canonical, Discrepant, and Intermediate.

    Classification per research document Section 3.1:
    - Canonical: High eCLIP AND high RBNS (Score_max >= 3.0)
      -> Standard binding events (explained by sequence affinity)
    - Discrepant: High eCLIP BUT low RBNS (Score_max < 1.5)
      -> Hypothesized modification-enabled binding
    - Intermediate: Everything else (1.5 <= Score_max < 3.0)

    Parameters
    ----------
    scored_df : pd.DataFrame
        DataFrame with 'score_max' column
    canonical_threshold : float
        Threshold for canonical classification (default: 3.0)
    discrepant_threshold : float
        Threshold for discrepant classification (default: 1.5)

    Returns
    -------
    tuple
        (canonical_df, discrepant_df, intermediate_df)
    """
    logger.info("Classifying peaks...")

    canonical = scored_df[scored_df['score_max'] >= canonical_threshold].copy()
    discrepant = scored_df[scored_df['score_max'] < discrepant_threshold].copy()
    intermediate = scored_df[
        (scored_df['score_max'] >= discrepant_threshold) &
        (scored_df['score_max'] < canonical_threshold)
    ].copy()

    logger.info(
        f"Classification: Canonical={len(canonical)}, "
        f"Discrepant={len(discrepant)}, Intermediate={len(intermediate)}"
    )

    return canonical, discrepant, intermediate


def scored_df_to_bed(
    scored_df: pd.DataFrame,
    output_path: str
) -> None:
    """
    Convert scored DataFrame back to BED format.

    Parses peak_id (format: chrom:start-end) back to BED columns.

    Parameters
    ----------
    scored_df : pd.DataFrame
        DataFrame with 'peak_id' column in format chrom:start-end
    output_path : str
        Path to output BED file
    """
    bed_lines = []
    for _, row in scored_df.iterrows():
        peak_id = row['peak_id']
        # Parse peak_id format: chrom:start-end or chrom:start-end(strand)
        if '(' in peak_id:
            coords, strand = peak_id.rsplit('(', 1)
            strand = strand.rstrip(')')
        else:
            coords = peak_id
            strand = '.'

        # Use rsplit to handle chromosome names with colons (e.g., chrUn_gl000220)
        # Split from the right to get the last colon which separates chrom from positions
        chrom, positions = coords.rsplit(':', 1)
        start, end = positions.split('-')

        # BED format: chrom, start, end, name, score, strand
        score = row.get('score_max', 0)
        bed_lines.append(f"{chrom}\t{start}\t{end}\t{peak_id}\t{score}\t{strand}")

    with open(output_path, 'w') as f:
        f.write('\n'.join(bed_lines))

    logger.info(f"Wrote {len(bed_lines)} peaks to {output_path}")


def run_peak_analysis(
    rbns_path: str,
    eclip_path: str,
    genome_path: str,
    chrom_sizes_path: str,
    output_dir: str,
    extension: int = PEAK_EXTENSION,
    min_enrichment: float = ECLIP_ENRICHMENT_MIN
) -> pd.DataFrame:
    """
    Run complete peak analysis pipeline.

    This is a convenience function that runs all steps:
    1. Load RBNS Z-scores
    2. Extend eCLIP peaks
    3. Filter by enrichment
    4. Extract sequences
    5. Score peaks
    6. Classify peaks

    Parameters
    ----------
    rbns_path : str
        Path to RBNS Z-scores CSV
    eclip_path : str
        Path to eCLIP narrowPeak BED
    genome_path : str
        Path to reference genome FASTA
    chrom_sizes_path : str
        Path to chromosome sizes file
    output_dir : str
        Output directory for intermediate and result files
    extension : int
        5' extension for peaks (default: 50)
    min_enrichment : float
        Minimum enrichment threshold (default: 2.0)

    Returns
    -------
    pd.DataFrame
        Scored peaks DataFrame with classifications
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Load RBNS Z-scores
    rbns_dict = load_rbns_zscores(rbns_path)

    # Step 2: Extend peaks
    peaks = extend_peaks_5prime(eclip_path, chrom_sizes_path, extension)

    # Step 3: Filter by enrichment
    peaks = filter_by_enrichment(peaks, min_enrichment)
    filtered_bed = str(output_dir / 'peaks_filtered.bed')
    peaks.saveas(filtered_bed)

    # Step 4: Extract sequences
    sequences = extract_sequences(filtered_bed, genome_path)

    # Step 5: Score peaks
    scored_data = []
    for peak_id, seq in sequences.items():
        score_max, score_sum = calculate_peak_scores(seq, rbns_dict)
        scored_data.append({
            'peak_id': peak_id,
            'sequence': seq,
            'score_max': score_max,
            'score_sum': score_sum
        })

    scored_df = pd.DataFrame(scored_data)
    scored_df.to_csv(output_dir / 'scored_peaks.csv', index=False)

    # Step 6: Classify peaks
    canonical, discrepant, intermediate = classify_peaks(scored_df)

    # Save classified peaks as BED files
    if len(canonical) > 0:
        scored_df_to_bed(canonical, str(output_dir / 'canonical_peaks.bed'))
    if len(discrepant) > 0:
        scored_df_to_bed(discrepant, str(output_dir / 'discrepant_peaks.bed'))
    if len(intermediate) > 0:
        scored_df_to_bed(intermediate, str(output_dir / 'intermediate_peaks.bed'))

    # Add category column
    scored_df['category'] = 'intermediate'
    scored_df.loc[scored_df['score_max'] >= CANONICAL_THRESHOLD, 'category'] = 'canonical'
    scored_df.loc[scored_df['score_max'] < DISCREPANT_THRESHOLD, 'category'] = 'discrepant'

    return scored_df


if __name__ == '__main__':
    # Basic test/example usage
    import sys
    logging.basicConfig(level=logging.INFO)

    if len(sys.argv) < 5:
        print("Usage: python peak_analysis.py <rbns.csv> <eclip.bed> <genome.fa> <chrom.sizes>")
        sys.exit(1)

    rbns_path, eclip_path, genome_path, chrom_sizes = sys.argv[1:5]

    result = run_peak_analysis(
        rbns_path, eclip_path, genome_path, chrom_sizes,
        output_dir='./test_output'
    )
    print(f"\nResults summary:")
    print(result['category'].value_counts())
