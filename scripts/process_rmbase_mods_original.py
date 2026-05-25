#!/usr/bin/env python3
"""
Process RMBase v3.0 Modification Data for RNAMod-RBNS Pipeline

Extracts and filters RNA modification sites from RMBase v3.0 files,
converting them to BED6 format for use with bedtools.

RMBase v3.0 Format (29 columns for human):
    1: Chrom, 2: Start, 3: End, 4: modID, 5: Score, 6: Strand,
    7: modType, 8: supportNum, 9-11: support/pub lists,
    12: cellList, 13: seqTypeList, 14-18: gene info, 19: Seq,
    20: motifScore, 21-29: additional annotations

Output BED6 format:
    chrom  start  end  modID  score  strand

Usage:
    # Process all modifications for a cell line
    python scripts/process_rmbase_mods.py --cell-line HepG2

    # Process specific modification
    python scripts/process_rmbase_mods.py --mod-type m6A --cell-line HepG2

    # Use all sites (no cell line filter)
    python scripts/process_rmbase_mods.py --mod-type Pseudo --all-sites

    # Process all modifications for both cell lines
    python scripts/process_rmbase_mods.py --setup-all
"""

import argparse
import gzip
import tarfile
import os
import sys
from pathlib import Path
from collections import defaultdict


# Mapping of modification types to filenames
MOD_FILES = {
    'm6A': 'human.hg38.m6A.result.col29.bed',
    'Pseudo': 'human.hg38.Pseudo.result.col29.bed',  # Pseudouridine (Ψ)
    'm5C': 'human.hg38.m5C.result.col29.bed',
    'ac4C': 'human.hg38.ac4C.result.col29.re.bed',
    'm1A': 'human.hg38.m1A.result.col29.bed',
    'm7G': 'human.hg38.m7G.result.col29.bed',
    'Nm': 'human.hg38.Nm.result.col29.bed',
}

# Archive files
MOD_ARCHIVES = {
    'm6A': 'hg38.m6A.tar.gz',
    'Pseudo': 'hg38.Pseudo.tar.gz',
    'm5C': 'hg38.m5C.tar.gz',
    'ac4C': 'hg38.ac4C.tar.gz',
    'm1A': 'hg38.m1A.tar.gz',
    'm7G': 'hg38.m7G.tar.gz',
    'Nm': 'hg38.Nm.tar.gz',
}

# Output filename mapping
OUTPUT_NAMES = {
    'm6A': 'm6A.bed',
    'Pseudo': 'pseudoU.bed',
    'm5C': 'm5C.bed',
    'ac4C': 'ac4C.bed',
}


def extract_archive(archive_path, output_dir):
    """Extract tar.gz archive if not already extracted."""
    archive_path = Path(archive_path)
    if not archive_path.exists():
        return None

    with tarfile.open(archive_path, 'r:gz') as tar:
        members = tar.getnames()
        for member in members:
            if member.endswith('.bed'):
                target = output_dir / member
                if not target.exists():
                    tar.extract(member, output_dir)
                return output_dir / member
    return None


def process_rmbase_file(input_path, output_path, cell_line=None, min_support=1,
                        min_motif_score=None):
    """
    Process RMBase file and extract sites.

    Args:
        input_path: Path to RMBase BED file
        output_path: Path to output BED6 file
        cell_line: Filter for specific cell line (None = all sites)
        min_support: Minimum number of supporting datasets
        min_motif_score: Minimum motif score (for m6A, 0-5)

    Returns:
        dict with statistics
    """
    input_path = Path(input_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    stats = {
        'total_lines': 0,
        'kept_lines': 0,
        'cell_line_matches': 0,
        'support_filter': 0,
        'motif_filter': 0,
    }

    print(f"Processing: {input_path.name}")
    print(f"  Cell line filter: {cell_line or 'None (all sites)'}")
    print(f"  Min support: {min_support}")
    print(f"  Min motif score: {min_motif_score or 'None'}")

    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            stats['total_lines'] += 1

            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue

            chrom = fields[0]
            start = fields[1]
            end = fields[2]
            mod_id = fields[3]
            score = fields[4]
            strand = fields[5]
            support_num = int(fields[7]) if fields[7].isdigit() else 1
            cell_list = fields[11]

            # Get motif score if available (column 20, index 19)
            motif_score = None
            if len(fields) > 19 and fields[19] != 'na':
                try:
                    motif_score = float(fields[19])
                except ValueError:
                    pass

            # Apply filters
            # Support filter
            if support_num < min_support:
                stats['support_filter'] += 1
                continue

            # Motif score filter
            if min_motif_score is not None and motif_score is not None:
                if motif_score < min_motif_score:
                    stats['motif_filter'] += 1
                    continue

            # Cell line filter
            if cell_line:
                # Check if cell line appears in the comma-separated list
                cell_types = [c.strip() for c in cell_list.split(',')]
                if cell_line not in cell_types:
                    continue
                stats['cell_line_matches'] += 1

            # Write BED6 format
            # Use motif score as BED score if available, otherwise support count
            bed_score = int(motif_score * 200) if motif_score else min(support_num * 100, 1000)
            outfile.write(f"{chrom}\t{start}\t{end}\t{mod_id}\t{bed_score}\t{strand}\n")
            stats['kept_lines'] += 1

    print(f"  Total lines: {stats['total_lines']:,}")
    print(f"  Kept lines: {stats['kept_lines']:,}")
    if cell_line:
        print(f"  Cell line matches: {stats['cell_line_matches']:,}")
    print(f"  Filtered (support): {stats['support_filter']:,}")
    print(f"  Filtered (motif): {stats['motif_filter']:,}")
    print(f"  Output: {output_path}")

    return stats


def setup_all_modifications(mods_dir, cell_lines=None):
    """
    Set up modification files for all cell lines.

    Args:
        mods_dir: Base directory containing RMBase files
        cell_lines: List of cell lines to process (default: ['K562', 'HepG2'])
    """
    mods_dir = Path(mods_dir)
    cell_lines = cell_lines or ['K562', 'HepG2']

    # Primary modifications for analysis
    primary_mods = ['m6A', 'Pseudo', 'm5C', 'ac4C']

    print("=" * 60)
    print("Setting up RMBase modification data")
    print("=" * 60)
    print(f"Modifications: {', '.join(primary_mods)}")
    print(f"Cell lines: {', '.join(cell_lines)}")
    print()

    # First, extract all archives
    print("Extracting archives...")
    for mod_type in primary_mods:
        archive = mods_dir / MOD_ARCHIVES.get(mod_type, '')
        if archive.exists():
            extract_archive(archive, mods_dir)
            print(f"  Extracted: {archive.name}")
    print()

    # Process for each cell line
    for cell_line in cell_lines:
        print(f"\n{'='*60}")
        print(f"Processing for {cell_line}")
        print('='*60)

        output_dir = mods_dir / cell_line
        output_dir.mkdir(exist_ok=True)

        for mod_type in primary_mods:
            input_file = mods_dir / MOD_FILES.get(mod_type, '')
            output_file = output_dir / OUTPUT_NAMES.get(mod_type, f'{mod_type}.bed')

            if not input_file.exists():
                print(f"\n  WARNING: {input_file.name} not found, skipping {mod_type}")
                continue

            print(f"\n--- {mod_type} ---")

            # For m6A, filter by cell line if available
            # For others (Ψ, m5C, ac4C), use all sites (no cell-specific data)
            if mod_type == 'm6A':
                # Check if this cell line has data
                with open(input_file, 'r') as f:
                    sample_line = f.readline()
                    has_cell_data = cell_line in sample_line

                # Read more lines to check
                cell_line_count = 0
                with open(input_file, 'r') as f:
                    for i, line in enumerate(f):
                        if cell_line in line:
                            cell_line_count += 1
                        if i > 10000:
                            break

                if cell_line_count > 0:
                    # Filter for this cell line
                    process_rmbase_file(input_file, output_file,
                                       cell_line=cell_line, min_support=1)
                else:
                    # No cell-specific data, use all sites
                    print(f"  NOTE: No {cell_line}-specific {mod_type} data found")
                    print(f"  Using all {mod_type} sites (may include other cell types)")
                    process_rmbase_file(input_file, output_file,
                                       cell_line=None, min_support=2)
            else:
                # For Ψ, m5C, ac4C - no cell-specific data available
                # Use all sites with support >= 1
                process_rmbase_file(input_file, output_file,
                                   cell_line=None, min_support=1)

    print("\n" + "=" * 60)
    print("Setup complete!")
    print("=" * 60)

    # Verify outputs
    print("\nVerifying outputs:")
    for cell_line in cell_lines:
        output_dir = mods_dir / cell_line
        print(f"\n{cell_line}/")
        for mod_type in primary_mods:
            output_file = output_dir / OUTPUT_NAMES.get(mod_type, f'{mod_type}.bed')
            if output_file.exists():
                with open(output_file, 'r') as f:
                    line_count = sum(1 for _ in f)
                print(f"  {output_file.name}: {line_count:,} sites")
            else:
                print(f"  {output_file.name}: NOT CREATED")


def main():
    parser = argparse.ArgumentParser(
        description='Process RMBase v3.0 modification data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Set up all modifications for K562 and HepG2
    python scripts/process_rmbase_mods.py --setup-all

    # Process single modification for specific cell line
    python scripts/process_rmbase_mods.py --mod-type m6A --cell-line HepG2

    # Process using all sites (no cell filter)
    python scripts/process_rmbase_mods.py --mod-type Pseudo --all-sites

    # Filter with minimum support
    python scripts/process_rmbase_mods.py --mod-type m6A --min-support 2
        """
    )

    parser.add_argument('--mods-dir', default='data/mods',
                        help='Directory containing RMBase files (default: data/mods)')
    parser.add_argument('--setup-all', action='store_true',
                        help='Set up all modifications for K562 and HepG2')
    parser.add_argument('--mod-type', choices=list(MOD_FILES.keys()),
                        help='Specific modification type to process')
    parser.add_argument('--cell-line',
                        help='Filter for specific cell line (e.g., HepG2, K562)')
    parser.add_argument('--all-sites', action='store_true',
                        help='Use all sites without cell line filtering')
    parser.add_argument('--output', '-o',
                        help='Output BED file path')
    parser.add_argument('--min-support', type=int, default=1,
                        help='Minimum number of supporting datasets (default: 1)')
    parser.add_argument('--min-motif-score', type=float,
                        help='Minimum motif score for m6A (0-5)')

    args = parser.parse_args()

    mods_dir = Path(args.mods_dir)

    if args.setup_all:
        setup_all_modifications(mods_dir)
        return

    if not args.mod_type:
        print("Error: --mod-type is required unless using --setup-all")
        sys.exit(1)

    # Find input file
    input_file = mods_dir / MOD_FILES[args.mod_type]

    # Extract archive if needed
    if not input_file.exists():
        archive = mods_dir / MOD_ARCHIVES.get(args.mod_type, '')
        if archive.exists():
            extracted = extract_archive(archive, mods_dir)
            if extracted:
                input_file = extracted

    if not input_file.exists():
        print(f"Error: Input file not found: {input_file}")
        sys.exit(1)

    # Determine output path
    if args.output:
        output_file = Path(args.output)
    elif args.cell_line:
        output_dir = mods_dir / args.cell_line
        output_file = output_dir / OUTPUT_NAMES.get(args.mod_type, f'{args.mod_type}.bed')
    else:
        output_file = mods_dir / f'{args.mod_type}_all.bed'

    # Process
    cell_line = None if args.all_sites else args.cell_line

    process_rmbase_file(
        input_file,
        output_file,
        cell_line=cell_line,
        min_support=args.min_support,
        min_motif_score=args.min_motif_score
    )


if __name__ == '__main__':
    main()
