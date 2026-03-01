#!/usr/bin/env python3
"""
Generate AlphaFold3 Server JSON input files for RBM22-pseudoU validation.

Creates 30 JSON files (up to 15 discrepant peaks × 2 conditions: unmodified /
pseudoU-modified) for upload to https://alphafoldserver.com.

Strategy:
  1. Load RBM22 discrepant peaks and their pre-extracted sequences
     (from results/RBM22/scored_peaks.csv — no genome FASTA required).
  2. Prioritise peaks that overlap pseudoU sites (from data/mods/K562/pseudoU.bed).
     Fill remainder with the most-discrepant non-pseudoU peaks if < 15 overlap.
  3. For each selected peak, locate the pseudoU site within the peak sequence and
     extract a 30 nt window centred on it. If a peak has no pseudoU overlap it is
     used unmodified-only (single job, not a pair).
  4. Load the full-length RBM22 protein sequence from data/proteins/RBM22.fasta.
  5. Write two JSON files per peak: one with pseudoU (modificationType = "PSU"),
     one without any modification.
  6. Write a manifest CSV and run basic validation.

AF3 Server JSON format:
  dialect       : "alphafoldserver"
  version       : 2
  sequences[].proteinChain.sequence
  sequences[].rnaSequence.sequence
  sequences[].rnaSequence.modifications[].modificationType   (CCD code "PSU")
  sequences[].rnaSequence.modifications[].basePosition       (1-indexed)

Usage:
  cd /path/to/RNAMod-RBNS
  python scripts/generate_af3_rbm22.py [--output results/RBM22/af3_inputs]
                                        [--n-peaks 15]
                                        [--seeds 1 2 3 4 5]
                                        [--window 15]
"""

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import pybedtools


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MOD_CCD  = "CCD_PSU"   # CCD code for pseudouridine (must include CCD_ prefix)
HALF_WIN = 15      # nucleotides on each side → 31 nt total centred on Ψ


# ---------------------------------------------------------------------------
# Step 1: Load RBM22 protein sequence from local FASTA
# ---------------------------------------------------------------------------

def load_fasta_sequence(fasta_path: Path) -> Tuple[str, str]:
    """
    Parse a single-entry FASTA file and return (name, sequence).
    Raises FileNotFoundError / RuntimeError on missing or malformed input.
    """
    if not fasta_path.exists():
        raise FileNotFoundError(f"Protein FASTA not found: {fasta_path}")
    text = fasta_path.read_text(encoding="utf-8")
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    if not lines or not lines[0].startswith(">"):
        raise RuntimeError(f"Invalid FASTA format in {fasta_path}")
    header   = lines[0]                        # >sp|Q9NW64|RBM22_HUMAN ...
    sequence = "".join(lines[1:])
    match    = re.search(r"\|[A-Z0-9]+\|(\S+)", header)
    name     = match.group(1) if match else fasta_path.stem
    print(f"Loaded protein from {fasta_path.name}: {name} ({len(sequence)} aa)")
    return name, sequence


# ---------------------------------------------------------------------------
# Step 2: Load discrepant peaks and select pseudoU-overlapping ones
# ---------------------------------------------------------------------------

def load_discrepant_peaks(results_dir: Path) -> pd.DataFrame:
    """
    Load scored_peaks.csv and return discrepant rows sorted by score_max (asc).
    The 'sequence' column already contains the 5′-extended eCLIP peak sequence.
    """
    scored_csv = results_dir / "scored_peaks.csv"
    if not scored_csv.exists():
        raise FileNotFoundError(f"Missing: {scored_csv}")
    df = pd.read_csv(scored_csv)
    discrepant = (
        df[df["category"] == "discrepant"]
        .sort_values("score_max")
        .reset_index(drop=True)
    )
    print(f"Discrepant peaks: {len(discrepant)}")
    return discrepant


def parse_peak_coords(peak_id: str) -> Tuple[str, int, int, str]:
    """
    Parse chrom, start, end, strand from a peak_id like:
      RBM22_K562_IDR::chr19:50804490-50804540(-)
    Returns (chrom, start, end, strand).
    """
    m = re.search(r"(chr[^:()]+):(\d+)-(\d+)\(([+-])\)", peak_id)
    if not m:
        raise ValueError(f"Cannot parse coords from peak_id: {peak_id!r}")
    return m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)


def select_peaks_with_pseudou(
    discrepant: pd.DataFrame,
    pseudou_bed_path: str,
    discrepant_bed_path: str,
    n_peaks: int = 15,
) -> pd.DataFrame:
    """
    Return up to n_peaks discrepant rows, prioritising those that overlap
    pseudoU sites.  Adds columns: has_pseudou (bool), pseudou_genome_pos (int|None).
    """
    pseudou = pybedtools.BedTool(pseudou_bed_path)
    disc_bed = pybedtools.BedTool(discrepant_bed_path)

    # Intersect: get discrepant peaks that overlap a pseudoU site
    # -wa: write original A record; -wb: write B record (pseudoU site)
    overlap = disc_bed.intersect(pseudou, wa=True, wb=True)

    # Build map: coord_key → first pseudoU genome position
    _coord_re = re.compile(r"(chr[^:()]+):(\d+-\d+)")

    def _coord_key(peak_id: str) -> str:
        pid = peak_id.split("(")[0]
        m = _coord_re.search(pid)
        return f"{m.group(1)}:{m.group(2)}" if m else pid

    pseudou_pos_map: Dict[str, int] = {}
    for feat in overlap:
        # feat.fields[3] is peak_id (from 6-col BED written by scored_df_to_bed)
        # feat.fields[6] is pseudoU chrom, [7] start, [8] end
        pid = feat.fields[3]
        key = _coord_key(pid)
        if key not in pseudou_pos_map:
            pseudou_pos_map[key] = int(feat.fields[7])  # 0-based start of pseudoU

    print(f"Discrepant peaks overlapping pseudoU: {len(pseudou_pos_map)}")

    # Annotate discrepant DataFrame
    discrepant = discrepant.copy()
    discrepant["_coord_key"] = discrepant["peak_id"].apply(_coord_key)
    discrepant["has_pseudou"] = discrepant["_coord_key"].isin(pseudou_pos_map)
    discrepant["pseudou_genome_pos"] = discrepant["_coord_key"].map(pseudou_pos_map)

    # Prioritise pseudoU-overlapping peaks
    has_psu = discrepant[discrepant["has_pseudou"]].head(n_peaks)
    remaining = n_peaks - len(has_psu)
    no_psu = discrepant[~discrepant["has_pseudou"]].head(remaining)

    selected = pd.concat([has_psu, no_psu]).reset_index(drop=True)
    print(
        f"Selected {len(selected)} peaks "
        f"({len(has_psu)} with pseudoU, {len(no_psu)} without)"
    )
    return selected


# ---------------------------------------------------------------------------
# Step 3: Locate pseudoU in the peak sequence and extract 30 nt window
# ---------------------------------------------------------------------------

RC_MAP = str.maketrans("ACGUacgu", "UGCAugca")


def reverse_complement_rna(seq: str) -> str:
    return seq.translate(RC_MAP)[::-1]


def extract_rna_window(
    row: pd.Series,
    half_window: int = HALF_WIN,
) -> Optional[Dict]:
    """
    Given a scored_peaks row with has_pseudou=True, return a dict containing:
      rna_sequence   : 30 nt (or shorter at chromosome edges) RNA string
      pseudou_rna_pos: 1-indexed position of Ψ within rna_sequence
      strand         : '+' or '-'
    Returns None if the pseudoU position cannot be mapped into the sequence.

    The peak 'sequence' column holds the 50 nt 5′-extended eCLIP peak sequence
    (already RNA, from scored_peaks.csv). We locate the pseudoU relative to the
    peak start and slice a ±half_window window from the raw peak sequence.
    """
    chrom, peak_start, peak_end, strand = parse_peak_coords(row["peak_id"])
    psu_genome = int(row["pseudou_genome_pos"])   # 0-based

    peak_seq = row["sequence"]    # already RNA (U not T) from extract_sequences()
    peak_len = len(peak_seq)

    if strand == "+":
        # peak_seq[0] corresponds to genome position peak_start
        psu_in_peak = psu_genome - peak_start      # 0-indexed offset in peak_seq
    else:
        # For minus strand the sequence was reverse-complemented during extraction.
        # peak_seq[0] corresponds to genome position peak_end - 1.
        psu_in_peak = peak_end - 1 - psu_genome   # 0-indexed offset in peak_seq

    if not (0 <= psu_in_peak < peak_len):
        return None   # pseudoU falls outside extracted peak sequence

    # Slice ±half_window around the pseudoU position
    win_start = max(0, psu_in_peak - half_window)
    win_end   = min(peak_len, psu_in_peak + half_window + 1)
    rna_window = peak_seq[win_start:win_end]

    # Position of pseudoU within the window (1-indexed for AF3)
    psu_in_window_0 = psu_in_peak - win_start
    psu_in_window_1 = psu_in_window_0 + 1

    # Validate: should be U at the pseudoU position
    base = rna_window[psu_in_window_0] if 0 <= psu_in_window_0 < len(rna_window) else "?"
    if base != "U":
        print(
            f"  WARNING: base at pseudoU site is '{base}' (expected U) "
            f"in {row['peak_id']}"
        )

    return {
        "rna_sequence":    rna_window,
        "rna_length":      len(rna_window),
        "pseudou_rna_pos": psu_in_window_1,
        "base_at_site":    base,
        "strand":          strand,
    }


# ---------------------------------------------------------------------------
# Step 4: Build AF3 server JSON
# ---------------------------------------------------------------------------

def build_af3_json(
    job_name: str,
    protein_seq: str,
    rna_seq: str,
    pseudou_pos: Optional[int],
    seeds: List[int],
) -> dict:
    """
    Return an AlphaFold Server-dialect JSON dict.
    pseudou_pos: 1-indexed position in rna_seq, or None for unmodified.

    AF3 Server dialect:
      top-level keys: name, modelSeeds, sequences, dialect, version
      protein entity : {"proteinChain": {"sequence": ..., "count": 1}}
      RNA entity     : {"rnaSequence": {"sequence": ..., "modifications": [...]}}
    """
    modifications = []
    if pseudou_pos is not None:
        modifications.append({
            "modificationType": MOD_CCD,   # "PSU"
            "basePosition": pseudou_pos
        })

    rna_entity: dict = {"sequence": rna_seq}
    if modifications:
        rna_entity["modifications"] = modifications

    return {
        "name": job_name,
        "modelSeeds": [str(s) for s in seeds],  # must be strings per AF3 spec
        "sequences": [
            {"proteinChain": {"sequence": protein_seq, "count": 1}},
            {"rnaSequence": rna_entity},
        ],
        "dialect": "alphafoldserver",
        "version": 1,
    }


# ---------------------------------------------------------------------------
# Step 5: Validate a single JSON dict
# ---------------------------------------------------------------------------

def validate_job(data: dict) -> List[str]:
    """Return a list of error strings (empty = valid)."""
    errors = []
    for field in ("name", "modelSeeds", "sequences", "dialect", "version"):
        if field not in data:
            errors.append(f"Missing field: {field}")
    if errors:
        return errors

    protein_seq = rna_seq = None
    modifications = []
    for entry in data["sequences"]:
        if "proteinChain" in entry:
            protein_seq = entry["proteinChain"].get("sequence", "")
        if "rnaSequence" in entry:
            rna_seq = entry["rnaSequence"].get("sequence", "")
            modifications = entry["rnaSequence"].get("modifications", [])

    if not protein_seq:
        errors.append("No protein sequence")
    else:
        invalid_aa = set(protein_seq) - set("ACDEFGHIKLMNPQRSTVWY")
        if invalid_aa:
            errors.append(f"Invalid amino acids: {invalid_aa}")

    if not rna_seq:
        errors.append("No RNA sequence")
    else:
        invalid_bases = set(rna_seq) - set("ACGU")
        if invalid_bases:
            errors.append(f"Invalid RNA bases: {invalid_bases}")
        for mod in modifications:
            pos = mod.get("basePosition")
            if pos is None:
                errors.append("Modification missing basePosition")
            elif not (1 <= pos <= len(rna_seq)):
                errors.append(f"Modification position {pos} out of range 1–{len(rna_seq)}")
            else:
                base = rna_seq[pos - 1]
                if base != "U":
                    errors.append(
                        f"PseudoU at position {pos} is '{base}' not 'U'"
                    )

    return errors


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate AlphaFold3 Server JSON files for RBM22-pseudoU validation"
    )
    parser.add_argument(
        "--output", default="results/RBM22/af3_inputs",
        help="Output directory for JSON files (default: results/RBM22/af3_inputs)"
    )
    parser.add_argument(
        "--n-peaks", type=int, default=15,
        help="Max number of discrepant peaks to model (default: 15)"
    )
    parser.add_argument(
        "--seeds", type=int, nargs="+", default=[1, 2, 3, 4, 5],
        help="Model seeds (default: 1 2 3 4 5)"
    )
    parser.add_argument(
        "--window", type=int, default=HALF_WIN,
        help=f"Nucleotides on each side of pseudoU (default: {HALF_WIN} → 31 nt total)"
    )
    parser.add_argument(
        "--results-dir", default="results/RBM22",
        help="RBM22 results directory (default: results/RBM22)"
    )
    parser.add_argument(
        "--mods-dir", default="data/mods/K562",
        help="Modification BED directory (default: data/mods/K562)"
    )
    parser.add_argument(
        "--protein-fasta", default="data/proteins/RBM22.fasta",
        help="Path to RBM22 protein FASTA (default: data/proteins/RBM22.fasta)"
    )
    args = parser.parse_args()

    results_dir  = Path(args.results_dir)
    mods_dir     = Path(args.mods_dir)
    output_dir   = Path(args.output)
    pseudou_bed  = str(mods_dir / "pseudoU.bed")
    disc_bed     = str(results_dir / "discrepant_peaks.bed")

    output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Load RBM22 protein sequence from local FASTA
    # ------------------------------------------------------------------
    protein_name, protein_seq = load_fasta_sequence(Path(args.protein_fasta))

    # ------------------------------------------------------------------
    # 2. Load and select discrepant peaks
    # ------------------------------------------------------------------
    discrepant = load_discrepant_peaks(results_dir)
    selected   = select_peaks_with_pseudou(
        discrepant, pseudou_bed, disc_bed, n_peaks=args.n_peaks
    )

    # ------------------------------------------------------------------
    # 3. Build all job dicts
    # ------------------------------------------------------------------
    all_jobs      = []   # accumulated job dicts for the single output JSON
    manifest_rows = []
    n_pairs  = 0
    n_single = 0
    n_errors = 0

    for i, row in selected.iterrows():
        peak_id  = row["peak_id"]
        # Safe job-name stem: replace special chars
        safe_id  = re.sub(r"[^A-Za-z0-9_\-]", "_", peak_id)
        safe_id  = re.sub(r"_+", "_", safe_id).strip("_")[:60]

        # Extract pseudoU window if this peak overlaps a Ψ site
        rna_info = None
        if row["has_pseudou"]:
            rna_info = extract_rna_window(row, half_window=args.window)
            if rna_info is None:
                print(f"  [{i+1}] SKIP: pseudoU outside peak sequence — {peak_id}")
                continue

        if rna_info is not None:
            rna_seq = rna_info["rna_sequence"]
            psu_pos = rna_info["pseudou_rna_pos"]
            base    = rna_info["base_at_site"]
        else:
            # No pseudoU overlap — use full peak sequence, no modification
            rna_seq = row["sequence"]
            psu_pos = None
            base    = "N/A"

        # --- Unmodified job ---
        name_unmod = f"RBM22_{safe_id}_unmod"
        job_unmod  = build_af3_json(name_unmod, protein_seq, rna_seq,
                                    pseudou_pos=None, seeds=args.seeds)
        errs_unmod = validate_job(job_unmod)
        all_jobs.append(job_unmod)

        # --- Modified job (only when pseudoU site is known) ---
        errs_mod = []
        name_mod = None
        if psu_pos is not None:
            name_mod = f"RBM22_{safe_id}_pseudoU"
            job_mod  = build_af3_json(name_mod, protein_seq, rna_seq,
                                      pseudou_pos=psu_pos, seeds=args.seeds)
            errs_mod = validate_job(job_mod)
            all_jobs.append(job_mod)
            n_pairs += 1
        else:
            n_single += 1

        # Report
        status     = "✓" if not errs_unmod and not errs_mod else "✗"
        pair_label = "(pair)" if psu_pos else "(unmod only)"
        print(f"  {status} [{i+1:02d}] {safe_id[:40]} {pair_label}")
        if errs_unmod:
            for e in errs_unmod:
                print(f"       UNMOD ERROR: {e}")
            n_errors += 1
        if errs_mod:
            for e in errs_mod:
                print(f"       MOD   ERROR: {e}")
            n_errors += 1

        manifest_rows.append({
            "peak_number":     i + 1,
            "peak_id":         peak_id,
            "has_pseudou":     row["has_pseudou"],
            "score_max":       row["score_max"],
            "rna_sequence":    rna_seq,
            "rna_length":      len(rna_seq),
            "pseudou_rna_pos": psu_pos,
            "base_at_site":    base,
            "protein_id":      protein_name,
            "protein_length":  len(protein_seq),
            "job_unmodified":  name_unmod,
            "job_modified":    name_mod if name_mod else "",
        })

    # ------------------------------------------------------------------
    # 4. Write single combined JSON (list of job dicts) + manifest
    # ------------------------------------------------------------------
    combined_path = output_dir / "RBM22_pseudoU_jobs.json"
    with open(combined_path, "w") as fh:
        json.dump(all_jobs, fh, indent=2)

    manifest_df   = pd.DataFrame(manifest_rows)
    manifest_path = output_dir / "manifest.csv"
    manifest_df.to_csv(manifest_path, index=False)

    total_jobs = n_pairs * 2 + n_single
    print(f"\n{'='*60}")
    print(f"Done.")
    print(f"  Peaks processed : {len(manifest_rows)}")
    print(f"  Paired jobs     : {n_pairs} pairs")
    print(f"  Unmod-only jobs : {n_single}")
    print(f"  Total jobs      : {total_jobs}")
    print(f"  Validation errors: {n_errors}")
    print(f"  Output JSON     : {combined_path}")
    print(f"  Manifest        : {manifest_path}")
    print(f"{'='*60}")

    if total_jobs > 30:
        print(
            f"\nWARNING: {total_jobs} jobs exceeds the AF3 server daily limit of 30."
            f"\nReduce --n-peaks (currently {args.n_peaks}) to stay within the limit."
        )

    if n_errors > 0:
        print(f"\n{n_errors} validation error(s) found — review output above before uploading.")
        sys.exit(1)
    else:
        print("\nAll JSONs valid. Ready to upload to https://alphafoldserver.com")


if __name__ == "__main__":
    main()
