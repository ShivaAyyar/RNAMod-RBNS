# Pipeline Debugging: Zero Modification Overlap

*Source: External review of the GitHub repository, Feb 18 2026*

The enrichment analysis shows no significant overlap between discrepant peaks and RNA modification sites, including in the positive controls IGF2BP1 and IGF2BP2. Below are the suggested causes and diagnostics, evaluated against the actual pipeline code.

---

## Suggested Issues

### 1. Coordinate System Mismatch

**Hypothesis:** eCLIP peaks and RMBase modification sites may use different coordinate systems (0-based vs 1-based), causing intersections to always miss.

**Diagnostic suggested:**
```python
def validate_coordinate_overlap(peaks_bed, mod_bed):
    peaks = pybedtools.BedTool(peaks_bed)
    mods = pybedtools.BedTool(mod_bed)

    exact = len(peaks.intersect(mods, u=True))

    peaks_window = peaks.slop(b=100, g=chrom_sizes)
    window = len(peaks_window.intersect(mods, u=True))

    print(f"Exact overlap: {exact}/{len(peaks)}")
    print(f"100bp window: {window}/{len(peaks)}")
    # If window >> exact: coordinates are correct but intersection is too stringent
    # If both ~0: coordinate system mismatch
```

**Bash checks:**
```bash
# Check chromosome naming consistency (chr1 vs 1)
head -5 data/eclip/IGF2BP1_K562.bed | cut -f1
head -5 data/mods/K562/m6A.bed | cut -f1

# Check for invalid coordinates (start >= end)
awk '$2 >= $3' data/mods/K562/m6A.bed

# Check strand distribution
cut -f6 data/eclip/IGF2BP1_K562.bed | sort | uniq -c
cut -f6 data/mods/K562/m6A.bed | sort | uniq -c
```

---

### 2. Peak Extension Logic May Be Incorrect

**Hypothesis:** The current 50nt 5' extension extends from the **peak start**, not the **crosslink summit**. The eCLIP narrowPeak column 10 encodes the summit offset from the start — the actual crosslink site is `start + summit_offset`.

**Current code (`peak_analysis.py`):**
```python
def extend_peaks_5prime(peaks_bed, chrom_sizes, extension=50):
    peaks = pybedtools.BedTool(peaks_bed)
    extended = peaks.slop(l=extension, r=0, s=True, g=chrom_sizes)
    return extended
```

**Suggested fix — extend from crosslink summit instead:**
```python
def extend_from_crosslink_site(peaks_bed, chrom_sizes, extension=50):
    """
    Extend from the crosslink site (peak summit), not peak start.
    eCLIP narrowPeak column 10 = summit offset from start.
    Actual crosslink = start + summit_offset
    """
    peaks = []
    with open(peaks_bed) as f:
        for line in f:
            fields = line.strip().split('\t')
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            strand = fields[5] if len(fields) > 5 else '.'

            summit_offset = int(fields[9]) if len(fields) > 9 else (end - start) // 2
            crosslink = start + summit_offset

            if strand == '-':
                new_start = crosslink
                new_end = min(crosslink + extension, chrom_sizes[chrom])
            else:
                new_start = max(0, crosslink - extension)
                new_end = crosslink

            peaks.append(f"{chrom}\t{new_start}\t{new_end}\t{fields[3]}\t{fields[4]}\t{strand}")

    return pybedtools.BedTool('\n'.join(peaks), from_string=True)
```

---

### 3. Modification Site Sparsity

**Hypothesis:** RMBase modification coverage may be insufficient in the genomic regions where eCLIP peaks occur (e.g., 3' UTRs for IGF2BP1/2).

**Site counts from current data:**

| Modification | Sites | Notes |
|---|---|---|
| m6A (K562) | 558,360 | Cumulative multi-study |
| m6A (HepG2) | 178,800 | HepG2-specific |
| Ψ | 5,705 | May be too sparse |
| m5C | 46,025 | |
| ac4C | 1,861 | Very sparse vs 5,885 TIA1 peaks |

**Diagnostic suggested:**
```python
def check_modification_coverage(eclip_bed, mod_bed):
    utr_regions = pybedtools.BedTool('gencode.v38.3UTR.bed')
    mods = pybedtools.BedTool(mod_bed)
    mods_in_utrs = mods.intersect(utr_regions, u=True)
    print(f"Total mods: {len(mods)}")
    print(f"Mods in 3'UTRs: {len(mods_in_utrs)} ({100*len(mods_in_utrs)/len(mods):.1f}%)")
    peaks = pybedtools.BedTool(eclip_bed)
    peaks_in_utrs = peaks.intersect(utr_regions, u=True)
    print(f"Total peaks: {len(peaks)}")
    print(f"Peaks in 3'UTRs: {len(peaks_in_utrs)} ({100*len(peaks_in_utrs)/len(peaks):.1f}%)")
```

---

### 4. RMBase Processing Filtering Too Aggressive

**Hypothesis:** The `process_rmbase_mods.py` support number filter or cell-line filter may be excluding legitimate sites.

**Bash validation:**
```bash
# Check cell-line labels present in raw RMBase m6A file
head -100 data/mods/raw/human.hg38.m6A.result.col29.bed | grep -E "K562|HepG2"

# Count K562 vs HepG2 specific sites
awk -F'\t' '$12 ~ /K562/' data/mods/raw/human.hg38.m6A.result.col29.bed | wc -l
awk -F'\t' '$12 ~ /HepG2/' data/mods/raw/human.hg38.m6A.result.col29.bed | wc -l
```

---

### 5. BED Format Issues

**Hypothesis:** The classified peak BED files (`discrepant_peaks.bed`, `canonical_peaks.bed`) may have formatting inconsistencies that prevent proper intersection.

**Bash checks:**
```bash
# Check chromosome naming
head -5 data/eclip/IGF2BP1_K562.bed | cut -f1
head -5 data/mods/K562/m6A.bed | cut -f1

# Check sorting
sort -k1,1 -k2,2n data/mods/K562/m6A.bed | head -5

# Check strand consistency
cut -f6 data/eclip/IGF2BP1_K562.bed | sort | uniq -c
cut -f6 data/mods/K562/m6A.bed | sort | uniq -c
```

---

## Comprehensive Diagnostic Script

```python
#!/usr/bin/env python3
"""
Diagnostic script to identify why enrichment analysis is failing.
Run from the project root on the HPC.
"""
import pybedtools
import pandas as pd
from pathlib import Path

def diagnose_overlap_failure(eclip_bed, mod_bed, genome, chrom_sizes):
    print("="*60)
    print("DIAGNOSTIC: Modification Overlap Analysis")
    print("="*60)

    peaks = pybedtools.BedTool(eclip_bed)
    mods = pybedtools.BedTool(mod_bed)

    print(f"\n1. Basic Counts:")
    print(f"   eCLIP peaks: {len(peaks)}")
    print(f"   Modification sites: {len(mods)}")

    print(f"\n2. Exact Overlap:")
    exact_overlap = peaks.intersect(mods, u=True)
    print(f"   {len(exact_overlap)}/{len(peaks)} ({100*len(exact_overlap)/len(peaks):.2f}%)")

    for window in [50, 100, 500]:
        print(f"\n3. Window-based (±{window}bp):")
        peaks_window = peaks.slop(b=window, g=chrom_sizes)
        window_overlap = peaks_window.intersect(mods, u=True)
        print(f"   {len(window_overlap)}/{len(peaks)} ({100*len(window_overlap)/len(peaks):.2f}%)")

    print(f"\n4. Chromosome Coverage:")
    peak_chroms = set([i.chrom for i in peaks])
    mod_chroms = set([i.chrom for i in mods])
    print(f"   Peak chroms (first 5): {sorted(peak_chroms)[:5]}")
    print(f"   Mod chroms  (first 5): {sorted(mod_chroms)[:5]}")
    print(f"   Shared: {len(peak_chroms & mod_chroms)}")
    print(f"   Peaks-only: {peak_chroms - mod_chroms}")

    print(f"\n5. Strand Analysis:")
    peak_strands = {}
    for i in peaks:
        peak_strands[i.strand] = peak_strands.get(i.strand, 0) + 1
    mod_strands = {}
    for i in mods:
        mod_strands[i.strand] = mod_strands.get(i.strand, 0) + 1
    print(f"   Peak strands: {peak_strands}")
    print(f"   Mod strands:  {mod_strands}")

    print(f"\n6. Same-strand Overlap:")
    same_strand = peaks.intersect(mods, u=True, s=True)
    print(f"   {len(same_strand)}/{len(peaks)} ({100*len(same_strand)/len(peaks):.2f}%)")

    print(f"\n7. Example Coordinates (chr1):")
    for p in [p for p in peaks if p.chrom == 'chr1'][:5]:
        print(f"   peak: {p.chrom}:{p.start}-{p.end} ({p.strand})")
    for m in [m for m in mods if m.chrom == 'chr1'][:5]:
        print(f"   mod:  {m.chrom}:{m.start}-{m.end} ({m.strand})")

    print("\n" + "="*60)
    print("DIAGNOSTIC COMPLETE")
    print("="*60)

if __name__ == '__main__':
    diagnose_overlap_failure(
        'data/eclip/IGF2BP1_K562.bed',
        'data/mods/K562/m6A.bed',
        'data/genome/hg38.fa',
        'data/genome/hg38.chrom.sizes'
    )
```

---

## Expected Results (If Pipeline is Working)

For IGF2BP1 (known m6A reader), if the pipeline is functioning correctly:
- 15–30% of discrepant peaks should overlap m6A sites within ±100bp
- Odds ratio > 2.0, p-value < 0.01

If exact overlap is 0–1% but window overlap is >10%, the intersection is too stringent.
If both exact and window overlap are ~0%, the issue is a coordinate system mismatch.

---

## Priority Order for Investigation

1. Run the diagnostic script on IGF2BP1 (positive control) to characterize the failure mode
2. Check chromosome naming consistency (`chr1` vs `1`)
3. Test window-based overlap (±100bp) vs exact intersection in `enrichment_analysis.py`
4. Validate that peak extension uses the summit position (column 10), not peak start
5. Check 3' UTR coverage overlap between eCLIP peaks and m6A sites
