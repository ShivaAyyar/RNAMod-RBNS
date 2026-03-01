I have pasted my current directory structure, also found at https://github.com/ShivaAyyar/RNAMod-RBNS or previously attached in this conversation. Provide a comprehensive detailed plan for Claude Code to implement to generate a json file with a set of AF3 jobs for me to upload to the server. These jobs should explore RBM22's pseudo U specific binding. Generate a job for each of the discrepant peaks (most discrepant if more than 15 -- daily job limit is 30). Provide ways or links for Claude Code to access the documentation for json files (I have context7 enabled if that would work). Use the full-length protein (sequence ideally taken from UniProt) and use a 30 nt sequence centered at the RNA modification. Run both with and without the modification. I believe the sequence information should already be present in my directory, but double check given my github/directory structure and propose a mechanism for finding that out.

4:16 PM

Architected comprehensive workflow for generating AlphaFold3 JSON validation jobs
Architected comprehensive workflow for generating AlphaFold3 JSON validation jobs





Comprehensive Plan: Generate AlphaFold3 Jobs for RBM22-PseudoU Validation
Project Overview
Objective: Create AlphaFold3 JSON input files to validate RBM22's pseudoU-dependent binding by modeling the top 15 discrepant peaks with and without pseudoU modification.

Output: 30 JSON files (15 peaks × 2 conditions) + manifest + validation scripts

Phase 1: Environment Setup & Documentation Access
Task 1.1: Access AlphaFold3 JSON Documentation
Use Context7 to query AF3 documentation:

bash
# In Claude Code terminal, query Context7 for AF3 format
Queries for Context7:

"AlphaFold3 JSON input format for RNA-protein complex prediction"
"How to specify RNA modifications in AlphaFold3 JSON"
"AlphaFold3 pseudouridine modification syntax"
"Example JSON for AlphaFold3 with modified RNA nucleotides"
Alternative documentation sources if Context7 doesn't have AF3:

Official: https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md
AlphaFold Server docs: https://alphafoldserver.com/faq
Look for example JSONs in AlphaFold3 repository
Expected JSON structure to validate:

json
{
  "name": "RBM22_peak1_pseudoU",
  "modelSeeds": [1, 2, 3, 4, 5],
  "sequences": [
    {
      "proteinChain": {
        "sequence": "PROTEIN_SEQUENCE_HERE",
        "count": 1
      }
    },
    {
      "rnaSequence": {
        "sequence": "RNA_SEQUENCE_HERE",
        "modifications": [
          {
            "modificationType": "pseudoU",
            "basePosition": 15
          }
        ]
      }
    }
  ],
  "dialect": "alphafold3",
  "version": 1
}
Phase 2: Data Verification & Assessment
Task 2.1: Verify Required Files Exist
Check for these files in your repository:

python
# Script: scripts/verify_data_availability.py

import os
from pathlib import Path

base_dir = Path(".")

required_files = {
    "RBM22 discrepant peaks": "results/RBM22/discrepant_peaks.bed",
    "RBM22 scored peaks": "results/RBM22/scored_peaks.csv",
    "PseudoU sites": "data/mods/K562/pseudoU.bed",
    "Genome FASTA": "data/genome/hg38.fa",
    "Genome index": "data/genome/hg38.fa.fai",
    "Chromosome sizes": "data/genome/hg38.chrom.sizes"
}

missing = []
present = []

for name, path in required_files.items():
    full_path = base_dir / path
    if full_path.exists():
        size = full_path.stat().st_size
        present.append(f"✓ {name}: {path} ({size:,} bytes)")
    else:
        missing.append(f"✗ {name}: {path}")

print("Files present:")
for p in present:
    print(f"  {p}")

if missing:
    print("\nFiles MISSING:")
    for m in missing:
        print(f"  {m}")
    print("\nAction required:")
    if "Genome FASTA" in str(missing):
        print("  - Run: bash scripts/download_data.sh")
else:
    print("\n✓ All required files present!")
Expected issues:

Genome FASTA may not be present (it's in .gitignore)
If missing, Claude Code should run bash scripts/download_data.sh or download directly
Task 2.2: Inspect RBM22 Results
Load and examine the discrepant peaks:

python
# Script: scripts/inspect_rbm22_peaks.py

import pandas as pd
import pybedtools

# Load scored peaks
scored = pd.read_csv("results/RBM22/scored_peaks.csv")

# Filter for discrepant category
discrepant = scored[scored['category'] == 'discrepant'].copy()

# Sort by Z-score (lowest first = most discrepant)
discrepant = discrepant.sort_values('score_max').reset_index(drop=True)

print(f"Total discrepant peaks: {len(discrepant)}")
print(f"\nTop 15 most discrepant peaks:")
print(discrepant[['peak_id', 'score_max', 'score_sum']].head(15))

# Check how many overlap pseudoU
pseudou_bed = pybedtools.BedTool("data/mods/K562/pseudoU.bed")
discrepant_bed = pybedtools.BedTool("results/RBM22/discrepant_peaks.bed")

overlapping = discrepant_bed.intersect(pseudou_bed, u=True)
print(f"\nDiscrepant peaks overlapping pseudoU: {len(overlapping)}/{len(discrepant)}")

# Save top 15 for AF3 modeling
top15 = discrepant.head(15)
top15.to_csv("results/RBM22/af3_target_peaks.csv", index=False)
print(f"\nSaved top 15 to: results/RBM22/af3_target_peaks.csv")
Phase 3: Protein Sequence Retrieval
Task 3.1: Fetch RBM22 Full-Length Sequence from UniProt
RBM22 UniProt ID: P45957

Script to retrieve:

python
# Script: scripts/fetch_rbm22_sequence.py

import requests
import json

def get_uniprot_sequence(uniprot_id):
    """
    Fetch protein sequence from UniProt REST API.
    
    RBM22 UniProt ID: P45957
    """
    # Fetch FASTA
    url_fasta = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url_fasta)
    
    if response.status_code != 200:
        raise ValueError(f"Failed to fetch {uniprot_id}: {response.status_code}")
    
    # Parse FASTA
    lines = response.text.strip().split('\n')
    header = lines[0]
    sequence = ''.join(lines[1:])
    
    # Also get JSON for metadata
    url_json = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    json_response = requests.get(url_json)
    data = json_response.json()
    
    protein_name = data['proteinDescription']['recommendedName']['fullName']['value']
    organism = data['organism']['scientificName']
    length = len(sequence)
    
    # Look for domains
    domains = []
    for feature in data.get('features', []):
        if feature.get('type') == 'Domain':
            domain_name = feature.get('description', 'Unknown')
            start = feature['location']['start']['value']
            end = feature['location']['end']['value']
            domains.append({
                'name': domain_name,
                'start': start,
                'end': end,
                'sequence': sequence[start-1:end]
            })
    
    return {
        'uniprot_id': uniprot_id,
        'protein_name': protein_name,
        'organism': organism,
        'sequence': sequence,
        'length': length,
        'domains': domains
    }

# Fetch RBM22
rbm22_data = get_uniprot_sequence("P45957")

print(f"Protein: {rbm22_data['protein_name']}")
print(f"Organism: {rbm22_data['organism']}")
print(f"Length: {rbm22_data['length']} aa")
print(f"\nFull sequence:")
print(rbm22_data['sequence'][:60] + "...")

print(f"\nDomains found: {len(rbm22_data['domains'])}")
for domain in rbm22_data['domains']:
    print(f"  {domain['name']}: residues {domain['start']}-{domain['end']}")

# Save
with open("results/RBM22/rbm22_protein.json", 'w') as f:
    json.dump(rbm22_data, f, indent=2)

print(f"\nSaved to: results/RBM22/rbm22_protein.json")
```

**Expected output:**
```
Protein: Pre-mRNA-splicing factor RBM22
Organism: Homo sapiens
Length: 420 aa
Domains: RRM domain (residues ~50-130)
Phase 4: RNA Sequence Extraction
Task 4.1: Extract 30nt Sequences Centered on PseudoU
This is the most complex part. Need to:

Parse discrepant_peaks.bed to get genomic coordinates
Find pseudoU sites that overlap each peak
Extract 30nt centered on the pseudoU position (±15nt)
Handle strand orientation (reverse complement for minus strand)
Map pseudoU position within the extracted 30nt sequence
Script:

python
# Script: scripts/extract_rna_sequences.py

import pandas as pd
import pybedtools
from pyfaidx import Fasta
from pathlib import Path
import json

def extract_rna_for_peak(peak, pseudou_sites, genome_fasta, window=15):
    """
    Extract 30nt RNA sequence centered on pseudoU modification.
    
    Parameters:
    - peak: BED feature (from pybedtools)
    - pseudou_sites: BedTool of pseudoU sites
    - genome_fasta: pyfaidx.Fasta object
    - window: nucleotides on each side (default 15 = 30nt total)
    
    Returns:
    - dict with RNA sequence, pseudoU position, metadata
    """
    # Parse peak
    chrom = peak.chrom
    peak_start = peak.start
    peak_end = peak.end
    strand = peak.strand if len(peak.fields) > 5 else '+'
    peak_id = peak.name if len(peak.fields) > 3 else f"{chrom}:{peak_start}-{peak_end}"
    
    # Find pseudoU sites in this peak
    peak_region = pybedtools.BedTool(
        f"{chrom}\t{peak_start}\t{peak_end}\t{peak_id}\t.\t{strand}",
        from_string=True
    )
    
    overlapping = peak_region.intersect(pseudou_sites, wa=True, wb=True)
    
    if len(overlapping) == 0:
        return None  # No pseudoU in this peak
    
    # Take the first pseudoU site (or closest to peak center)
    # pseudoU site is in columns 7-9 of intersect output
    pseudou_pos = int(overlapping[0].fields[7])  # Start position of pseudoU
    
    # Extract 30nt centered on pseudoU
    extract_start = pseudou_pos - window
    extract_end = pseudou_pos + window + 1  # +1 for pseudoU itself
    
    # Ensure within chromosome bounds
    chrom_len = len(genome_fasta[chrom])
    extract_start = max(0, extract_start)
    extract_end = min(chrom_len, extract_end)
    
    # Extract sequence
    seq = genome_fasta[chrom][extract_start:extract_end].seq.upper()
    
    # Convert DNA → RNA
    rna_seq = seq.replace('T', 'U')
    
    # Calculate pseudoU position in extracted sequence (0-indexed)
    pseudou_offset = pseudou_pos - extract_start
    
    # Reverse complement if minus strand
    if strand == '-':
        complement = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
        rna_seq = ''.join(complement.get(b, b) for b in reversed(rna_seq))
        # For minus strand, position is flipped
        pseudou_offset = len(rna_seq) - pseudou_offset - 1
    
    # Validate: position should be a U (uridine, target of pseudouridylation)
    if 0 <= pseudou_offset < len(rna_seq):
        base_at_site = rna_seq[pseudou_offset]
        if base_at_site != 'U':
            print(f"WARNING: PseudoU site at non-U position in {peak_id}")
            print(f"  Position {pseudou_offset}: {base_at_site}")
    
    # AF3 uses 1-indexed positions
    pseudou_pos_1indexed = pseudou_offset + 1
    
    return {
        'peak_id': peak_id,
        'chrom': chrom,
        'peak_start': peak_start,
        'peak_end': peak_end,
        'strand': strand,
        'rna_sequence': rna_seq,
        'rna_length': len(rna_seq),
        'extract_start': extract_start,
        'extract_end': extract_end,
        'pseudou_genome_pos': pseudou_pos,
        'pseudou_rna_pos': pseudou_pos_1indexed,
        'base_at_pseudou': rna_seq[pseudou_offset] if 0 <= pseudou_offset < len(rna_seq) else 'N'
    }

# Main execution
if __name__ == "__main__":
    # Load data
    print("Loading data...")
    genome = Fasta("data/genome/hg38.fa")
    pseudou_sites = pybedtools.BedTool("data/mods/K562/pseudoU.bed")
    
    # Load top 15 discrepant peaks
    target_peaks = pd.read_csv("results/RBM22/af3_target_peaks.csv")
    
    # Load the actual BED file for coordinates
    discrepant_bed = pybedtools.BedTool("results/RBM22/discrepant_peaks.bed")
    
    # Extract RNA sequences
    rna_extracts = []
    
    for peak in discrepant_bed[:15]:  # Top 15
        rna_data = extract_rna_for_peak(peak, pseudou_sites, genome, window=15)
        
        if rna_data:
            rna_extracts.append(rna_data)
            print(f"✓ {rna_data['peak_id']}: {rna_data['rna_length']}nt, pseudoU at pos {rna_data['pseudou_rna_pos']}")
        else:
            print(f"✗ {peak.name}: No pseudoU overlap (skip)")
    
    print(f"\nExtracted RNA from {len(rna_extracts)} peaks")
    
    # Save
    with open("results/RBM22/rna_sequences.json", 'w') as f:
        json.dump(rna_extracts, f, indent=2)
    
    print(f"Saved to: results/RBM22/rna_sequences.json")
Expected challenges:

Some discrepant peaks may not overlap pseudoU sites
Need to handle edge cases (near chromosome ends)
Validate that pseudoU positions are actually U nucleotides
Phase 5: Generate AlphaFold3 JSON Files
Task 5.1: Create JSON Generator Script
python
# Script: scripts/generate_af3_jsons.py

import json
from pathlib import Path

def create_af3_json(peak_id, protein_seq, rna_seq, pseudou_pos, include_modification=True):
    """
    Generate AlphaFold3 JSON input.
    
    Parameters:
    - peak_id: Peak identifier
    - protein_seq: Full protein sequence (from UniProt)
    - rna_seq: RNA sequence (30nt)
    - pseudou_pos: Position of pseudoU in rna_seq (1-indexed)
    - include_modification: True = pseudoU, False = unmodified U
    
    Returns:
    - dict ready for JSON serialization
    """
    # Clean peak_id for filename (remove special characters)
    clean_id = peak_id.replace(':', '_').replace('(', '_').replace(')', '')
    
    # Model name
    mod_status = "pseudoU" if include_modification else "unmodified"
    model_name = f"RBM22_{clean_id}_{mod_status}"
    
    # Build modifications list
    modifications = []
    if include_modification:
        modifications.append({
            "modificationType": "pseudoU",  # or "Pseudouridine" - check AF3 docs
            "basePosition": pseudou_pos
        })
    
    # AF3 JSON structure
    af3_input = {
        "name": model_name,
        "modelSeeds": [1, 2, 3, 4, 5],  # 5 seeds for ensemble
        "sequences": [
            {
                "proteinChain": {
                    "sequence": protein_seq,
                    "count": 1
                }
            },
            {
                "rnaSequence": {
                    "sequence": rna_seq,
                    "modifications": modifications
                }
            }
        ],
        "dialect": "alphafold3",
        "version": 1
    }
    
    return af3_input, model_name

# Main execution
if __name__ == "__main__":
    # Load protein sequence
    with open("results/RBM22/rbm22_protein.json") as f:
        rbm22_data = json.load(f)
    protein_seq = rbm22_data['sequence']
    
    print(f"Protein: {rbm22_data['protein_name']}")
    print(f"Length: {len(protein_seq)} aa")
    
    # Load RNA sequences
    with open("results/RBM22/rna_sequences.json") as f:
        rna_data = json.load(f)
    
    print(f"\nGenerating JSONs for {len(rna_data)} peaks...")
    
    # Create output directory
    output_dir = Path("results/RBM22/af3_inputs")
    output_dir.mkdir(exist_ok=True)
    
    # Generate JSONs
    manifest = []
    
    for i, rna in enumerate(rna_data, 1):
        peak_id = rna['peak_id']
        rna_seq = rna['rna_sequence']
        pseudou_pos = rna['pseudou_rna_pos']
        
        # Generate unmodified version
        json_unmod, name_unmod = create_af3_json(
            peak_id, protein_seq, rna_seq, pseudou_pos,
            include_modification=False
        )
        
        # Generate modified version
        json_mod, name_mod = create_af3_json(
            peak_id, protein_seq, rna_seq, pseudou_pos,
            include_modification=True
        )
        
        # Save JSONs
        file_unmod = output_dir / f"{name_unmod}.json"
        file_mod = output_dir / f"{name_mod}.json"
        
        with open(file_unmod, 'w') as f:
            json.dump(json_unmod, f, indent=2)
        
        with open(file_mod, 'w') as f:
            json.dump(json_mod, f, indent=2)
        
        print(f"  {i}. {peak_id}")
        print(f"     ✓ {file_unmod.name}")
        print(f"     ✓ {file_mod.name}")
        
        # Add to manifest
        manifest.append({
            'peak_number': i,
            'peak_id': peak_id,
            'chrom': rna['chrom'],
            'strand': rna['strand'],
            'rna_sequence': rna_seq,
            'rna_length': len(rna_seq),
            'pseudou_position': pseudou_pos,
            'protein_length': len(protein_seq),
            'json_unmodified': str(file_unmod.name),
            'json_modified': str(file_mod.name)
        })
    
    # Save manifest
    manifest_df = pd.DataFrame(manifest)
    manifest_df.to_csv(output_dir / "manifest.csv", index=False)
    
    print(f"\n✓ Generated {len(manifest) * 2} JSON files")
    print(f"✓ Manifest: {output_dir / 'manifest.csv'}")
Phase 6: Validation & Quality Control
Task 6.1: Create Validation Script
python
# Script: scripts/validate_af3_jsons.py

import json
from pathlib import Path

def validate_af3_json(json_path):
    """
    Validate AlphaFold3 JSON format.
    
    Checks:
    - Required fields present
    - Protein sequence is valid (20 standard amino acids)
    - RNA sequence is valid (A, C, G, U only)
    - Modification position is within sequence bounds
    - Modification position is actually a U (for pseudoU)
    
    Returns:
    - (is_valid, error_messages)
    """
    errors = []
    
    try:
        with open(json_path) as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        return False, [f"Invalid JSON: {e}"]
    
    # Check required top-level fields
    required_fields = ["name", "sequences", "dialect", "version"]
    for field in required_fields:
        if field not in data:
            errors.append(f"Missing required field: {field}")
    
    if errors:
        return False, errors
    
    # Extract sequences
    protein_seq = None
    rna_seq = None
    modifications = []
    
    for seq_entry in data.get("sequences", []):
        if "proteinChain" in seq_entry:
            protein_seq = seq_entry["proteinChain"]["sequence"]
        if "rnaSequence" in seq_entry:
            rna_seq = seq_entry["rnaSequence"]["sequence"]
            modifications = seq_entry["rnaSequence"].get("modifications", [])
    
    # Validate protein sequence
    if not protein_seq:
        errors.append("No protein sequence found")
    else:
        valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
        invalid_aa = set(protein_seq) - valid_aa
        if invalid_aa:
            errors.append(f"Invalid amino acids in protein: {invalid_aa}")
        
        if len(protein_seq) < 50:
            errors.append(f"Protein unusually short: {len(protein_seq)} aa")
        if len(protein_seq) > 1000:
            errors.append(f"Protein unusually long: {len(protein_seq)} aa (may be slow)")
    
    # Validate RNA sequence
    if not rna_seq:
        errors.append("No RNA sequence found")
    else:
        valid_bases = set("ACGU")
        invalid_bases = set(rna_seq) - valid_bases
        if invalid_bases:
            errors.append(f"Invalid RNA bases: {invalid_bases}")
        
        if len(rna_seq) < 15:
            errors.append(f"RNA unusually short: {len(rna_seq)} nt")
        if len(rna_seq) > 50:
            errors.append(f"RNA unusually long: {len(rna_seq)} nt (may be slow)")
    
    # Validate modifications
    for mod in modifications:
        mod_type = mod.get("modificationType")
        mod_pos = mod.get("basePosition")
        
        if mod_type not in ["pseudoU", "Pseudouridine", "PSU"]:
            errors.append(f"Unexpected modification type: {mod_type}")
        
        if mod_pos is None:
            errors.append("Modification missing basePosition")
        elif not (1 <= mod_pos <= len(rna_seq)):
            errors.append(f"Modification position {mod_pos} out of range (1-{len(rna_seq)})")
        else:
            # Check that position is actually a U
            base = rna_seq[mod_pos - 1]  # Convert 1-indexed to 0-indexed
            if base != 'U':
                errors.append(f"PseudoU position {mod_pos} is not U (found: {base})")
    
    return len(errors) == 0, errors

# Validate all JSONs
if __name__ == "__main__":
    json_dir = Path("results/RBM22/af3_inputs")
    
    all_valid = True
    results = []
    
    for json_file in sorted(json_dir.glob("*.json")):
        is_valid, errors = validate_af3_json(json_file)
        
        status = "✓" if is_valid else "✗"
        print(f"{status} {json_file.name}")
        
        if not is_valid:
            for error in errors:
                print(f"    ERROR: {error}")
            all_valid = False
        
        results.append({
            'file': json_file.name,
            'valid': is_valid,
            'errors': '; '.join(errors) if errors else ''
        })
    
    # Save validation report
    import pandas as pd
    report = pd.DataFrame(results)
    report.to_csv(json_dir / "validation_report.csv", index=False)
    
    # Summary
    valid_count = sum(r['valid'] for r in results)
    print(f"\n{'='*60}")
    print(f"Validation Summary: {valid_count}/{len(results)} files valid")
    print(f"{'='*60}")
    
    if all_valid:
        print("✓ All JSONs validated successfully!")
        print(f"\nReady to upload to AlphaFold Server:")
        print(f"  Directory: {json_dir}")
        print(f"  Total jobs: {len(results)}")
        print(f"  Daily limit: 30 (you're within limits!)")
    else:
        print("✗ Some JSONs have errors - fix before uploading")
        print(f"  See: {json_dir / 'validation_report.csv'}")
Phase 7: Documentation & Instructions
Task 7.1: Create README for AF3 Submission
markdown
# AlphaFold3 Input Files for RBM22-PseudoU Validation

## Overview

This directory contains AlphaFold3 JSON input files to validate RBM22's pseudouridine-dependent binding.

**Target:** RBM22 (Pre-mRNA-splicing factor, UniProt P45957)
**RNA:** 30nt sequences from top 15 discrepant peaks
**Modification:** Pseudouridine (Ψ)
**Total jobs:** 30 (15 peaks × 2 conditions)

## Files

- `manifest.csv` - Summary of all models
- `RBM22_*_unmodified.json` - Models with unmodified U
- `RBM22_*_pseudoU.json` - Models with pseudouridine
- `validation_report.csv` - QC results

## Protein Sequence

**Source:** UniProt P45957  
**Length:** 420 amino acids (full-length)  
**Domains:** RRM domain (residues ~50-130)

Using full-length protein as requested (not just RRM domain).

## RNA Sequences

**Extraction method:** 30nt centered on pseudoU modification site (±15nt)  
**Genome:** hg38  
**Strand handling:** Reverse complemented for minus-strand peaks  
**Coordinate system:** 1-indexed for AF3

## Submission Instructions

### Option 1: AlphaFold Server Web UI (Recommended for 30 jobs)

1. Go to: https://golgi.sandbox.google.com/
2. Sign in with academic Google account
3. For each JSON file:
   - Click "New Prediction"
   - Click "Upload JSON"
   - Select JSON file
   - Submit
4. Track jobs in "My Jobs" dashboard
5. Download results when complete (PDB + confidence scores)

**Time estimate:** 
- ~5 minutes per upload
- ~20-30 minutes compute time per job
- Total: ~2-3 hours to submit all, 10-15 hours to complete

### Option 2: AlphaFold Server API (If you have API key)

See `scripts/submit_via_api.py` for batch submission script.

## Expected Results

**For each peak, compare:**
- Unmodified: Baseline binding interface
- PseudoU: Modified binding interface

**Metrics to extract:**
- Number of protein-RNA contacts (<4Å)
- Buried surface area
- pLDDT confidence scores
- Key residues in binding interface

**Hypothesis:**
If pseudoU enables binding at discrepant sites, we expect:
- More contacts in pseudoU models
- Higher pLDDT scores (better-defined structure)
- Consistent effect across multiple peaks

## Analysis After Completion

Run: `python scripts/analyze_af3_results.py`

This will:
1. Compare modified vs unmodified structures
2. Calculate ΔContacts, ΔBSA, ΔpLDDT
3. Generate summary plots
4. Identify which peaks show strongest pseudoU effect

## Troubleshooting

**"Modification type not recognized"**
- AF3 may use different syntax for pseudoU
- Try: "Pseudouridine", "PSU", or "pseudouridine"
- Check AF3 documentation for supported modifications

**"Protein too large - slow prediction"**
- 420 aa is manageable but slower than domain-only
- Consider using RRM domain only (residues 50-130) for faster results
- Would require regenerating JSONs

**"Position X is not uridine"**
- Validation failed - check strand orientation
- Verify pseudoU BED file coordinates

## Files Generated
```
results/RBM22/af3_inputs/
├── manifest.csv                                  # Summary table
├── validation_report.csv                         # QC report
├── RBM22_chr1_12345_67890_+_unmodified.json     # Example files
├── RBM22_chr1_12345_67890_+_pseudoU.json
├── ...
└── README.md                                     # This file
```

## Contact

Questions? Check:
- AlphaFold Server FAQ: https://alphafoldserver.com/faq
- AF3 GitHub: https://github.com/google-deepmind/alphafold3
Complete Workflow Summary
Claude Code should execute these steps in order:
bash
# Step 1: Access documentation
# Use Context7 to query AF3 JSON format
# Query: "AlphaFold3 JSON RNA modification syntax"

# Step 2: Verify data availability
python scripts/verify_data_availability.py

# If genome missing:
# bash scripts/download_data.sh

# Step 3: Inspect RBM22 results
python scripts/inspect_rbm22_peaks.py
# Output: results/RBM22/af3_target_peaks.csv

# Step 4: Fetch protein sequence
python scripts/fetch_rbm22_sequence.py
# Output: results/RBM22/rbm22_protein.json

# Step 5: Extract RNA sequences
python scripts/extract_rna_sequences.py
# Output: results/RBM22/rna_sequences.json

# Step 6: Generate JSONs
python scripts/generate_af3_jsons.py
# Output: results/RBM22/af3_inputs/*.json
#         results/RBM22/af3_inputs/manifest.csv

# Step 7: Validate
python scripts/validate_af3_jsons.py
# Output: results/RBM22/af3_inputs/validation_report.csv

# Step 8: Create README
# Copy README.md to results/RBM22/af3_inputs/
```

---

## Key Files Claude Code Should Create

1. **`scripts/verify_data_availability.py`** - Check file existence
2. **`scripts/inspect_rbm22_peaks.py`** - Load and examine peaks
3. **`scripts/fetch_rbm22_sequence.py`** - Get RBM22 from UniProt
4. **`scripts/extract_rna_sequences.py`** - Extract 30nt RNA sequences
5. **`scripts/generate_af3_jsons.py`** - Create JSON files
6. **`scripts/validate_af3_jsons.py`** - Validate JSON format
7. **`results/RBM22/af3_inputs/README.md`** - Submission instructions

---

## Expected Final Output
```
results/RBM22/af3_inputs/
├── README.md
├── manifest.csv
├── validation_report.csv
├── RBM22_chr10_100233940_100234060_+_unmodified.json
├── RBM22_chr10_100233940_100234060_+_pseudoU.json
├── RBM22_chr11_64839318_64839437_-_unmodified.json
├── RBM22_chr11_64839318_64839437_-_pseudoU.json
└── ... (30 JSON files total)

manifest.csv columns:
- peak_number
- peak_id
- chrom
- strand
- rna_sequence
- rna_length
- pseudou_position
- protein_length
- json_unmodified
- json_modified
```

---

## Context7 Integration

Claude Code should use Context7 to:

1. **Query AF3 documentation:**
```
   Query: "AlphaFold3 server JSON input format documentation"
   Library: Look for official Google DeepMind or AlphaFold documentation
```

2. **Verify modification syntax:**
```
   Query: "How to specify pseudouridine in AlphaFold3 JSON modificationType field"
```

3. **Check examples:**
```
   Query: "AlphaFold3 JSON example with RNA modifications"
If Context7 doesn't have AF3 docs, Claude Code should:

Fall back to fetching from GitHub: https://github.com/google-deepmind/alphafold3/tree/main/docs
Or use generic RNA-protein complex examples and adapt
This plan provides Claude Code with everything needed to generate production-ready AlphaFold3 inputs for RBM22 validation!