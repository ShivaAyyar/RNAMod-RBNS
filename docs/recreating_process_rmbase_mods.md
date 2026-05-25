# Recreating `scripts/process_rmbase_mods.py`

This guide walks you through rebuilding the RMBase modification processing script from scratch. The script reads large raw RNA modification data files from RMBase v3.0, filters them by cell line and quality thresholds, and writes clean BED6 files used later in the pipeline.

---

## What This Script Does (Big Picture)

RMBase provides RNA modification sites (m6A, pseudoU, m5C, ac4C, etc.) as large tab-delimited files with 29 columns. This script:

1. **Extracts** those files from `.tar.gz` archives (if needed)
2. **Reads** each line, parses the 29 columns
3. **Filters** rows by cell line, minimum supporting datasets, and optionally motif score
4. **Writes** a clean 6-column BED file: `chrom`, `start`, `end`, `modID`, `score`, `strand`

The output BED files are placed in `data/mods/K562/` and `data/mods/HepG2/`.

---

## File Structure Overview

```
imports
constants (3 dictionaries)
extract_archive()
process_rmbase_file()
setup_all_modifications()
main()
```

Build in this order — each function uses the ones above it.

---

## Step 1: Imports

```python
import argparse      # for building the command-line interface
import gzip          # not used directly here but useful for .gz files
import tarfile       # for extracting .tar.gz archives
import os            # general OS utilities (not heavily used)
import sys           # for sys.exit() when errors occur
from pathlib import Path           # modern, clean way to handle file paths
from collections import defaultdict  # not actually used in the final code — can omit
```

**Key concept — `pathlib.Path`:** Instead of joining strings like `"data/mods/" + "K562"`, you write `Path("data/mods") / "K562"`. It handles slashes correctly on all operating systems and has useful methods like `.exists()`, `.mkdir()`, `.parent`, `.name`, `.stem`.

---

## Step 2: Module-Level Constants (Three Dictionaries)

These are defined at the top of the file, outside any function, so all functions can use them.

### `MOD_FILES`

Maps each modification short name to its raw filename inside the archive.

```python
MOD_FILES = {
    'm6A':   'human.hg38.m6A.result.col29.bed',
    'Pseudo': 'human.hg38.Pseudo.result.col29.bed',
    'm5C':   'human.hg38.m5C.result.col29.bed',
    'ac4C':  'human.hg38.ac4C.result.col29.re.bed',
    'm1A':   'human.hg38.m1A.result.col29.bed',
    'm7G':   'human.hg38.m7G.result.col29.bed',
    'Nm':    'human.hg38.Nm.result.col29.bed',
}
```

### `MOD_ARCHIVES`

Maps each modification to its `.tar.gz` archive filename.

```python
MOD_ARCHIVES = {
    'm6A':   'hg38.m6A.tar.gz',
    'Pseudo': 'hg38.Pseudo.tar.gz',
    'm5C':   'hg38.m5C.tar.gz',
    'ac4C':  'hg38.ac4C.tar.gz',
    'm1A':   'hg38.m1A.tar.gz',
    'm7G':   'hg38.m7G.tar.gz',
    'Nm':    'hg38.Nm.tar.gz',
}
```

### `OUTPUT_NAMES`

Maps each modification to the clean output filename used downstream. Only the 4 primary modifications are listed here — they are the ones the pipeline actually uses.

```python
OUTPUT_NAMES = {
    'm6A':   'm6A.bed',
    'Pseudo': 'pseudoU.bed',
    'm5C':   'm5C.bed',
    'ac4C':  'ac4C.bed',
}
```

**Why dictionaries?** It avoids long `if/elif` chains. You can look up the filename for any modification type with one line: `MOD_FILES['m6A']`.

---

## Step 3: `extract_archive(archive_path, output_dir)`

### Goal
Open a `.tar.gz` file, find the `.bed` file inside it, extract it to `output_dir` if not already there, and return the path to the extracted file.

### Inputs
| Parameter | Type | Description |
|---|---|---|
| `archive_path` | `str` or `Path` | Path to the `.tar.gz` file |
| `output_dir` | `Path` | Directory to extract into |

### Output
Returns a `Path` pointing to the extracted `.bed` file, or `None` if the archive doesn't exist or contains no `.bed` files.

### Logic
1. Convert `archive_path` to a `Path` object
2. If the file doesn't exist, return `None`
3. Open the archive with `tarfile.open(archive_path, 'r:gz')`
4. Get the list of filenames inside: `tar.getnames()`
5. Loop over them, find the first one ending in `.'bed'`
6. Build the target path: `output_dir / member`
7. If that path doesn't already exist, extract it: `tar.extract(member, output_dir)`
8. Return the target path

### Key Python Patterns
```python
# Opening a tar.gz file
import tarfile
with tarfile.open('file.tar.gz', 'r:gz') as tar:
    names = tar.getnames()        # list of filenames inside the archive
    tar.extract('file.bed', '/output/dir/')  # extract one file

# Check if string ends with something
if filename.endswith('.bed'):
    ...

# Check if a file already exists
path = Path('/some/dir/file.bed')
if not path.exists():
    ...
```

---

## Step 4: `process_rmbase_file(input_path, output_path, cell_line=None, min_support=1, min_motif_score=None)`

### Goal
This is the core function. It reads a raw 29-column RMBase file line by line, applies filters, and writes a clean 6-column BED file.

### Inputs
| Parameter | Type | Default | Description |
|---|---|---|---|
| `input_path` | `str` or `Path` | required | Raw RMBase `.bed` file to read |
| `output_path` | `str` or `Path` | required | Where to write the filtered BED6 output |
| `cell_line` | `str` or `None` | `None` | If provided (e.g. `'K562'`), only keep rows where this cell line appears in column 12 |
| `min_support` | `int` | `1` | Minimum value in column 8 (number of supporting experiments); skip rows below this |
| `min_motif_score` | `float` or `None` | `None` | If provided, skip rows where column 20 is below this value |

### Output
Returns a `dict` with counts: `total_lines`, `kept_lines`, `cell_line_matches`, `support_filter`, `motif_filter`.

Also writes the filtered BED6 file to `output_path`.

### The RMBase Column Layout (0-indexed)
| Index | Column name | Description |
|---|---|---|
| 0 | chrom | chromosome (e.g. `chr1`) |
| 1 | start | 0-based start coordinate |
| 2 | end | end coordinate |
| 3 | modID | modification identifier |
| 4 | score | original RMBase score |
| 5 | strand | `+` or `-` |
| 7 | supportNum | integer count of supporting datasets |
| 11 | cellList | comma-separated list of cell lines (e.g. `K562,HepG2,HeLa`) |
| 19 | motifScore | float 0–5 confidence score (may be `'na'`) |

### Logic — Walk Through Each Line

```
for each line in input file:
    increment total_lines counter
    split the line on tabs → list of fields
    skip if fewer than 12 fields (malformed line)

    extract: chrom, start, end, mod_id, score, strand (fields 0–5)
    extract: support_num from fields[7] — convert to int if it's a digit, else default to 1
    extract: cell_list from fields[11]

    try to extract motif_score from fields[19]:
        - skip if fields[19] == 'na'
        - try converting to float, catch ValueError

    FILTER 1 — support:
        if support_num < min_support: increment support_filter, skip

    FILTER 2 — motif score:
        if min_motif_score is set AND motif_score was found:
            if motif_score < min_motif_score: increment motif_filter, skip

    FILTER 3 — cell line:
        if cell_line is specified:
            split cell_list on comma, strip whitespace → list of cell types
            if cell_line is NOT in that list: skip (don't increment, just continue)
            if it IS in the list: increment cell_line_matches

    PASSED all filters — write BED6 line:
        compute bed_score:
            if motif_score exists: bed_score = int(motif_score * 200)
            else: bed_score = min(support_num * 100, 1000)
        write: chrom \t start \t end \t mod_id \t bed_score \t strand \n
        increment kept_lines
```

### Key Python Patterns
```python
# Open two files at once (read and write)
with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
    for line in infile:
        fields = line.strip().split('\t')

# Safe integer conversion
support_num = int(fields[7]) if fields[7].isdigit() else 1

# Safe float conversion with try/except
try:
    motif_score = float(fields[19])
except ValueError:
    motif_score = None

# Split a comma-separated list and strip whitespace
cell_types = [c.strip() for c in cell_list.split(',')]

# Cap a value at a maximum
bed_score = min(support_num * 100, 1000)

# Create parent directories automatically
from pathlib import Path
Path(output_path).parent.mkdir(parents=True, exist_ok=True)

# Format large numbers with commas in f-strings
print(f"  Total lines: {stats['total_lines']:,}")
```

---

## Step 5: `setup_all_modifications(mods_dir, cell_lines=None)`

### Goal
A high-level convenience function that loops over all 4 primary modifications and both cell lines, calling `extract_archive()` and `process_rmbase_file()` for each combination. This is what runs when you pass `--setup-all` on the command line.

### Inputs
| Parameter | Type | Default | Description |
|---|---|---|---|
| `mods_dir` | `str` or `Path` | required | Base directory (`data/mods/`) containing archives and raw files |
| `cell_lines` | `list` or `None` | `None` → `['K562', 'HepG2']` | Cell lines to process |

### Output
No return value. Creates output BED files in `mods_dir/{cell_line}/` and prints a verification summary showing how many sites were written for each file.

### Logic
```
set primary_mods = ['m6A', 'Pseudo', 'm5C', 'ac4C']
if cell_lines is None, default to ['K562', 'HepG2']

PHASE 1 — extract archives:
    for each mod_type in primary_mods:
        build archive path from MOD_ARCHIVES[mod_type]
        if archive exists: call extract_archive()

PHASE 2 — process for each cell line:
    for each cell_line in cell_lines:
        create output directory: mods_dir / cell_line
        for each mod_type in primary_mods:
            build input_file path from MOD_FILES[mod_type]
            build output_file path from OUTPUT_NAMES[mod_type]
            if input_file doesn't exist: print warning, skip

            SPECIAL CASE for m6A only:
                scan first 10,000 lines to count occurrences of cell_line
                if cell_line is found (count > 0):
                    call process_rmbase_file() with cell_line filter
                else:
                    print note about missing cell data
                    call process_rmbase_file() with cell_line=None, min_support=2

            FOR Pseudo, m5C, ac4C:
                no cell-specific data exists in RMBase for these
                call process_rmbase_file() with cell_line=None, min_support=1

PHASE 3 — verify outputs:
    for each cell line and mod type:
        if output file exists: count lines and print
        else: print 'NOT CREATED'
```

### Key Python Patterns
```python
# Default mutable argument pattern (use None, then set inside)
def my_func(items=None):
    items = items or ['default1', 'default2']

# Count lines in a file
with open(output_file, 'r') as f:
    line_count = sum(1 for _ in f)

# Check a large file for a string, stopping early
count = 0
with open(input_file, 'r') as f:
    for i, line in enumerate(f):
        if 'K562' in line:
            count += 1
        if i > 10000:
            break
```

---

## Step 6: `main()`

### Goal
Defines the command-line interface using `argparse`. Parses arguments and routes to the right function.

### Logic
```
create ArgumentParser with description and epilog (usage examples)

add arguments:
    --mods-dir      str, default='data/mods'
    --setup-all     flag (store_true)
    --mod-type      choice from MOD_FILES keys
    --cell-line     str
    --all-sites     flag (store_true)
    --output / -o   str
    --min-support   int, default=1
    --min-motif-score  float

parse args

if --setup-all:
    call setup_all_modifications(mods_dir)
    return

if no --mod-type given:
    print error and exit

build input_file path from MOD_FILES[args.mod_type]

if input_file doesn't exist:
    try to extract from archive
    if still doesn't exist: print error and exit

determine output_file:
    if --output given: use it
    elif --cell-line given: mods_dir / cell_line / OUTPUT_NAMES[mod_type]
    else: mods_dir / '{mod_type}_all.bed'

resolve cell_line:
    if --all-sites: cell_line = None
    else: cell_line = args.cell_line

call process_rmbase_file(input_file, output_file, cell_line, min_support, min_motif_score)
```

### Key Python Pattern — `argparse`
```python
import argparse, sys

parser = argparse.ArgumentParser(
    description='Short description shown in --help',
    formatter_class=argparse.RawDescriptionHelpFormatter,  # preserves newlines in epilog
    epilog="""
Examples:
    python script.py --mod-type m6A --cell-line K562
    """
)

# Required choice argument
parser.add_argument('--mod-type', choices=['m6A', 'Pseudo', 'm5C'])

# Flag (True if present, False if absent)
parser.add_argument('--setup-all', action='store_true')

# Optional with type conversion and default
parser.add_argument('--min-support', type=int, default=1)

# Short and long form
parser.add_argument('--output', '-o', help='Output file path')

args = parser.parse_args()

# Access values:
print(args.mod_type)
print(args.setup_all)  # True or False
print(args.min_support)

# Exit with error
sys.exit(1)
```

---

## Recommended Build Order

1. Write the 3 dictionaries (`MOD_FILES`, `MOD_ARCHIVES`, `OUTPUT_NAMES`)
2. Write `extract_archive()` — test it manually on one archive
3. Write `process_rmbase_file()` — this is the most important function; test on a small file first
4. Write `setup_all_modifications()` — just calls the above two in loops
5. Write `main()` — wire everything to the CLI

---

## Validation

To check your version produces the same output as the original:

```bash
# Run original
python scripts/process_rmbase_mods_original.py --mod-type m6A --cell-line K562 --output /tmp/original_m6A.bed

# Run yours
python scripts/process_rmbase_mods.py --mod-type m6A --cell-line K562 --output /tmp/new_m6A.bed

# Compare line counts
wc -l /tmp/original_m6A.bed /tmp/new_m6A.bed

# Compare content (should produce no output if identical)
diff /tmp/original_m6A.bed /tmp/new_m6A.bed
```

---

## Common Mistakes to Watch For

| Mistake | Effect | Fix |
|---|---|---|
| Using `split('\t')` without `.strip()` on the line first | Last field includes a `\n`, causing comparison failures | Always call `line.strip()` before splitting |
| Not handling `fields[7] == '.'` or non-digit support values | Crash on `int()` conversion | Use `int(fields[7]) if fields[7].isdigit() else 1` |
| Forgetting `output_path.parent.mkdir(parents=True, exist_ok=True)` | Crash if output directory doesn't exist | Add this before opening the output file |
| Cell line check using `in` on the raw string instead of the split list | `'K562'` matching inside `'K562_ENCODE'` unexpectedly | Split on comma first, then check membership in the list |
