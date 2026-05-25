# Recreating `scripts/process_rbns_enrichment.py`

This script converts raw RBNS enrichment files (R-values from ENCODE) into Z-scores, which are the format expected by the main pipeline. It is heavier on data manipulation than the previous script, making use of **pandas** and **numpy**.

---

## What This Script Does (Big Picture)

ENCODE provides RBNS data as tab-separated files where each row is a k-mer (e.g. `GCAUG`) and each column is a protein concentration. The values are R-values (fold-enrichment ratios, where 1.0 = no enrichment).

This script:
1. **Parses** the raw TSV into a pandas DataFrame
2. **Selects** R-values from the highest concentration column (or the max across all)
3. **Normalizes** them to Z-scores: `Z = (R - mean_R) / std_R`
4. **Writes** a clean CSV: `kmer, z_score, r_value`, sorted by Z-score descending
5. **Optionally validates** existing output files

The output CSVs are what `src/peak_analysis.py` loads to score eCLIP peaks.

---

## File Structure Overview

```
imports + logging setup
parse_encode_enrichment_tsv()
convert_rvalues_to_zscores()
process_single_file()
process_all_files()
validate_zscore_file()
main()
```

---

## Key New Concepts: pandas and numpy

These two libraries are used throughout. Here is what you need to know:

```python
import pandas as pd
import numpy as np

# A DataFrame is like a spreadsheet — rows and columns
df = pd.read_csv('file.csv')         # read a CSV into a DataFrame
df = pd.read_csv('file.tsv', sep='\t')  # read a tab-separated file

df.columns          # list of column names
df.columns[0]       # first column name
df['col']           # a single column (a Series — like a list with an index)
df[['col1','col2']] # multiple columns (a smaller DataFrame)

# String operations on a column
df['kmer'] = df['kmer'].str.strip()   # strip whitespace from every value
df['kmer'] = df['kmer'].str.upper()   # uppercase every value

# Convert a column to numeric (non-numeric becomes NaN)
df['col'] = pd.to_numeric(df['col'], errors='coerce')

# Max across multiple columns, row-by-row
df[['a', 'b', 'c']].max(axis=1)   # returns a Series of per-row maximums

# Boolean mask — True where values are NOT NaN
mask = ~df['col'].isna()
df['col'][mask]     # only the non-NaN values

# Statistics
series.mean()       # arithmetic mean
series.std()        # standard deviation

# Create a new DataFrame from scratch
result = pd.DataFrame({'col1': [1, 2], 'col2': [3, 4]})

# Drop rows with any NaN
df = df.dropna()

# Sort by a column
df = df.sort_values('z_score', ascending=False)

# Reset the row index after sorting
df = df.reset_index(drop=True)

# Save to CSV (index=False means don't write row numbers)
df.to_csv('output.csv', index=False)

# Access the first row
df.iloc[0]['kmer']
```

---

## Step 1: Imports and Logging Setup

```python
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import sys
import logging
```

**Logging** is used here instead of plain `print()`. It automatically adds a timestamp and severity level to every message, which is useful for longer-running scripts.

```python
# Set up at the top of the file, outside all functions
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Then use like this anywhere in the file:
logger.info("Processing file...")        # normal message
logger.warning("Missing data found")    # warning
logger.error("File not found")          # error
```

The format string `%(asctime)s`, `%(levelname)s`, and `%(message)s` are placeholders that get filled in automatically — you never touch them directly.

---

## Step 2: `parse_encode_enrichment_tsv(tsv_path)`

### Goal
Read the raw ENCODE enrichment TSV into a pandas DataFrame and return the column names needed by the next function.

### Input
| Parameter | Type | Description |
|---|---|---|
| `tsv_path` | `str` or `Path` | Path to the ENCODE enrichment TSV file |

### Output
Returns a **tuple** of three values: `(kmer_col_name, concentration_columns, dataframe)`

| Return value | Type | Description |
|---|---|---|
| `kmer_col` | `str` | Name of the first column (e.g. `'IGF2BP1'`) |
| `conc_cols` | `list of str` | Names of the concentration columns (e.g. `['0', '5', '20', '80', '320', '1300']`) |
| `df` | `DataFrame` | Full table, cleaned |

### Input File Format
```
IGF2BP1    0    5    20    80    320    1300
AAAAA      0.98 1.00 1.02  1.05  1.08   1.10
AAAAC      0.97 0.99 1.01  1.04  1.07   1.09
...
```
The first row is a header. The first column is the k-mer. All other columns are protein concentrations (nM).

### Logic
```
read the TSV into a DataFrame with pd.read_csv(..., sep='\t', header=0)
extract kmer_col = df.columns[0]       (first column name)
extract conc_cols = list(df.columns[1:])  (all other column names)
clean the kmer column: strip whitespace and uppercase
for each concentration column: convert to numeric (use errors='coerce' for safety)
return kmer_col, conc_cols, df
```

---

## Step 3: `convert_rvalues_to_zscores(df, kmer_col, conc_cols, use_max_conc=True)`

### Goal
This is the core math function. It takes R-values from the DataFrame, computes Z-scores, and returns a clean sorted result table.

### Inputs
| Parameter | Type | Default | Description |
|---|---|---|---|
| `df` | `DataFrame` | required | Parsed table from `parse_encode_enrichment_tsv` |
| `kmer_col` | `str` | required | Name of the k-mer column |
| `conc_cols` | `list` | required | List of concentration column names |
| `use_max_conc` | `bool` | `True` | If True, use only the last (highest) concentration column; if False, take the max R-value across all concentration columns per k-mer |

### Output
Returns a **tuple**: `(result_df, mean_r, std_r, conc_used)`

| Return value | Type | Description |
|---|---|---|
| `result` | `DataFrame` | Columns: `kmer`, `z_score`, `r_value`; sorted Z-score descending |
| `mean_r` | `float` | Mean R-value used in normalization |
| `std_r` | `float` | Standard deviation used in normalization |
| `conc_used` | `str` | Which concentration was used (for logging) |

### The Z-score Formula
```
Z = (R - mean_R) / std_R
```
Where `mean_R` and `std_R` are computed only from non-NaN values in the selected R-value column.

### Logic
```
if use_max_conc is True:
    r_values = df[conc_cols[-1]]    ← last column = highest concentration
    conc_used = conc_cols[-1]
else:
    r_values = df[conc_cols].max(axis=1)   ← row-wise max across all columns
    conc_used = "max_across_all"

compute valid_mask = rows where r_values is NOT NaN
valid_r = r_values filtered by valid_mask

mean_r = valid_r.mean()
std_r = valid_r.std()

if std_r == 0 or std_r is NaN:
    log a warning
    z_scores = Series of 0.0, same length as r_values
else:
    z_scores = (r_values - mean_r) / std_r   ← pandas handles this element-wise

build result DataFrame:
    column 'kmer'    = df[kmer_col]
    column 'z_score' = z_scores
    column 'r_value' = r_values

drop any rows with NaN values
sort by 'z_score' descending
reset the index

return result, mean_r, std_r, conc_used
```

### Key Pattern — Checking for NaN or Zero std
```python
import numpy as np

if std_r == 0 or np.isnan(std_r):
    # Can't divide by zero — return flat Z-scores
    z_scores = pd.Series([0.0] * len(r_values))
else:
    z_scores = (r_values - mean_r) / std_r
```

---

## Step 4: `process_single_file(input_path, output_path, use_max_conc=True)`

### Goal
Orchestrates the full pipeline for one file: parse → convert → save → return stats.

### Inputs
| Parameter | Type | Default | Description |
|---|---|---|---|
| `input_path` | `str` or `Path` | required | TSV file to read |
| `output_path` | `str` or `Path` | required | CSV file to write |
| `use_max_conc` | `bool` | `True` | Passed through to `convert_rvalues_to_zscores` |

### Output
Returns a `dict` of statistics, or `None` if the input file doesn't exist.

```python
{
    'input_file': 'IGF2BP1_enrichment.tsv',
    'output_file': 'IGF2BP1_zscores.csv',
    'n_kmers': 1024,
    'concentration_used': '1300',
    'mean_r': 1.02,
    'std_r': 0.08,
    'z_min': -3.4,
    'z_max': 12.1,
    'top_kmer': 'GCAUG',
    'top_zscore': 12.1
}
```

### Logic
```
convert input_path and output_path to Path objects
if input_path doesn't exist: log error, return None
log that we are processing this file

call parse_encode_enrichment_tsv(input_path)  → kmer_col, conc_cols, df
log: number of k-mers and concentration columns found

call convert_rvalues_to_zscores(df, kmer_col, conc_cols, use_max_conc)  → result, mean_r, std_r, conc_used

create output directory if it doesn't exist
save result DataFrame to CSV with result.to_csv(output_path, index=False)

build and return stats dict
    use result.iloc[0]['kmer'] for top_kmer  (first row = highest Z-score)
    guard against empty result with:  if len(result) > 0 else 'N/A'
```

### Key Patterns
```python
# Accessing the first row of a sorted DataFrame
top_kmer = result.iloc[0]['kmer'] if len(result) > 0 else 'N/A'
top_z    = result.iloc[0]['z_score'] if len(result) > 0 else 0

# Logging with formatted floats
logger.info(f"  Z-score range: {z_min:.2f} to {z_max:.2f}")
```

---

## Step 5: `process_all_files(input_dir, output_dir=None, use_max_conc=True)`

### Goal
Finds all `*_enrichment.tsv` files in a directory and calls `process_single_file()` on each one.

### Inputs
| Parameter | Type | Default | Description |
|---|---|---|---|
| `input_dir` | `str` or `Path` | required | Directory to search for `*_enrichment.tsv` files |
| `output_dir` | `str`, `Path`, or `None` | `None` → same as `input_dir` | Where to write the output CSVs |
| `use_max_conc` | `bool` | `True` | Passed through |

### Output
Returns a `list` of stats dicts (one per successfully processed file). Empty list if no files found.

### Deriving Output Filename from Input
Given `IGF2BP1_enrichment.tsv`, the output should be `IGF2BP1_zscores.csv`.

```python
rbp_name = tsv_path.stem.replace('_enrichment', '')
# tsv_path.stem gives 'IGF2BP1_enrichment' (filename without extension)
# .replace('_enrichment', '') gives 'IGF2BP1'
output_path = output_dir / f"{rbp_name}_zscores.csv"
```

### Logic
```
convert paths to Path objects
if output_dir is None: output_dir = input_dir
find all files matching pattern "*_enrichment.tsv" using Path.glob()
if none found: log warning, return []

initialize: all_stats = [], success_count = 0, fail_count = 0

for each tsv_path (sorted alphabetically):
    derive rbp_name and output_path
    call process_single_file() → stats
    if stats is not None: append to all_stats, success_count += 1
    else: fail_count += 1

log summary: successful, failed counts
log top k-mer for each RBP (sorted by input filename)

return all_stats
```

### Key Pattern — Glob
```python
# Find all files matching a pattern in a directory
files = list(input_dir.glob("*_enrichment.tsv"))
# '*' matches any characters; this finds files ending in '_enrichment.tsv'

# Sort alphabetically
for tsv_path in sorted(files):
    ...

# Sort a list of dicts by a key
sorted(all_stats, key=lambda x: x['input_file'])
```

---

## Step 6: `validate_zscore_file(csv_path)`

### Goal
Reads a previously-written Z-score CSV and checks it has the right columns, isn't empty, has consistent k-mer lengths, and has no NaN Z-scores. Useful for sanity-checking after batch processing.

### Input
| Parameter | Type | Description |
|---|---|---|
| `csv_path` | `str` or `Path` | Path to a `*_zscores.csv` file |

### Output
Returns `True` if all checks pass, `False` otherwise. Also logs details via `logger`.

### Logic
```
wrap everything in try/except Exception to catch any read errors

load CSV with pd.read_csv()
check that columns {'kmer', 'z_score', 'r_value'} are all present
    use: required_cols.issubset(df.columns)
    if not: log error, return False

check that len(df) > 0 (not empty)
    if empty: log error, return False

check k-mer length consistency:
    df['kmer'].str.len().unique() → should be length-1 array (all same length)
    if more than one unique length: log warning (not a hard failure)

check for NaN Z-scores:
    nan_count = df['z_score'].isna().sum()
    if > 0: log warning

log summary and return True

if exception: log error, return False
```

### Key Pattern — Set Subset Check
```python
required_cols = {'kmer', 'z_score', 'r_value'}   # a Python set
actual_cols = set(df.columns)

if not required_cols.issubset(actual_cols):
    # not all required columns are present
    ...
```

---

## Step 7: `main()`

### Goal
Defines the CLI. This script has **three modes** controlled by which arguments are passed:

| Mode | Trigger | What it does |
|---|---|---|
| Validate | `--validate DIR` | Validates all `*_zscores.csv` files in DIR |
| Single file | `--single INPUT OUTPUT` | Processes one TSV → one CSV |
| Batch (default) | no special flag | Processes all TSVs in `--input-dir` |

### Arguments
| Flag | Type | Default | Description |
|---|---|---|---|
| `--input-dir` | `str` | `'data/rbns'` | Directory with `*_enrichment.tsv` files |
| `--output-dir` | `str` | `None` | Where to write output CSVs |
| `--single INPUT OUTPUT` | 2 strings | — | Process one file |
| `--validate DIR` | `str` | — | Validate output files in DIR |
| `--use-max-across-concentrations` | flag | False | Use max across all concs instead of just the last |

### The `--single` Argument Pattern
`nargs=2` means this argument expects exactly 2 values after it:
```python
parser.add_argument('--single', nargs=2, metavar=('INPUT', 'OUTPUT'))
# Usage: --single data/rbns/IGF2BP1_enrichment.tsv data/rbns/IGF2BP1_zscores.csv
# args.single will be ['data/rbns/IGF2BP1_enrichment.tsv', 'data/rbns/IGF2BP1_zscores.csv']
input_path, output_path = args.single   # unpack the list
```

### The `use_max_conc` Inversion
The flag is named `--use-max-across-concentrations` (use max across all cols), but the internal variable is `use_max_conc` (use only the highest col). They are opposites:
```python
use_max_conc = not args.use_max_across_concentrations
# Flag absent (default) → use_max_conc = True  → use last column
# Flag present          → use_max_conc = False → use max across all columns
```

### Logic
```
create ArgumentParser, add all arguments, call parse_args()

compute use_max_conc = not args.use_max_across_concentrations

if args.validate:
    find all *_zscores.csv files in the validate directory
    if none found: log error, sys.exit(1)
    call validate_zscore_file() on each
    sys.exit(0) if all valid, else sys.exit(1)

elif args.single:
    unpack input_path, output_path = args.single
    call process_single_file()
    sys.exit(0) if stats returned, else sys.exit(1)

else (batch mode):
    log header info
    call process_all_files(args.input_dir, args.output_dir, use_max_conc)
    if no files processed: log error, sys.exit(1)
    log next-step hint
```

### Key Pattern — `sys.exit()`
```python
import sys
sys.exit(0)   # exit successfully (0 = success in Unix convention)
sys.exit(1)   # exit with error code (non-zero = something went wrong)
```
This matters if other scripts or shell pipelines check the exit code of this script.

---

## Recommended Build Order

1. Write imports and logging setup
2. Write `parse_encode_enrichment_tsv()` — test it by printing the DataFrame on a real file
3. Write `convert_rvalues_to_zscores()` — test the math on a small example
4. Write `process_single_file()` — test end-to-end on one TSV
5. Write `process_all_files()` — test on the `data/rbns/` directory
6. Write `validate_zscore_file()` — test on one of your outputs
7. Write `main()` — wire everything to the CLI

---

## Validation

```bash
# Run original on one file
python scripts/process_rbns_enrichment_original.py --single data/rbns/IGF2BP1_enrichment.tsv /tmp/original_zscores.csv

# Run yours
python scripts/process_rbns_enrichment.py --single data/rbns/IGF2BP1_enrichment.tsv /tmp/new_zscores.csv

# Compare (should be identical)
diff /tmp/original_zscores.csv /tmp/new_zscores.csv
```

---

## Common Mistakes to Watch For

| Mistake | Effect | Fix |
|---|---|---|
| Using `df.columns[1:]` without wrapping in `list()` | Returns a pandas Index object, not a plain list — works in most places but can cause subtle issues | Use `list(df.columns[1:])` |
| Forgetting `~` before `.isna()` | `valid_mask = r_values.isna()` keeps NaN rows instead of discarding them | Use `~r_values.isna()` for "not NaN" |
| Computing `std()` on a Series with all-identical values | Returns 0.0 → division by zero | Guard with `if std_r == 0 or np.isnan(std_r)` |
| Using `df['kmer'].str.len().unique()` and assuming it's a Python list | It returns a numpy array — check `len(...)` or use `!= 1` directly | `if len(df['kmer'].str.len().unique()) != 1` |
| Forgetting `index=False` in `to_csv()` | Adds an extra unnamed row-number column to the output CSV | Always use `result.to_csv(output_path, index=False)` |
