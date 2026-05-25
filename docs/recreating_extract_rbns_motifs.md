# Recreating `scripts/extract_rbns_motifs.py`

This script extracts RNA-binding protein (RBP) motif data from the supplementary Excel file of Dominguez et al. (2018). It is simpler than the previous scripts — no CLI, no filtering loops — but introduces **reading Excel files with pandas** and navigating irregular spreadsheet layouts.

> **Important context:** This script only extracts the top ~5–10 motifs per RBP from the paper's supplementary table. It does NOT produce comprehensive Z-scores for all 1024 5-mers. The output is a partial fallback for when the full RBNS pipeline hasn't been run.

---

## What This Script Does (Big Picture)

1. Opens `data/rbns/mmc4.xlsx` (Supplementary Table S3 from Dominguez et al. 2018)
2. Reads two different sheets in that file:
   - Sheet 1: top enriched motifs per RBP (each in the format `UUUUU_1_15.38`)
   - Sheet `logo_5mers.prop_in_logo`: k-mer proportions from motif logos
3. Prints a summary of which RBPs were found
4. Writes one CSV file per RBP to `data/rbns/` in the format `{RBP}_zscores.csv`

---

## File Structure Overview

```
imports
constants (DATA_DIR, OUTPUT_DIR, TARGET_RBPS)
parse_motif_string()
extract_top_motifs_from_table_s3()
extract_logo_proportions()
write_csv_for_rbp()
main()
```

---

## Key New Concept: Reading Excel Files with pandas

```python
import pandas as pd

# Read first sheet (index 0), with the first row as header
df = pd.read_excel('file.xlsx', sheet_name=0, header=0)

# Read a sheet by name, with no header (you handle the layout manually)
df = pd.read_excel('file.xlsx', sheet_name='logo_5mers.prop_in_logo', header=None)

# Access a specific cell by row and column index (0-based)
value = df.iloc[row_index, col_index]

# Iterate over all rows as (index, row) pairs
for _, row in df.iterrows():
    val = row['ColumnName']    # access by column name
    val = row.get('ColName')   # same, but returns None if column doesn't exist

# Check if a value is NaN (missing)
import pandas as pd
pd.isna(value)    # True if NaN, None, or np.nan
pd.notna(value)   # True if it has a real value

# Enumerate a pandas Series (like enumerate on a list)
for i, val in enumerate(df.iloc[1]):   # iterate over row 1
    print(i, val)
```

---

## Step 1: Imports and Constants

```python
import pandas as pd
import os
import re   # imported but not used — you can omit it
```

> `os.path.join()` and `os.path.exists()` are used here instead of `pathlib.Path`. Either works — the original was written this way. You can use either in your version.

### Constants

```python
DATA_DIR = "data/rbns"
OUTPUT_DIR = "data/rbns"

TARGET_RBPS = [
    'IGF2BP1', 'IGF2BP2', 'HNRNPC', 'TIA1', 'HNRNPK', 'PCBP2', 'RBFOX2',
    'PTBP1', 'TARDBP', 'QKI', 'SRSF1', 'SRSF9', 'RBM22', 'TRA2A',
    'HNRNPL', 'LIN28B', 'FUS', 'MATR3', 'HNRNPA1', 'HNRNPM',
    'NONO', 'U2AF2', 'EWSR1'
]
```

`TARGET_RBPS` is a plain Python list of strings. It's used to filter — only keep rows in the Excel sheet that correspond to one of these 23 proteins.

---

## Step 2: `parse_motif_string(motif_str)`

### Goal
Parse a single motif string like `'UUUUU_1_15.38'` into its three components.

### Input
| Parameter | Type | Description |
|---|---|---|
| `motif_str` | `str`, `float` (NaN), or `None` | Raw cell value from the Excel sheet |

### Output
Returns a **tuple** of three values: `(kmer, logo_num, zscore)`

| Value | Type | Example | Description |
|---|---|---|---|
| `kmer` | `str` | `'UUUUU'` | The 5-mer sequence |
| `logo_num` | `int` | `1` | Which logo this motif belongs to |
| `zscore` | `float` | `15.38` | Enrichment Z-score |

Returns `(None, None, None)` if the input is NaN or doesn't have at least 3 underscore-separated parts.

### The motif string format
```
UUUUU_1_15.38
 ↑     ↑  ↑
kmer  logo  Z-score
```

### Logic
```
if motif_str is NaN (check with pd.isna()): return None, None, None

convert to string with str(), then split on '_'
if fewer than 3 parts: return None, None, None

kmer    = parts[0]             (e.g. 'UUUUU')
logo_num = int(parts[1])       (e.g. 1)
zscore  = float(parts[2])      (e.g. 15.38)

return kmer, logo_num, zscore
```

### Key Pattern
```python
# Safe split when you're not sure how many parts there are
parts = str(motif_str).split('_')
if len(parts) >= 3:
    kmer = parts[0]
    logo_num = int(parts[1])
    zscore = float(parts[2])
```

---

## Step 3: `extract_top_motifs_from_table_s3()`

### Goal
Read the main sheet of `mmc4.xlsx` and extract the top motifs for each of the target RBPs.

### Input
None (uses `DATA_DIR` constant).

### Output
Returns a **dict**: `{RBP_name: [(kmer, zscore), ...]}`

```python
{
    'IGF2BP1': [('GCATG', 12.5), ('CATGT', 8.3), ...],
    'HNRNPC':  [('TTTTT', 15.2), ...],
    ...
}
```

An empty dict `{}` is returned if the Excel file doesn't exist.

### The Excel sheet layout (Sheet 1)
The first row is a header. The relevant columns are:
- `'RBP'` — the protein name
- `'Motif 5mer_logonum_stepwiseRminus1'` — first motif column
- Several `'Unnamed: 6'`, `'Unnamed: 7'`, etc. columns — additional motifs

### Logic
```
build xlsx_path from DATA_DIR + 'mmc4.xlsx'
if file doesn't exist: print error, return {}

read sheet 0 with header: df = pd.read_excel(xlsx_path, sheet_name=0, header=0)

identify motif columns:
    motif_cols = all column names where 'Motif' or 'Unnamed' appears in the name
    use a list comprehension with: [col for col in df.columns if 'Motif' in str(col) or 'Unnamed' in str(col)]

initialize results = {}

for each row in df (using df.iterrows()):
    rbp = row['RBP']
    if rbp is not in TARGET_RBPS: skip

    motifs = []
    for each col in motif_cols:
        call parse_motif_string(row.get(col))  → kmer, logo_num, zscore
        if kmer is not None and zscore is not None:
            convert kmer from RNA to DNA: kmer_dna = kmer.replace('U', 'T')
            append (kmer_dna, zscore) to motifs

    if motifs is not empty:
        results[rbp] = motifs

return results
```

### Key Patterns
```python
# Filter column names by substring
motif_cols = [col for col in df.columns if 'Motif' in str(col) or 'Unnamed' in str(col)]

# Safe row access — returns None if column doesn't exist
value = row.get('ColumnName')

# Convert RNA to DNA (replace U with T)
kmer_dna = kmer.replace('U', 'T')
```

---

## Step 4: `extract_logo_proportions()`

### Goal
Read the `logo_5mers.prop_in_logo` sheet from `mmc4.xlsx`, which has an irregular layout. Extract k-mer proportion data for each RBP.

### Input
None (uses `DATA_DIR` constant).

### Output
Returns a **dict**: `{RBP_name: [(kmer, proportion), ...]}`

```python
{
    'IGF2BP1': [('GCATG', 0.85), ('CATGT', 0.72), ...],
    ...
}
```

An empty dict `{}` is returned if the file doesn't exist.

### The Irregular Sheet Layout
This sheet has no standard header. The structure is:

```
Row 0:  NAME         Stepwise_R-1   [blank]  [blank]   NAME    Stepwise_R-1  ...
Row 1:  A1CF         NaN            BOLL     NaN       CELF1   NaN           ...
Row 2+: GCAUG        0.85           UUUUU    0.91      AUGCA   0.77          ...
        CATGU        0.72           UUUUA    0.88      ...
```

RBP names are in **row 1** (index 1), at **even-numbered columns** (0, 2, 4, ...).
K-mers are in **rows 2+**, in the same column as the RBP name.
Proportions are in the column immediately **to the right** of the k-mer column.

Because of this layout, you must read it with `header=None` and navigate by row/column index.

### Logic
```
build xlsx_path
if not exists: return {}

read with no header: df = pd.read_excel(xlsx_path, sheet_name='logo_5mers.prop_in_logo', header=None)

initialize results = {}
initialize rbp_cols = {}   (maps column_index → RBP_name)

PHASE 1 — find RBP column positions from row 1:
    rbp_row = df.iloc[1]    ← row at index 1
    for each (i, val) in enumerate(rbp_row):
        if val is not NaN AND val is in TARGET_RBPS:
            rbp_cols[i] = val    ← record which column this RBP is in

PHASE 2 — extract k-mers and proportions:
    for each col_idx, rbp in rbp_cols.items():
        motifs = []
        for row_idx in range(2, len(df)):     ← data starts at row 2
            kmer = df.iloc[row_idx, col_idx]           ← k-mer in this column
            prop = df.iloc[row_idx, col_idx + 1]       ← proportion in next column

            if kmer is not NaN AND prop is not NaN:
                kmer_str = str(kmer).upper()
                try:
                    prop_float = float(prop)
                    motifs.append((kmer_str, prop_float))
                except:
                    pass   ← skip non-numeric proportions

        if motifs:
            results[rbp] = motifs

return results
```

### Key Patterns
```python
# Read Excel with no header (you navigate by index)
df = pd.read_excel('file.xlsx', sheet_name='SheetName', header=None)

# Access a specific cell
cell_value = df.iloc[row_index, col_index]

# Get a whole row as a Series
row = df.iloc[1]   # row at index 1

# Enumerate a pandas Series
for i, val in enumerate(row):
    if pd.notna(val):
        print(i, val)

# Safe float conversion with bare except
try:
    prop_float = float(prop)
except:
    pass   # skip if it can't be converted
```

---

## Step 5: `write_csv_for_rbp(rbp, motifs, output_dir)`

### Goal
Write one CSV file for a single RBP containing its motifs and Z-scores.

### Inputs
| Parameter | Type | Description |
|---|---|---|
| `rbp` | `str` | RBP name (e.g. `'IGF2BP1'`) |
| `motifs` | `list of (str, float)` | List of `(kmer, zscore)` tuples |
| `output_dir` | `str` | Directory to write into |

### Output
No return value. Creates `{output_dir}/{rbp}_zscores.csv` with header `kmer,z_score`.

### Output file format
```
kmer,z_score
GCATG,12.5000
CATGT,8.3000
```

### Logic
```
build output_path = os.path.join(output_dir, f'{rbp}_zscores.csv')

open output_path for writing
write header line: 'kmer,z_score\n'
for each (kmer, zscore) in motifs:
    write: f'{kmer},{zscore:.4f}\n'     ← 4 decimal places

print confirmation message
```

### Key Pattern
```python
# f-string with fixed decimal places
f'{zscore:.4f}'     # always 4 decimal places: 15.38 → '15.3800'
```

---

## Step 6: `main()`

### Goal
Orchestrate everything: call the two extraction functions, print a detailed summary, and write the CSV files.

### Logic
```
print header banner

call extract_top_motifs_from_table_s3()  → top_motifs dict

print how many of the 23 target RBPs were found
for each rbp in TARGET_RBPS:
    if rbp is in top_motifs:
        print: RBP name, count of motifs, top kmer, top Z-score
    else:
        print: RBP name + "NOT FOUND in supplementary data"

call extract_logo_proportions()  → logo_props dict
print how many RBPs have logo data

print summary section:
    list any missing RBPs
    print instructions for getting comprehensive data (ENCODE, RBNS pipeline, etc.)

if top_motifs is not empty:
    create output directory: os.makedirs(OUTPUT_DIR, exist_ok=True)
    for each rbp, motifs in top_motifs.items():
        call write_csv_for_rbp(rbp, motifs, OUTPUT_DIR)
    print note about partial data
else:
    print that mmc4.xlsx was not found
```

### Key Patterns
```python
# Create a directory (no error if it already exists)
os.makedirs('data/rbns', exist_ok=True)

# Iterate over a dict's keys and values
for rbp, motifs in top_motifs.items():
    ...

# Safe access with a default
top_kmer = top_motifs[rbp][0][0] if top_motifs[rbp] else 'N/A'
top_z    = top_motifs[rbp][0][1] if top_motifs[rbp] else 0

# Fixed-width string formatting in f-strings (useful for aligned columns)
print(f"  {rbp:12} - {n_motifs:2} motifs (top: {top_kmer} Z={top_z:.2f})")
#        ^^^^                                                ^^^^
#   left-aligned, 12 chars wide                      2 decimal places
```

### The guard at the bottom
```python
if __name__ == '__main__':
    main()
```
This means `main()` only runs when you execute this file directly (`python scripts/extract_rbns_motifs.py`). If another script imports from this file, `main()` won't automatically run. All your other scripts should have this too.

---

## Recommended Build Order

1. Write the constants (`DATA_DIR`, `OUTPUT_DIR`, `TARGET_RBPS`)
2. Write `parse_motif_string()` — test it on a few example strings
3. Write `write_csv_for_rbp()` — simple file write, test it standalone
4. Write `extract_top_motifs_from_table_s3()` — test by printing the returned dict
5. Write `extract_logo_proportions()` — trickier due to the irregular layout; test by checking a known RBP
6. Write `main()` — assemble and test end-to-end

---

## Validation

```bash
# Run original
python scripts/extract_rbns_motifs_original.py

# Run yours
python scripts/extract_rbns_motifs.py

# Compare one output file
diff data/rbns/IGF2BP1_zscores.csv /tmp/backup_IGF2BP1_zscores.csv
```

---

## Common Mistakes to Watch For

| Mistake | Effect | Fix |
|---|---|---|
| Not wrapping `pd.read_excel()` in an existence check | Crashes with FileNotFoundError if mmc4.xlsx isn't present | Check `os.path.exists(xlsx_path)` first |
| Forgetting `str(col)` when checking column names | Crashes on numeric column names (pandas auto-names unnamed columns like `'Unnamed: 6'`) | Always wrap in `str()`: `if 'Motif' in str(col)` |
| Using `row['col']` instead of `row.get('col')` | Raises `KeyError` if the column doesn't exist | Use `.get()` for safety when accessing motif columns |
| Off-by-one in the logo proportions sheet | Gets the wrong k-mers or proportions | Data starts at row index 2; proportion is at `col_idx + 1` |
| Not converting `U` → `T` | K-mers remain as RNA (e.g. `UUUUU`) but downstream code uses DNA (e.g. `TTTTT`) | Always apply `.replace('U', 'T')` after extracting the k-mer |
| Forgetting `index=False` on `to_csv` | Not applicable here (using manual `write()`) — but relevant if you refactor to use pandas for output | Use `df.to_csv(path, index=False)` if switching to pandas |
