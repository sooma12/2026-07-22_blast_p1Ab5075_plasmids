"""
Extract plasmid accessions from TableS2_plasmids_curated_for_publication_310823.xlsx
(sheet: Final_plasmid_set, column: "Plasmid accession") into a plain text file,
one accession per line -- ready for fetch_plasmids.py.
"""

import pandas as pd

INPUT_XLSX = "input/TableS2_plasmids_curated_for_publication_310823.xlsx"
SHEET_NAME = "Final_plasmid_set"
ACCESSION_COLUMN = "Plasmid accession"
OUTPUT_FILE = "input/accessions.txt"

print('script is working, starting to read excel file')

df = pd.read_excel(INPUT_XLSX, sheet_name=SHEET_NAME)

print(f"Loaded sheet '{SHEET_NAME}' with {len(df)} rows and columns:")
print(list(df.columns))

if ACCESSION_COLUMN not in df.columns:
    raise ValueError(
        f"Column '{ACCESSION_COLUMN}' not found. "
        f"Available columns: {list(df.columns)}"
    )

accessions = df[ACCESSION_COLUMN].dropna().astype(str).str.strip()

# drop any blank strings and obvious placeholder rows
accessions = accessions[accessions != ""]

# check for duplicates -- worth knowing about, not necessarily an error
n_dupes = accessions.duplicated().sum()
if n_dupes:
    print(f"Warning: {n_dupes} duplicate accessions found (keeping all, deduplicate if unwanted).")

with open(OUTPUT_FILE, "w") as f:
    f.write("\n".join(accessions))

print(f"\nWrote {len(accessions)} accessions to {OUTPUT_FILE}")
print("First few:", accessions.head(5).tolist())
