"""
Extract plasmid accessions from TableS2_plasmids_curated_for_publication_310823.xlsx
(sheet: Final_plasmid_set, column: "Plasmid accession").

Writes two outputs:
  1. accessions.txt              -- DEDUPLICATED list, one accession per line.
                                     Use this for fetching (avoids the
                                     epost/efetch batch shortfall caused by
                                     posting the same accession twice in one
                                     batch, which NCBI's history server
                                     silently collapses).
  2. accession_mapping.csv       -- every original row from the table, with
                                     its accession, so you can re-join back
                                     to full metadata (and to every duplicate
                                     row) after fetching/BLASTing.
"""

import pandas as pd

INPUT_XLSX = "input/TableS2_plasmids_curated_for_publication_310823.xlsx"
SHEET_NAME = "Final_plasmid_set"
ACCESSION_COLUMN = "Plasmid accession"
ACCESSIONS_OUTPUT = "input/accessions.txt"
MAPPING_OUTPUT = "input/accession_mapping.csv"

print("Beginning plasmid accession extraction script")

df = pd.read_excel(INPUT_XLSX, sheet_name=SHEET_NAME)

print(f"Loaded sheet '{SHEET_NAME}' with {len(df)} rows and columns:")
print(list(df.columns))

if ACCESSION_COLUMN not in df.columns:
    raise ValueError(
        f"Column '{ACCESSION_COLUMN}' not found. "
        f"Available columns: {list(df.columns)}"
    )

# Clean up the accession column but keep every row (including duplicates)
# in the mapping table -- we only deduplicate the list we send to NCBI.
df[ACCESSION_COLUMN] = df[ACCESSION_COLUMN].astype(str).str.strip()
df = df[df[ACCESSION_COLUMN].notna() & (df[ACCESSION_COLUMN] != "") & (df[ACCESSION_COLUMN] != "nan")]

print(f"\n{len(df)} rows have a non-blank accession.")

# --- Report duplicates so you can eyeball whether they're a real biological
#     duplicate (e.g. resequenced/reassembled) or a table entry issue ---
dup_mask = df[ACCESSION_COLUMN].duplicated(keep=False)
if dup_mask.any():
    dup_accessions = sorted(df.loc[dup_mask, ACCESSION_COLUMN].unique())
    print(f"\n{len(dup_accessions)} accession(s) appear on more than one row:")
    for acc in dup_accessions:
        rows = df.index[df[ACCESSION_COLUMN] == acc].tolist()
        print(f"  {acc}: rows {rows}")
    print("(Worth a quick manual check -- see accession_mapping.csv for full row detail.)")

# --- Write the full mapping table (all rows preserved, duplicates included) ---
df.to_csv(MAPPING_OUTPUT, index_label="original_row_index")
print(f"\nWrote full row -> accession mapping to {MAPPING_OUTPUT}")

# --- Write the deduplicated accession list for fetching ---
unique_accessions = df[ACCESSION_COLUMN].drop_duplicates()
with open(ACCESSIONS_OUTPUT, "w") as f:
    f.write("\n".join(unique_accessions))

print(f"Wrote {len(unique_accessions)} unique accessions to {ACCESSIONS_OUTPUT} "
      f"(deduplicated from {len(df)} rows)")
print("First few:", unique_accessions.head(5).tolist())