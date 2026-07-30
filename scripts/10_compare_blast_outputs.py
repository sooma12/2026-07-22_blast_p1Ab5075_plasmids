import pandas as pd

# Input files
tblastn_file = "identity_matrix.csv"
blastn_file = "blastn_identity_matrix.csv"

# Read spreadsheets
tblastn = pd.read_csv(tblastn_file)
blastn = pd.read_csv(blastn_file)

# Columns that identify a row
key_cols = [
    "original_row_index",
    "Plasmid accession",
    "Plasmid name",
    "unique plasmid numerical",
    "strain name",
    "unique strain numerical",
    "Biosample",
    "Data source",
    "ST (clone)",
    "Rep_family",
    "Type",
    "Exact rep sequence",
    "Rep (Bertini scheme)",
    "Region",
    "Country",
]

# Rename tblastn columns to match blastn
rename_dict = {
    "AKA33575.1": "ABUW_4004",
    "AKA33576.1": "ABUW_4005",
    "AKA33577.1": "ABUW_4006",
    "AKA33578.1": "ABUW_4007",
    "AKA33651.1": "ABUW_4094",
}

tblastn = tblastn.rename(columns=rename_dict)

gene_cols = [
    "ABUW_4004",
    "ABUW_4005",
    "ABUW_4006",
    "ABUW_4007",
    "ABUW_4094",
]

# Merge the two tables
merged = tblastn.merge(
    blastn,
    on=key_cols,
    how="outer",
    suffixes=("_tblastn", "_blastn"),
    indicator=True,
)

print(merged.columns.tolist()) # TODO
exit()
# Report rows present in only one file
only_left = merged[merged["_merge"] == "left_only"]
only_right = merged[merged["_merge"] == "right_only"]

print(f"Rows only in tblastn: {len(only_left)}")
print(f"Rows only in blastn: {len(only_right)}")

# Compare gene presence/absence
comparison = pd.DataFrame(index=merged.index)

for gene in gene_cols:
    comparison[gene] = (
        merged[f"{gene}_tblastn"] == merged[f"{gene}_blastn"]
    )

comparison["all_match"] = comparison.all(axis=1)

mismatches = merged.loc[~comparison["all_match"]].copy()

# Add columns showing which genes differ
for gene in gene_cols:
    mismatches[f"{gene}_match"] = comparison.loc[mismatches.index, gene]

print(f"Rows with gene mismatches: {len(mismatches)}")

# Save results
mismatches.to_excel("blast_comparison_mismatches.xlsx", index=False)

print("Mismatch report written to blast_comparison_mismatches.xlsx")