import pandas as pd
import matplotlib.pyplot as plt

# Read files
blastn = pd.read_csv("blastn_identity_matrix.csv")
tblastn = pd.read_csv("identity_matrix.csv")

# Metadata columns used to identify rows
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

# Rename tblastn columns
rename = {
    "AKA33575.1_pident": "ABUW_4004_pident",
    "AKA33576.1_pident": "ABUW_4005_pident",
    "AKA33577.1_pident": "ABUW_4006_pident",
    "AKA33578.1_pident": "ABUW_4007_pident",
    "AKA33651.1_pident": "ABUW_4094_pident",

    "AKA33575.1_qcovs": "ABUW_4004_qcovs",
    "AKA33576.1_qcovs": "ABUW_4005_qcovs",
    "AKA33577.1_qcovs": "ABUW_4006_qcovs",
    "AKA33578.1_qcovs": "ABUW_4007_qcovs",
    "AKA33651.1_qcovs": "ABUW_4094_qcovs",
}

tblastn = tblastn.rename(columns=rename)

merged = tblastn.merge(
    blastn,
    on=key_cols,
    suffixes=("_tblastn", "_blastn"),
    validate="one_to_one"
)
comparison_rows = []

genes = [
    "ABUW_4004",
    "ABUW_4005",
    "ABUW_4006",
    "ABUW_4007",
    "ABUW_4094",
]

for _, row in merged.iterrows():
    for gene in genes:

        comparison_rows.append({
            "original_row_index": row["original_row_index"],
            "Plasmid accession": row["Plasmid accession"],
            "Plasmid name": row["Plasmid name"],
            "strain name": row["strain name"],

            "Gene": gene,

            "tblastn_pident": row[f"{gene}_pident_tblastn"],
            "blastn_pident": row[f"{gene}_pident_blastn"],
            "delta_pident":
                row[f"{gene}_pident_blastn"] -
                row[f"{gene}_pident_tblastn"],

            "tblastn_qcovs": row[f"{gene}_qcovs_tblastn"],
            "blastn_qcovs": row[f"{gene}_qcovs_blastn"],
            "delta_qcovs":
                row[f"{gene}_qcovs_blastn"] -
                row[f"{gene}_qcovs_tblastn"],
        })



comparison = pd.DataFrame(comparison_rows)
comparison["abs_delta_pident"] = comparison["delta_pident"].abs()
comparison["abs_delta_qcovs"] = comparison["delta_qcovs"].abs()

comparison = comparison.sort_values(
    ["abs_delta_pident", "abs_delta_qcovs"],
    ascending=False,
)
comparison.to_csv("blast_vs_tblastn_values.csv", index=False)

print(f"Mean |Δ pident|: {comparison['abs_delta_pident'].mean():.3f}")
print(f"Max  |Δ pident|: {comparison['abs_delta_pident'].max():.3f}")

print(f"Mean |Δ qcovs|: {comparison['abs_delta_qcovs'].mean():.3f}")
print(f"Max  |Δ qcovs|: {comparison['abs_delta_qcovs'].max():.3f}")
