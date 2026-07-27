"""
Parse tblastn results into:
  1. best_hits.tsv          -- one row per (gene, plasmid) pair: the best-
                                scoring HSP, UNFILTERED, for inspecting the
                                identity/coverage distribution before
                                deciding thresholds.
  2. presence_absence.csv   -- gene x plasmid binary matrix, using the
                                thresholds set below, merged back onto the
                                full original 839-row plasmid table (via
                                accession_mapping.csv) so every original
                                row is represented.
  3. identity_matrix.csv    -- same shape as (2), but with raw %identity
                                values instead of 1/0 (blank where no hit
                                passed the e-value cutoff at all), useful
                                for judging borderline cases.

Re-run just the matrix-building section (bottom) after changing
PIDENT_THRESHOLD / QCOV_THRESHOLD -- no need to re-parse BLAST output.
"""

import pandas as pd
import re

# ----------------------------------------------------------------
# INPUTS -- adjust paths as needed
# ----------------------------------------------------------------
BLAST_RESULTS = "output/tblastn_results.tsv"
ACCESSION_MAPPING = "input/accession_mapping.csv"   # from the extraction step

# ----------------------------------------------------------------
# THRESHOLDS -- the only two knobs you should need to turn
# ----------------------------------------------------------------
PIDENT_THRESHOLD = 80    # minimum % identity to call a gene "present"
QCOV_THRESHOLD = 80      # minimum % query coverage (qcovs) to call "present"

# ----------------------------------------------------------------
# Column names, matching the -outfmt string used in run_tblastn.sbatch
# ----------------------------------------------------------------
COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qlen", "slen", "qcovs", "qcovhsp", "stitle",
]

df = pd.read_csv(BLAST_RESULTS, sep="\t", names=COLS)
print(f"Loaded {len(df)} raw HSP rows from {BLAST_RESULTS}")
print(f"Genes: {sorted(df['qseqid'].unique())}")
print(f"Unique plasmid subjects hit: {df['sseqid'].nunique()}")

# --- Clean up sseqid: BLAST db built with -parse_seqids gives IDs like
#     "gb|CP012345.1|" -- extract the plain accession ---
def clean_accession(sseqid):
    m = re.search(r"\|([A-Za-z0-9_]+\.\d+)\|", sseqid)
    if m:
        return m.group(1)
    return sseqid  # fallback: leave as-is if format is unexpected

df["accession"] = df["sseqid"].apply(clean_accession)

# --- Collapse multiple HSPs per (gene, plasmid) pair to the single best
#     hit, by bitscore (standard convention for "best hit") ---
df_sorted = df.sort_values("bitscore", ascending=False)
best_hits = df_sorted.drop_duplicates(subset=["qseqid", "accession"], keep="first")
best_hits = best_hits.sort_values(["qseqid", "accession"]).reset_index(drop=True)

best_hits.to_csv("best_hits.tsv", sep="\t", index=False)
print(f"\nWrote {len(best_hits)} best-hit rows (one per gene x plasmid pair) to best_hits.tsv")

# --- Quick distribution summary to help pick thresholds ---
print("\n=== %identity distribution per gene (best hits only) ===")
print(best_hits.groupby("qseqid")["pident"].describe()[["min", "25%", "50%", "75%", "max", "count"]])

print("\n=== %query coverage (qcovs) distribution per gene (best hits only) ===")
print(best_hits.groupby("qseqid")["qcovs"].describe()[["min", "25%", "50%", "75%", "max", "count"]])

# ==================================================================
# From here down: build the presence/absence and %identity matrices
# using the thresholds set above. Re-run just this section after
# tweaking PIDENT_THRESHOLD / QCOV_THRESHOLD.
# ==================================================================

passing = best_hits[
    (best_hits["pident"] >= PIDENT_THRESHOLD) & (best_hits["qcovs"] >= QCOV_THRESHOLD)
]
print(f"\n{len(passing)}/{len(best_hits)} best-hit pairs pass "
      f"pident>={PIDENT_THRESHOLD} and qcovs>={QCOV_THRESHOLD}")

# Presence/absence: 1 if a passing hit exists for that (gene, accession) pair
presence = passing.pivot_table(
    index="accession", columns="qseqid", values="pident", aggfunc="max"
)
presence_binary = presence.notna().astype(int)

# %identity matrix: best hit's %identity regardless of pass/fail, blank if
# no hit at all passed the e-value cutoff (i.e. wasn't in best_hits at all)
identity_matrix = best_hits.pivot_table(
    index="accession", columns="qseqid", values="pident", aggfunc="max"
)

# --- Merge onto the full original table (all 839 rows, incl. duplicates) ---
mapping = pd.read_csv(ACCESSION_MAPPING)

acc_col = next((c for c in mapping.columns if "accession" in c.lower()), None)
if acc_col is None:
    raise ValueError(f"Could not find accession column in {ACCESSION_MAPPING}")

# base-accession join to tolerate version-suffix mismatches, same approach
# as the earlier verification scripts
mapping["_base_acc"] = mapping[acc_col].astype(str).str.split(".").str[0]
presence_binary_reset = presence_binary.reset_index()
presence_binary_reset["_base_acc"] = presence_binary_reset["accession"].str.split(".").str[0]

identity_matrix_reset = identity_matrix.reset_index()
identity_matrix_reset["_base_acc"] = identity_matrix_reset["accession"].str.split(".").str[0]
identity_matrix_reset = identity_matrix_reset.add_suffix("_pident")
identity_matrix_reset = identity_matrix_reset.rename(columns={"_base_acc_pident": "_base_acc"})

final_presence = mapping.merge(
    presence_binary_reset.drop(columns=["accession"]), on="_base_acc", how="left"
)
gene_cols = [c for c in presence_binary_reset.columns if c not in ("accession", "_base_acc")]
final_presence[gene_cols] = final_presence[gene_cols].fillna(0).astype(int)

final_identity = mapping.merge(
    identity_matrix_reset.drop(columns=["accession_pident"], errors="ignore"),
    on="_base_acc", how="left"
)

final_presence.drop(columns=["_base_acc"]).to_csv("presence_absence.csv", index=False)
final_identity.drop(columns=["_base_acc"]).to_csv("identity_matrix.csv", index=False)

print(f"\nWrote presence_absence.csv and identity_matrix.csv "
      f"({len(final_presence)} rows, matching original table row count)")
print(f"\nGene presence counts (current thresholds: pident>={PIDENT_THRESHOLD}, qcovs>={QCOV_THRESHOLD}):")
print(final_presence[gene_cols].sum())