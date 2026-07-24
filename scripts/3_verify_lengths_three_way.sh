#!/bin/bash
# Three-way verification of plasmid sequence lengths:
#   1. Actual length of the sequence in your downloaded FASTA
#   2. Length NCBI reports for that accession (via esummary)
#   3. Length reported in the paper's table (accession_mapping.csv,
#      auto-detects a column containing "length", e.g. "Length (bp)")
#
# Flags any accession where these three don't all agree.
#
# Usage:
#   ./verify_lengths_three_way.sh accessions.txt all_plasmids.fasta accession_mapping.csv

ACCESSION_FILE="${1:-accessions.txt}"
FASTA_FILE="${2:-all_plasmids.fasta}"
MAPPING_CSV="${3:-accession_mapping.csv}"
WORK_DIR="./length_check_work"
BATCH_SIZE=25
REPORT_FILE="three_way_length_report.tsv"

mkdir -p "$WORK_DIR"

if ! command -v esummary >/dev/null 2>&1; then
    echo "ERROR: esummary not found on PATH. Is EDirect installed and sourced?"
    exit 1
fi

# --- Step 1: actual lengths from local FASTA ---
echo "Computing actual sequence lengths from $FASTA_FILE..."
awk '
    /^>/ {
        if (acc != "") print acc "\t" len
        split(substr($0, 2), a, " ")
        acc = a[1]
        len = 0
        next
    }
    {
        len += length($0)
    }
    END {
        if (acc != "") print acc "\t" len
    }
' "$FASTA_FILE" > "$WORK_DIR/actual_lengths.tsv"
echo "  Found $(wc -l < "$WORK_DIR/actual_lengths.tsv") sequences."

# --- Step 2: expected lengths from NCBI esummary, in batches ---
echo "Querying NCBI for expected lengths..."
> "$WORK_DIR/expected_lengths.tsv"
rm -f "$WORK_DIR"/lc_batch_*
split -l "$BATCH_SIZE" "$ACCESSION_FILE" "$WORK_DIR/lc_batch_"

# --- Step 2: expected lengths from NCBI esummary, in batches, with retries ---
echo "Querying NCBI for expected lengths..."
> "$WORK_DIR/expected_lengths.tsv"
rm -f "$WORK_DIR"/lc_batch_*
split -l "$BATCH_SIZE" "$ACCESSION_FILE" "$WORK_DIR/lc_batch_"

MAX_RETRIES=4
SLEEP_BASE=2

batch_num=0
n_batches=$(ls "$WORK_DIR"/lc_batch_* 2>/dev/null | wc -l)
for batch_file in "$WORK_DIR"/lc_batch_*; do
    batch_num=$((batch_num + 1))
    n_expected_in_batch=$(wc -l < "$batch_file")
    id_list=$(paste -sd, "$batch_file")

    attempt=1
    success=0
    batch_out="$WORK_DIR/esummary_batch_result.tsv"

    while [ "$attempt" -le "$MAX_RETRIES" ]; do
        echo "  Batch $batch_num/$n_batches, attempt $attempt..."
        esummary -db nuccore -id "$id_list" 2>"$WORK_DIR/esummary_err.log" \
            | xtract -pattern DocumentSummary -element Caption,Slen \
            > "$batch_out"

        if [ -s "$WORK_DIR/esummary_err.log" ]; then
            echo "    esummary stderr:"
            sed 's/^/      /' "$WORK_DIR/esummary_err.log"
        fi

        n_got=$(wc -l < "$batch_out")
        if [ "$n_got" -eq "$n_expected_in_batch" ]; then
            cat "$batch_out" >> "$WORK_DIR/expected_lengths.tsv"
            echo "    OK: got $n_got/$n_expected_in_batch summaries"
            success=1
            break
        else
            echo "    Partial/failed: got $n_got/$n_expected_in_batch -- retrying"
        fi

        sleep_time=$((SLEEP_BASE * (2 ** (attempt - 1))))
        sleep "$sleep_time"
        attempt=$((attempt + 1))
    done

    if [ "$success" -ne 1 ]; then
        echo "  Batch $batch_num failed after $MAX_RETRIES attempts -- falling back to one-at-a-time"
        while IFS= read -r single_acc; do
            [ -z "$single_acc" ] && continue
            single_out="$WORK_DIR/esummary_single.tsv"
            esummary -db nuccore -id "$single_acc" 2>"$WORK_DIR/esummary_single_err.log" \
                | xtract -pattern DocumentSummary -element Caption,Slen \
                > "$single_out"
            if [ -s "$single_out" ]; then
                cat "$single_out" >> "$WORK_DIR/expected_lengths.tsv"
            else
                echo "    STILL FAILED for accession: $single_acc"
                sed 's/^/      /' "$WORK_DIR/esummary_single_err.log"
            fi
            sleep 0.34
        done < "$batch_file"
    fi

    rm -f "$batch_file"
    sleep 0.4
done
echo "  Got expected lengths for $(wc -l < "$WORK_DIR/expected_lengths.tsv") accessions."

# --- Step 3: three-way comparison, including the paper's table ---
echo "Comparing FASTA vs NCBI vs paper's table..."

python3 << PYEOF
import pandas as pd

def base(acc):
    return str(acc).split(".")[0]

# --- actual (from FASTA) ---
actual = {}
with open("$WORK_DIR/actual_lengths.tsv") as f:
    for line in f:
        acc, length = line.strip().split("\t")
        actual[base(acc)] = (acc, int(length))

# --- expected (from NCBI) ---
expected = {}
with open("$WORK_DIR/expected_lengths.tsv") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) != 2:
            continue
        acc, length = parts
        expected[base(acc)] = (acc, int(length))

# --- table (from the paper's Excel, via accession_mapping.csv) ---
mapping = pd.read_csv("$MAPPING_CSV")

# auto-detect the accession column and a length column
acc_col = None
for c in mapping.columns:
    if "accession" in c.lower():
        acc_col = c
        break
if acc_col is None:
    raise ValueError(f"Could not find an accession column in $MAPPING_CSV. Columns: {list(mapping.columns)}")

length_col = None
for c in mapping.columns:
    if "length" in c.lower():
        length_col = c
        break
if length_col is None:
    raise ValueError(f"Could not find a length column in $MAPPING_CSV. Columns: {list(mapping.columns)}")

print(f"Using accession column: '{acc_col}', length column: '{length_col}'")

table_lengths = {}
table_conflicts = {}
for _, row in mapping.iterrows():
    acc = row[acc_col]
    if pd.isna(acc):
        continue
    b = base(acc)
    try:
        length_val = int(row[length_col])
    except (ValueError, TypeError):
        continue
    if b in table_lengths and table_lengths[b][1] != length_val:
        table_conflicts.setdefault(b, set()).add(table_lengths[b][1])
        table_conflicts[b].add(length_val)
    table_lengths[b] = (acc, length_val)

all_bases = set(actual) | set(expected) | set(table_lengths)
rows = []
n_all_match = 0

for b in sorted(all_bases):
    a = actual.get(b)
    e = expected.get(b)
    t = table_lengths.get(b)

    a_len = a[1] if a else None
    e_len = e[1] if e else None
    t_len = t[1] if t else None

    lengths_present = [x for x in (a_len, e_len, t_len) if x is not None]
    all_agree = len(set(lengths_present)) <= 1 and len(lengths_present) == 3

    if all_agree:
        n_all_match += 1
        continue  # only report problems

    status_parts = []
    if a_len is None:
        status_parts.append("MISSING_FROM_FASTA")
    if e_len is None:
        status_parts.append("MISSING_FROM_NCBI_SUMMARY")
    if t_len is None:
        status_parts.append("MISSING_FROM_TABLE")
    if b in table_conflicts:
        status_parts.append("TABLE_HAS_CONFLICTING_VALUES")
    if a_len is not None and e_len is not None and a_len != e_len:
        status_parts.append("FASTA_VS_NCBI_MISMATCH")
    if a_len is not None and t_len is not None and a_len != t_len:
        status_parts.append("FASTA_VS_TABLE_MISMATCH")
    if e_len is not None and t_len is not None and e_len != t_len:
        status_parts.append("NCBI_VS_TABLE_MISMATCH")

    rows.append((b, ";".join(status_parts), a_len, e_len, t_len))

with open("$REPORT_FILE", "w") as f:
    f.write("accession\tissues\tfasta_length\tncbi_length\ttable_length\n")
    for r in rows:
        f.write("\t".join("" if x is None else str(x) for x in r) + "\n")

print(f"\n{n_all_match} accession(s) match across all three sources (FASTA, NCBI, table).")
print(f"{len(rows)} accession(s) have a discrepancy -- see $REPORT_FILE")
PYEOF