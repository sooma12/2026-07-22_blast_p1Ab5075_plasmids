#!/bin/bash
# Fetch plasmid nucleotide sequences from NCBI by accession, using EDirect.
#
# v3 changes:
#   - Uses a local working directory (./fetch_work by default) instead of
#     /tmp, created explicitly with mkdir -p. Some environments have a
#     read-only, non-persistent, or restricted /tmp, which can cause
#     silent failures.
#   - Errors are printed directly to the terminal (not just logged to a
#     file), so failures are visible immediately.
#   - set -e is intentionally NOT used (we handle per-batch failures
#     ourselves); but every step now checks its own exit status explicitly
#     so failures can't pass silently.
#
# Usage:
#   bash scripts/2_fetch_plasmid_seqs.sh input/accessions.txt data/all_plasmids.fasta [work_dir]
#
# Safe to re-run: appends to the output file, so you can run this against
# a "missing_accessions.txt" or "failed_accessions.txt" file to top up an
# existing FASTA.

ACCESSION_FILE="${1:-accessions.txt}"
OUTPUT_FASTA="${2:-all_plasmids.fasta}"
WORK_DIR="${3:-./fetch_work}"
BATCH_SIZE=25
MAX_RETRIES=4
SLEEP_BASE=2      # seconds; doubles each retry (2, 4, 8, 16)

if [ ! -f "$ACCESSION_FILE" ]; then
    echo "Accession file not found: $ACCESSION_FILE"
    exit 1
fi

# --- Set up working directory explicitly, and confirm it's actually writable ---
mkdir -p "$WORK_DIR"
if [ $? -ne 0 ]; then
    echo "ERROR: could not create working directory: $WORK_DIR"
    exit 1
fi
if [ ! -w "$WORK_DIR" ]; then
    echo "ERROR: working directory is not writable: $WORK_DIR"
    exit 1
fi
echo "Using working directory: $WORK_DIR"

# Check that epost/efetch are actually on PATH before we start
if ! command -v epost >/dev/null 2>&1 || ! command -v efetch >/dev/null 2>&1; then
    echo "ERROR: epost/efetch not found on PATH. Is EDirect installed and sourced?"
    echo "  (e.g. run: source ~/.bashrc  -- or re-check the EDirect install instructions)"
    exit 1
fi

# Optional: export NCBI_API_KEY="your_key_here"

touch "$OUTPUT_FASTA"
FAILED_LOG="$WORK_DIR/failed_accessions.txt"
> "$FAILED_LOG"

total=$(wc -l < "$ACCESSION_FILE")
echo "Fetching $total accessions in batches of $BATCH_SIZE, up to $MAX_RETRIES retries per batch..."

rm -f "$WORK_DIR"/acc_batch_*
split -l "$BATCH_SIZE" "$ACCESSION_FILE" "$WORK_DIR/acc_batch_"
if [ $? -ne 0 ]; then
    echo "ERROR: 'split' command failed -- check that $WORK_DIR is writable and has space."
    exit 1
fi

n_batches=$(ls "$WORK_DIR"/acc_batch_* 2>/dev/null | wc -l)
if [ "$n_batches" -eq 0 ]; then
    echo "ERROR: no batch files were created in $WORK_DIR -- 'split' likely failed silently."
    exit 1
fi

batch_num=0
for batch_file in "$WORK_DIR"/acc_batch_*; do
    batch_num=$((batch_num + 1))
    n_expected=$(wc -l < "$batch_file")

    attempt=1
    success=0
    result_file="$WORK_DIR/batch_result.fasta"
    err_file="$WORK_DIR/batch_err.log"

    while [ "$attempt" -le "$MAX_RETRIES" ]; do
        echo "Batch $batch_num/$n_batches, attempt $attempt..."

        epost -db nuccore -input "$batch_file" -format acc 2>"$WORK_DIR/epost_err.log" \
            | efetch -format fasta > "$result_file" 2>"$err_file"
        exit_code=$?

        if [ -s "$WORK_DIR/epost_err.log" ]; then
            echo "  epost stderr:"
            sed 's/^/    /' "$WORK_DIR/epost_err.log"
        fi
        if [ -s "$err_file" ]; then
            echo "  efetch stderr:"
            sed 's/^/    /' "$err_file"
        fi

        if [ "$exit_code" -eq 0 ]; then
            n_got=$(grep -c "^>" "$result_file" 2>/dev/null || echo 0)

            if [ "$n_got" -eq "$n_expected" ]; then
                cat "$result_file" >> "$OUTPUT_FASTA"
                echo "  OK: got $n_got/$n_expected sequences"
                success=1
                break
            else
                echo "  Partial result: got $n_got/$n_expected -- retrying batch"
            fi
        else
            echo "  Command exited with code $exit_code -- retrying"
        fi

        sleep_time=$((SLEEP_BASE * (2 ** (attempt - 1))))
        echo "  waiting ${sleep_time}s before retry..."
        sleep "$sleep_time"
        attempt=$((attempt + 1))
    done

    if [ "$success" -ne 1 ]; then
        echo "  Batch $batch_num FAILED after $MAX_RETRIES attempts -- logging accessions"
        cat "$batch_file" >> "$FAILED_LOG"
        # Keep any partial result from the last attempt rather than discard it
        if [ -s "$result_file" ]; then
            cat "$result_file" >> "$OUTPUT_FASTA"
        fi
    fi

    rm -f "$batch_file"
    sleep 0.5
done

n_seqs=$(grep -c "^>" "$OUTPUT_FASTA" 2>/dev/null || echo 0)
echo ""
echo "Total sequences now in $OUTPUT_FASTA: $n_seqs (requested $total unique accessions)"

if [ -s "$FAILED_LOG" ]; then
    n_failed_lines=$(wc -l < "$FAILED_LOG")
    echo "$n_failed_lines accessions in failed batches -- see $FAILED_LOG"
    echo "Re-run this script pointed at $FAILED_LOG to retry just those."
fi