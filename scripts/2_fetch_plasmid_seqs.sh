#!/bin/bash
# Fetch plasmid nucleotide sequences from NCBI by accession, using EDirect.
# Improved version: smaller batches + automatic retry on transient failures
# (dropped SSL connections, timeouts, etc).
#
# Usage:
#   ./fetch_plasmids_edirect_v2.sh accessions.txt all_plasmids.fasta
#
# Safe to re-run: appends to the output file, so you can run this against
# a "missing_accessions.txt" file to top up an existing FASTA.

set -uo pipefail   # NOTE: not using -e here, since we want to handle
                   # per-batch failures ourselves without killing the script

ACCESSION_FILE="${1:-accessions.txt}"
OUTPUT_FASTA="${2:-all_plasmids.fasta}"
BATCH_SIZE=25
MAX_RETRIES=4
SLEEP_BASE=2      # seconds; doubles each retry (2, 4, 8, 16)

if [ ! -f "$ACCESSION_FILE" ]; then
    echo "Accession file not found: $ACCESSION_FILE"
    exit 1
fi

# Optional: export NCBI_API_KEY="your_key_here"  # once the login issue is fixed

touch "$OUTPUT_FASTA"
> failed_accessions.txt

total=$(wc -l < "$ACCESSION_FILE")
echo "Fetching $total accessions in batches of $BATCH_SIZE, up to $MAX_RETRIES retries per batch..."

rm -f /tmp/acc_batch_*
split -l "$BATCH_SIZE" "$ACCESSION_FILE" /tmp/acc_batch_

batch_num=0
n_batches=$(ls /tmp/acc_batch_* | wc -l)

for batch_file in /tmp/acc_batch_*; do
    batch_num=$((batch_num + 1))
    n_expected=$(wc -l < "$batch_file")

    attempt=1
    success=0

    while [ "$attempt" -le "$MAX_RETRIES" ]; do
        echo "Batch $batch_num/$n_batches, attempt $attempt..."

        # Fetch into a temp file first so we can verify before appending
        if epost -db nuccore -input "$batch_file" -format acc 2>/tmp/epost_err.log \
            | efetch -format fasta > /tmp/batch_result.fasta 2>/tmp/efetch_err.log; then

            n_got=$(grep -c "^>" /tmp/batch_result.fasta || true)

            if [ "$n_got" -eq "$n_expected" ]; then
                cat /tmp/batch_result.fasta >> "$OUTPUT_FASTA"
                echo "  OK: got $n_got/$n_expected sequences"
                success=1
                break
            else
                echo "  Partial result: got $n_got/$n_expected -- retrying batch"
            fi
        else
            echo "  Command failed (see /tmp/efetch_err.log) -- retrying"
        fi

        sleep_time=$((SLEEP_BASE * (2 ** (attempt - 1))))
        echo "  waiting ${sleep_time}s before retry..."
        sleep "$sleep_time"
        attempt=$((attempt + 1))
    done

    if [ "$success" -ne 1 ]; then
        echo "  Batch $batch_num FAILED after $MAX_RETRIES attempts -- logging accessions"
        cat "$batch_file" >> failed_accessions.txt
        # Still keep any partial result we got on the last attempt, better than nothing
        if [ -s /tmp/batch_result.fasta ]; then
            cat /tmp/batch_result.fasta >> "$OUTPUT_FASTA"
        fi
    fi

    rm -f "$batch_file"
    sleep 0.5   # small courtesy pause between batches regardless of outcome
done

n_seqs=$(grep -c "^>" "$OUTPUT_FASTA" || true)
echo ""
echo "Total sequences now in $OUTPUT_FASTA: $n_seqs"

if [ -s failed_accessions.txt ]; then
    n_failed_lines=$(wc -l < failed_accessions.txt)
    echo "$n_failed_lines accessions in failed batches -- see failed_accessions.txt"
    echo "Re-run this script pointed at failed_accessions.txt to retry just those."
fi
