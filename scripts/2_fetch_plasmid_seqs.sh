#!/bin/bash
# Fetch plasmid nucleotide sequences from NCBI by accession, using EDirect.
#
# v4 changes (see v3 for the working-directory / visible-error fixes,
# both kept here):
#   - Switched from `epost -format acc | efetch` to `efetch -id acc1,acc2,...`
#     directly. The epost/esearch pipeline translates accessions to UIDs
#     via NCBI's [ACCN] field search, which was silently dropping a
#     consistent subset of accessions per batch (visible as retmax being
#     smaller than the batch size, before any fetch even happened) --
#     particularly affecting WGS-contig-style accessions. Calling
#     `efetch -id` directly fetches by accession with no translation step,
#     avoiding that failure mode and cutting out a network round-trip.
#
# Usage:
#   ./2_fetch_plasmids_seqs.sh [accessions.txt] data/all_plasmids.fasta [work_dir]
#
# Safe to re-run: appends to the output file, so you can run this against
# a "failed_accessions.txt" file to top up an existing FASTA.

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

if ! command -v efetch >/dev/null 2>&1; then
    echo "ERROR: efetch not found on PATH. Is EDirect installed and sourced?"
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

    # Build comma-separated accession list for this batch
    id_list=$(paste -sd, "$batch_file")

    attempt=1
    success=0
    result_file="$WORK_DIR/batch_result.fasta"
    err_file="$WORK_DIR/batch_err.log"

    while [ "$attempt" -le "$MAX_RETRIES" ]; do
        echo "Batch $batch_num/$n_batches, attempt $attempt..."

        efetch -db nuccore -id "$id_list" -format fasta > "$result_file" 2>"$err_file"
        exit_code=$?

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
        echo "  Batch $batch_num FAILED after $MAX_RETRIES attempts."
        echo "  Falling back to fetching this batch one accession at a time..."

        # Per-accession fallback: identifies exactly which accession(s)
        # in the batch are the problem, and still recovers the good ones.
        while IFS= read -r single_acc; do
            [ -z "$single_acc" ] && continue
            single_out="$WORK_DIR/single_result.fasta"
            efetch -db nuccore -id "$single_acc" -format fasta > "$single_out" 2>"$WORK_DIR/single_err.log"
            n_single=$(grep -c "^>" "$single_out" 2>/dev/null || echo 0)

            if [ "$n_single" -ge 1 ]; then
                cat "$single_out" >> "$OUTPUT_FASTA"
            else
                echo "    FAILED: $single_acc"
                sed 's/^/      /' "$WORK_DIR/single_err.log"
                echo "$single_acc" >> "$FAILED_LOG"
            fi
            sleep 0.34
        done < "$batch_file"
    fi

    rm -f "$batch_file"
    sleep 0.5
done

n_seqs=$(grep -c "^>" "$OUTPUT_FASTA" 2>/dev/null || echo 0)
echo ""
echo "Total sequences now in $OUTPUT_FASTA: $n_seqs (requested $total unique accessions)"

if [ -s "$FAILED_LOG" ]; then
    n_failed_lines=$(wc -l < "$FAILED_LOG")
    echo "$n_failed_lines accessions genuinely failed -- see $FAILED_LOG"
    echo "These are true problem accessions (not batch artifacts) -- worth checking individually,"
    echo "e.g. searching the accession directly at https://www.ncbi.nlm.nih.gov/nuccore/"
fi