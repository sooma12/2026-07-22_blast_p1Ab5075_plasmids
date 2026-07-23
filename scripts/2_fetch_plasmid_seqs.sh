#!/bin/bash
# Fetch plasmid nucleotide sequences from NCBI by accession, using EDirect.
#
# Requires EDirect installed:
#   sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
# (adds `epost`/`efetch` to your PATH)
#
# Usage:
#   ./fetch_plasmids_edirect.sh accessions.txt all_plasmids.fasta

set -euo pipefail

ACCESSION_FILE="${1:-accessions.txt}"
OUTPUT_FASTA="${2:-all_plasmids.fasta}"
BATCH_SIZE=200

if [ ! -f "$ACCESSION_FILE" ]; then
    echo "Accession file not found: $ACCESSION_FILE"
    exit 1
fi

# Optional: set your NCBI API key to raise the rate limit (3/sec -> 10/sec)
# export NCBI_API_KEY="your_key_here"

> "$OUTPUT_FASTA"           # truncate/create output file
> failed_accessions.txt     # truncate/create failure log

total=$(wc -l < "$ACCESSION_FILE")
echo "Fetching $total accessions in batches of $BATCH_SIZE..."

i=0
split -l "$BATCH_SIZE" "$ACCESSION_FILE" /tmp/acc_batch_

for batch_file in /tmp/acc_batch_*; do
    i=$((i + BATCH_SIZE))
    echo "  Fetching batch (up to record $i of $total)..."

    if epost -db nuccore -input "$batch_file" -format acc \
        | efetch -format fasta >> "$OUTPUT_FASTA" 2>/tmp/efetch_err.log; then
        :
    else
        echo "  Batch failed -- see /tmp/efetch_err.log; logging accessions"
        cat "$batch_file" >> failed_accessions.txt
    fi

    rm -f "$batch_file"
    sleep 0.4
done

n_seqs=$(grep -c "^>" "$OUTPUT_FASTA" || true)
echo ""
echo "Wrote $n_seqs sequences to $OUTPUT_FASTA"

if [ -s failed_accessions.txt ]; then
    echo "Some accessions failed -- see failed_accessions.txt"
fi

# Sanity check: compare requested vs retrieved count
if [ "$n_seqs" -ne "$total" ]; then
    echo "WARNING: requested $total accessions but got $n_seqs sequences."
    echo "This can happen with withdrawn/replaced accessions -- check manually."
fi