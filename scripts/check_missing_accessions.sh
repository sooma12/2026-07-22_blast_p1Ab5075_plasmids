#!/bin/bash
# Compare accessions.txt against what's actually in all_plasmids.fasta,
# and write the missing ones to a new file for a resume/retry pass.
#
# Usage:
#   ./check_missing_accessions.sh accessions.txt all_plasmids.fasta missing_accessions.txt

set -euo pipefail

ACCESSION_FILE="${1:-accessions.txt}"
FASTA_FILE="${2:-all_plasmids.fasta}"
MISSING_FILE="${3:-missing_accessions.txt}"

# Extract accessions actually present in the fasta (first token after '>',
# strip any version-less/version-mismatched comparison by matching on
# the accession prefix before checking).
grep "^>" "$FASTA_FILE" | sed 's/^>//' | awk '{print $1}' | sort -u > /tmp/present_ids.txt

sort -u "$ACCESSION_FILE" > /tmp/requested_ids.txt

# Compare ignoring version suffix mismatches (e.g. requested CP012345.1,
# but NCBI returned CP012345.2) -- strip trailing ".N" for comparison,
# then map back to the original requested accession.
awk -F'.' '{print $1}' /tmp/present_ids.txt | sort -u > /tmp/present_noversion.txt

> "$MISSING_FILE"
n_present=0
n_missing=0
while IFS= read -r acc; do
    base="${acc%%.*}"
    if grep -qx "$base" /tmp/present_noversion.txt; then
        n_present=$((n_present + 1))
    else
        echo "$acc" >> "$MISSING_FILE"
        n_missing=$((n_missing + 1))
    fi
done < /tmp/requested_ids.txt

total=$(wc -l < /tmp/requested_ids.txt)
echo "Requested: $total"
echo "Found in FASTA: $n_present"
echo "Missing: $n_missing (written to $MISSING_FILE)"