#!/bin/bash
# Build a nucleotide BLAST database from the plasmid sequences, for use
# with tblastn (protein query vs. translated nucleotide database).
#
# Usage:
#   ./make_blast_db.sh all_plasmids.fasta plasmid_db

set -euo pipefail

INPUT_FASTA="${1:-all_plasmids.fasta}"
DB_NAME="${2:-plasmid_db}"
DB_DIR="$(dirname "$DB_NAME")"

if [ ! -f "$INPUT_FASTA" ]; then
    echo "Input FASTA not found: $INPUT_FASTA"
    exit 1
fi

if ! command -v makeblastdb >/dev/null 2>&1; then
    echo "ERROR: makeblastdb not found on PATH. Is BLAST+ installed?"
    echo "  (e.g. conda install -c bioconda blast, or module load blast on an HPC system)"
    exit 1
fi

[ -n "$DB_DIR" ] && [ "$DB_DIR" != "." ] && mkdir -p "$DB_DIR"

n_seqs=$(grep -c "^>" "$INPUT_FASTA" || echo 0)
echo "Building nucleotide BLAST database from $n_seqs sequences in $INPUT_FASTA..."

makeblastdb \
    -in "$INPUT_FASTA" \
    -dbtype nucl \
    -out "$DB_NAME" \
    -parse_seqids \
    -title "Plasmid database (${n_seqs} sequences)"

echo ""
echo "Database built: $DB_NAME"
echo ""
echo "Summary:"
blastdbcmd -db "$DB_NAME" -info