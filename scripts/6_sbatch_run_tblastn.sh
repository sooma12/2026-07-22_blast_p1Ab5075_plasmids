#!/bin/bash
#SBATCH --job-name=tblastn_plasmid_search
#SBATCH --output=logs/tblastn_%j.out
#SBATCH --error=logs/tblastn_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --partition=short

# ------------------------------------------------------------------
# Runs tblastn for all query genes (protein) against the plasmid
# nucleotide database, with NO identity/coverage filtering applied --
# outputs every hit above the e-value cutoff so thresholds can be
# chosen later by inspecting the full results.
# ------------------------------------------------------------------

set -euo pipefail

mkdir -p logs output

QUERY_FASTA="data/genes_protein.fasta"
DB_NAME="plasmid_db"
OUTPUT_TSV="output/tblastn_results.tsv"
EVALUE="1e-10"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

# Load BLAST+ if your cluster uses environment modules -- adjust module
# name/version as needed for your system, or comment out if BLAST+ is
# already on PATH (e.g. via conda environment).
module load blast+ 2>/dev/null || module load blast 2>/dev/null || true

if ! command -v tblastn >/dev/null 2>&1; then
    echo "ERROR: tblastn not found on PATH. Check module load / conda env activation."
    exit 1
fi

echo "Running tblastn with $THREADS threads..."
echo "Query: $QUERY_FASTA"
echo "Database: $DB_NAME"

# Column choices: includes qcovs/qcovhsp (BLAST's own built-in coverage
# calculations) so you don't need to compute coverage manually later,
# plus stitle for a human-readable description of each hit subject.
tblastn \
    -query "$QUERY_FASTA" \
    -db "$DB_NAME" \
    -out "$OUTPUT_TSV" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp stitle" \
    -evalue "$EVALUE" \
    -max_target_seqs 1000 \
    -num_threads "$THREADS"

echo "Done. Results written to $OUTPUT_TSV"
echo "Row count: $(wc -l < "$OUTPUT_TSV")"
