#!/bin/bash
# Fetch protein sequences for a list of locus tags from NCBI.
#
# Primary method: search the protein database directly by [Locus Tag] field.
# Fallback method: if that fails, search the Gene database instead, then
# follow the link from gene -> protein (handles cases where the locus tag
# is annotated in Gene but the direct protein-db field search misses it).
#
# Usage:
#   ./fetch_genes_by_locus_tag.sh locus_tags.txt genes_protein.fasta

LOCUS_TAG_FILE="${1:-locus_tags.txt}"
OUTPUT_FASTA="${2:-genes_protein.fasta}"
WORK_DIR="./gene_fetch_work"

mkdir -p "$WORK_DIR"

if [ ! -f "$LOCUS_TAG_FILE" ]; then
    echo "Locus tag file not found: $LOCUS_TAG_FILE"
    exit 1
fi

for cmd in esearch efetch elink; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: $cmd not found on PATH. Is EDirect installed and sourced?"
        exit 1
    fi
done

> "$OUTPUT_FASTA"
FAILED_LOG="$WORK_DIR/failed_locus_tags.txt"
> "$FAILED_LOG"
REPORT_LOG="$WORK_DIR/fetch_report.tsv"
echo -e "locus_tag\tmethod\tstatus\tprotein_accession" > "$REPORT_LOG"

total=$(wc -l < "$LOCUS_TAG_FILE")
echo "Fetching $total locus tags..."

count=0
while IFS= read -r tag <&3; do
    [ -z "$tag" ] && continue
    count=$((count + 1))
    echo "[$count/$total] $tag"

    # --- Primary: direct protein-db search by Locus Tag field ---
    n_hits=$(esearch -db protein -query "${tag}[Locus Tag]" </dev/null 2>"$WORK_DIR/err.log" \
             | xtract -pattern ENTREZ_DIRECT -element Count)
    n_hits=${n_hits:-0}

    if [ "$n_hits" -eq 1 ]; then
        result=$(esearch -db protein -query "${tag}[Locus Tag]" </dev/null 2>>"$WORK_DIR/err.log" \
                 | efetch -format fasta 2>>"$WORK_DIR/err.log")
        acc=$(echo "$result" | head -1 | sed 's/^>//' | awk '{print $1}')
        echo "$result" >> "$OUTPUT_FASTA"
        echo -e "${tag}\tprotein_locus_tag\tOK\t${acc}" >> "$REPORT_LOG"
        echo "  OK (protein db direct): $acc"

    elif [ "$n_hits" -gt 1 ]; then
        echo "  WARNING: $n_hits protein hits for this locus tag -- fetching all, please review"
        result=$(esearch -db protein -query "${tag}[Locus Tag]" </dev/null 2>>"$WORK_DIR/err.log" \
                 | efetch -format fasta 2>>"$WORK_DIR/err.log")
        echo "$result" >> "$OUTPUT_FASTA"
        echo -e "${tag}\tprotein_locus_tag\tMULTIPLE_HITS(${n_hits})\tsee_fasta" >> "$REPORT_LOG"

    else
        # --- Fallback: Gene db -> elink -> protein ---
        echo "  No direct protein hit -- trying Gene db fallback..."
        gene_uid=$(esearch -db gene -query "${tag}[Locus Tag]" </dev/null 2>"$WORK_DIR/err.log" \
                   | efetch -format uid 2>>"$WORK_DIR/err.log" | head -1)

        if [ -n "$gene_uid" ]; then
            result=$(esearch -db gene -query "${tag}[Locus Tag]" </dev/null 2>>"$WORK_DIR/err.log" \
                     | elink -target protein 2>>"$WORK_DIR/err.log" \
                     | efetch -format fasta 2>>"$WORK_DIR/err.log")

            if [ -n "$result" ]; then
                acc=$(echo "$result" | head -1 | sed 's/^>//' | awk '{print $1}')
                echo "$result" >> "$OUTPUT_FASTA"
                echo -e "${tag}\tgene_elink_protein\tOK\t${acc}" >> "$REPORT_LOG"
                echo "  OK (via Gene db + elink): $acc"
            else
                echo "  FAILED: found Gene record but no linked protein"
                echo "$tag" >> "$FAILED_LOG"
                echo -e "${tag}\tgene_elink_protein\tFAILED_NO_PROTEIN_LINK\t" >> "$REPORT_LOG"
            fi
        else
            echo "  FAILED: no hits in protein or gene databases"
            echo "$tag" >> "$FAILED_LOG"
            echo -e "${tag}\tnone\tFAILED_NOT_FOUND\t" >> "$REPORT_LOG"
        fi
    fi

    sleep 0.34
done 3< "$LOCUS_TAG_FILE"

n_seqs=$(grep -c "^>" "$OUTPUT_FASTA" || echo 0)
echo ""
echo "Wrote $n_seqs sequences to $OUTPUT_FASTA"
echo "Full report: $REPORT_LOG"

if [ -s "$FAILED_LOG" ]; then
    echo "Some locus tags failed entirely -- see $FAILED_LOG"
fi