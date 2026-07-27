# 2026-07-22_blast_p1Ab5075_plasmids

## Inputs:

Search query genes:
ABUW_4004
ABUW_4005
ABUW_4006
ABUW_4007
ABUW_4094

Search subject plasmids:

From TableS2_plasmids_curated_for_publication_310823.xlsx

## Extract plasmid accessions

Use `1_extract_plasmid_accessions.py`

Writes accessions.txt containing a deduplicated list of plasmid accessions

Writes accession_mapping.csv, which preserves all original rows from the excel file input.  Important because this file contained some duplicates which failed during plasmid fetching.

## Download sequences

For the plasmid subject sequences... start with nucleic acid fastas.

Activate conda env: `conda activate /projects/geisingerlab/conda_env/blast_corr/`

Request a compute node... then, command to fetch plasmid fastas: `bash scripts/2_fetch_plasmid_seqs.sh input/accessions.txt data/all_plasmids.fasta`

**Note, plasmid sequences were fetched from NCBI on July 24th, 2026**

Verify fasta sequence lengths via `3_verify_lengths_three_ways.sh`

This produces a .tsv file noting mismatches between the downloaded fasta and NCBI's sequence metadata, or between the downloaded fasta and the lengths provided in the input Excel table.

10 sequences had mismatches:
CP080453 is missing from NCBI metadata search, but I double checked with `esummary -db nuccore -id CP080453` and matched the result.  OK
The following mismatched length with the input Excel table, but I suspect these were manually entered in the table and contained typos (e.g. 110967 vs. 119067):
CP033245
CP033871
CP050429
CP050434
CP050908
CP051868
CU468232

Two had unambiguous length mismatches: CM009085 and CP058731 both had the downloaded FASTA length match the NCBI metadata, but they were shorter than the length given in the input Excel table.  Possible version differences?

For the query genes (above), get both nucleic acid and protein sequences

Fetched protein queries with script 4

## BLAST

Ran scripts 5 and 6 to make blast database and run tblastn

Generated presence/absence matrix with script 7_build_presence_matrix.py  (Parses raw `tblastn` results into a best-hit table and a gene x plasmid presence/absence matrix.)

### What it does

1. **Loads raw BLAST output** (`results/tblastn_results.tsv`, tab-separated,
   no header), using the column order matching the `-outfmt` string from
   `run_tblastn.sbatch`.

2. **Cleans up subject IDs.** The BLAST database was built with
   `-parse_seqids`, so subject IDs look like `gb|CP012345.1|` — these are
   parsed down to a plain accession (`CP012345.1`) for matching against
   the plasmid table.

3. **Collapses multiple HSPs to one best hit per (gene, plasmid) pair.**
   A single gene can align to a single plasmid in more than one local
   alignment (HSP). For each `(qseqid, accession)` pair, only the
   highest-`bitscore` HSP is kept — this is the standard convention for
   "best hit" in BLAST-based presence/absence calls.

4. **Writes `best_hits.tsv`** — one row per gene x plasmid pair,
   **completely unfiltered** (only constrained by the e-value cutoff
   already applied during the BLAST search itself). This is meant for
   manual inspection: sort by `pident` or `qcovs` to see where real hits
   drop off into noise, before committing to thresholds.

5. **Prints summary statistics** (min / 25% / median / 75% / max) of
   `%identity` and `%query coverage (qcovs)` per gene, so you can eyeball
   the distribution directly in the terminal.

6. **Applies thresholds and builds two matrices:**
   - `presence_absence.csv` — binary (1 = present, 0 = absent) gene x
     plasmid matrix, using a hit as "present" only if it passes **both**
     `PIDENT_THRESHOLD` and `QCOV_THRESHOLD` (set as plain variables near
     the top of the script — no need to touch the logic below to adjust
     them).
   - `identity_matrix.csv` — same shape, but with the best hit's raw
     `%identity` value instead of 1/0 (blank if no hit at all was found
     for that pair, even below threshold). Useful for judging borderline
     calls that the binary matrix would otherwise hide.

7. **Merges both matrices back onto the full original plasmid table**
   (`accession_mapping.csv`, from the earlier Excel-extraction step), so
   the final output has one row per original table entry — including any
   duplicate accessions — rather than being limited to only the unique
   accessions that were actually BLASTed. The join tolerates version-suffix
   mismatches (e.g. `.1` vs `.2`) by matching on the accession's base ID.

### Inputs required

| File | Description |
|---|---|
| `results/tblastn_results.tsv` | Raw `tblastn` output from `run_tblastn.sbatch` |
| `accession_mapping.csv` | Full plasmid table with accessions, from the Excel-extraction step |

### Outputs

| File | Description |
|---|---|
| `best_hits.tsv` | One row per gene x plasmid pair, best HSP only, unfiltered |
| `presence_absence.csv` | Binary gene presence/absence, merged onto full original table |
| `identity_matrix.csv` | Raw %identity per gene x plasmid, merged onto full original table |

### Adjusting thresholds

Edit these two lines near the top of the script, then re-run — no need to
redo the BLAST search itself, since `best_hits.tsv` is generated
independently of the thresholds:

```python
PIDENT_THRESHOLD = 80    # minimum % identity to call a gene "present"
QCOV_THRESHOLD = 80      # minimum % query coverage to call "present"
```
