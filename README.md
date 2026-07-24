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