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

For the query genes (above), get both nucleic acid and protein sequences