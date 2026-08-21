#!/bin/bash

# Converts a single fastGEAR output directory into an importation status
# file that maskAlignFastGearImportationFile.py can read.
# Usage: bash fastgear_import.sh <fastgear_output_dir>

DIR=$1

if [ ! -d "$DIR" ]; then
    echo "ERROR: directory not found -> $DIR"
    exit 1
fi

echo "Processing: $DIR"

# tidy up the lineage file into a proper tsv, not strictly needed
# downstream but useful to keep around for checking lineage calls
if [ -f "$DIR/lineage_information.txt" ]; then
    sed 's/  */ /g' "$DIR/lineage_information.txt" \
        | tr ' ' '\t' \
        > "$DIR/lineage_information.tsv"
else
    echo "WARNING: lineage_information.txt missing in $DIR"
fi

printf "Node\tBeg\tEnd\n" > "$DIR/FastGear.importation_status.txt"

# recent recombination events (recipient lineage, start, end)
if [ -f "$DIR/recombinations_recent.txt" ]; then
    grep -v "RECOMB\|Start" "$DIR/recombinations_recent.txt" \
        | sed 's/  */ /g' \
        | cut -f1,2,6 -d' ' \
        | awk '{print $3,$1,$2}' \
        | tr ' ' '\t' \
        | sort | uniq \
        >> "$DIR/FastGear.importation_status.txt"
else
    echo "WARNING: recombinations_recent.txt missing in $DIR"
fi

# ancestral recombination events
if [ -f "$DIR/recombinations_ancestral.txt" ]; then
    grep -v "LINEAGE" "$DIR/recombinations_ancestral.txt" \
        | sed 's/  */ /g' \
        | tr ' ' '\t' \
        | awk '{print $1,$2,$3,$4}' \
        >> "$DIR/FastGear.importation_status.txt"
else
    echo "WARNING: recombinations_ancestral.txt missing in $DIR"
fi

echo "DONE -> $DIR/FastGear.importation_status.txt"
