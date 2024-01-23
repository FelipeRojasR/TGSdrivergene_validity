#!/bin/bash

#### dndscv
# Define global variables
base_path="PATH/"

cd "$base_path/data/original_files" || exit 1

# Targeted genome sequencing - BCAST
for filename in ./BCAST_*.tsv ;
do
  Rscript "$base_path/scripts/dndscv/dndscv_complete.R" -i "$filename" -p "BCAST" -o "$base_path/results/dndscv/${filename%.tsv}"
done
