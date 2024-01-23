#!/bin/bash

# Activate the Conda environment
# conda activate oncodrivefml

#### OncodriveFML
# Define global variables
base_path="PATH/"

Make files suitable for OncodriveFML
cd "$base_path/data/original_files" || exit 1

for filename in ./*.tsv :
  do
	  Rscript "$base_path/scripts/oncodrivefml/modification_files.R" "$filename" "$base_path/data/oncodrivefml/${filename%*.tsv}.tsv"
  done


cd "$base_path/data/oncodrivefml" || exit 1

# Whole exome sequencing - It starts with WE_
for filename in ./WE_*.tsv ; do
  oncodrivefml -i "$filename" -e "$base_path/scripts/oncodrivefml/cds.tsv" -t coding -s wes -c "$base_path/scripts/oncodrivefml/oncodrivefml_configure_updated.conf" --seed 1234 --debug -o "$base_path/results/oncodrivefml/${filename%}"
done

# Targeted genome sequencing - MSK
for filename in ./MSK_*.tsv ; do
  oncodrivefml -i "$filename" -e "$base_path/scripts/oncodrivefml/cds.tsv" -t coding -s targeted -c "$base_path/scripts/oncodrivefml/oncodrivefml_configure_updated.conf" --seed 1234 --debug -o "$base_path/results/oncodrivefml/${filename%}"
done

# Targeted genome sequencing
for filename in ./BCAST_*.tsv ; do
  oncodrivefml -i "$filename" -e "$base_path/scripts/oncodrivefml/cds.tsv" -t coding -s targeted -c "$base_path/scripts/oncodrivefml/oncodrivefml_configure_updated.conf" --seed 1234 --debug -o "$base_path/results/oncodrivefml/${filename%}"
done