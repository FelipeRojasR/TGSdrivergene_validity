#!/bin/bash
export LC_ALL=C.UTF-8
export LANG=C.UTF-8

# conda activate oncodriveclustl

#### OncodriveCLUSTL
base_path="PATH/"
oncodriveclustl_folder_path="PATH/oncodriveclustl"

# Modify input file
cd "$base_path/data/original_files" || exit 1

for filename in ./*.tsv :
  do
	  Rscript "$base_path/scripts/oncodriveclustl/modification_files.R" "$filename" "$base_path/data/oncodriveclustl/${filename%*.tsv}.tsv"
  done


cd "$base_path/data/oncodriveclustl" || exit 1

# Run oncodriveCLUSTL
for filename in ./*.tsv ;
  do       
    oncodriveclustl -i "$filename" -r "$oncodriveclustl_folder_path/cds.hg19.regions" --signature-calculation region_normalized --concatenate -g hg19 -n 1000 -c 30 --seed 1234 --simulation-mode region_restricted -o "$base_path/results/oncodriveclustl/${filename}"
  done
  