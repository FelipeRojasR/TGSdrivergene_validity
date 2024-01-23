#!/bin/bash

# conda activate 2020plus
# export PATH=$PATH:PATH/2020plus

### 20/20+
base_path="PATH/"

cd "$base_path/data/original_files" || exit 1

# Prepare files
for filename in ./*.tsv.gz ;
  do
    zcat "$filename" | python3 "$base_path/scripts/official_scripts/tsv2maf.py" -o "$base_path/data/plus2020/${filename%.tsv.gz}.maf" -
  done


# Correct maf file
cd "$base_path/data/plus2020/" || exit 1

for filename in ./*.maf ;
  do
    Rscript "$base_path/scripts/official_scripts/plus2020/maf_corrected_2020plus.R" "$filename" "$base_path/data/plus2020/${filename%*.maf}.txt"
  done


# Run 20/20+
export PATH=$PATH:PATH/2020plus

for filename in *.txt;
  do
    path=$(realpath "$filename")
    bash "$base_path/scripts/official_scripts/plus2020/snakemake_script.sh" "$path" "$base_path/results/plus2020/${filename%}"
  done
