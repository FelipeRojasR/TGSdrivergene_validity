#!/bin/bash

# Define global variables
base_path="PATH/"

# Enter folder with official data
cd "$base_path/data/original_files" || exit 1

# Transform tsv to maf file
for filename in ./*.tsv.gz ;
do
   zcat "$filename" | python3 "$base_path/scripts/tsv2maf.py" -o "$base_path/data/mutsig/${filename%.tsv.gz}.maf" -
done

# Modify maf file to create maf file that can be used for MutSigCV
cd "$base_path/data/mutsig/" || exit 1
for filename in ./*.maf ;
do
	Rscript "$base_path/scripts/mutsig/maf_corrected.R" "$filename" "${filename%*.maf}".maf
done

# Run Mutsig in all files
for filename in ./*.maf ;
do
    bash PATH/mutsigcv/MutSigCV_1/run_MutSigCV.sh PATH/programs/matlab_runtime/v901 "$filename" PATH/mutsig2cv/exome_completeTCGA_full192.coverage.txt PATH/mutsigcv/gene.covariates.txt "$base_path/results/mutsig/${filename%}" PATH/mutsigcv/mutation_type_dictionary_file.txt PATH/mutsigcv/chr_files_hg19
done
