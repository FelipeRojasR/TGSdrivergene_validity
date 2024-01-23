#!/bin/bash

export LANGUAGE=en_US.UTF-8
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8

bash -c ''


# DriverML
base_path="PATH/"


# Modify input file
cd "$base_path/data/original_files" || exit 1
for filename in ./*.tsv :
  do
	  Rscript "$base_path/scripts/driverml/modification_files_driverml.R" "$filename" "$base_path/data/driverml/${filename%*.tsv}.txt"
  done


cd "$base_path/data/driverml" || exit 1
for filename in ./*.txt :
  do
	  bash "$base_path/scripts/driverml/driverml_run.sh" "$(realpath $filename)" "$base_path/results/driverml/${filename%.*}_driverml"
  done
