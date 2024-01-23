#!/usr/bin/env Rscript

library(dplyr)
library(stringr)
library(readr)

args = commandArgs(trailingOnly=TRUE)
warning()
print("WARNING: THE INPUT DATASETS MUST INCLUDE ONLY MUTATIONAL DATASETS")

input <- readr::read_tsv(args[1])
cancer_types <- as.character(args[2])

activedriver_data_modification <- function(dataset, cancer_types){
  dataset_activeDriver <- dataset %>% 
    dplyr::filter(variant_type == "SNV" & !is.na(hgvsp_short)) %>% # Select only the single nucleotide variants of the dataset
    dplyr::mutate(hgvsp_short = gsub("p.", "", hgvsp_short)) %>%
    mutate(wt_residue = gsub("[0-9].*", "", hgvsp_short)) %>% # New column extracting only CHARACTERS before the number of the hgvsp_short string.
    mutate(position = as.numeric(gsub("[^[:digit:]]", "", hgvsp_short))) %>% # New column extracting only the NUMBER of the hgvsp_short string.
    mutate(mut_residue = sub("[:A-Z].*[0-9]", "", hgvsp_short)) %>% # New column extracting the CHARACTERS after the number of the hgvsp_short string.
    dplyr::filter(str_detect(mut_residue, "^[A-Za-z]+$") & str_detect(wt_residue, "^[A-Za-z]+$")) %>% # Only letters, excludes symbols like *.
    mutate(cancer_type = as.character(cancer_types)) %>% # "pancancer"
    dplyr::select(gene_symbol,
                  cancer_type, 
                  TCGA_ID, 
                  position, 
                  wt_residue,
                  mut_residue) %>%
    dplyr::filter(mut_residue != "fs") %>%
    dplyr::rename(
      gene = gene_symbol,
      sample_id = TCGA_ID
    )
  return(dataset_activeDriver)
}

export <- activedriver_data_modification(input, cancer_types)

write.table(export, file = args[3], 
            append = TRUE, sep = "\t", row.names=FALSE, quote=FALSE)

