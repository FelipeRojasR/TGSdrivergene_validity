#!/usr/bin/env Rscript

library(dplyr)
library(stringr)
library(readr)

args = commandArgs(trailingOnly=TRUE)
warning()
print("WARNING: THE INPUT DATASETS MUST INCLUDE ONLY MUTATIONAL DATASETS")

dataset_input <- read_tsv(args[1])

oncodrive_tuning_TCGA <- function(df) { # df = TCGA_mutational_BRCA
  x <- df %>%
    dplyr::select(chrom, start_position, ref, alt, TCGA_ID) %>%
    dplyr::rename(
      Chromosome = chrom,
      Position = start_position,
      Ref = ref,
      Alt = alt,
      Sample = TCGA_ID
    )
  
  return(x)
}

for_output <- oncodrive_tuning_TCGA(dataset_input)

# writes a TSV file 
write.table(for_output, file=args[2], sep = "\t", row.names=FALSE, quote=FALSE)
