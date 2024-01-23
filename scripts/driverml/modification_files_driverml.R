#!/usr/bin/env Rscript

library(dplyr)
library(stringr)
library(readr)

args = commandArgs(trailingOnly=TRUE)
warning()

dataset_input <- read_tsv(args[1])

# DriverML
driverml_data_creation <- function(dataset){
  dataset_driverml <- dataset %>%
    dplyr::rename(
      Hugo_Symbol = gene_symbol,
      Chromosome = chrom,
      Start_Position = start_position,
      Variant_Classification = variant_classification,
      Variant_Type = variant_type,
      Reference_Allele = ref,
      Tumor_Seq_Allele2 = alt,
      Tumor_Sample_Barcode = TCGA_ID
    ) %>%
    dplyr::select(Hugo_Symbol, Chromosome, Start_Position, 
                  Variant_Classification, Variant_Type,
                  Reference_Allele, Tumor_Seq_Allele2, 
                  Tumor_Sample_Barcode)
  return(dataset_driverml)
}

for_output <- driverml_data_creation(dataset_input)

# writes a TXT file 
write.table(for_output, file=args[2], append = TRUE, sep = "\t", row.names=FALSE, quote=FALSE)
