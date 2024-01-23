#!/usr/bin/env Rscript

library(dplyr)

args = commandArgs(trailingOnly=TRUE)

# Testing if there are at least one argument: return an error if there is non
if (length(args) != 2) {
    stop("Two arguments should be provided, input.maf and output.txt", call.=FALSE)
}

df <- read.delim(args[1])

df <- df %>%
  dplyr::rename(
    Gene = Hugo_Symbol,
    Tumor_Allele = Tumor_Seq_Allele2,
    Tumor_Sample = Tumor_Sample_Barcode
  ) %>%
  dplyr::select(-c(NCBI_Build))
 
df$Chromosome <- paste("chr", df$Chromosome, sep="")  

write.table(df, file=args[2], sep = "\t", row.names=FALSE, quote=FALSE)
