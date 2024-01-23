#!/usr/bin/env Rscript

library(dplyr)
library(readr)

args = commandArgs(trailingOnly=TRUE)

# Testing if there are at least one argument: return an error if there is non
if (length(args) == 0) {
    stop("At least one argument should be provided, input.maf and output.maf", call.=FALSE)
} else if (length(args) == 1) {
    # default output file
    args[2] = "out.txt"
} else if (length(args) >= 3) {
    stop("Input.maf and output.maf arguments should be provided", call.=FALSE)
}

df = read_tsv(args[1])

mutsig_TCGA_fix <- function(df){
  df_fix <- df %>% dplyr::rename(
    gene = Hugo_Symbol,
    patient = Tumor_Sample_Barcode
  ) %>%
  dplyr::mutate(patient = gsub("-", "_", patient))

  return(df_fix)
}

df_fixed <- mutsig_TCGA_fix(df)

write.table(df_fixed, file=args[2], sep = "\t", row.names=FALSE, quote=FALSE)
