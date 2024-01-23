#!/usr/bin/env Rscript

library(dplyr)
library(ActiveDriver)
library(stringr)
library(readr)

args = commandArgs(trailingOnly=TRUE)
warning()
print("WARNING: THE INPUT DATASETS MUST INCLUDE ONLY MUTATIONAL DATASETS")

# Testing if there are at least one argument: return an error if there is non
if (length(args) != 2) {
  stop("Only two arguments should be provided, input.txt and output.txt", call.=FALSE)
} else if (length(args) == 2) {
  # default output file
  print("correct number of arguments")
}


set.seed(5984)

## Preparation of input data.
phospho_table <- load("PATH/activedriver/psite_table.rsav") # psite_table object created
psite_table <- psite_table %>%
  dplyr::mutate(gene = gsub("_.*", "", gene))
cols <- c("gene", "residue", "kinase", "pmid")
psite_table[cols] <- lapply(psite_table[cols], factor)


disorder_ens70 <- load("PATH/activedriver/ens70_protein_seqs_disorder.fa.rsav") # seqs_disorder object created
for (i in length(attr(seqs_disorder, "names"))) {
  attr(seqs_disorder, "names") <- lapply(
    attr(seqs_disorder, "names"), gsub, pattern = "_.*", replacement = "" 
  )
}


sequences <- load("PATH/activedriver/ens70_protein_seqs.fa.rsav") # seqs object created
for (i in length(attr(seqs, "names"))) {
  attr(seqs, "names") <- lapply(
    attr(seqs, "names"), gsub, pattern = "_.*", replacement = "" 
  )
}


### THIS LINE IS IMPORTANT
seqs_disorder <- seqs_disorder[-(which(seqs == "Sequence unavailable"))]
seqs <- seqs[-(which(seqs == "Sequence unavailable"))]

mutational_dataset = read.delim(args[1])

activedriver_results <- ActiveDriver(seqs, seqs_disorder, mutational_dataset, psite_table, simplified = TRUE)

activedriver_results_output <- activedriver_results$all_gene_based_fdr

write.table(activedriver_results_output, file=args[2], sep = "\t", row.names=FALSE, quote=FALSE)
