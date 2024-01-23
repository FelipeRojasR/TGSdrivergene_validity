#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(dndscv)
library(optparse)

# DNDSCV

##################################################################################################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Mutational dataset", metavar="character"),
  make_option(c("-p", "--panel"), type="character", default=NULL, help="Which gene panel? is whole-exome then WE", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output file of dndscv", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
##################################################################################################


input_dataset <- read_tsv(opt$input)

# Import panels
BCAST_gene_panel <- as_tibble(read.table("PATH/BCAST_genepanel.txt", header = TRUE, stringsAsFactors = FALSE))

MSK_Impact_extended_cols_fix <- function(dataset){
  dataset %>%
    dplyr::rename(
      hugo_symbol = `Hugo Symbol`,
      entrez_gene_id = `Entrez Gene ID`,
      GRCh37_isoform = `GRCh37 Isoform`,
      GRCh37_RefSeq = `GRCh37 RefSeq`,
      GRCh38_isoform = `GRCh38 Isoform`,
      GRCh38_RefSeq = `GRCh38 RefSeq`,
      ocurrences_within_resources = `# of occurrence within resources (Column D-J)`,
      OncoKB_annotated = `OncoKB Annotated`,
      is_oncogene = `Is Oncogene`,
      is_tumor_suppressor_gene = `Is Tumor Suppressor Gene`,
      MSK_IMPACT = `MSK-IMPACT`,
      MSK_HEME = `MSK-HEME`,
      foundation_one = `FOUNDATION ONE`,
      foundation_one_heme = `FOUNDATION ONE HEME`,
      vogelstein = Vogelstein,
      sanger_CGC = `SANGER CGC(05/30/2017)`,
      gene_aliases = `Gene Aliases`
    )
}

MSK_Impact_panel <- read_tsv("PATH/MSK_Impact_cancerGeneList.tsv") %>%
  MSK_Impact_extended_cols_fix(.) %>% 
  dplyr::filter(MSK_IMPACT == "Yes") %>%
  dplyr::select(hugo_symbol, entrez_gene_id)

data(refcds_hg19, package = "dndscv")

MSK_panel_for_dndscv <- MSK_Impact_panel %>%
  dplyr::select(hugo_symbol) %>%
  dplyr::filter(hugo_symbol %in% unique(gr_genes$names))


dndscv_data_creation_running <- function(dataset, target = FALSE, gene_panel = NULL){
  complete18000_genepanel <- dataset %>%
    dplyr::select(gene_symbol) %>%
    dplyr::filter(gene_symbol %in% unique(gr_genes$names) )
  
  dndscv_ready_dataset <- as.data.frame(dataset %>%
                                          dplyr::filter(!is.na(gene_symbol)) %>%
                                          dplyr::filter(gene_symbol %in% complete18000_genepanel$gene_symbol) %>%
                                          dplyr::select(TCGA_ID,
                                                        chrom,
                                                        start_position, 
                                                        ref,
                                                        alt) %>%
                                          dplyr::rename(
                                            sampleID = TCGA_ID,
                                            chr = chrom,
                                            pos = start_position,
                                            mut = alt
                                          )
  )
  
  if (target == FALSE) {
    complete_dndscv <- dndscv(dndscv_ready_dataset, 
                              refdb = "hg19",
                              cv = "hg19", 
                              max_muts_per_gene_per_sample = Inf, 
                              max_coding_muts_per_sample = Inf)
    complete_dndscv <- complete_dndscv$sel_cv %>%
      dplyr::select(gene_name,
                    pglobal_cv,
                    qglobal_cv)
    return(complete_dndscv)
    
  } else if (target == TRUE) {
    target_dndscv <- dndscv(dndscv_ready_dataset,
                            gene_list = gene_panel,
                            refdb = "hg19",
                            cv = "hg19", 
                            max_muts_per_gene_per_sample = Inf, 
                            max_coding_muts_per_sample = Inf)
    target_dndscv <- target_dndscv$sel_cv %>%
      dplyr::select(gene_name,
                    pglobal_cv,
                    qglobal_cv)
    return(target_dndscv)
  }
}


if(opt$panel == "BCAST"){
  output_results <- dndscv_data_creation_running(input_dataset, target = TRUE, gene_panel = BCAST_gene_panel$gene_symbol)
} else if(opt$panel == "MSK"){
  output_results <- dndscv_data_creation_running(input_dataset, target = TRUE, gene_panel = MSK_panel_for_dndscv$hugo_symbol)
} else if(opt$panel == "WE"){
  output_results <- dndscv_data_creation_running(input_dataset, target = FALSE)
}

write_tsv(output_results, file = opt$output)
