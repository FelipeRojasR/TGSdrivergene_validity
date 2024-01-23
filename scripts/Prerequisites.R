library(maftools)
library(readr)
library(dplyr)
library(purrr)

# Clinical dataset - This includes the clinical data of the TCGA database
TCGA_clinical <- read_table("PATH/clinical_PANCAN_patient_with_followup.tsv") %>%
  dplyr::select(bcr_patient_uuid, bcr_patient_barcode, acronym) %>%
  mutate(acronym = case_when(acronym == "LUAD" | acronym == "LUSC" ~ "LC", TRUE ~ acronym))


# TCGA mutational dataset
TCGA_mutational <- read.delim(file = "PATH/mc3.v0.2.8.PUBLIC.maf",
                              sep = '', header = TRUE, stringsAsFactors = FALSE)

transform_labels_Matched_Norm_Sample_Barcode <- function(){
  TCGA_mutational_Norm_SampleBarcode_fix <- gsub("-", ":", TCGA$Matched_Norm_Sample_Barcode)
  n <- 3
  pat <- paste0('^([^:]+(?::[^:]+){',n-1,'}).*')
  
  TCGA_mutational_Norm_SampleBarcode_fix <- sub(pat, '\\1', TCGA_mutational_Norm_SampleBarcode_fix)
  TCGA_mutational_Norm_SampleBarcode_fix <- gsub(":", "-", TCGA_mutational_Norm_SampleBarcode_fix)
  return(TCGA_mutational_Norm_SampleBarcode_fix)
}

TCGA_mutational$Matched_Norm_Sample_Barcode_fix <- transform_labels_Matched_Norm_Sample_Barcode()

TCGA_mutational[TCGA_mutational == "."] = NA
TCGA_mutational[TCGA_mutational == "SNP"] = "SNV"


TCGA_mutational <- TCGA_mutational %>% dplyr::rename(TCGA_ID = Matched_Norm_Sample_Barcode_fix,
                                                     chrom = Chromosome, start_position = Start_Position, end_position = End_Position,
                                                     ref = Reference_Allele, alt = Tumor_Seq_Allele2,
                                                     gene_symbol = Hugo_Symbol, Entrez_Gene_Id = Entrez_Gene_Id,
                                                     impact = IMPACT,
                                                     variant_classification = Variant_Classification, variant_type = Variant_Type,
                                                     hgvsc = HGVSc, hgvsp = HGVSp, hgvsp_short = HGVSp_Short
) %>% dplyr::select(TCGA_ID, chrom, start_position, end_position, ref, alt, gene_symbol, Entrez_Gene_Id, impact, variant_type, variant_classification, hgvsc, hgvsp, hgvsp_short, NCBI_Build) %>%
  mutate(chrom = as.character(chrom), start_position = as.numeric(start_position), end_position = as.numeric(end_position),
         ref = as.character(ref), alt = as.character(alt),
         gene_symbol = as.character(gene_symbol), Entrez_Gene_Id = as.numeric(Entrez_Gene_Id),
         impact = as.character(impact),
         variant_type = as.character(variant_type), variant_classification = as.character(variant_classification),
         hgvsc = as.character(hgvsc), hgvsp = as.character(hgvsp), hgvsp_short = as.character(hgvsp_short))


# BCAST gene panel
BCAST_gene_panel <- as_tibble(read.table("PATH/BCAST_genepanel.txt", header = TRUE, stringsAsFactors = FALSE))


# MSK-Impact gene panel
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

# Cancer types - Here we include all the cancer types to be analyzed
cancer_types <- c("BRCA", "LC", "GBM", "OV", "UCEC", "KIRC", "HNSC", "LGG", "THCA", "PRAD", "SKCM", "COAD", "STAD", "BLCA", "ACC")

# Function to filter clinical IDs by cancer type
filter_clinical_ids <- function(data, cancer_types) {
  map(cancer_types, ~ data %>% filter(acronym == .x))
}

filtered_ids <- filter_clinical_ids(TCGA_clinical, cancer_types)
names(filtered_ids) <- paste(cancer_types, "_IDS", sep = "")


# Function to filter mutational data by gene panel
filter_mutational_data_by_gene_panel <- function(data, gene_panel) {
  if (!is.null(gene_panel)) {
    data %>% filter(gene_symbol %in% gene_panel)
  } else {
    data
  }
}

# Filter mutational data by gene panel
WE_complete_TCGA_mutational <- filter_mutational_data_by_gene_panel(TCGA_mutational, NULL)
BCAST_complete_TCGA_mutational <- filter_mutational_data_by_gene_panel(TCGA_mutational, BCAST_gene_panel$gene_symbol)
MSK_complete_TCGA_mutational <- filter_mutational_data_by_gene_panel(TCGA_mutational, MSK_Impact_panel$hugo_symbol)


mutational_data_list <- list()

for (acronym in cancer_types) {
  acronym_name <- paste(acronym, "_IDS", sep = "")
  we_mutational <- filter_mutational_data_by_gene_panel(TCGA_mutational, NULL) %>% dplyr::filter(TCGA_ID %in% filtered_ids[[acronym_name]]$bcr_patient_barcode)
  msk_mutational <- filter_mutational_data_by_gene_panel(TCGA_mutational, MSK_Impact_panel$hugo_symbol) %>% dplyr::filter(TCGA_ID %in% filtered_ids[[acronym_name]]$bcr_patient_barcode)
  mutational_data_list[[paste("WE", acronym, "TCGA_mutational", sep = "_")]] <- we_mutational
  mutational_data_list[[paste("MSK", acronym, "TCGA_mutational", sep = "_")]] <- msk_mutational
}

# TCGA and BCAST BRCA have to be done outside the loop to avoid addint too many unnecessary objects
mutational_data_list[[paste("BCAST", "BRCA", "TCGA_mutational", sep = "_")]] <- filter_mutational_data_by_gene_panel(TCGA_mutational, BCAST_gene_panel$gene_symbol) %>% 
  dplyr::filter(TCGA_ID %in% filtered_ids[["BRCA_IDS"]]$bcr_patient_barcode)
mutational_data_list[[paste("WE", "COMPLETE", "TCGA_mutational", sep = "_")]] <- filter_mutational_data_by_gene_panel(TCGA_mutational, NULL) 
mutational_data_list[[paste("BCAST", "COMPLETE", "TCGA_mutational", sep = "_")]] <- filter_mutational_data_by_gene_panel(TCGA_mutational, BCAST_gene_panel$gene_symbol) 
mutational_data_list[[paste("MSK", "COMPLETE", "TCGA_mutational", sep = "_")]] <- filter_mutational_data_by_gene_panel(TCGA_mutational, MSK_Impact_panel$hugo_symbol) 

for(i in seq_along(mutational_data_list)){
  object_name <- names(mutational_data_list[i])
  file_name <- paste0("PATH/data/original_files/", object_name, ".tsv")
  write.table(mutational_data_list[[i]], file = file_name, sep = "\t", row.names = FALSE)
}

# gzip -k *.tsv &&
# rm(TCGA_mutational)
# gc()