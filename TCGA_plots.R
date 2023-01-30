library(ggplot2)

maf_dataset <- read.maf(maf = "/PATH/mc3.v0.2.8.PUBLIC.maf") 
maf_dataset_clinical <- readr::read_table("/PATH/clinical_PANCAN_patient_with_followup.tsv")

maf_dataset@data %>%
  select(Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, SYMBOL) %>%
  mutate(bcr_patient_barcode = stringr::str_extract(Tumor_Sample_Barcode, '^([^-]+-[^-]+-[^-]+)')) %>%
  inner_join(., BRCA, by = "bcr_patient_barcode") %>%
  # filter(gender == "FEMALE") %>%
  group_by(bcr_patient_barcode) %>% distinct(bcr_patient_barcode) %>%
  summarise(n())

clinical <- maf_dataset_clinical %>%
  select(bcr_patient_uuid, bcr_patient_barcode, acronym) %>%
  mutate(acronym = case_when(acronym == "LUAD" | acronym == "LUSC" ~ "LC", TRUE ~ acronym))

#################################
#################################

supplementary_1_A <- maf_dataset@data %>%
  select(Tumor_Sample_Barcode, SYMBOL, Hugo_Symbol, 
         Entrez_Gene_Id, Chromosome, Start_Position, End_Position) %>%
  mutate(bcr_patient_barcode = stringr::str_extract(Tumor_Sample_Barcode, '^([^-]+-[^-]+-[^-]+)')) %>%
  inner_join(., clinical, by = "bcr_patient_barcode") %>%
  mutate(labeller = case_when(acronym == "BRCA" | acronym == "LC" ~ 1, TRUE ~ 0)) %>%
  group_by(bcr_patient_barcode, acronym, labeller) %>%
  summarise(n = n()) %>%
  ggplot(mapping = aes(x = reorder(acronym, -n), y = log10(n), fill = as.factor(labeller))) + 
  geom_boxplot(alpha = 0.8, show.legend = F) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  scale_fill_manual(values=c("#eeeeee","#444444")) + 
  xlab("Cancer Type") + ylab("Log10(Mutation count)")

supplementary_1_B <- maf_dataset@data %>%
  select(Tumor_Sample_Barcode, SYMBOL, Hugo_Symbol, 
         Entrez_Gene_Id, Chromosome, Start_Position, End_Position) %>%
  mutate(bcr_patient_barcode = stringr::str_extract(Tumor_Sample_Barcode, '^([^-]+-[^-]+-[^-]+)')) %>%
  inner_join(., clinical, by = "bcr_patient_barcode") %>%
  mutate(labeller = case_when(acronym == "BRCA" | acronym == "LC" ~ 1, TRUE ~ 0)) %>%
  group_by(acronym, labeller) %>% summarise(n = n_distinct(bcr_patient_barcode)) %>%
  ggplot(mapping = aes(x = reorder(acronym, -n), y = n, fill = as.factor(labeller))) + 
  geom_col(show.legend = F) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  scale_fill_manual(values=c("#bcbcbc","#444444")) + 
  xlab("Cancer Type") + ylab("Number of samples") + 
  scale_y_continuous(expand = c(0,0), limits = c(0,1050))

p = list(supplementary_1_A, supplementary_1_B) %>%
  map(~.x + labs(x = NULL))
grid_plot <- grid.arrange(grobs = p,
             ncol = 2, nrow = 1, bottom = "Cancer type")

