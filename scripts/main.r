
## Install or extract necessary packages.
list.of.packages <- c("purrr", "dplyr", "readr", "tidyr", "ggplot2",
                      "forcats", "ggpubr", "egg", "gridExtra", "ggh4x",
                      "grid", "scales")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0){install.packages(new.packages, dep=TRUE)}

for(package.i in list.of.packages){suppressPackageStartupMessages(
  library(package.i, character.only = TRUE))}
options(digits=2)

##########
# Importation of B-CAST and MSK-IMPACT gene panels 
BCAST_gene_panel <- as_tibble(read.table("PATH/data/panels/bcast_panel_v2.txt", 
                              header = TRUE, stringsAsFactors = FALSE))

MSK_Impact_panel <- as_tibble(read.table("PATH/data/panels/MSK_Impact_cancerGeneList.tsv",
                                         header = TRUE, stringsAsFactors = FALSE))

panel_genes <- list(
  `B-CAST` = BCAST_gene_panel$gene_symbol,
  `MSK-IMPACT` = MSK_Impact_panel$hugo_symbol)


##########
# Inclusion of the Cancer gene census tier 1 gene list.
# downloaded from https://cancer.sanger.ac.uk/cosmic/census?tier=1
CGC_tier1 <- read_tsv("PATH/data/CGC_tier1_18122023.tsv"); colnames(CGC_tier1) <- gsub(" ", "_",colnames(CGC_tier1), fixed = TRUE); 
CGC_tier1 <- CGC_tier1 %>%
  dplyr::rename("gene_symbol" = "Gene_Symbol", 
         "entrez_gene_id" = "Entrez_GeneId", 
         "tumor_type" = "Tumour_Types(Somatic)") %>%
  dplyr::select(gene_symbol, entrez_gene_id, tumor_type)

CGC_tier1_BCAST <- CGC_tier1 %>% filter(gene_symbol %in% BCAST_gene_panel$gene_symbol)
CGC_tier1_MSK <- CGC_tier1 %>% filter(gene_symbol %in% MSK_Impact_panel$hugo_symbol)


##########
# Expected variables for each of the seven methods that were implemented
method_col_mapping <- list(
  MutSigCV = list(),
  ActiveDriver = list(),
  dNdScv = list(gene = "gene_name", p = "pglobal_cv"),
  DriverML = list(p = "p_no_N_2s"),
  OncodriveCLUSTL = list(gene = "SYMBOL", p = "P_ANALYTICAL"),
  OncodriveFML = list(gene = "SYMBOL", p = "P_VALUE"),
  `20/20+` = list(p = "combined p-value")
)


##########
# Function to compute the estimates and confidence interval for 2x2 metrics
ratio <- function (x, n) {
  out <- tryCatch(
    exp = {
      b <- binom.test(x, n)
      list(estimate = unname(b$estimate), lb = b$conf.int[1], ub = b$conf.int[2])
      },
    error = function(e){
      list(estimate = NA, lb = NA, ub = NA)
      }
    )
  return(out)
}


##########
# Read all files into a single object
read_driver_results <- function (filename, method, panel=NULL) {
  results <-
    read_tsv(filename, show_col_types = FALSE) %>%
    dplyr::rename(!!!method_col_mapping[[method]]) %>%
    dplyr::select(gene, p)
  
  if (!is.null(panel)) {
    results <-
      results %>%
      dplyr::filter(gene %in% panel_genes[[panel]])
  }
  
  results
}

driver_results <- driver_result_files %>%
  pmap_df(
    function (filename, method, panel, cancer_type) {
      if (is.na(panel)) {
        read_driver_results(filename, method) %>% mutate(method = method, cancer_type = cancer_type, panel = NA)
      } else {
        read_driver_results(filename, method, panel) %>% mutate(method = method, cancer_type = cancer_type, panel = panel)
      }
    }
  )


driver_results <- driver_results %>%
  filter(!is.na(panel)) %>%
  left_join(
    driver_results %>%
      dplyr::filter(is.na(panel)) %>%
      dplyr::select(-panel),
    by=c("gene", "method", "cancer_type"),
    suffix=c("_panel", "_wex")
  )

driver_results$panel <- factor(driver_results$panel, levels = c("MSK-IMPACT", "B-CAST"))
driver_results$method <- factor(driver_results$method, 
                                levels = c("20/20+", "ActiveDriver", 
                                           "dNdScv", "DriverML", 
                                           "MutSigCV", "OncodriveCLUSTL", "OncodriveFML"))

# Results with significance threshold of FDR < 0.05
driver_results0.05 <- driver_results %>%
  filter(!is.na(p_panel), !is.na(p_wex)) %>%
  group_by(method, cancer_type, panel) %>%
  mutate(
    sig_wex = p.adjust(p_wex, "BH") < 0.05,
    sig_panel = p.adjust(p_panel, "BH") < 0.05,
    in_CGC = case_when(panel == "B-CAST" & gene %in% CGC_tier1_BCAST$gene_symbol ~ 1,
                       panel == "MSK-IMPACT" & gene %in% CGC_tier1_MSK$gene_symbol ~ 1,
                       TRUE ~ 0)
    )

# Results with significance threshold of FDR < 0.01
# This is the object that will be primarly used in the paper.
driver_results <- driver_results %>%
  dplyr::filter(!is.na(p_panel), !is.na(p_wex)) %>%
  group_by(method, cancer_type, panel) %>%
  mutate(
    sig_wex = p.adjust(p_wex, "BH") < 0.01,
    sig_panel = p.adjust(p_panel, "BH") < 0.01,
    in_CGC = case_when(panel == "B-CAST" & gene %in% CGC_tier1_BCAST$gene_symbol ~ 1,
                       panel == "MSK-IMPACT" & gene %in% CGC_tier1_MSK$gene_symbol ~ 1,
                       TRUE ~ 0)
    )


################
# Scatterplot of methods
figure_1_A <- driver_results %>%
  dplyr::filter(panel == "MSK-IMPACT") %>%
  ggplot(mapping = aes(x=p_wex, y=p_panel)) +
  geom_point(size=1, alpha = 0.15) +
  geom_abline(linetype="dashed", alpha = 0.5, size = 0.5) +
  facet_grid(rows = vars(method), cols = vars(cancer_type)) +
  theme_light() + theme(panel.grid = element_blank()) + 
  scale_y_continuous(name ="", breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_x_continuous(name ="", breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) + 
  theme(strip.text = element_text(colour = 'black', face = "bold"))


figure_1_B <- driver_results %>%
  dplyr::filter(panel == "B-CAST") %>%
  ggplot(mapping = aes(x=p_wex, y=p_panel)) +
  geom_point(size=1, alpha = 0.15) +
  geom_abline(linetype="dashed", alpha = 0.5, size = 0.5) +
  facet_grid(rows = vars(method), cols = vars(cancer_type)) +
  theme_light() + theme(panel.grid = element_blank()) + 
  scale_y_continuous(name ="", breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_x_continuous(name ="", breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) + 
  theme(strip.text = element_text(colour = 'black', face = "bold"))

figure_1 <- ggpubr::ggarrange(figure_1_A, figure_1_B, labels = c("A", "B"), ncol=2,
                              widths = c(7,1)) %>%
  annotate_figure(., left = textGrob("Targeted genome sequencing p-values", rot = 90, vjust = 2, gp = gpar(cex = 1)),
                  bottom = textGrob("Whole exome sequencing p-values", gp = gpar(cex = 1), vjust = -1.5))


################
# Forest plot of the main results
figure_2_A <- driver_results %>%
  dplyr::filter(panel == "MSK-IMPACT") %>%
  dplyr::filter(!is.na(p_panel), !is.na(p_wex)) %>%
  group_by(method, cancer_type, panel) %>%
  summarise(
    Sensitivity = as_tibble(ratio(sum(sig_wex & sig_panel), sum(sig_wex))),
    `False discovery rate` = as_tibble(ratio(sum(!sig_wex & sig_panel), sum(sig_panel)))
  ) %>%
  pivot_longer(c("Sensitivity", "False discovery rate")) %>%
  tidyr::unpack(value, names_sep="_") %>%

  ggplot(aes(y = fct_rev(method), x=value_estimate)) +
  geom_point(alpha = 0.7, size=2, show.legend = FALSE) +
  theme_bw() +
  geom_segment(aes(y=method, yend=method, x=value_lb, xend=value_ub), size=0.5, alpha = 0.7, show.legend = FALSE) +
  facet_grid(vars(name), vars(cancer_type, panel), scales="free_x") +
  scale_y_discrete(name ="") + scale_x_continuous(name ="", breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1))

  
figure_2_B <- driver_results %>%
  dplyr::filter(panel == "B-CAST") %>%
  dplyr::filter(!is.na(p_panel), !is.na(p_wex)) %>%
  group_by(method, cancer_type, panel) %>%
  summarise(
    Sensitivity = as_tibble(ratio(sum(sig_wex & sig_panel), sum(sig_wex))),
    `False discovery rate` = as_tibble(ratio(sum(!sig_wex & sig_panel), sum(sig_panel)))
  ) %>%
  pivot_longer(c("Sensitivity", "False discovery rate")) %>% # "specificity" was removed
  tidyr::unpack(value, names_sep="_") %>%

  ggplot(aes(y = fct_rev(method), x=value_estimate)) +
  geom_point(alpha = 0.7, size=2, show.legend = FALSE) +
  theme_bw() +
  geom_segment(aes(y=method, yend=method, x=value_lb, xend=value_ub), size=0.5, alpha = 0.7, show.legend = FALSE) +
  facet_grid(vars(name), vars(cancer_type, panel), scales="free_x") +
  scale_y_discrete(name ="") + scale_x_continuous(name ="", breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1))

figure_2 <- ggpubr::ggarrange(figure_2_A, figure_2_B, labels = c("A", "B"), ncol=2,
                              widths = c(6,1)) %>%
  annotate_figure(., left = textGrob("Methods", rot = 90, vjust = 2, gp = gpar(cex = 1)),
                  bottom = textGrob("Estimate", gp = gpar(cex = 1), vjust = -1.5))


################
# Figure for the comparison with CGC
figure_3_A_1 <- driver_results %>%
  group_by(method, cancer_type, panel) %>%
  summarise(sum = sum(sig_panel)) %>%
  ggplot(aes(x = fct_rev(method), y = sum)) +
  geom_boxplot() + theme_bw() + 
  scale_x_discrete(name = "Methods") + 
  scale_y_continuous(name ="Number of driver genes (FDR < 0.01)") + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

figure_3_A_2 <- driver_results %>%
  group_by(method, cancer_type, panel) %>%
  summarise(
    # Goldstandard
    Sensitivity = as_tibble(ratio(sum(in_CGC & sig_panel), sum(in_CGC))), # `Sensitivity CGC`
    `False discovery rate` = as_tibble(ratio(sum(!in_CGC & sig_panel), sum(sig_panel))) # `False discovery rate CGC`
  ) %>%
  pivot_longer(c("Sensitivity", "False discovery rate")) %>%
  tidyr::unpack(value, names_sep="_") %>%
  mutate(name = factor(name, levels = c("False discovery rate", "Sensitivity"))) %>%
  
  ggplot(aes(x = fct_rev(method), y = value_estimate)) +
  geom_boxplot() + theme_bw() +
  facet_grid(cols = vars(name)) + 
  scale_x_discrete(name = "Methods") + 
  scale_y_continuous(name ="Estimate", 
                     limits = c(0,1)) + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

figure_3_A <- ggpubr::ggarrange(figure_3_A_1, figure_3_A_2, nrow = 1, labels = c("A", "B"), widths = c(1,2))


figure_3_B_1 <- driver_results %>%
  dplyr::filter(panel == "MSK-IMPACT") %>%
  group_by(method, cancer_type, panel) %>%
  summarise(
    # Goldstandard
    Sensitivity = as_tibble(ratio(sum(in_CGC & sig_panel), sum(in_CGC))), # `Sensitivity CGC`
    `False discovery rate` = as_tibble(ratio(sum(!in_CGC & sig_panel), sum(sig_panel))) # `False discovery rate CGC`
  ) %>%
  pivot_longer(c("Sensitivity", "False discovery rate")) %>%
  tidyr::unpack(value, names_sep="_") %>%
  
  ggplot(aes(y = fct_rev(method), x = value_estimate)) +
  geom_point(position = position_dodge2(width = 0.7), alpha = 0.7, size=2, show.legend = FALSE) +
  theme_bw() +
  geom_linerange(aes(xmin = value_lb, xmax = value_ub), position = position_dodge2(width = 0.7), size=0.5, alpha = 0.7) + 
  facet_grid(vars(name), vars(cancer_type, panel), scales="free_x") +
  scale_y_discrete(name ="Methods") + scale_x_continuous(name ="", limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1))

figure_3_B_2 <- driver_results %>%
  dplyr::filter(panel == "B-CAST") %>%
  group_by(method, cancer_type, panel) %>%
  summarise(
    # Goldstandard
    Sensitivity = as_tibble(ratio(sum(in_CGC & sig_panel), sum(in_CGC))), # `Sensitivity CGC`
    `False discovery rate` = as_tibble(ratio(sum(!in_CGC & sig_panel), sum(sig_panel))) # `False discovery rate CGC`
  ) %>%
  pivot_longer(c("Sensitivity", "False discovery rate")) %>%
  tidyr::unpack(value, names_sep="_") %>%
  
  ggplot(aes(y = fct_rev(method), x = value_estimate)) +
  geom_point(position = position_dodge2(width = 0.7), alpha = 0.7, size=2, show.legend = FALSE) +
  theme_bw() +
  geom_linerange(aes(xmin = value_lb, xmax = value_ub), position = position_dodge2(width = 0.7), size=0.5, alpha = 0.7) + 
  facet_grid(vars(name), vars(cancer_type, panel), scales="free_x") +
  scale_y_discrete(name ="") + scale_x_continuous(name ="", limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1))

figure_3_B <- ggpubr::ggarrange(figure_3_B_1, figure_3_B_2, labels = c("C", "D"), ncol=2, widths = c(6,1)) %>%
  annotate_figure(., bottom = textGrob("Estimate", gp = gpar(cex = 1), vjust = -1.5))

figure_3 <- ggpubr::ggarrange(figure_3_A, figure_3_B, nrow = 2, widths = c(1,1), heights = c(1,1.3))


################
# Inclusion of TCGA baseline characteristics of cancer datasets
# mc3.v0.2.8.PUBLIC.maf file was obtained from teh TCGA webside as described in the manuscript
library(maftools)
maf_dataset <- read.maf(maf = "PATH/mc3.v0.2.8.PUBLIC.maf") 
maf_dataset_clinical <- readr::read_table("PATH/clinical_PANCAN_patient_with_followup.tsv")

clinical <- maf_dataset_clinical %>%
  dplyr::select(bcr_patient_uuid, bcr_patient_barcode, acronym) %>%
  mutate(acronym = case_when(acronym == "LUAD" | acronym == "LUSC" ~ "LC", TRUE ~ acronym))

# All 14 cancer types that were included
cancer_types <- c("BRCA", "LC", "GBM", "OV", "UCEC", "KIRC", "HNSC", "LGG", "THCA", "PRAD", "SKCM", "COAD", "STAD", "BLCA")

supplementary_1_A <- maf_dataset@data %>%
  dplyr::select(Tumor_Sample_Barcode, SYMBOL, Hugo_Symbol, 
         Entrez_Gene_Id, Chromosome, Start_Position, End_Position) %>%
  mutate(bcr_patient_barcode = stringr::str_extract(Tumor_Sample_Barcode, '^([^-]+-[^-]+-[^-]+)')) %>%
  inner_join(., clinical, by = "bcr_patient_barcode") %>%
  mutate(labeller = case_when(acronym %in% cancer_types ~ 1, TRUE ~ 0)) %>%
  dplyr::group_by(bcr_patient_barcode, acronym, labeller) %>%
  summarise(n = n()) %>%
  ggplot(mapping = aes(x = reorder(acronym, -n), y = log10(n), fill = as.factor(labeller))) + 
  geom_boxplot(alpha = 0.8, show.legend = F) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  scale_fill_manual(values=c("#eeeeee","#444444")) + 
  xlab("Cancer Type") + ylab("Log10(Mutation count)")

supplementary_1_B <- maf_dataset@data %>%
  dplyr::select(Tumor_Sample_Barcode, SYMBOL, Hugo_Symbol, 
         Entrez_Gene_Id, Chromosome, Start_Position, End_Position) %>%
  mutate(bcr_patient_barcode = stringr::str_extract(Tumor_Sample_Barcode, '^([^-]+-[^-]+-[^-]+)')) %>%
  inner_join(., clinical, by = "bcr_patient_barcode") %>%
  mutate(labeller = case_when(acronym %in% cancer_types ~ 1, TRUE ~ 0)) %>%
  group_by(acronym, labeller) %>% summarise(n = n_distinct(bcr_patient_barcode)) %>%
  ggplot(mapping = aes(x = reorder(acronym, -n), y = n, fill = as.factor(labeller))) + 
  geom_col(show.legend = F) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  scale_fill_manual(values=c("#bcbcbc","#444444")) + 
  xlab("Cancer Type") + ylab("Number of samples") + 
  scale_y_continuous(expand = c(0,0), limits = c(0,1050))

p = list(supplementary_1_A, supplementary_1_B) %>%
  map(~.x + labs(x = NULL))

# Number of cases and average mutational load of each dataset
supplementary_figure_1 <- grid.arrange(grobs = p,
                                       ncol = 2, nrow = 1, bottom = "Cancer type")




################
# Forest plot using a FDR thershold of 0.05
supplementary_figure_2_A <- driver_results0.05 %>%
  filter(panel == "MSK-IMPACT") %>%
  filter(!is.na(p_panel), !is.na(p_wex)) %>%
  group_by(method, cancer_type, panel) %>%
  summarise(
    Sensitivity = as_tibble(ratio(sum(sig_wex & sig_panel), sum(sig_wex))),
    `False discovery rate` = as_tibble(ratio(sum(!sig_wex & sig_panel), sum(sig_panel)))
  ) %>%
  pivot_longer(c("Sensitivity", "False discovery rate")) %>% 
  tidyr::unpack(value, names_sep="_") %>%

  ggplot(aes(y = fct_rev(method), x=value_estimate)) +
  geom_point(alpha = 0.7, size=2, show.legend = FALSE) +
  theme_bw() +
  geom_segment(aes(y=method, yend=method, x=value_lb, xend=value_ub), size=0.5, alpha = 0.7, show.legend = FALSE) +
  facet_grid(vars(name), vars(cancer_type, panel), scales="free_x") +
  scale_y_discrete(name ="") + scale_x_continuous(name ="", breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1))

supplementary_figure_2_B <- driver_results0.05 %>%
  filter(panel == "B-CAST") %>%
  filter(!is.na(p_panel), !is.na(p_wex)) %>%
  group_by(method, cancer_type, panel) %>%
  summarise(
    Sensitivity = as_tibble(ratio(sum(sig_wex & sig_panel), sum(sig_wex))),
    `False discovery rate` = as_tibble(ratio(sum(!sig_wex & sig_panel), sum(sig_panel)))
  ) %>%
  pivot_longer(c("Sensitivity", "False discovery rate")) %>% 
  tidyr::unpack(value, names_sep="_") %>%

  ggplot(aes(y = fct_rev(method), x=value_estimate)) +
  geom_point(alpha = 0.7, size=2, show.legend = FALSE) +
  theme_bw() +
  geom_segment(aes(y=method, yend=method, x=value_lb, xend=value_ub), size=0.5, alpha = 0.7, show.legend = FALSE) +
  facet_grid(vars(name), vars(cancer_type, panel), scales="free_x") +
  scale_y_discrete(name ="") + scale_x_continuous(name ="", breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1))

supplementary_figure_2 <- ggpubr::ggarrange(supplementary_figure_2_A, supplementary_figure_2_B, labels = c("A", "B"),
                                            ncol=2, widths = c(6,1)) %>%
  annotate_figure(., left = textGrob("Methods", rot = 90, vjust = 2, gp = gpar(cex = 1)),
                  bottom = textGrob("Estimate", gp = gpar(cex = 1), vjust = -1.5))


################
# Comparison of gene length across methods
gene_length_dataset <- rbind(
  dplyr::select(gene_symbol, gene_length, .data = BCAST_gene_panel), 
  dplyr::select(hugo_symbol, gene_length, .data = MSK_Impact_panel) %>%
    dplyr::rename(gene_symbol = hugo_symbol)
) %>%
  distinct() %>%
  dplyr::rename(gene = gene_symbol)


supplementary_figure_3 <- left_join(driver_results, gene_length_dataset, by = "gene") %>%
  mutate(ocurrence = case_when(sig_wex == TRUE & sig_panel == TRUE ~ "Both datasets",
                               sig_wex == TRUE & sig_panel == FALSE ~ "Whole-exome",
                               sig_wex == FALSE & sig_panel == TRUE ~ "Targeted", 
                               TRUE ~ NA_character_)) %>%
  filter(!is.na(ocurrence) & ocurrence != "Both datasets") %>%
  dplyr::select(-c(p_panel, p_wex, sig_wex, sig_panel)) %>%
  filter(!(method %in% c("ActiveDriver", "20/20+"))) %>%
  mutate(gene_length_kb = gene_length/1000) %>%
  
  ggplot(mapping = aes(x = ocurrence, y = gene_length_kb)) +
  facet_grid(cols = vars(method)) + 
  geom_boxplot(position = position_dodge(width = 1), alpha = 0.7) + 
  theme_bw() + scale_y_sqrt(limits = c(0,3000)) +
  scale_fill_manual(values = c("#424242", "#696969", "#bdbdbd")) + 
  stat_compare_means(label = "p.format", method = "wilcox.test") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Concordance between datasets",
       y = "Gene length (Kbp)")


################
# Comparison of average mutational load per cancer type across methods
gene_mutational_frequency <- maf_dataset@data %>%
  dplyr::select(Tumor_Sample_Barcode, SYMBOL, Hugo_Symbol, 
         Entrez_Gene_Id, Chromosome, Start_Position, End_Position) %>%
  mutate(bcr_patient_barcode = stringr::str_extract(Tumor_Sample_Barcode, '^([^-]+-[^-]+-[^-]+)')) %>%
  inner_join(., clinical, by = "bcr_patient_barcode") %>%
  mutate(labeller = case_when(acronym %in% cancer_types ~ 1, TRUE ~ 0)) %>% 
  filter(acronym %in% cancer_types) %>%
  group_by(Hugo_Symbol, acronym) %>%
  summarise(mutational_frequency = n()) %>%
  dplyr::rename(gene = Hugo_Symbol, 
                cancer_type = acronym)

supplementary_figure_4 <- left_join(driver_results, gene_mutational_frequency, join_by("gene", "cancer_type")) %>%
  mutate(ocurrence = case_when(sig_wex == TRUE & sig_panel == TRUE ~ "Both datasets",
                               sig_wex == TRUE & sig_panel == FALSE ~ "Whole-exome",
                               sig_wex == FALSE & sig_panel == TRUE ~ "Targeted",
                               TRUE ~ NA_character_)) %>%
  dplyr::filter(!is.na(ocurrence) & ocurrence != "Both datasets") %>%
  dplyr::select(-c(p_panel, p_wex, sig_wex, sig_panel)) %>%
  filter(!(method %in% c("ActiveDriver", "20/20+"))) %>%

  ggplot(mapping = aes(x = ocurrence, y = mutational_frequency)) +
  facet_grid(cols = vars(method)) +
  geom_boxplot(position = position_dodge(width = 1), alpha = 0.7) +
  theme_bw() + stat_compare_means(label = "p.format", method = "wilcox.test") + 
  scale_y_sqrt() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x = "Concordance between datasets",
       y = "Gene-specific mutational load")


################
# Supplementary table 2 - Mean Absolute Error
sup_table2 <- driver_results %>%
  filter(!is.na(p_panel), !is.na(p_wex)) %>%  ##!!!!!
  group_by(method, cancer_type, panel) %>%
  summarise(
    mean_abs_error = mean(abs(p_panel - p_wex))
  ) %>%
  tidyr::pivot_longer(mean_abs_error) %>%
  tidyr::unpack(value, names_sep="_")


sup_table2 <- sup_table2 %>% unite("merged", c(cancer_type, panel), remove = TRUE) %>%
  dplyr::select(-name) %>% pivot_wider(names_from = merged, values_from = value) %>%
  mutate_all(., scientific)

################
# Supplementary table 3 - Point estimate and 95% CI for sensitivity and false discovery rate (FDR < 0.01)
sup_table3 <- driver_results %>%
  group_by(method, cancer_type, panel) %>%
  summarise(
    sensitivity = as_tibble(ratio(sum(sig_wex & sig_panel), sum(sig_wex))),
    fdr = as_tibble(ratio(sum(!sig_wex & sig_panel), sum(sig_panel)))
  ) %>%
  pivot_longer(c("sensitivity", "fdr")) %>%
  tidyr::unpack(value, names_sep="_") %>%
  mutate(to_print = paste0(round(value_estimate, 2), " (", round(value_lb,2), " - ", round(value_ub,2), ")")) %>%
  dplyr::select(-c(value_estimate, value_lb, value_ub))

sup_table3 <- sup_table3 %>% unite("merged", c(cancer_type, panel), remove = TRUE) %>%
  pivot_wider(names_from = "merged", values_from = "to_print") %>%
  arrange(name)


################
# Supplementary table 4 - Point estimate and 95% CI for sensitivity and false discovery rate (FDR < 0.05)
sup_table4 <- driver_results0.05 %>%
  summarise(
    sensitivity = as_tibble(ratio(sum(sig_wex & sig_panel), sum(sig_wex))),
    fdr = as_tibble(ratio(sum(!sig_wex & sig_panel), sum(sig_panel)))
  ) %>%
  pivot_longer(c("sensitivity", "fdr")) %>% # "specificity" was removed
  unpack(value, names_sep="_") %>%
  mutate(to_print = paste0(round(value_estimate, 2), " (", round(value_lb,2), " - ", round(value_ub,2), ")")) %>%
  dplyr::select(-c(value_estimate, value_lb, value_ub))

sup_table4 <- sup_table4 %>% unite("merged", c(cancer_type, panel), remove = TRUE) %>%
  pivot_wider(names_from = "merged", values_from = "to_print") %>%
  arrange(name)


################
# Supplementary table 5 - Point estimate and 95% CI for sensitivity and false discovery rate using CGC as reference (FDR < 0.01)
sup_table5 <- driver_results %>%
  summarise(
    # Goldstandard
    Sensitivity = as_tibble(ratio(sum(in_CGC & sig_panel), sum(in_CGC))), # `Sensitivity CGC`
    `False discovery rate` = as_tibble(ratio(sum(!in_CGC & sig_panel), sum(sig_panel))) # `False discovery rate CGC`
  ) %>%
  pivot_longer(c("Sensitivity", "False discovery rate")) %>%
  tidyr::unpack(value, names_sep="_") %>%
  mutate(to_print = paste0(round(value_estimate, 2), " (", round(value_lb,2), " - ", round(value_ub,2), ")")) %>%
  dplyr::select(-c(value_estimate, value_lb, value_ub))


sup_table5 <- sup_table5 %>% unite("merged", c(cancer_type, panel), remove = TRUE) %>%
  pivot_wider(names_from = "merged", values_from = "to_print") %>%
  arrange(name)


# DO NOT RUN
# This is how the gene length was obtained from BioMart

# final.genesBCAST <- annotations %>% dplyr::filter(ensembl_gene_id %in% BCAST_gene_panel$ensembl_gene_id)
# final.genesBCAST <- final.genesBCAST[order(match(final.genesBCAST$ensembl_gene_id, BCAST_gene_panel$ensembl_gene_id)),]; rownames(final.genesBCAST) <-NULL
# final.genesBCAST <- final.genesBCAST %>%
#   select(ensembl_gene_id, gene_length) %>%
#   distinct() %>%
#   group_by(ensembl_gene_id) %>%
#   filter(gene_length == max(gene_length))
# BCAST_gene_panel_new <- left_join(BCAST_gene_panel, final.genesBCAST, by = "ensembl_gene_id") %>%
#   filter(gene_symbol %in% BCAST_gene_panel$gene_symbol)
# 
# 
# # MSK_Impact_panel$
# final.genesMSK <- annotations %>% dplyr::filter(entrezgene_id %in% MSK_Impact_panel$entrez_gene_id)
# final.genesMSK <- final.genesMSK[order(match(final.genesMSK$entrezgene_id, MSK_Impact_panel$entrez_gene_id)),]; rownames(final.genesMSK) <-NULL
# final.genesMSK <- final.genesMSK %>%
#   select(entrezgene_id, gene_length) %>%
#   distinct() %>%
#   group_by(entrezgene_id) %>%
#   filter(gene_length == max(gene_length))
# MSK_Impact_panel_new <- left_join(MSK_Impact_panel, final.genesMSK, by = c("entrez_gene_id" = "entrezgene_id")) %>%
#   distinct() %>%
#   filter(hugo_symbol %in% MSK_Impact_panel$hugo_symbol)