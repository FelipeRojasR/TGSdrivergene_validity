# Methods for driver detection
These are the list of commands used to generate all the results in the manuscript.

#### [20/20+](https://2020plus.readthedocs.io/en/latest/)
```
snakemake -s Snakefile pretrained_predict -p --cores 20 \
	--config mutations = full_path/INPUT.txt output_dir = full_path/OUTPUT trained_classifier = 2020plus_100k.Rdata
```

#### [ActiveDriver](http://www.baderlab.org/Software/ActiveDriver)
```
library(dplyr)
library(ActiveDriver)
library(stringr)

set.seed(5984)

## Preparation of input data.
phospho_table <- load("PATH/psite_table.rsav") # psite_table probided in the package
psite_table <- psite_table %>% dplyr::mutate(gene = gsub("_.*", "", gene))
cols <- c("gene", "residue", "kinase", "pmid")
psite_table[cols] <- lapply(psite_table[cols], factor)

disorder_ens70 <- load("PATH/ens70_protein_seqs_disorder.fa.rsav") # seqs_disorder probided in the original package
for (i in length(attr(seqs_disorder, "names"))) {
  attr(seqs_disorder, "names") <- lapply(
    attr(seqs_disorder, "names"), gsub, pattern = "_.*", replacement = ""
  )
}

sequences <- load("PATH/ens70_protein_seqs.fa.rsav") # seqs object probided in the original package
for (i in length(attr(seqs, "names"))) {
  attr(seqs, "names") <- lapply(
    attr(seqs, "names"), gsub, pattern = "_.*", replacement = ""
  )
}

activedriver_results <- ActiveDriver(seqs, seqs_disorder, INPUT_mutational_dataset, psite_table)
activedriver_results_output <- activedriver_results$all_gene_based_fdr ## Only p values were selected
```


#### [dNdScv](https://github.com/im3sanger/dndscv)
```
dndscv(dndscv_ready_dataset,
       refdb = "hg19",
       cv = "hg19",
       gene_list = gene_panel, # Only use in case of targeted datasets, for whole exome use NA
       max_muts_per_gene_per_sample = Inf,
       max_coding_muts_per_sample = Inf)
```

#### [DriverML](https://github.com/HelloYiHan/DriverML)
```
bash PATH/run.driverml.sh -w PATH/DriverML-master -i INPUT -f hg19.fa -r hg19.fa -m 10 -o OUTPUT
```

#### [MutSigCV](https://software.broadinstitute.org/cancer/cga/mutsig)
```
bash run_MutSigCV.sh matlab_runtime/v901 INPUT exome_completeTCGA_full192.coverage.txt gene.covariates.txt OUTPUT mutation_type_dictionary_file.txt chr_files_hg19
```

#### [MutSigCV_TS]
Refer to the scripts/mutsig2cv_adaptation.Rmd 


#### [OncodriveCLUSTL](https://bitbucket.org/bbglab/oncodriveclustl/src/master/)
```
oncodriveclustl -i INPUT -r cds.hg19.regions --signature-calculation region_normalized --concatenate -g hg19 -n 1000 -c 30 --seed 1234 --simulation-mode region_restricted -o OUTPUT --qqplot
```

#### [OncodriveFML](https://bitbucket.org/bbglab/oncodrivefml/src/master/)
```
oncodrivefml -i INPUT -e cds.tsv -t coding -s wes -c oncodrivefml_configure_updated.conf --seed 1234 --debug -o OUTPUT
# The updated oncodrivefml_configure_updated.conf can be find in the dependencies folder 
```

### License
GNU general public licence v.3 from 29 June 2007.