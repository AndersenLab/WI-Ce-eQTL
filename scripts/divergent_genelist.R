
#########################
####### bash ############
#########################

zcat ../raw_data/lee2020.divergent_regions_strain.bed.gz > ../processed_data/div_regions_st402.bed

zcat ../raw_data/c_elegans.PRJNA13758.WS276.canonical_geneset.gtf.gz | grep "protein_coding\|pseudogene" | awk '{print $1,$3,$4,$5,$10,$12}' | grep "transcript"  | sed 's/"//g' | sed 's/;//g' | awk '{print $1,$3,$4,$6}' OFS="\\t" > ../processed_data/WS276_tx_bed.bed

bedtools intersect -wao -a ../processed_data/WS276_tx_bed.bed -b ../processed_data/div_regions_st402.bed > ../processed_data/WS276_tx_divergent_402st.bed


#########################
########## R ############
#########################

library(tidyverse)
 

iso_ref <- read.delim("../raw_data/iso_ref_20210226.tsv")


c_elegans.PRJNA13758.WS276.geneIDs <- read.csv("../raw_data/c_elegans.PRJNA13758.WS276.geneIDs.txt.gz", header=FALSE, stringsAsFactors=FALSE, na.strings=c("","NA")) %>% 
  dplyr::select(-V1)


WS276_t2g_all <- data.table::fread("../raw_data/WS276_t2g_all.tsv")

samples609 <- read.delim("../raw_data/kallisto_path_609samples.tsv", stringsAsFactors=FALSE)






WS276_tx_divergent_402st <- data.table::fread("../processed_data/WS276_tx_divergent_402st.bed")

colnames(WS276_tx_divergent_402st) <- c("Chr","start","end","transcript", "div_Chr","div_start","div_end","strain","over")

gene_div_list <- WS276_tx_divergent_402st %>% 
  dplyr::filter(over>0) %>% 
  dplyr::distinct(transcript,strain)  %>% 
  dplyr::left_join(WS276_t2g_all) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(gene_div="divergent") %>% 
  dplyr::left_join(iso_ref) %>% 
  na.omit() %>% 
  dplyr::filter(strain %in% samples609$condition)

 
write.table(gene_div_list, paste("../processed_data/iso",length(unique(gene_div_list$isotype)),"_div_g",length(unique(gene_div_list$ens_gene)),"t",length(unique(gene_div_list$transcript)),"_list.tsv",sep=""), sep = "\t", row.names = F, quote = FALSE)


