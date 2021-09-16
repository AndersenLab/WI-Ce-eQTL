

#### sleuth normalize data in quest #######
library(tidyverse)
library(sleuth)

fig_dir <- "../processed_data/"


tx2gene <- data.table::fread("../raw_data/WS276_t2g_all.tsv")


cele_Dfam_157 <- read.table("../raw_data/cele_Dfam_157_singleline_list.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)

te_list <- cele_Dfam_157 %>% 
  dplyr::rename(transcript=V1) %>% 
  dplyr::mutate(ens_gene=transcript,ext_gene=transcript,biotype="TE") %>% 
  dplyr::select(transcript, ens_gene, ext_gene, biotype)

t2g <- bind_rows(te_list,tx2gene) %>%
  dplyr::rename(target_id=transcript) 

samples609 <- read.delim("../raw_data/kallisto_path_609samples.tsv", stringsAsFactors=FALSE)

samples<-samples609 %>%
  dplyr::filter(!condition =="JU2800")



####### 561 samples with filtered transcripts  : protein-coding and pseudogene #######

dist_higher_50 <- read.delim(paste("../processed_data/all561_allgene__intraDist_passQ50_561.tsv",sep = ""), stringsAsFactors=FALSE)

count5_2 <- read.delim("../processed_data/all608_filtered_26043_transcript_list.tsv", stringsAsFactors=FALSE)


s2c <- samples  %>%
  dplyr::mutate(samplename2=paste(condition,id,sep = "_")) %>%
  dplyr::filter(samplename2 %in% dist_higher_50$sample) %>% 
  dplyr::mutate(path=paste(filepath,"/",sample,sep = "")) %>%
  dplyr::select(-filepath,-sample) %>%
  dplyr::rename(sample=samplename) 



datasett<-paste("all",nrow(s2c),sep = "") 



mp <- t2g %>% 
  dplyr::distinct() %>% 
  dplyr::filter(target_id %in% count5_2$transcript)  %>% 
  dplyr::filter((biotype %in% c("mRNA","pseudogenic_transcript")) ) 


so561_mp_transcript <- sleuth_prep(s2c, 
                                   num_cores=16, 
                                   target_mapping = t2g,
                                   filter_target_id = mp$target_id)

# Raw filtered transcript counts

transcript_table_mp_raw <- kallisto_table(so561_mp_transcript, use_filtered = TRUE, normalized = FALSE) %>%
  left_join(t2g, by = "target_id") %>%
  dplyr::rename(transcript=target_id) 


## estimate sizefacotr using genes not in the divergent regions
 

iso207_div_g4949t7329_list <- data.table::fread("../processed_data/iso207_div_g4949t7329_list.tsv") %>% 
  dplyr::filter(strain %in% s2c$condition)


sleuth_matrix <- sleuth_to_matrix(so561_mp_transcript, 'obs_raw', 'tpm') %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var="transcript") %>% 
  dplyr::filter(transcript %in% mp$target_id) %>% 
  dplyr::filter(!transcript %in% iso207_div_g4949t7329_list$transcript)%>% 
  tibble::column_to_rownames(var="transcript") %>% 
  as.matrix()

sizefactor <- sleuth::norm_factors(sleuth_matrix) %>% 
  as.data.frame()%>% 
  tibble::rownames_to_column(var="sample") %>% 
  dplyr::rename(sizefactor=".")


# norm  transcript tpm

transcript_table_mp_fd <- transcript_table_mp_raw %>% 
  dplyr::left_join(sizefactor) %>% 
  dplyr::rename(raw_tpm=tpm) %>% 
  dplyr::mutate(tpm=raw_tpm/sizefactor)



write.table(transcript_table_mp_fd, paste(fig_dir,datasett ,"mp_fd_transcript_norm.tsv",sep = ""), sep = "\t", row.names = F,quote = FALSE)



###### GWAS data of protein-coding and pseudogene ########

iso_ref <- read.delim("../raw_data/iso_ref_20210226.tsv")

 
sleuth_transcript_gwas <- transcript_table_mp_fd %>% 
  dplyr::mutate(lg2tpm=log2(tpm+0.5)) %>%
  dplyr::group_by(condition,transcript) %>%
  dplyr::mutate(lgmtpm=mean(lg2tpm)) %>%
  dplyr::distinct(condition,transcript,lgmtpm) 


sleuth_gwas_fd_div <- sleuth_transcript_gwas %>%
  dplyr::select(strain=condition,transcript,lgmtpm) %>%
  dplyr::left_join(iso207_div_g4949t7329_list) %>%  
  dplyr::filter(is.na(gene_div)) %>% 
  dplyr::group_by(transcript) %>% 
  dplyr::add_count(name= "num_divdata" ) %>% 
  dplyr::filter(num_divdata>=100) %>% 
  dplyr::distinct(strain,transcript,lgmtpm) %>% 
  dplyr::filter(!dplyr::n_distinct(lgmtpm) == 1 ) %>% 
  dplyr::left_join(iso_ref) %>% 
  dplyr::distinct(isotype,transcript,lgmtpm) %>% 
  dplyr::rename(strain=isotype) %>%
  tidyr::spread(transcript,lgmtpm)

 

write.table(sleuth_gwas_fd_div, paste(fig_dir,datasett,"_mp_tx_fddiv_gwas.tsv",sep=""), sep = "\t", row.names = F, quote = FALSE)
 






