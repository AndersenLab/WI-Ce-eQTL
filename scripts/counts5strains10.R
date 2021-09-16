
library(tidyverse)
library(sleuth)

 
tx2gene <- read.delim("../raw_data/WS276_t2g_all.tsv", stringsAsFactors=FALSE)


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


datasett<-paste("all",nrow(samples),sep = "")


s2c <- samples %>% 
  dplyr::mutate(path=paste(filepath,"/",sample,sep = "")) %>%
  dplyr::select(-filepath,-sample) %>%
  dplyr::rename(sample=samplename) 

####### to filter transcripts ######

 
mpt <-  bind_rows(te_list,tx2gene)  %>% dplyr::filter(biotype %in% c("mRNA","pseudogenic_transcript")) 


#transcript
so609_full_transcript <- sleuth_prep(s2c, 
                                     read_bootstrap_tpm = TRUE, 
                                     extra_bootstrap_summary = TRUE,
                                     num_cores=16, 
                                     target_mapping = t2g,
                                     filter_target_id = mpt$transcript)


 

# norm transcript counts
transcript_table_mpt <- kallisto_table(so609_full_transcript, use_filtered = TRUE, normalized = TRUE) 

 

transcript_table <- transcript_table_mpt %>%
  dplyr::rename(transcript=target_id) 


count_st <- transcript_table %>% 
  dplyr::group_by(transcript,condition) %>% 
  dplyr::count()


count5 <- transcript_table %>% 
  dplyr::mutate(countlevel=ifelse(est_counts>=5,"Yes","No")) %>% 
  dplyr::group_by(transcript,condition,countlevel) %>% 
  dplyr::count() %>% 
  ungroup() %>% 
  dplyr::filter(countlevel=="Yes") %>% 
  dplyr::rename(c5=n) %>% 
  left_join(count_st) %>% 
  dplyr::mutate(countlevel_st=ifelse(c5>=n,"Yes","No"))

 

count5_1 <- count5  %>% 
  dplyr::group_by(transcript,countlevel_st) %>% 
  dplyr::count() %>%
  dplyr::filter(countlevel_st=="Yes") %>% 
  dplyr::arrange(desc(n))


count5_2 <- count5_1 %>%  
  dplyr::filter(n>=10) %>% 
  dplyr::arrange(desc(n)) 

nt <- length(unique(count5_2$transcript))

 

write.table(count5_2, paste("../processed_data/", datasett ,"_filtered_",nt,"_transcript_list.tsv",sep = ""), sep = "\t",quote = FALSE)


