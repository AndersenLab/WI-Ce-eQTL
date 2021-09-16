
library(tidyverse)
library(tximport)
library(sleuth)

### normalized gene by non-div size factor ####

count5_2 <- read.delim("../processed_data/all608_filtered_26043_transcript_list.tsv", stringsAsFactors=FALSE)

fig_dir <- "../processed_data/"

WS276_t2g_all <- data.table::fread("../raw_data/WS276_t2g_all.tsv")

tx2gene <- WS276_t2g_all %>% 
  dplyr::filter(transcript %in% count5_2$transcript)


samples609 <- read.delim("../raw_data/kallisto_path_609samples.tsv", stringsAsFactors=FALSE) %>%
  dplyr::mutate(samplename2=paste(condition,id,sep = "_"))


dist_higher_50 <- read.delim(paste("../processed_data/all561_allgene__intraDist_passQ50_561.tsv",sep = ""), stringsAsFactors=FALSE)


samples<-samples609 %>%
  dplyr::filter(!condition =="JU2800") %>%
  dplyr::filter(samplename2 %in% dist_higher_50$sample) 


datasett<-paste("all",nrow(samples),sep = "") 



iso207_div_g4949t7329_list <- data.table::fread("../processed_data/iso207_div_g4949t7329_list.tsv") %>% 
  dplyr::filter(strain %in% samples$condition)

niso=length(unique(iso207_div_g4949t7329_list$isotype))
nst=length(unique(iso207_div_g4949t7329_list$strain))
ntx=length(unique(iso207_div_g4949t7329_list$transcript))
ng=length(unique(iso207_div_g4949t7329_list$ens_gene))

### use in age estimation
write.table(iso207_div_g4949t7329_list, paste(fig_dir,datasett ,"isost",niso,"_divg",ng,"t",ntx,"_list.tsv",sep = ""), sep = "\t", row.names = F,quote = FALSE)


# import using tximport
meta<-samples

files <- file.path(meta$filepath, meta$sample,"abundance.h5")
names(files) <- meta$samplename

txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene)




sleuth_matrix <- txi.kallisto$abundance %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var="ens_gene") %>% 
  dplyr::filter(ens_gene %in% tx2gene$ens_gene) %>% 
  dplyr::filter(!ens_gene %in% iso207_div_g4949t7329_list$ens_gene)%>% 
  tibble::column_to_rownames(var="ens_gene") %>% 
  as.matrix()

sizefactor <- sleuth::norm_factors(sleuth_matrix) %>% 
  as.data.frame()%>% 
  tibble::rownames_to_column(var="sample") %>% 
  dplyr::rename(sizefactor=".")

norm_gene <- txi.kallisto$abundance %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var="ens_gene") %>%  
  tidyr::gather(sample,tpm,-ens_gene) %>% 
  dplyr::left_join(sizefactor) %>% 
  na.omit() %>% 
  dplyr::rename(raw_tpm=tpm) %>% 
  dplyr::mutate(tpm=raw_tpm/sizefactor)


write.table(norm_gene, paste(fig_dir,datasett ,"mp_fd_gene_norm.tsv",sep = ""), sep = "\t", row.names = F,quote = FALSE)



