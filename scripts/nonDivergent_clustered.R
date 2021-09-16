
library(tidyverse)
library(tximport)
library(DESeq2)
library(harrietr)


args <- commandArgs(trailingOnly = TRUE)


count5_2 <- read.delim("../processed_data/all608_filtered_26043_transcript_list.tsv", stringsAsFactors=FALSE)


iso207_div_g4949t7329_list <- data.table::fread("../processed_data/iso207_div_g4949t7329_list.tsv")

fig_dir <- "../processed_data/"

 
WS276_t2g_all <- data.table::fread("../raw_data/WS276_t2g_all.tsv")



tx2gene <- WS276_t2g_all %>% 
  dplyr::filter(transcript %in% count5_2$transcript)



samples609 <- read.delim("../raw_data/kallisto_path_609samples.tsv", stringsAsFactors=FALSE) %>%
  dplyr::mutate(samplename2=paste(condition,id,sep = "_"))


###### function for kallisto tximport and deseq2 minimal filtering #####
txi_kallisto_dds <- function(meta){
  
  files <- file.path(meta$filepath, meta$sample,"abundance.h5")
  names(files) <- meta$samplename
  
  txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene)
  
  
  
  
  dds <- DESeqDataSetFromTximport(txi.kallisto,
                                  colData = meta,
                                  design = ~ condition)
  
 
  ## transform
  vsd <- vst(dds, blind=FALSE, fitType = "local")
  
 
  return(vsd)
  
}



###### function for Euc distance with transformed data #####
vsd_dist <- function(dds_vsd){
  
  
  ## remove genes in divergent regions for each strain
  
  
  dds_vsd_fd <- assay(dds_vsd) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "ens_gene") %>% 
    tidyr::gather(sample,values,-ens_gene) %>% 
    dplyr::mutate(strain=sub("(.*)(_GZ.*)(GZ)(p.*)(_2.*)","\\1",sample)) %>% 
    left_join(iso207_div_g4949t7329_list) %>% 
    dplyr::filter(is.na(gene_div)) %>% 
    dplyr::distinct(sample,values,ens_gene)  %>% 
    tidyr::spread(sample,values) %>% 
    tibble::column_to_rownames(var = "ens_gene") %>% 
    as.matrix()
  
   
  
  sampleDists <- dist(t(dds_vsd_fd))
  
 
  sampleDistMatrix <- as.matrix( sampleDists )
  
 
  ## reformat distrance matrix 
  dist_melt <- melt_dist(sampleDistMatrix) %>% 
    dplyr::mutate(iso1=sub("(.*)(_GZpool.*)","\\1",iso1),
                  iso2=sub("(.*)(_GZpool.*)","\\1",iso2))
  
  return(dist_melt)
  
}


###### function for pick samples with intra-strain distance < q.50 inter-strain distance ####
q50_inter_dist <- function(dist_df) {
  
  AB1 <- dist_df %>%
    dplyr::mutate(v1=sub("(*)(_GZ.*)","\\1",iso1),
                  v2=sub("(*)(_GZ.*)","\\1",iso2)) %>%
    dplyr::mutate(v3=ifelse(v1==v2,"intra-strain","inter-strain")) 
  
  
  dist_higher_25_list=list()
  
  for(i in unique(AB1$v1)) {
    
    AB2 <-AB1 %>%
      dplyr::filter(v1==i | v2==i )
    
    AB2_inter <- AB2 %>% dplyr::filter(!v1==v2)
    
    q50.dist = as.numeric(quantile(AB2_inter$dist, probs = 0.50, na.rm = TRUE)[1])
    
    
    AB3 <- AB2 %>%
      dplyr::filter(v1==v2) 
    
    AB4<- AB3 %>%
      dplyr::filter(dist<q50.dist) 
    
    
    
    if (nrow(AB3)==nrow(AB4)){
      AB5=AB3
    } else if (nrow(AB4)==0) {
      AB5=AB4
    } else {
      AB5=dplyr::filter(AB3,dist==min(dist))
    }
    
    
    
    dist_higher_25_list[[i]] <- AB5
    
    
  }
  
  
  dist_higher_25 <- bind_rows(dist_higher_25_list) %>%
    dplyr::select(iso1,iso2) %>%
    gather(iso,sample) 
  
  return(dist_higher_25)
}





 

try( dist_higher_25 <- read.delim(args[1], stringsAsFactors=FALSE) ) 


if((exists('dist_higher_25') && is.data.frame(get('dist_higher_25')))) {
  
  samples<-samples609 %>%
    dplyr::filter(!condition =="JU2800")  %>%
    dplyr::filter(samplename2 %in% dist_higher_25$sample) 
} else {
  
  samples<-samples609 %>%
    dplyr::filter(!condition =="JU2800")
}




datasett<-paste("all",nrow(samples),sep = "") 
 

print(datasett)


# dds vst
vsd <- txi_kallisto_dds(meta=samples)

# dist reformat
dist_melt <- vsd_dist(dds_vsd=vsd)


#save(vsd,dist_melt, file = paste(fig_dir,datasett ,"_allgene_vst_dist.RData",sep = ""))

# filter by q50
dist_intraDist_passQ50 <- q50_inter_dist(dist_df=dist_melt)

np <- length(unique(dist_intraDist_passQ50$sample))

 

write.table(dist_intraDist_passQ50, paste(fig_dir,datasett ,"_allgene__intraDist_passQ50_",np,".tsv",sep = ""), sep = "\t",quote = FALSE)


