
#### age ####
library(RAPToR)
library(tidyverse)

fig_dir <- "../processed_data/"

isost206_divg4949t7329 <- data.table::fread("../processed_data/all561isost206_divg4949t7329_list.tsv")


 
#561
nc561<-  data.table::fread("../processed_data/all561mp_fd_gene_norm.tsv")


#ref2
r_YA2 <- prepare_refdata("Cel_YA_2", "wormRef", n.inter = 2000)
 

##### non div gene ######

sample561_df_nondiv <- nc561  %>% 
  dplyr::select(sample,ens_gene,tpm) %>%
  dplyr::filter(!ens_gene %in% isost206_divg4949t7329$ens_gene) %>% 
  dplyr::mutate(lgt=log2(tpm+0.5)) %>%
  dplyr::select(-tpm) %>%
  spread(sample,lgt) %>% 
  column_to_rownames(var="ens_gene")


sample561_mt_nondiv <- data.matrix(sample561_df_nondiv) 

 
# age #
ae_sample561_nondiv_YA2 <- ae(samp = sample561_mt_nondiv,                        
                       refdata = r_YA2$interpGE,            
                       ref.time_series = r_YA2$time.series)





'''
nb.genes
refdata            15710
samp               13637
intersect.genes    10489
Bootstrap set size is 3496


'''


save( ae_sample561_nondiv_YA2,  file = paste(fig_dir,"age_esti_561samples.RData",sep = ""))



