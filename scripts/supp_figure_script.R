
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/"))
setwd("../")




### load R packages####

library(tidyverse)
library(ggthemes) 

 

### ggplot theme ####
theme_cust <- theme_bw() + 
  theme(plot.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size=12,  color = "black"),
        legend.title =  ggplot2::element_text(size=12,  color = "black"),
        axis.title =  ggplot2::element_text(size=12,  color = "black"),
        axis.text =  ggplot2::element_text(size=12,  color = "black"),
        strip.text = ggplot2::element_text(size=12, vjust = 1,  color = "black"),
        strip.background = ggplot2::element_blank(), 
        panel.grid = ggplot2::element_blank(),
        text = ggplot2::element_text(family="Helvetica"))



ancestry.colours <- c(  "darkorange1",   "lightskyblue2", 
                      "springgreen4", "lightpink2",  "deepskyblue4", 
                      "yellow3",  "yellow4",  
                      'black','red2', 'cornflowerblue', 'magenta', 'darkolivegreen4', 
                       'tan4', 'darkblue', 'yellowgreen', "tan1",
                        'wheat4', '#DDAD4B', 'chartreuse','seagreen1',
                      'moccasin', 'mediumvioletred', 'cadetblue1',"darkolivegreen1" ,"#7CE3D8",
                      "gainsboro","#E69F00","#009E73", "#F0E442", "sienna4", "#0072B2", 
                      "mediumpurple4","#D55E00", "burlywood3","gray51","#CC79A7","gray19", 'indianred1',"firebrick") 

rodon_color <- c( "#aadce0", "springgreen4","#ffe6b7",
                  "#1e466e","gold2","red","magenta",
                  "#376795","lightpink2","#528fad")

##### 

##############################
##### Supplementary Figs #####
##############################
 
 
#### Figure S2  ####   
##  
#### wormCat enrichment of all genes with eQTL

#input

data_figS2  <- data.table::fread("processed_data/Supp_Data_3.tsv") %>% 
  dplyr::filter(Input=="All genes with eQTL")

data_figS2$level= gsub("category", "Category ", data_figS2$level)
# PLOT


fig_S2 <- ggplot(data_figS2,
                 aes(x=-log10(Bonferroni),y=fct_reorder(Category, -log10(Bonferroni)),size=RGS,color=GeneRatio)) +
  geom_point()  + 
  theme_cust +
  theme(axis.title.y =  element_blank(),
        legend.position = c(-0.8,0.73))+      
  scale_color_gradient(high = "#D7263D", low = "#0072B2") +
  labs(x =expression(-log[10](adjusted~italic(p))),
       color = "Gene ratio",
       size = "Input gene counts")+ 
  guides(size = guide_legend(order = 2),color = guide_colourbar(order = 1)) +
  facet_grid( level~.,scales = "free",space = "free")  

ggsave(fig_S2, filename = paste("figures/p1_Supp_Fig2_eQTLenrich.png",sep = ""), units = "mm",height = 180, width = 170)


#### Figure S3  ####   
##  
#### A histogram showing the distribution of linkage disequilibrium (LD) 

#input
data_figS3 <- data.table::fread("processed_data/FileS13_eQTL_LD.tsv")

# PLOT
fig_S3 <- ggplot(data_figS3,aes(LD)) + 
  geom_histogram(#aes(fill=chroms),
                 alpha=0.8,binwidth = 0.02,color="black",fill="gray69") + 
  theme_cust +
  theme( legend.position = "bottom",
         legend.background = element_rect(),
         legend.title = element_blank())+
  ylab("Number of pairwise LD calculation") + 
  xlab( expression(paste("Pairwise LD (", italic(r)^2,") among eQTL" ))) +
 # guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  facet_grid(chroms~ld_betw,scales = "free_y")  

ggsave(fig_S3, filename = paste("figures/p1_Supp_Fig3_eQTL_LD.png",sep = ""), units = "mm",height = 170, width = 140)
 
#### Figure S4 ####
#### wormCat enrichment of genes with distant eQTL in each hotspot
  
#input
data_figS4  <- data.table::fread("processed_data/Supp_Data_3.tsv") %>% 
  dplyr::filter(Input=="Genes with distant eQTL in each merged hotspot") %>% 
  dplyr::mutate( cM1=Merged_hotspot_cM,
                Bonferroni_log10=-log10(Bonferroni)) %>% 
  dplyr::filter(!grepl("Unknown",Category)) %>% 
  tidyr::separate(cM1, into=c("cM2","cM3"),sep=",") %>% 
  dplyr::mutate( cM2=as.numeric(cM2)) %>% 
  dplyr::arrange(hotspot_Chr,cM2) %>% 
  dplyr::group_by(hotspot_Chr,cM2) %>% 
  dplyr::mutate(group_no = dplyr::group_indices( )) %>% 
  dplyr::mutate(cM_marker=paste(hotspot_Chr,Merged_hotspot_cM,sep=":" ))

data_figS4$level= gsub("category", "Category ", data_figS4$level)


# PLOT

fig_S4 <- ggplot(data_figS4,
                  aes(x=fct_reorder(cM_marker, group_no),y=Category, color=Bonferroni_log10)) +
  geom_point() + 
  theme_bw() + 
  theme(plot.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size=12,  color = "black"),
        legend.title =  ggplot2::element_text(size=12,  color = "black"),
        axis.title =  ggplot2::element_text(size=12,  color = "black"),
        axis.text =  ggplot2::element_text(size=12,  color = "black"),
        strip.text = ggplot2::element_text(size=12, vjust = 1,  color = "black"),
        strip.background = ggplot2::element_blank(), 
        text = ggplot2::element_text(family="Helvetica"),
        axis.text.x = element_text(size=12,angle = 90, hjust =  1,vjust = 0.5),
        axis.title.y = element_blank(),
        legend.position = "bottom" ) +      
  scale_color_gradient(high = "#D7263D", low = "#0072B2") +
  labs(x ="Hotspots (Chromosome:cM)",
       color = expression(-log[10](adjusted~italic(p)))) + 
  facet_grid( level~.,scales = "free",space = "free")  

ggsave(fig_S4, filename = paste("figures/p1_Supp_Fig4_HSenrich.png",sep = ""), units = "mm",height = 225, width = 170)




#### Figure S5 ####
 
# plot 
data_figS5 <- data.table::fread("processed_data/FileS14_hil2_enrich_finemap.tsv")
 
fine_tf_raw <- data_figS5 %>% 
  dplyr::mutate(finemap_variant=as.numeric(sub("(.*)(_)(.*)","\\3",finemarker))) %>% 
  dplyr::mutate(TF=ifelse(TF %in% c("modERN_TF","wormcat_TF"),"TF",
                          ifelse(TF=="wormcat_chroma_cof","cof",NA)) ) %>% 
  dplyr::mutate(impact=dplyr::case_when(
    impact=="hil-2" & BLOSUM <0  ~ "hil-2",
    BLOSUM <0 ~ "HIGH",
    BLOSUM >=0 ~ "LOW",
    TRUE ~ "Intergenic")) %>% 
  dplyr::mutate(cate2=dplyr::case_when(
    Category=="Development: somatic"  ~ "Development\nsomatic",
    Category=="Proteolysis proteasome: E3: F box"  ~ "Proteolysis proteasome\nE3: F box",
    Category=="Transcription factor: homeodomain" ~ "Transcription factor\nhomeodomain",
    TRUE ~ "Stress response\nheat")) %>% 
  dplyr::filter(Merged_hotspots_cM=="30.5, 31, 31.5, 32, 32.5" & grepl("IV",QTL))%>% 
  dplyr::mutate(transcript=paste(GeneName,transcript,sep = "\n"))
 

fine_tf_peak_raw <- fine_tf_raw %>% 
 # dplyr::filter(Category==cat) %>% 
  dplyr::distinct(transcript,QTL) %>% 
  dplyr::mutate(eQTL_peak=as.numeric(sub("(.*)(_)(.*)","\\3",QTL)))%>% 
  na.omit() %>% 
  dplyr::mutate(eQTL_peak= as.numeric(eQTL_peak)) 


hil2_fine_plt_list=list()

for(cat in unique(fine_tf_raw$Category)) {
  
  fine_tf <- fine_tf_raw   %>% dplyr::filter(Category==cat)
  
  fine_tf_peak  <- fine_tf   %>% 
    dplyr::distinct(transcript,QTL) %>% 
    dplyr::mutate(eQTL_peak=as.numeric(sub("(.*)(_)(.*)","\\3",QTL)))%>% 
    na.omit() %>% 
    dplyr::mutate(eQTL_peak= as.numeric(eQTL_peak)) 
  
  
  fine_tf_plt1 <- ggplot() +
    geom_point(data=subset(fine_tf, impact == "Intergenic"  ),aes(x=finemap_variant/1e6,y=finemap_log10p ),size=0.1, color="gray80") + 
    geom_point(data=subset(fine_tf, impact == "LOW"),aes(x=finemap_variant/1e6,y=finemap_log10p ),size=0.1, color="gray50") + 
    geom_point(data=subset(fine_tf, impact == "HIGH"),aes(x=finemap_variant/1e6,y=finemap_log10p ),size=0.1, color="orange")+
    geom_point(data=subset(fine_tf, !impact %in% c("Intergenic","LOW","HIGH")),
               aes(x=finemap_variant/1e6,y=finemap_log10p, color=impact  ),size=2,shape=18 )+
    geom_point(data = fine_tf_peak,aes(x=eQTL_peak/1e6,y = 0 ), size = 1 , alpha = 1, shape = 25, stroke = 0.5,color="black",fill="plum4") +
    #  facet_wrap(Category~transcript, scales = "free",ncol=4) + 
    facet_grid(cate2~transcript, scales = "free" ) + 
    theme_cust+
    theme( legend.text = ggplot2::element_text(size=12, color = "black",face = "italic"),
           plot.title = ggplot2::element_text(size=12,  color = "black"),
           legend.position = "none",
           #   axis.ticks.x = element_blank(),
           #   axis.text.x = element_blank(),
           panel.spacing = unit(0.01,"line"),
           plot.margin = unit(c(0, 2, 0, 2), "mm"),
           strip.text.y = ggplot2::element_text(size=11, vjust = 1,  color = "black" ),
           strip.text.x = ggplot2::element_text(size=11, vjust = 1,  color = "black",face = "italic")) +
    scale_color_manual(values = c( "red", "deepskyblue4", "blue",  "purple")) +
    #  scale_shape_manual(values = c("TF"=18,"cof"=15), guide = "none") +
    labs(x="Genomic position (Mb)",
         color="Candidate genes",
         y=expression(-log[10](italic(p))))  +
    scale_x_continuous(breaks=c(4,6,8,10,12,14,16) )+
    scale_y_continuous(breaks=c(0,5,10,15,20,25) )
  
  hil2_fine_plt_list[[cat]] <- fine_tf_plt1
}


fig_S5cd <- cowplot::plot_grid(hil2_fine_plt_list$`Stress response: heat`, hil2_fine_plt_list$`Transcription factor: homeodomain`,  
                               labels = c('', 'c'  ), 
                               label_size = 12, 
                               label_fontfamily="Helvetica",
                               rel_widths  = c(2.5,3),
                               
                               align = "h",
                               axis = "tb",
                               nrow = 1)


fig_S5 <- cowplot::plot_grid( hil2_fine_plt_list$`Proteolysis proteasome: E3: F box`, fig_S5cd, 
  labels = c('a', 'b'   ), 
  label_size = 12, 
  label_fontfamily="Helvetica",
  #  rel_widths  = c(2,3),
  #  rel_heights = 
  axis = "lr",
  # align = "v",
  
  nrow = 2)

ggsave(fig_S5, filename = paste("figures/p1_Supp_Fig5_hil2.png",sep = ""), units = "mm",height = 130, width = 170)

 
#### Figure S6 ####
data_figS6 <- data.table::fread("processed_data/FileS10_ABZ_GWA.tsv")  

data_figS6$dataset2 = factor( data_figS6$dataset,levels = c("202 strains","167 strains"))

fig_S6   <-   ggplot2::ggplot(data_figS6) +
  ggplot2::aes(x = POS/1e6, y = log10p) +
  ggplot2::scale_color_manual(values = c("0" = "black", 
                                         "1" = "red",
                                         "2" = "hotpink3",
                                         "3" = "darkorange1")) +
  ggplot2::scale_shape_manual(values = c("0" = 20, 
                                         "1" = 20,
                                         "2" = 20,
                                         "3" = 25))+
  ggplot2::scale_size_manual(values = c("0" = 0.5, 
                                        "1" = 0.5,
                                        "2" = 0.5,
                                        "3" = 1.5))+
  ggplot2::geom_hline(ggplot2::aes(yintercept = BF),
                      color = "gray", 
                      alpha = .75,  
                      size = 0.5) +
  ggplot2::geom_hline(ggplot2::aes(yintercept = EIGEN_CUTOFF),
                      color = "gray", 
                      alpha = .75,  
                      size = 0.5,
                      linetype = 2) +
  ggplot2::geom_point( ggplot2::aes(color= factor(EIGEN_SIG),
                                    shape= factor(EIGEN_SIG),
                                    size = factor(EIGEN_SIG)) ) +
  ggplot2::facet_grid( dataset2 ~ CHROM, scales = "free"  ) +
  theme_cust+   
  ggplot2::theme(legend.position = "none") +
  ggplot2::labs(x = "Genomic position (Mb)",
                y = expression(-log[10](italic(p)))) + 
  scale_y_continuous(expand = c(0, 1.5) ) 

ggsave(fig_S6, filename = paste("figures/p1_Supp_Fig6_abz_manh.png",sep = ""), units = "mm",height = 90, width = 170)





 

#### Figure S7 ####

# RIAILs eQTL map  
data_figS7 <- data.table::fread("processed_data/FileS15_RIAILs_eQTL.tsv")   

data_figS7$gene_Chr[data_figS7$gene_Chr=="MtDNA"] <- "M"
data_figS7$RILe_chr[data_figS7$RILe_chr=="MtDNA"] <- "M"

data_figS7$Chr_pos<- factor(data_figS7$gene_Chr,levels = c("M","X","V","IV","III","II","I"))
data_figS7$eChr_pos<- factor(data_figS7$RILe_chr,levels = c("I","II","III","IV","V","X","M"))


fig_S7a <- ggplot()  + 
  geom_point(data=subset(data_figS7,is.na(Detected_in_WIeQTL) ), aes(x=RILe_peak/1E6,y=gene_start,color=RILeQTL_classification),size=0.25,alpha=0.5)+ 
   geom_point(data=subset(data_figS7,Detected_in_WIeQTL=="Yes"), aes(x=RILe_peak/1E6,y=gene_start ,fill=RILeQTL_classification),color="black", alpha=0.5,shape=25) +
  scale_color_manual(values = c("Distant eQTL"="plum4","Local eQTL"="gold2")) +
  scale_fill_manual(values = c("Distant eQTL"="plum4","Local eQTL"="gold2")) +
  facet_grid(cols=vars(eChr_pos), rows=vars(Chr_pos), 
              scales = "free",  switch="both"
            ) +
  theme_cust +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(color = "grey", fill = NA, size = 1),
        plot.margin = unit(c(2, 2, 0, 4), "mm"),
        panel.spacing = unit(0.01,"line"),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position = "none")  +
  ylab("Gene position (Mb)") + 
  xlab("eQTL position (Mb)") +
  labs(fill="eQTL classification")
#fig_S10a

Sgt = ggplot_gtable(ggplot_build(fig_S7a))

 
Sgt$heights[7] = 0.3*Sgt$heights[9]


fig_S7a_gt <- ggplotify::as.ggplot(Sgt) 


# RIAILs hotspots

 
RILhotspot_pos <- data_figS7 %>% 
  dplyr::filter(!is.na(RILHotspot)) %>% 
  dplyr::mutate( Merged_RILHotspots_cM=ifelse(is.na(Merged_RILHotspots_cM) ,RILhotspot_cM,Merged_RILHotspots_cM),
                 merged_RILHotspots_QTL_count=ifelse(is.na(merged_RILHotspots_QTL_count) , RILHotspot_QTL_count,merged_RILHotspots_QTL_count))  %>% 
  dplyr::distinct(RILhotspot_Chr, RILhotspot_cM, RILHotspot,  Merged_RILHotspots_cM,merged_RILHotspots_QTL_count,Hotspots_Detected_in_WIeQTL ) %>% 
  dplyr::group_by(RILhotspot_Chr,Merged_RILHotspots_cM) %>% 
  dplyr::mutate(merged_Hotspot_center=mean(RILhotspot_cM))  %>% 
  dplyr::ungroup() %>% 
  dplyr::distinct(RILhotspot_Chr, Merged_RILHotspots_cM,merged_Hotspot_center, RILHotspot, merged_RILHotspots_QTL_count,Hotspots_Detected_in_WIeQTL )  %>% 
  dplyr::mutate(hotPOS = ifelse( Hotspots_Detected_in_WIeQTL=="Yes" , "overlap",RILHotspot)) %>% 
  dplyr::mutate(hotPOS=ifelse(is.na(hotPOS),"no",hotPOS))


fig_S7b <- ggplot(RILhotspot_pos,aes(x=merged_Hotspot_center,y=merged_RILHotspots_QTL_count)) +
  geom_bar(stat='identity',
           aes(color=hotPOS),size=0.5) +
  facet_grid(.~RILhotspot_Chr,scales = "free_x" ) +
  scale_color_manual(values = c( "gray69", "blue","red")) +
  geom_hline(yintercept=6,  color="black",alpha=0.8,linetype="dashed") + 
  ylab("Number of\ndistant eQTL") +
  xlab("Hotspot position (cM)")  + 
  theme_cust  +
  theme(legend.position = "none",
        #  axis.text.x = element_blank(),
          panel.spacing = unit( 0.5,"line"),
        axis.text.y =  ggplot2::element_text(size=12,  color = "black",angle = 90 ,vjust = 1,hjust=0.50),
        axis.text.x =  ggplot2::element_text(size=10,  color = "black" )) +
  scale_y_continuous(breaks=seq(0, 300, 150), limits = c(0,300),expand = c(0, 0))  +
  scale_x_continuous(breaks =  c(0, 15, 30, 45), limits = c(0,NA) ) 
 
 
fig_S7 <- cowplot::plot_grid(fig_S7a_gt, fig_S7b,  
                             labels = c('a', 'b'  ), 
                             label_size = 12, 
                             label_fontfamily="Helvetica",
                             rel_heights = c(3,1.5),
                            
                             align = "v",
                            #   axis = "lr",
                             nrow = 2)

ggsave(fig_S7, filename = paste("figures/p1_Supp_Fig7_RIAILs_eQTL.png",sep = ""), units = "mm",height = 160, width = 120)



#### Figure S8 ####

data_fig2a <- data.table::fread("processed_data/FileS4_eQTLmap.tsv")  

data_figS8 <- data.table::fread("processed_data/FileS16_eqtlStudies13_Wi207eQTL_common.tsv")  

data_figS8_count <- data_figS8 %>% 
  dplyr::distinct( ens_gene,transcript, eQTL_Chr,eQTL_peak,study) %>% 
  dplyr::group_by(ens_gene,transcript, eQTL_Chr,eQTL_peak ) %>% 
  dplyr::count(name = "study_count")  %>% 
  dplyr::left_join(data_fig2a) %>% 
  dplyr::mutate(study_count=study_count+1)


data_figS8_count$transcript_Chr[data_figS8_count$transcript_Chr=="MtDNA"] <- "M"
data_figS8_count$eQTL_Chr[data_figS8_count$eQTL_Chr=="MtDNA"] <- "M"

data_figS8_count$Chr_pos<- factor(data_figS8_count$transcript_Chr,levels = c("M","X","V","IV","III","II","I"))
data_figS8_count$eChr_pos<- factor(data_figS8_count$eQTL_Chr,levels = c("I","II","III","IV","V","X","M"))


fig_S8 <- ggplot(data=data_figS8_count)  + 
  geom_point(aes(x=eQTL_peak/1E6,y=transcript_start,color=eQTL_classification,size=study_count),alpha=0.5) +
  scale_color_manual(values = c("Distant eQTL"="plum4","Local eQTL"="gold2"), name="eQTL classification") +
  facet_grid(cols=vars(eChr_pos), rows=vars(Chr_pos), 
             scales = "free", 
             switch="both") +
  theme_cust +
  theme( panel.border = element_rect(color = "grey", fill = NA, size = 1),
        panel.spacing = unit(0.01,"line"),
        axis.text=element_blank(),
        axis.ticks=element_blank())  +
  ylab("Transcript position (Mb)") + 
  xlab("eQTL position (Mb)") +
  scale_size_continuous(breaks=seq(2, 14, 3),name="Number of detections\nin 14 conditions\nacross nine studies" )


ggsave(fig_S8, filename = paste("figures/p1_Supp_Fig8_eQTLacrossStudies.png",sep = ""), units = "mm",height = 130, width = 170)



#########