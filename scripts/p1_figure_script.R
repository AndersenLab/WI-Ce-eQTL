
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


#####

##########################################
#           Figure 1                     #
#           workflow                     #
##########################################
#### Figure 1a  #####
#### Origins  
#input

data_fig1a <- read.csv("processed_data/FileS1_distribution.csv", stringsAsFactors=FALSE)


# load world map 
world <- map_data("world") %>% 
  dplyr::filter(!region=="Antarctica")

strain_to_plot <- data_fig1a %>%
  dplyr::mutate(strain2=" ")


fig_1a <- ggplot()+ 
  geom_map(data=world, map=world,
           aes(x=long, y=lat, map_id=region),
           color="gray51", fill="white", 
           size=0.5 ) +
  geom_point(data = strain_to_plot, 
             aes(x=as.numeric(longitude), 
                 y=as.numeric(latitude)), 
             color = "black",
             shape = 21,   
             size = 0.1,
             alpha=0.1) + 
  ggrepel::geom_label_repel(data=strain_to_plot, 
                            aes(x=as.numeric(longitude), 
                                y=as.numeric(latitude),
                                label = strain2),
                            fill="#EC9E13",
                            box.padding = 0.005,
                            label.padding = 0.15, 
                            max.overlaps = Inf,
                            point.padding = 0.1,
                            segment.alpha =0.5,
                            segment.size=0.1,
                            size=0.1) + 
  theme_map() +
  theme(legend.position ="bottom",
        plot.margin = unit(c(0, 0, 0,0), "mm"),
        panel.background = element_rect(fill = 'white' , color = NA),
        plot.background = element_rect(fill = 'white' , color = NA),
        legend.text = element_text(size=12,  color = "black"),
        legend.title =  element_text(size=12,  color = "black")) 


#### Figure 1b  #####
#### Genetic relatedness   


load("processed_data/FileS2_tree.RData")

fig_1b <- ggtree::ggtree(FileS2, layout="rectangular", branch.length="rate", size = 0.3,color="gray51") +
  ggtree::geom_tippoint( aes(fill=label,  color=label),  
                         shape=21, 
                         size= 0.3) +
  coord_flip() + 
  scale_x_reverse() + 
  theme(legend.position =  "none",
        text=element_text(family="Helvetica") ,
        plot.margin = unit(c(0, 2, 0, 5), "mm")
  ) + 
  scale_fill_manual(values = c("RNAseq"="#EC9E13" )) +
  scale_color_manual(values = c("RNAseq"="#EC9E13" )) 




#### Figure 1c #####
#### sample age 
#input

data_fig1c <- data.table::fread("processed_data/FileS3_age.tsv")

# PLOT
fig_1c <- ggplot(data_fig1c,
                 aes(x=fct_reorder(strain, age),
                     y=age)) + 
  geom_point(size= 0.5,color="#EC9E13") + 
  geom_errorbar(aes(ymin=(age-age_sd), 
                    ymax=(age+age_sd)),
                size= 0.2, width = 0.2) +
  ylab("Estimated ages after\nhatching (hour)") +
  xlab("Strain") +
  theme_cust + 
  theme( axis.text.x = element_blank() ,
         axis.ticks.x = element_blank() ) 



## cowplot 
fig_1abc <- cowplot::plot_grid(fig_1a, fig_1b, fig_1c,
                              labels = c('a', 'b','c'), 
                              label_size = 12, 
                              label_fontfamily="Helvetica",
                            #  align = "v",
                              nrow = 3)

ggsave(fig_1abc, filename = paste("figures/p1_Fig_1_abc.png",sep = ""),  units = "mm",height = 160, width = 100)
 


##########################################
#           Figure 2                     #
#            eQTL map                    #
##########################################

#### Figure 2a ####
#eQTL map

data_fig2a <- data.table::fread("processed_data/FileS4_eQTLmap.tsv")  

data_fig2a$transcript_Chr[data_fig2a$transcript_Chr=="MtDNA"] <- "M"
data_fig2a$eQTL_Chr[data_fig2a$eQTL_Chr=="MtDNA"] <- "M"

data_fig2a$Chr_pos<- factor(data_fig2a$transcript_Chr,levels = c("M","X","V","IV","III","II","I"))
data_fig2a$eChr_pos<- factor(data_fig2a$eQTL_Chr,levels = c("I","II","III","IV","V","X","M"))


fig_2a <- ggplot(data=data_fig2a)  + 
  geom_point(aes(x=eQTL_peak/1E6,y=transcript_start,color=eQTL_classification),size=0.25,alpha=0.5) +
  scale_color_manual(values = c("Distant eQTL"="plum4","Local eQTL"="gold2")) +
  facet_grid(cols=vars(eChr_pos), rows=vars(Chr_pos), 
             scales = "free", 
             switch="both") +
  theme_cust +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(color = "grey", fill = NA, size = 1),
        plot.margin = unit(c(2, 2, 0, 0), "mm"),
        panel.spacing = unit(0.01,"line"),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position = "none")  +
  ylab("Transcript position (Mb)") + 
  xlab("eQTL position (Mb)") +
  labs(fill="eQTL classification")
 

gt = ggplot_gtable(ggplot_build(fig_2a))

gt$widths[18] = 0.3*gt$widths[10]

gt$heights[7] = 0.3*gt$heights[9]

 
fig_2a_gt <- ggplotify::as.ggplot(gt) 




#### Figure 2b #######
#input

data_fig2b <- data.table::fread("processed_data/FileS5_eQTLperTx.tsv") 

# PLOT

fig_2b <- ggplot(data_fig2b,aes(x=factor(n_dis),fill=both )) +
  geom_histogram(stat="count" ,color="plum4" ) +
  scale_fill_manual(values = c("single"="white","both"="gold2"))+
  theme_cust +
  theme(legend.position = "none",
        axis.text.y = element_text(size=12,angle = 90, hjust =  0.5,vjust = 1))+
  ylab("Number of traits") + 
  xlab("Number of distant eQTL per trait")   +
  geom_text(aes(label = distant_both,y=eQTL_count+100),color="gray6" )+
  scale_y_continuous(breaks=seq(0, 2500, 1000) )






#### Figure 2c ######
#### Heritability for 25,849 transcript expression traits  ###
data_fig2c <- data.table::fread("processed_data/FileS6_expression_H2h2.tsv")

 
heri_eQTL <- data_fig2c  %>% 
  dplyr::select(transcript, h2,H2) %>% 
  dplyr::mutate(eQTL = ifelse(transcript %in% data_fig2a$transcript,"Traits with eQTL", "Traits without eQTL"))  


# PLOT


fig_2c <- ggplot() + 
  geom_point(data=subset(heri_eQTL,eQTL=="Traits without eQTL"),aes(y=h2,x=H2), size=0.05,alpha=0.5,color="black" ) +
  geom_point(data=subset(heri_eQTL,eQTL=="Traits with eQTL"),aes(y=h2,x=H2), size=0.05,alpha=0.3,color="orange" ) +
  theme_cust + 
  theme(plot.margin = unit(c(2, 3, 0, 0), "mm")) +
  xlab(expression(italic(H^2))) +
  ylab(expression(italic(h^2))) +
  geom_abline(intercept=0,slope=1,colour="black",linetype=2) +
  scale_y_continuous(breaks=c(0,0.5,  1),expand = c(0, 0), limits = c(0,1))  +
  scale_x_continuous(breaks=c(0, 0.5, 1),expand = c(0, 0), limits = c(0,1))
 
  


 
#### Figure 2d  ####
##  Variance explained by h2 
 
var_h2 <- data_fig2a %>%
  dplyr::select(-biotype) %>% 
  dplyr::left_join(data_fig2c) %>% 
  dplyr::select(transcript, h2,H2,var_exp,eQTL_classification,Thresholds)

var_h2$eQTL_classification2 <- factor(var_h2$eQTL_classification,levels = c("Local eQTL","Distant eQTL"))


fig_2d <- ggplot(var_h2,aes(x=h2,y=var_exp ,color=eQTL_classification2)) + 
  geom_point(size=0.25,alpha=0.5)  +
  scale_color_manual(values = c("Distant eQTL"="plum4","Local eQTL"="gold2")) +
  theme_cust + 
  facet_grid(Thresholds~.) +
  theme(plot.margin = unit(c(2, 3, 1, 0), "mm"),
        panel.spacing = unit(1,"line"),
        legend.position = "none" )+
  ylab("VE") + 
  xlab(expression(paste( italic(h^2)))) +
  geom_abline(intercept=0,slope=1,colour="black",linetype=2) +
  scale_y_continuous(breaks=c(0,0.5, 1),expand = c(0, 0), limits = c(0,1))  +
  scale_x_continuous(breaks=c(0, 0.5, 1),expand = c(0, 0), limits = c(0,1))

 
#### Figure 2e ####
##Comparison of variance explained between detected local and distant eQTL ###

 
#input
 
var_h2_c <- data_fig2a %>%
  dplyr::left_join(data_fig2c) %>% 
  dplyr::select(transcript, h2,H2,var_exp,eQTL_classification,Thresholds)


var_mean <- var_h2_c %>% 
  dplyr::select(transcript ,var_exp,eQTL_classification)  %>% 
  dplyr::mutate(eQTL_classification=ifelse(eQTL_classification=="Local eQTL","Local","Distant")) %>% 
  dplyr::group_by(eQTL_classification) %>% 
  dplyr::mutate(varexp_mean=mean(var_exp)) %>% 
  dplyr::mutate(id=row_number(),
                mvar =ifelse(id==1,varexp_mean,NA))
var_mean$eQTL_classification2 <- factor(var_mean$eQTL_classification,levels = c("Local","Distant"))

#ggpubr::compare_means(var_exp ~ eQTL_classification,  data = var_mean)

# PLOT


fig_2e <- ggplot(var_mean,aes(x=eQTL_classification2,y=var_exp,color= eQTL_classification2))+
  geom_violin(width=1.1 ) +
  geom_boxplot(width=0.08, color="black",fill="grey", alpha=0.2,outlier.shape = NA,size=0.2)+ 
  geom_point(aes(x=eQTL_classification2,y=mvar),size=0.5,color="red") +
    scale_color_manual(values = c("Distant"="plum4","Local"="gold2")) +
  theme_cust +
  theme(legend.position = "none")+
  xlab("eQTL")+
  ylab("VE") +
  scale_y_continuous(breaks=seq(0, 1.2, 0.5), limits=c(0, 1.2)) +
  ggpubr::stat_compare_means(comparisons = list(c("Local","Distant")), label = "p.signif",label.y = c(1.01), method = "wilcox.test" )

  
# cowplot

   

fig_2ab <- cowplot::plot_grid(fig_2a_gt, fig_2b ,
                               labels = c('a', 'b'), 
                               rel_heights = c(2,1),
                               label_size = 12, 
                               label_fontfamily="Helvetica",
                               axis = "lr",
                              # align = "h",
                               nrow = 2)
 



fig_2cde <- cowplot::plot_grid(fig_2c,fig_2d, fig_2e,
                              labels = c('c','d', 'e'), 
                              rel_heights  = c(1,2,1),
                              label_size = 12, 
                              label_fontfamily="Helvetica",
                              axis = "lr",
                              align = "v",
                              nrow = 3)


fig_2 <- cowplot::plot_grid(fig_2ab, fig_2cde, 
                            rel_widths = c(2,1),
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            axis = "tb",
                            nrow =1)

ggsave(fig_2, filename = paste("figures/p1_Fig_2.png",sep = ""), units = "mm",height = 170, width = 170)






##########################################
#           Figure 3                    #
#           hotspots                    #
#                                       #   
##########################################
####  Figure 3  ###

 
hotspot_pos <- data_fig2a %>% 
  dplyr::distinct(hotspot_Chr, hotspot_cM, Hotspot, Hotspot_QTL_count) %>% na.omit()


fig_3 <- ggplot(hotspot_pos,aes(x=hotspot_cM,y=Hotspot_QTL_count)) +
  geom_bar(stat='identity',aes(fill=Hotspot)) +
  facet_grid(.~hotspot_Chr,scales = "free_x" ) +
  scale_fill_manual(values = c( "black", "red")) +
  geom_hline(yintercept=12,  color="gray69",alpha=0.8) + 
  ylab("Number of distant eQTL") +
  xlab("Hotspot position (cM)")  + 
  theme_cust  +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        panel.spacing = unit(0.1,"line"),
        axis.text.y =  ggplot2::element_text(size=12,  color = "black",angle = 90 ,vjust = 1,hjust=0.50)) +
  scale_y_continuous(breaks=seq(0, 200, 100), limits = c(0,200),expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0) ) 


ggsave(fig_3, filename = paste("figures/p1_Fig_3.png",sep = ""), units = "mm",height = 60, width = 170)


##########################################
#           Figure 4                    #
#           mediation  ben1            #
##########################################

######## figure 4a mediation   ben1 ########

data_fig4a <- data.table::fread("processed_data/FileS7_ben1_mediation.tsv")  

data_fig4a_sig <- data_fig4a %>% dplyr::filter(!q99_mediator == "Other genes")  
data_fig4a_other <- data_fig4a %>% dplyr::filter(q99_mediator == "Other genes")  

 
q99_med = quantile(data_fig4a$mediation_estimate, probs = 0.99)[[1]]


fig_4a <- ggplot() +
  geom_point(data=data_fig4a_other, aes( x=eQTL/1e6,
                                         y=mediation_estimate ), color =  "gray80" ,size=0.3) +
  geom_point(data=data_fig4a_sig, aes( x=eQTL/1e6,
                                       y=mediation_estimate, color = q99_mediator ),size=1 ) +
  geom_hline( yintercept = q99_med , color = "grey") +
  scale_size_continuous(range = c(2, 0))+
  labs(x = "Genomic position (Mb)", 
       y = "Mediation estimate",
       color = "Mediator gene" ) +
  theme_cust + 
  theme(legend.text = ggplot2::element_text(size=12,  color = "black",face = "italic" ),
        legend.margin = margin(),
        legend.spacing = unit(0.04, "cm")) +
  scale_color_manual(values=c( "lightskyblue2", "orange", "burlywood3", 
                               'cornflowerblue', "deepskyblue4", "mediumpurple4", 
                               "yellow4", "sienna4",'yellowgreen', 'chartreuse','lightpink'))+ 
  scale_y_continuous(expand = c(0, 0),limits = c(0,0.165))



######## figure 4b cor ########



data_fig4b <- data.table::fread("processed_data/FileS8_ben1_corr.tsv")

regressed_phe_all<-data_fig4b

ymax <- max(regressed_phe_all$original_pheno,regressed_phe_all$regressed_pheno)*1.05

ymin_raw <- min(regressed_phe_all$original_pheno,regressed_phe_all$regressed_pheno)

ymin <- ifelse(ymin_raw>0,ymin_raw*0.95,ymin_raw*1.05)


#ori
corr<-cor.test(regressed_phe_all$original_pheno, regressed_phe_all$original_exp, method = "pearson")

estim<-format(corr$estimate,digits=2, nsmall = 2)

pvalue<-format(ifelse(corr$p.value==0,-log(2.2e-16,base=10),-log(corr$p.value,base=10)),digits=2, nsmall = 2)

corr_plt <- ggplot(regressed_phe_all,aes(y=original_pheno,x=original_exp))+ 
  geom_point(size=0.3) +
  xlab(expression(paste(italic('ben-1'), " expression"))) +
  ylab("Animal length") +
  theme_cust +
  ylim(ymin,ymax) +
  ggtitle(paste("ρ",estim,sep = " : ")) +
  theme(plot.title = ggplot2::element_text(size=12,  color = "black",hjust = 0,vjust = 1.5))

#corr_plt


#reged  

corr_reg<-cor.test(regressed_phe_all$regressed_pheno, regressed_phe_all$original_exp, method = "pearson")

estim_reg <- format(corr_reg$estimate,digits=2, nsmall = 2)


corr_reg_plt <- ggplot(regressed_phe_all,aes(y=regressed_pheno,x=original_exp))+ 
  geom_point(size=0.3)+
  xlab(expression(paste(italic('ben-1'), " expression"))) +
  ylab("Animal length\n(regressed)") +
  theme_cust +
  ylim(ymin,ymax) +
  ggtitle(paste("ρ",estim_reg,sep = " : ")) +
  theme(plot.title = ggplot2::element_text(size=12,  color = "black",hjust = 0,vjust = 1.5))



fig_4b <- cowplot::plot_grid(corr_plt, corr_reg_plt, 
                              label_fontfamily="Helvetica",
                              align = "v",
                              axis = "tblr",
                              nrow = 2)


######## figure 4c reg mapping ########

sig_mediator <- data_fig4a %>%  
  dplyr::filter(!q99_mediator == "Other genes") %>% 
  dplyr::select(transcript=mediator_transcript,mediator_gene)


data_fig4c <- data.table::fread("processed_data/FileS9_reg_mapping167.tsv")  

data_figS9 <- data.table::fread("processed_data/FileS10_ABZ_GWA.tsv")  

oriSig <-  subset(data_figS9, dataset == "167 strains" & marker == "III_3539640")$log10p

data_fig4c_sig <- data_fig4c %>% 
  dplyr::filter(transcript %in% sig_mediator$transcript) %>% 
  dplyr::left_join(sig_mediator)


data_fig4c_other <- data_fig4c %>% dplyr::filter(!transcript %in% sig_mediator$transcript) 

fig_4c <- ggplot() +
  geom_point(data = data_fig4c_other,aes(x=transcript_start/1e6,y=log10p), size=0.3,color="gray80" )+
  geom_point(data = data_fig4c_sig,  aes(x=transcript_start/1e6,y=log10p, color=mediator_gene), size=1 ) +
  scale_color_manual(values=c( "lightskyblue2", "orange", "burlywood3", 
                               'cornflowerblue', "deepskyblue4", "mediumpurple4", 
                               "yellow4", "sienna4",'yellowgreen', 'chartreuse','lightpink' )) +
  theme_cust +
  ggplot2::theme(legend.position = "none",
                 axis.ticks.x = element_blank(),
                 panel.spacing = unit(0.04,"line"),
                 axis.text.x = element_blank()) +
  ylab(expression(-log[10](italic(p)))) + 
  xlab("Gene position (Mb)") +
  facet_grid(.~transcript_Chr,scales="free") +
  geom_hline(  yintercept = as.numeric(oriSig) ,
               color = "orange")   




######## cow figure4 ####
fig_4ab <- cowplot::plot_grid(fig_4a, fig_4b,  
                              labels = c('', 'b'), 
                              rel_widths = c(1.5, 1),
                              label_size = 12, 
                              label_fontfamily="Helvetica",
                              axis = "tb",
                              nrow = 1)

fig_4 <- cowplot::plot_grid(fig_4ab, fig_4c,
                            labels = c('a',  'c'), 
                            rel_heights = c(1.5,0.6),
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            axis = "lr",
                            nrow = 2)

ggsave(fig_4, filename = paste("figures/p1_Fig_4.png",sep = ""), units = "mm",height = 120, width = 140)






##########################################
#           Figure 5                    #
#       mediation seven traits          #
##########################################
#input
data_fig5_gwa <- data.table::fread("processed_data/FileS11_7traits_GWA.tsv")  
data_fig5_med <- data.table::fread("processed_data/FileS12_7traits_MED.tsv")  


##fig_5a
fig_5a_manh <- manhaplot("telomere_resids") 

fig_5a_med <- medaplot("telomere_resids")

fig_5a <- cowplot::plot_grid(fig_5a_manh, fig_5a_med,
                             rel_heights = c(1 ,0.8),
                             label_size = 12, 
                             label_fontfamily="Helvetica",
                             # axis = "lr",
                             align = "v",
                             nrow = 2)

##fig_5b
fig_5b_manh <- manhaplot("arsenic_pc1")+
  theme(axis.title.y = element_blank()) + 
  scale_y_continuous(expand = c(0, 0.4),breaks = seq(0,6,2) )

fig_5b_med <- medaplot("arsenic_pc1")+
  theme(axis.title.y = element_blank())

fig_5b <- cowplot::plot_grid(fig_5b_manh, fig_5b_med,
                             rel_heights = c(1 ,0.8),
                             label_size = 12, 
                             label_fontfamily="Helvetica",
                             # axis = "lr",
                             align = "v",
                             nrow = 2)

##fig_5c
fig_5c_manh <- manhaplot("zinc_norm_EXT") + 
  scale_y_continuous(expand = c(0, 0.2),breaks = seq(0,6,2) )+
  theme(axis.title.y = element_blank()) 

fig_5c_med <- medaplot("zinc_norm_EXT") +
  theme(axis.title.y = element_blank())  

fig_5c <- cowplot::plot_grid(fig_5c_manh, fig_5c_med,
                             rel_heights = c(1 ,0.8),
                             label_size = 12, 
                             label_fontfamily="Helvetica",
                             # axis = "lr",
                             align = "v",
                             nrow = 2)


##fig_5d
fig_5d_manh <- manhaplot("etoposide_median.TOF")+ 
  scale_y_continuous(expand = c(0, 0.4),breaks = seq(0,8,4) ) 

fig_5d_med <- medaplot("etoposide_median.TOF") 

fig_5d <- cowplot::plot_grid(fig_5d_manh, fig_5d_med,
                             rel_heights = c(1 ,0.8),
                             label_size = 12, 
                             label_fontfamily="Helvetica",
                             # axis = "lr",
                             align = "v",
                             nrow = 2)


##fig_5e
fig_5e_manh <- manhaplot("propionicL1survival") + 
  scale_y_continuous(expand = c(0, 0.4),breaks = seq(0,8,4) ) +
  theme(axis.title.y = element_blank())

fig_5e_med <- medaplot("propionicL1survival")  + 
  scale_y_continuous( breaks = seq(0,0.2,0.1) ) +
  theme(axis.title.y = element_blank())


fig_5e <- cowplot::plot_grid(fig_5e_manh, fig_5e_med,
                             rel_heights = c(1 ,0.8),
                             label_size = 12, 
                             label_fontfamily="Helvetica",
                             # axis = "lr",
                             align = "v",
                             nrow = 2)



##fig_5f
fig_5f_manh <- manhaplot("abamectin_norm.n") + 
  scale_y_continuous(expand = c(0, 0.4),breaks = seq(0,15,5) )  

fig_5f_med <- medaplot_aba("abamectin_norm.n")  + 
  scale_y_continuous( breaks = seq(0,0.22,0.1) )  


fig_5f <- cowplot::plot_grid(fig_5f_manh, fig_5f_med,
                             rel_heights = c(1 ,1.8),
                             label_size = 12, 
                             label_fontfamily="Helvetica",
                             # axis = "lr",
                             align = "v",
                             nrow = 2)


##fig_5g
fig_5g_manh <- manhaplot("broods") + 
  scale_y_continuous(expand = c(0, 0.4),breaks = seq(0,15,5) ) +
  theme(axis.title.y = element_blank())

fig_5g_med <- medaplot_bro("broods")  + 
  scale_y_continuous( breaks = seq(0,0.25,0.1) ) +
  theme(axis.title.y = element_blank())


fig_5g <- cowplot::plot_grid(fig_5g_manh, fig_5g_med,
                             rel_heights = c(1 ,1.8),
                             label_size = 12, 
                             label_fontfamily="Helvetica",
                             # axis = "lr",
                             align = "v",
                             nrow = 2)





##cow fig5

fig_manh_med_abc <- cowplot::plot_grid(fig_5a,fig_5b,fig_5c,
                                       labels = c('a', 'b', 'c'), 
                                       rel_widths = c(2.3 ,2,1.2),
                                       label_size = 12, 
                                       label_fontfamily="Helvetica",
                                       nrow = 1)




fig_manh_med_de <- cowplot::plot_grid(fig_5d ,fig_5e,
                                      labels = c('d', 'e'), 
                                      label_size = 12, 
                                      label_fontfamily="Helvetica",
                                      nrow = 1)


fig_manh_med_fg <- cowplot::plot_grid(fig_5f ,fig_5g,
                                      labels = c('f', 'g'), 
                                      rel_widths = c(1,1.1),
                                      label_size = 12, 
                                      label_fontfamily="Helvetica",
                                      axis = "tb",
                                      nrow = 1)


fig_5  <- cowplot::plot_grid(fig_manh_med_abc,
                             fig_manh_med_de,
                             fig_manh_med_fg,
                             rel_heights = c(1,1,1.8),
                             label_size = 12, 
                             label_fontfamily="Helvetica",
                             nrow = 3)




ggsave(fig_5, filename = paste("figures/p1_Fig_5.png",sep = ""), units = "mm",height = 180, width = 160)




##############################
##### Supplementary Figs #####
##############################
 

#### Figure S2  ####               

#### Expression and genetic distance  


#input
data_figS2 <- data.table::fread("processed_data/FileS13_distance_exp_genet.tsv") 


# PLOT
dist_exp_gene <- data_figS2 %>% 
  dplyr::mutate(domain = ifelse(domain=="All","Whole genome",domain),
                n_variants =  format(n_variants,big.mark=",",scientific=FALSE),
                n_tx =  format(n_tx,big.mark=",",scientific=FALSE)) %>% 
  dplyr::mutate(Domain = paste0(domain,"\nTranscripts: ",n_tx, "\nGenetic variants: ", n_variants)) %>% 
  dplyr::mutate(Domain = ifelse(Domain=="Divergent\nTranscripts: 13,432\nGenetic variants: 584,731",
                                "Non-swept\nTranscripts: 13,432\nGenetic variants: 584,731",Domain))

dist_exp_gene$Domain2 <- factor(dist_exp_gene$Domain,levels = c("Whole genome\nTranscripts: 22,268\nGenetic variants: 851,105",
                                                                "Tip\nTranscripts:  1,308\nGenetic variants:  64,010",
                                                                "Arm\nTranscripts:  7,135\nGenetic variants: 404,638",
                                                                "Center\nTranscripts: 13,813\nGenetic variants: 381,740",
                                                                "Swept\nTranscripts:  8,836\nGenetic variants: 266,374",
                                                                "Non-swept\nTranscripts: 13,432\nGenetic variants: 584,731"))


fig_S2 <- ggplot2::ggplot(dist_exp_gene,
                          ggplot2::aes(x=dist_gene,y=dist_exp,color=genotype)) + 
  ggplot2::geom_point(size=1,alpha=0.5 ) +
  scale_color_manual(values = c("swept"="darkorange1","divergent"="deepskyblue4")) +
  theme_cust +
  theme(legend.position = "bottom") +
  xlab("Genetic distance") +
  ylab("Expression distance")   +
  labs(color="Genotype")+
  facet_wrap(.~Domain2, nrow = 3,scales = "free") +
  geom_text(aes(label = cv_dist_exp,y=cv_dist_exp_y, x=cv_dist_exp_x+5 ),color="gray6" )+
  geom_text(aes(label = cv_dist_gene,y=cv_dist_gene_y, x=cv_dist_gene_x-5 ),color="gray6" )



ggsave(fig_S2, filename = paste("figures/p1_Supp_Fig2_distance_geexp.png",sep = ""),  units = "mm",height = 225, width = 170)


#### Figure S3  ####   
##  
#### wormCat enrichment of all genes with eQTL

#input
data_figS3  <- data.table::fread("processed_data/FileS14_all_eQTL_wormCat.tsv") %>% 
  dplyr::filter(!grepl("Unknown",Category))

data_figS3$level= gsub("category", "Category ", data_figS3$level)
# PLOT


fig_S3 <- ggplot(data_figS3,
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

ggsave(fig_S3, filename = paste("figures/p1_Supp_Fig3_eQTLenrich.png",sep = ""), units = "mm",height = 180, width = 170)


#### Figure S4  ####   
##  
#### A histogram showing the distribution of linkage disequilibrium (LD) 

#input
data_figS4 <- data.table::fread("processed_data/FileS15_eQTL_LD.tsv")

# PLOT
fig_S4 <- ggplot(data_figS4,aes(LD)) + 
  geom_histogram(aes(fill=chroms),alpha=0.8,binwidth = 0.02,color="black") + 
  theme_cust +
  theme( legend.position = "bottom",
         legend.background = element_rect(),
         legend.title = element_blank())+
  ylab("Number of LD values") + 
  xlab("Pairwise LD among eQTL") +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  facet_grid(.~ld_betw)  

ggsave(fig_S4, filename = paste("figures/p1_Supp_Fig4_eQTL_LD.png",sep = ""), units = "mm",height = 100, width = 140)


#### Figure S5 ####
#### TajimaD

#input
data_figS5  <- data.table::fread("processed_data/FileS16_theta_pi_tjd.tsv")


SNP_0Miss_TjD_0.5cM_bin_median <- data_figS5 %>% 
  dplyr::group_by(chrom,cM,pop_stats) %>% 
  dplyr::mutate(median_value=median(value)) %>% 
  dplyr::select(-start,-end,-value) %>% 
  dplyr::distinct()


# PLOT
plt_theta <- ggplot() +
  geom_point(data=subset(data_figS5,pop_stats=="theta" & DiseQTL_Hotspot=="No"),aes(x=as.numeric(cM),y=as.numeric(value) ), size=0.1,alpha=0.5,color="gray70")+
  geom_point(data=subset(data_figS5,pop_stats=="theta" & DiseQTL_Hotspot=="Yes"),aes(x=as.numeric(cM),y=as.numeric(value) ), size=0.1,alpha=0.5,color="red") +
  geom_line(data=subset(SNP_0Miss_TjD_0.5cM_bin_median,pop_stats=="theta"  ),
            aes(x=as.numeric(cM),y=as.numeric(median_value) ), size= 0.5,color="black") +
  facet_grid(.~chrom,scales="free") + 
  labs(y=expression(paste("Watterson’s theta (", italic(θ),")")),
       x="Genomic position (cM)")+
  theme_cust  +
  theme(legend.position = "none",
        axis.text.x = element_blank())  

plt_pi <- ggplot() +
  geom_point(data=subset(data_figS5,pop_stats=="pi" & DiseQTL_Hotspot=="No"),aes(x=as.numeric(cM),y=as.numeric(value)),
             size=0.1,alpha=0.5,color="gray70")+
  geom_point(data=subset(data_figS5,pop_stats=="pi" & DiseQTL_Hotspot=="Yes"),aes(x=as.numeric(cM),y=as.numeric(value)),
             size=0.1,alpha=0.5,color="red") +
  geom_line(data=subset(SNP_0Miss_TjD_0.5cM_bin_median,pop_stats=="pi"  ),
            aes(x=as.numeric(cM),y=as.numeric(median_value) ), size= 0.5,color="black")+
  facet_grid(.~chrom,scales="free") + 
  labs(y="Pi",
       x="Genomic position (cM)")+
  theme_cust  +
  theme(legend.position = "none",
        axis.text.x = element_blank()) 

plt_td <- ggplot() +
  geom_point(data=subset(data_figS5,pop_stats=="td" & DiseQTL_Hotspot=="No"),aes(x=as.numeric(cM),y=as.numeric(value) ), size=0.1,alpha=0.5,color="gray70")+
  geom_point(data=subset(data_figS5,pop_stats=="td" & DiseQTL_Hotspot=="Yes"),aes(x=as.numeric(cM),y=as.numeric(value) ), size=0.1,alpha=0.5,color="red") +
  geom_line(data=subset(SNP_0Miss_TjD_0.5cM_bin_median,pop_stats=="td"  ),
            aes(x=as.numeric(cM),y=as.numeric(median_value) ), size= 0.5,color="black") +
  facet_grid(.~chrom,scales="free") + 
  labs(y=expression(paste("Tajima’s ", italic(D))),
       x="Genomic position (cM)")+
  theme_cust  +
  theme(legend.position = "none",
        axis.text.x = element_blank()) 


 
  
fig_S5 <- cowplot::plot_grid(plt_theta, plt_pi, plt_td, 
                              labels = c('a', 'b',"c" ), 
                               label_size = 12, 
                              label_fontfamily="Helvetica",
                              axis = "lr",
                              nrow = 3)

ggsave(fig_S5, filename = paste("figures/p1_Supp_Fig5_TJD.png",sep = ""), units = "mm",height = 150, width = 160)


#### Figure S6 ####
#### wormCat enrichment of genes with distant eQTL in each hotspot
 

#input
data_figS6  <- data.table::fread("processed_data/FileS17_hotspot_gene_wormCat.tsv") %>% 
  dplyr::mutate(chr=sub("(.*)(:)(.*)","\\1",cM_marker),
                cM=sub("(.*)(:)(.*)","\\3",cM_marker)) %>% 
  dplyr::filter(!grepl("Unknown",Category)) %>% 
  dplyr::mutate( cM=as.numeric(cM)) %>% 
  dplyr::arrange(chr,cM) %>% 
  dplyr::group_by(chr,cM) %>% 
  dplyr::mutate(group_no = dplyr::group_indices( ))

data_figS6$level= gsub("category", "Category ", data_figS6$level)


# PLOT

fig_S6 <- ggplot(data_figS6,
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

ggsave(fig_S6, filename = paste("figures/p1_Supp_Fig6_HSenrich.png",sep = ""), units = "mm",height = 225, width = 170)




#### Figure S7 ####

#### TFs in hotspots 
 
#input
data_figS7  <- data.table::fread("processed_data/FileS18_hotspot_TFs.tsv") 



data_fig2a <- data.table::fread("processed_data/FileS4_eQTLmap.tsv")  

de_hotspot <- data_fig2a  %>% 
  dplyr::filter(Hotspot=="Yes") %>% 
  dplyr::select(cM_Chr=hotspot_Chr,  cM=hotspot_cM, eQTL_Hotspot=Hotspot ) %>% 
  dplyr::distinct() %>% 
  na.omit()



tf_bin <- data_figS7 %>% 
  dplyr::mutate(type=ifelse(grepl("TF",source),"Transcription factor","Chromatin cofactor")) %>% 
  dplyr::group_by(type,cM_Chr,cM) %>% 
  dplyr::count(name="count") %>% 
  dplyr::left_join(de_hotspot) %>% 
  dplyr::mutate(eQTL_Hotspot=ifelse(is.na(eQTL_Hotspot),"No",eQTL_Hotspot))

# PLOT

fig_S7 <- ggplot(tf_bin,aes(x=cM,y=count)) +
  geom_bar(stat='identity',aes(fill=eQTL_Hotspot)) +
  facet_grid(type ~cM_Chr,scales = "free") +
  scale_fill_manual(values = c( "black", "red")) +
  ylab("Counts") +
  xlab("cM ")  + 
  theme_cust  +
  theme(legend.position = "none",
        axis.text.x = element_blank())   

 

ggsave(fig_S7, filename = paste("figures/p1_Supp_Fig7_dshotspot_tf_cof.png",sep = ""), units = "mm",height = 120, width = 170)


#### Figure S8 ####
#fine mappings 

 
for( i in 1:length(grep("FileS19_hotspotFine_TFs",list.files("processed_data/"), value = T))) {
  
#  i=13
  
  fine_file <-(grep("FileS19_hotspotFine_TFs",list.files("processed_data/"), value = T)[i])
  
  cMmarker <- sub("(FileS19_hotspotFine_TFs_)(.*)(.tsv)","\\2",fine_file)
  
  data_figS8 <- data.table::fread(paste0("processed_data/",fine_file))  
  
  
  tf_cof <- unique(dplyr::filter(data_figS8,!TF=="")$fine_ext_gene)
  
  
  for (tfc in tf_cof) {
    
   # tfc<- "hil-2"
    
    data_figS8_spe <- data_figS8 %>% 
      dplyr::filter((!TF=="") & fine_ext_gene==tfc) %>% 
      dplyr::distinct(transcript)
   
  
  fine_tf_raw <- data_figS8 %>% 
    dplyr::filter(transcript %in% data_figS8_spe$transcript) %>% 
    dplyr::mutate(finemap_variant=as.numeric(sub("(.*)(_)(.*)","\\3",finemarker))) %>% 
    dplyr::filter(cM_marker==cMmarker ) %>% 
    dplyr::mutate(TF=ifelse(TF %in% c("modERN_TF","wormcat_TF"),"TF",
                            ifelse(TF=="wormcat_chroma_cof","cof",NA))) %>% 
  dplyr::mutate(impact=ifelse(impact %in% c("Intergenic", "LOW",tfc),impact,"HIGH"))
 
  
  fine_tf_peak_raw <- fine_tf_raw %>% 
    dplyr::distinct(transcript,QTL) %>% 
    dplyr::mutate(eQTL_peak=as.numeric(sub("(.*)(_)(.*)","\\3",QTL)))%>% 
    na.omit() %>% 
    dplyr::mutate(eQTL_peak= as.numeric(eQTL_peak))
  
  ntrait <- length(unique(fine_tf_peak_raw$transcript))
  
 
  
  
  if(ntrait>24){
    
    
    half_trait <- sample( fine_tf_peak_raw$transcript, ntrait/2)
    
    
    fine_tf_1 <- fine_tf_raw %>% 
      dplyr::filter(transcript %in% half_trait )
    
    fine_tf_peak_1  <- fine_tf_peak_raw %>% 
      dplyr::filter(transcript %in% half_trait )
    
    
    
    
    fine_tf_plt_1 <- ggplot() +
      geom_point(data=subset(fine_tf_1, impact == "Intergenic"),aes(x=finemap_variant/1e6,y=finemap_log10p ),size=0.1, color="gray80") + 
      geom_point(data=subset(fine_tf_1, impact == "LOW"),aes(x=finemap_variant/1e6,y=finemap_log10p ),size=0.1, color="gray50") + 
      geom_point(data=subset(fine_tf_1, impact == "HIGH"),aes(x=finemap_variant/1e6,y=finemap_log10p ),size=0.1, color="orange")+
      geom_point(data=subset(fine_tf_1, !impact %in% c("Intergenic","LOW","HIGH")),aes(x=finemap_variant/1e6,y=finemap_log10p, color=impact ),size=2,shape=18 )+
      geom_point(data = fine_tf_peak_1,aes(x=eQTL_peak/1e6,y = 0 ), size = 1 , alpha = 1, shape = 25, stroke = 0.5,color="black",fill="plum4") +
      facet_wrap(.~transcript, scales = "free",ncol=4) + 
      theme_cust+
      theme( legend.text = ggplot2::element_text(size=12, color = "black",face = "italic"),
             legend.position = "bottom",
             plot.title = ggplot2::element_text(size=12,  color = "black"),
             axis.ticks.x = element_blank(),
             axis.text.x = element_blank(),
             panel.spacing = unit(0.01,"line"),
             plot.margin = unit(c(0, 2, 0, 2), "mm"),
             strip.text = ggplot2::element_text(size=11, vjust = 1,  color = "black",face = "italic")) +
      scale_color_manual(values = c( "red", "deepskyblue4", "blue",  "purple")) +
      labs(x="Genomic position (Mb)",
           color="Candidate genes",
           y=expression(-log[10](italic(p)))) +
      ggtitle( paste("Chromosome",sub("(.*)(_)(.*)","\\1:",cMmarker),sub("(.*)(_)(.*)","\\3",cMmarker),"cM" ))
    
    
 
    
    # plot_height <- ifelse((ceiling(ntrait/4))*40<200, ((ceiling(ntrait/4))*40+10),210)
    
    ggsave(fine_tf_plt_1, filename = paste("figures/p1_Supp_Fig8_",cMmarker,"_",tfc,"_1.png",sep = ""), units = "mm",height = 210, width = 175)
    
    
    
    
    fine_tf_2 <- fine_tf_raw %>% 
      dplyr::filter(!transcript %in% half_trait )
    
    fine_tf_peak_2  <- fine_tf_peak_raw %>% 
      dplyr::filter(!transcript %in% half_trait )
    
    
    
    
    fine_tf_plt_2 <- ggplot() +
      geom_point(data=subset(fine_tf_2, impact == "Intergenic"),aes(x=finemap_variant/1e6,y=finemap_log10p ),size=0.1, color="gray80") + 
      geom_point(data=subset(fine_tf_2, impact == "LOW"),aes(x=finemap_variant/1e6,y=finemap_log10p ),size=0.1, color="gray50") + 
      geom_point(data=subset(fine_tf_2, impact == "HIGH"),aes(x=finemap_variant/1e6,y=finemap_log10p ),size=0.1, color="orange")+
      geom_point(data=subset(fine_tf_2, !impact %in% c("Intergenic","LOW","HIGH")),
                 aes(x=finemap_variant/1e6,y=finemap_log10p, color=impact  ),size=2,shape=18)+
      geom_point(data = fine_tf_peak_2,aes(x=eQTL_peak/1e6,y = 0 ), size = 1 , alpha = 1, shape = 25, stroke = 0.5,color="black",fill="plum4") +
      facet_wrap(.~transcript, scales = "free",ncol=4) + 
      theme_cust+
      theme( legend.text = ggplot2::element_text(size=12, color = "black",face = "italic"),
             plot.title = ggplot2::element_text(size=12,  color = "black"),
             legend.position = "bottom",
             axis.ticks.x = element_blank(),
             axis.text.x = element_blank(),
             panel.spacing = unit(0.01,"line"),
             plot.margin = unit(c(0, 2, 0, 2), "mm"),
             strip.text = ggplot2::element_text(size=11, vjust = 1,  color = "black",face = "italic")) +
      scale_color_manual(values = c( "red", "deepskyblue4", "blue",  "purple")) +
      #  scale_shape_manual(values = c("TF"=18,"cof"=16)) +
      labs(x="Genomic position (Mb)",
           color="Candidate genes",
           y=expression(-log[10](italic(p)))) +
      ggtitle( paste("Chromosome",sub("(.*)(_)(.*)","\\1:",cMmarker),sub("(.*)(_)(.*)","\\3",cMmarker),"cM" ))
    
    
  
    
    
    #plot_height <- ifelse((ceiling(ntrait/4))*40<200, ((ceiling(ntrait/4))*40+10),210)
    
     ggsave(fine_tf_plt_2, filename = paste("figures/p1_Supp_Fig8_",cMmarker,"_",tfc,"_2.png",sep = ""), units = "mm",height = 210, width = 175)
    
    
  } else {
    
    
    fine_tf <- fine_tf_raw  
    
    fine_tf_peak  <- fine_tf_peak_raw  
    
    fine_tf_plt <- ggplot() +
      geom_point(data=subset(fine_tf, impact == "Intergenic"),aes(x=finemap_variant/1e6,y=finemap_log10p ),size=0.1, color="gray80") + 
      geom_point(data=subset(fine_tf, impact == "LOW"),aes(x=finemap_variant/1e6,y=finemap_log10p ),size=0.1, color="gray50") + 
      geom_point(data=subset(fine_tf, impact == "HIGH"),aes(x=finemap_variant/1e6,y=finemap_log10p ),size=0.1, color="orange")+
      # geom_segment( data = subset(fine_tf, !impact %in% c("Intergenic","LOW","HIGH")), arrow = arrow(length = unit(5, "points")), 
      #              aes(x = finemap_variant/1e6,xend = finemap_variant/1e6, y = finemap_log10p, yend = finemap_log10p - 0.5 , color=impact) ) + 
      geom_point(data=subset(fine_tf, !impact %in% c("Intergenic","LOW","HIGH")),
                 aes(x=finemap_variant/1e6,y=finemap_log10p, color=impact  ),size=2,shape=18 )+
      geom_point(data = fine_tf_peak,aes(x=eQTL_peak/1e6,y = 0 ), size = 1 , alpha = 1, shape = 25, stroke = 0.5,color="black",fill="plum4") +
      facet_wrap(.~transcript, scales = "free",ncol=4) + 
      theme_cust+
      theme( legend.text = ggplot2::element_text(size=12, color = "black",face = "italic"),
             plot.title = ggplot2::element_text(size=12,  color = "black"),
             legend.position = "bottom",
             axis.ticks.x = element_blank(),
             axis.text.x = element_blank(),
             panel.spacing = unit(0.01,"line"),
             plot.margin = unit(c(0, 2, 0, 2), "mm"),
             strip.text = ggplot2::element_text(size=11, vjust = 1,  color = "black",face = "italic")) +
      scale_color_manual(values = c( "red", "deepskyblue4", "blue",  "purple")) +
      #  scale_shape_manual(values = c("TF"=18,"cof"=15), guide = "none") +
      labs(x="Genomic position (Mb)",
           color="Candidate genes",
           y=expression(-log[10](italic(p)))) +
      ggtitle( paste("Chromosome",sub("(.*)(_)(.*)","\\1:",cMmarker),sub("(.*)(_)(.*)","\\3",cMmarker),"cM" ))
    
    
    
    
     plot_height <- ifelse((ceiling(ntrait/4))*40<200, ((ceiling(ntrait/4))*40+10),210)
    
    ggsave(fine_tf_plt, filename = paste("figures/p1_Supp_Fig8_",cMmarker,"_",tfc,".png",sep = ""), units = "mm",height = plot_height, width = 175)
    
 
    
  }
  
  
  
  
  
  }
  
  
}


 





#### Figure S9 ####
data_figS9 <- data.table::fread("processed_data/FileS10_ABZ_GWA.tsv")  

data_figS9$dataset2 = factor( data_figS9$dataset,levels = c("202 strains","167 strains"))

fig_S9   <-   ggplot2::ggplot(data_figS9) +
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
  ggplot2::geom_point( ggplot2::aes(color= factor(EIGEN_SIG),shape= factor(EIGEN_SIG),size = factor(EIGEN_SIG)) ) +
  ggplot2::facet_grid( dataset2 ~ CHROM, scales = "free"  ) +
  theme_cust+   
  ggplot2::theme(legend.position = "none") +
  ggplot2::labs(x = "Genomic position (Mb)",
                y = expression(-log[10](italic(p)))) + 
  scale_y_continuous(expand = c(0, 1.5) ) 

ggsave(fig_S9, filename = paste("figures/p1_Supp_Fig9_abz_manh.png",sep = ""), units = "mm",height = 90, width = 170)







 
 


