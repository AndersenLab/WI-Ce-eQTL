
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

##########################################
#           Figure 1                     #
#           workflow                     #
##########################################
#### Figure 1a  #####
#### Origins  
#input
WI_geo_20220111 <- data.table::fread("raw_data/WI_geo_20220111.txt") %>% 
  dplyr::select( isotype, geo)

data_fig1a <- read.csv("processed_data/FileS1_distribution.csv", stringsAsFactors=FALSE)  %>% 
  dplyr::left_join(WI_geo_20220111)

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
                                label = strain2,
                                fill = geo),
                            box.padding = 0.005,
                            label.padding = 0.15, 
                            max.overlaps = Inf,
                            point.padding = 0.1,
                            segment.alpha =0.5,
                            segment.size=0.1,
                            size=0.1) + 
  theme_map() +
  theme(  legend.position ="none",
        plot.margin = unit(c(0, 0, 0,0), "mm"),
        panel.background = element_rect(fill = 'white' , color = NA),
        plot.background = element_rect(fill = 'white' , color = NA),
        legend.text = element_text(size=12,  color = "black"),
        legend.title =  element_text(size=12,  color = "black")) +
   scale_fill_manual(values=rodon_color)  


#### Figure 1c  #####

#### sample age 
#input

data_fig1c <- data.table::fread("processed_data/FileS2_age.tsv") %>% 
  dplyr::left_join(data_fig1a)

# PLOT
fig_1c <- ggplot(data_fig1c,
                 aes(x=fct_reorder(strain, age),
                     y=age,color=geo)) + 
  geom_point(size= 0.3 ) + 
  geom_errorbar(aes(ymin=(age-age_sd), 
                    ymax=(age+age_sd)),
                size= 0.2, width = 0.2) +
  ylab("Estimated ages after\nhatching (hour)") +
  xlab("Strain") +
  theme_cust + 
  theme( axis.text.x = element_blank() ,
         legend.position = "none",
         axis.ticks.x = element_blank() )+
  scale_color_manual(values = rodon_color)

 
 

## cowplot 

fig_1ac <- cowplot::plot_grid(fig_1a,  fig_1c,
                              label_size = 12, 
                              label_fontfamily="Helvetica",
                              nrow = 2)

ggsave(fig_1ac, filename = paste("figures/p1_Fig_1ac.png",sep = ""),  units = "mm",height = 100, width = 100)


#### Figure 1d #####
#### Genetic relatedness   
load("processed_data/FileS3_trees.RData")
 
fig_1d <- ggtree::ggtree(tree_nondiv_genetic, layout="fan", branch.length="rate", size = 0.3,aes(color=label)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "transparent",colour = NA),  
      #  panel.grid.minor = theme_blank(), 
      #  panel.grid.major = theme_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA),
        plot.margin = unit(c(0, 0, 0, 0), "mm") ) +
  scale_color_manual(values=rodon_color)


ggsave(fig_1d, filename = paste("figures/p1_Fig_1d_gtree.png",sep = ""),  units = "mm",height = 200, width = 200)

#### Figure 1e #####
#### Exp relatedness  

 fig_1e <-  ggtree::ggtree(tree_ptc_ph_exp, 
                                    layout="fan", 
                                    branch.length="rate", 
                                    size = 0.3,aes(color=label)) +
  theme(legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "mm") ) +
  scale_color_manual(values=rodon_color)

ggsave(fig_1e, filename = paste("figures/p1_Fig_1e_etree.png",sep = ""),  units = "mm",height = 60, width = 60)

 

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
        legend.position = "bottom",
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        legend.spacing.x = unit(0.001, 'cm') ,
        axis.ticks=element_blank())  +
  ylab("Transcript position (Mb)") + 
  xlab("eQTL position (Mb)")# +
#  labs(color="eQTL classification")
 

gt = ggplot_gtable(ggplot_build(fig_2a))

gt$widths[18] = 0.3*gt$widths[10]

gt$heights[7] = 0.3*gt$heights[9]

 
fig_2a_gt <- ggplotify::as.ggplot(gt) 

#fig_2a_gt


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
  geom_point(data=subset(heri_eQTL,eQTL=="Traits without eQTL"),aes(y=h2,x=H2,color=eQTL), size=0.05,alpha=0.5#,color="black" 
             ) +
  geom_point(data=subset(heri_eQTL,eQTL=="Traits with eQTL"),aes(y=h2,x=H2,color=eQTL), size=0.05,alpha=0.3#,color="orange" 
             ) +
  theme_cust + 
  scale_color_manual(values = c("Traits with eQTL"="orange","Traits without eQTL"="black") )+
  theme(plot.margin = unit(c(2, 3, 0, 0), "mm"),
        legend.position = "bottom",
        legend.margin = margin(0,0,0,0),
        legend.spacing.x = unit(0.001, 'cm') ,
        legend.direction="vertical",
        legend.box.margin=margin(-13,-10,0,-10),
        legend.title = element_blank()) +
  xlab(expression(italic(H^2))) +
  ylab(expression(italic(h^2))) +
  geom_abline(intercept=0,slope=1,colour="black",linetype=2) +
  scale_y_continuous(breaks=c(0,0.5,  1),expand = c(0, 0), limits = c(0,1))  +
  scale_x_continuous(breaks=c(0, 0.5, 1),expand = c(0, 0), limits = c(0,1))
 
#fig_2c


 
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
  ylab("Variance explained (VE)") + 
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
                              labels = c('b', 'e'), 
                              rel_heights = c(2,1),
                              label_size = 12, 
                              label_fontfamily="Helvetica",
                              axis = "lr",
                              # align = "h",
                              nrow = 2)




fig_2cde <- cowplot::plot_grid(fig_2c,fig_2d, fig_2e,
                               labels = c('a','c', 'd'), 
                               rel_heights  = c(1.3,2,1),
                               label_size = 12, 
                               label_fontfamily="Helvetica",
                               axis = "lr",
                               align = "v",
                               nrow = 3)


fig_2 <- cowplot::plot_grid(fig_2cde, fig_2ab, 
                            rel_widths = c(1,2),
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            axis = "tb",
                            nrow =1)

ggsave(fig_2, filename = paste("figures/p1_Fig_2.pdf",sep = ""), units = "mm",height = 170, width = 170)




##########################################
#           Figure 3                    #
#           hotspots                    #
#                                       #   
##########################################
####  fig3 hotspot ####

 
hotspot_pos <- data_fig2a %>% 
  dplyr::filter(!is.na(Hotspot)) %>% 
  dplyr::mutate( Merged_hotspots_cM=ifelse(is.na(Merged_hotspots_cM) ,hotspot_cM,Merged_hotspots_cM),
                merged_Hotspot_QTL_count=ifelse(is.na(merged_Hotspot_QTL_count) , Hotspot_QTL_count,merged_Hotspot_QTL_count)) %>% 
  dplyr::distinct(hotspot_Chr, hotspot_cM, Hotspot,  Merged_hotspots_cM,merged_Hotspot_QTL_count ) %>% 
  dplyr::group_by(hotspot_Chr,Merged_hotspots_cM) %>% 
  dplyr::mutate(merged_Hotspot_center=mean(hotspot_cM))  %>% 
  dplyr::ungroup() %>% 
  dplyr::distinct(hotspot_Chr, Merged_hotspots_cM,merged_Hotspot_center, Hotspot, merged_Hotspot_QTL_count )  %>% 
   dplyr::mutate(hotPOS = ifelse( grepl(",", Merged_hotspots_cM) , "merged",ifelse(Hotspot=="Yes","hotspot","no")))


fig_3 <- ggplot(hotspot_pos,aes(x=merged_Hotspot_center,y=merged_Hotspot_QTL_count)) +
  geom_bar(stat='identity',
           aes(color=hotPOS),size=0.5) +
  facet_grid(.~hotspot_Chr,scales = "free_x" ) +
  scale_color_manual(values = c("red", "purple", "gray69")) +
  geom_hline(yintercept=12,  color="black",alpha=0.8,linetype="dashed") + 
  ylab("Number of distant eQTL") +
  xlab("Hotspot position (cM)")  + 
  theme_cust  +
  theme(legend.position = "none",
      #  axis.text.x = element_blank(),
       # panel.spacing = unit(0.1,"line"),
        axis.text.y =  ggplot2::element_text(size=12,  color = "black",angle = 90 ,vjust = 1,hjust=0.50)) +
  scale_y_continuous(breaks=seq(0, 200, 100), limits = c(0,200),expand = c(0, 0))  +
   scale_x_continuous(breaks =  c(0, 15, 30, 45) ) 


ggsave(fig_3, filename = paste("figures/p1_Fig_3.pdf",sep = ""), units = "mm",height = 60, width = 170)


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
       y = "Mediation\nestimate",
       color = "Mediator gene" ) +
  scale_color_manual(values=c( "lightskyblue2", "orange", "burlywood3", 
                               'cornflowerblue', "deepskyblue4", "mediumpurple4", 
                               "yellow4", "sienna4",'yellowgreen', 'chartreuse','lightpink'),
                       guide = guide_legend(
                       direction = "vertical",
                       ncol=3,
                       title.hjust = 0.5
                     ))+ 
  scale_y_continuous(expand = c(0, 0),limits = c(0,0.165))+
  theme_cust + 
  theme(legend.margin = margin(0,0,0,0),
        legend.text = element_text(size=12,  color = "black",face = "italic" ,  
                                   margin = margin(r = 0.1,l = 0.1, t = 2,b = 2,unit = "pt")),
        legend.spacing.x = unit(0.001, 'cm') ,
        plot.margin = unit(c(1, 1, 1, 1), "mm")) 

######## figure 4b cor ########



data_fig4b <- data.table::fread("processed_data/FileS8_ben1_corr.tsv")

regressed_phe_all<- data_fig4b %>% 
  dplyr::mutate(variant=ifelse(is.na(variantPOS),-1,1),
                regressed_exp = residuals(lm(original_exp ~ variant)))  

regressed_phe_all$ben1_prediction[is.na(regressed_phe_all$ben1_prediction)] <- "A"

 
ymax <- max(regressed_phe_all$original_pheno,regressed_phe_all$regressed_pheno)*1.05

ymin_raw <- min(regressed_phe_all$original_pheno,regressed_phe_all$regressed_pheno)

ymin <- ifelse(ymin_raw>0,ymin_raw*0.95,ymin_raw*1.05)


#ori
corr<-cor.test(regressed_phe_all$original_pheno, regressed_phe_all$original_exp, method = "pearson")

estim<-format(corr$estimate,digits=2, nsmall = 2)

cor_pvalue <- format(corr$p.value ,digits=2, nsmall = 2)

 

#reged  

corr_reg<-cor.test(regressed_phe_all$original_pheno, regressed_phe_all$regressed_exp, method = "pearson")

estim_reg <- format(corr_reg$estimate,digits=2, nsmall = 2)

cor_pvalue_reg <- format(corr_reg$p.value ,digits=2, nsmall = 2)

#

regressed_phe_1 <- regressed_phe_all %>% 
  dplyr::select(strain, pheno= original_pheno, exp = original_exp,ben1_prediction ) %>% 
  dplyr::mutate(type="Raw",
                pearson_cor=ifelse(strain=="AB1",estim,NA),
                pearson_p=ifelse(strain=="AB1",cor_pvalue,NA))
 
 
regressed_phe_2 <- regressed_phe_all %>% 
  dplyr::select(strain, pheno= original_pheno, exp = regressed_exp ,ben1_prediction) %>% 
  dplyr::mutate(type="Regressed",pearson_cor=ifelse(strain=="AB1",estim_reg,NA),
                pearson_p=ifelse(strain=="AB1",cor_pvalue_reg,NA))

regressed_phe_cor <- dplyr::bind_rows(regressed_phe_1,regressed_phe_2)  %>% 
  dplyr::mutate(pp="p",rr="r")
 

fig_4b <- ggplot() +
  geom_point(data=subset(regressed_phe_cor, ben1_prediction=="A"), aes(y=pheno,x=exp), color="gray75",size=0.3)+
  geom_point(data=subset(regressed_phe_cor, ben1_prediction!="A"), aes(y=pheno,x=exp,color= ben1_prediction) ,size=1 )+
  #geom_text(data= subset(regressed_phe_cor, strain=="AB1" & type=="Raw"), aes(label = paste("r",pearson_cor,sep = " : "), y=-136, x=3 ),color="gray6" ,parse = TRUE) +
 # geom_text(data= subset(regressed_phe_cor, strain=="AB1" & type=="Regressed"), aes(label = paste("r",pearson_cor,sep = " : "), y=-136, x=-0.3 ),color="gray6" )+ 
   geom_text(data= subset(regressed_phe_cor, strain=="AB1" & type=="Raw"), aes(label = paste0(": ",pearson_cor), y=-136, x=3 ),color="gray6"  ) +
   geom_text(data= subset(regressed_phe_cor, strain=="AB1" & type=="Regressed"), aes(label = paste0(": ",pearson_cor), y=-136, x=-0.3 ),color="gray6" ) + 
  geom_text(data= subset(regressed_phe_cor, strain=="AB1" & type=="Raw"), aes(label = rr, y=-136, x=2.3,fontface=3 ),color="gray6"  ) +
  geom_text(data= subset(regressed_phe_cor, strain=="AB1" & type=="Regressed"), aes(label = rr, y=-136, x=-1,fontface=3 ),color="gray6"  )+ 
  geom_text(data= subset(regressed_phe_cor, strain=="AB1"& type=="Raw"), aes(label = paste0(": ",pearson_p ), y=-136, x=4.8 ),color="gray6" )+ 
  geom_text(data= subset(regressed_phe_cor, strain=="AB1"& type=="Raw"), aes(label = pp, y=-136, x=4 ,fontface=3 ),color="gray6"  )+ 
  geom_text(data= subset(regressed_phe_cor, strain=="AB1"& type=="Regressed"), aes(label = paste0(": ",pearson_p ), y=-136, x=1.4 ),color="gray6" )+ 
  geom_text(data= subset(regressed_phe_cor, strain=="AB1"& type=="Regressed"), aes(label = pp, y=-136, x=0.8,fontface=3 ),color="gray6"  )+ 
  facet_grid(.~ type, scales = "free" )+
  xlab(expression(paste(italic('ben-1'), " expression"))) +
  ylab("Animal length") +
  labs(color=expression(paste("Variation at ", italic("ben-1"))))+
  ylim(ymin,ymax) +
  theme_cust +
  theme(plot.title = ggplot2::element_text(size=12,  color = "black",hjust = 0,vjust = 1.5),
        legend.position = "right",        
        plot.margin = unit(c(1, 1, 1, 4.5), "mm"))+
  scale_color_manual(values=c("Deletion"="springgreen4", "Frameshift"="red","Missense"="magenta",
                              "Stop gained"="#1e466e","Inversion"="gold2" ),
                     guide = guide_legend(
                       direction = "horizontal",nrow=5,
                       title.position = "top", title.hjust = 0.5
                     ))  
 


#fig_4b
 
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
                 plot.margin = unit(c(1, 1, 1, 7.5), "mm"),
                 axis.text.x = element_blank()) +
  ylab(expression(-log[10](italic(p)))) + 
  xlab("Gene position (Mb)") +
  facet_grid(.~transcript_Chr,scales="free") +
  geom_hline(  yintercept = as.numeric(oriSig) ,
               color = "orange")   




######## cow figure4 ####
 


fig_4 <- cowplot::plot_grid(fig_4a,fig_4b, fig_4c,
                            labels = c('a', 'b', 'c'), 
                             rel_heights = c(1,1.3,1),
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            axis = "lr",
                             nrow = 3)

ggsave(fig_4, filename = paste("figures/p1_Fig_4.pdf",sep = ""), units = "mm",height = 140,      width = 140)

#ggsave(fig_4, filename = paste("figures/p1_Fig_4.png",sep = ""), units = "mm",height = 140,      width = 140)





##########################################
#           Figure 5                    #
#       mediation seven traits          #
##########################################
#input
data_fig5_gwa <- data.table::fread("processed_data/FileS11_7traits_GWA.tsv")  
data_fig5_med <- data.table::fread("processed_data/FileS12_7traits_MED.tsv")  


### fig_5a ####
fig_5a_manh <- manhaplot("telomere_resids") + 
  theme(plot.title =  ggplot2::element_text(size=12,  color = "black")) + 
  ggtitle("Telomere length")

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
  scale_y_continuous(expand = c(0, 0.4),breaks = seq(0,6,2) )+ 
  theme(plot.title =  ggplot2::element_text(size=12,  color = "black")) + 
  ggtitle("Response to arsenic")

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
  theme(axis.title.y = element_blank()) + 
  theme(plot.title =  ggplot2::element_text(size=12,  color = "black")) + 
  ggtitle("Response to zinc")

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
  scale_y_continuous(expand = c(0, 0.4),breaks = seq(0,8,4) ) + 
  theme(plot.title =  ggplot2::element_text(size=12,  color = "black")) + 
  ggtitle("Response to etoposide")

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
  theme(axis.title.y = element_blank())+ 
  theme(plot.title =  ggplot2::element_text(size=12,  color = "black")) + 
  ggtitle("Response to propionate")

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
  scale_y_continuous(expand = c(0, 0.4),breaks = seq(0,15,5) ) + 
  theme(plot.title =  ggplot2::element_text(size=12,  color = "black")) + 
  ggtitle("Response to abamectin")

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
  theme(axis.title.y = element_blank(),plot.title =  ggplot2::element_text(size=12,  color = "black")) + 
  ggtitle("Lifetime fecundity")

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
                                       rel_widths = c(2.3 ,2.6,1.8),
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




ggsave(fig_5, filename = paste("figures/p1_Fig_5.pdf",sep = ""), units = "mm",height = 190, width = 160)




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