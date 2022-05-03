
manhaplot <- function (trait){
  
 
  
  data_manh_raw <- data_fig5_gwa %>% 
    dplyr::filter(organism_level_trait==trait)
  
  
  data_manh_QTLchr <- data_manh_raw %>% dplyr::filter(!is.na(startPOS))
  
  data_manh <- data_manh_raw %>% dplyr::filter(CHROM %in% data_manh_QTLchr$CHROM)
  
  
  figS_manh  <-   ggplot2::ggplot(data_manh) +
    ggplot2::aes(x = POS/1e6, y = log10p) +
    ggplot2::scale_color_manual(values = c("0" = "black", 
                                           "1" = "red",
                                           "2" = "hotpink3",
                                           "3" = "darkorange1")) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = BF),
                        color = "gray", 
                        alpha = .75,  
                        size = 0.5) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = EIGEN_CUTOFF),
                        color = "gray", 
                        alpha = .75,  
                        size = 0.5,
                        linetype = 2) +
    ggplot2::geom_point( ggplot2::aes(color= factor(EIGEN_SIG)),size=0.5 ) +
    ggplot2::facet_grid( . ~ CHROM, scales = "free"  ) +
    theme_cust + 
    theme( 
      axis.text.x = element_blank(),
      panel.spacing = unit(0.01,"line"),
    #  plot.margin = unit(c(0, 0, 2, 0), "mm"),
      axis.title.x = element_blank())+   
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = "Genomic position (Mb)",
                  y = expression(-log[10](italic(p)))) + 
    scale_y_continuous(expand = c(0, 0.8) ) 
  
  
  return(figS_manh)
  
  
}

medaplot <- function (trait){
  
  data_med <- data_fig5_med %>% 
    dplyr::filter(organism_level_trait==trait)
  
  
  q99_med = quantile(data_med$mediation_estimate, probs = 0.99)[[1]]
   
  data_manh_QTLchr <- data_fig5_gwa %>% 
    dplyr::filter(organism_level_trait==trait) %>% 
    dplyr::filter(!is.na(startPOS))
  
  
  newrows_size <- data.frame(organism_level_trait=NA,
                             QTL_Chr=c("I","II","III","IV","V","X"),
                             eQTL=0,
                             startPOS=0,
                             endPOS = c(14972282, 15173999, 13829314, 17450860, 20914693, 17748731),
                             mediation_estimate=0) %>% 
    dplyr::filter(QTL_Chr %in% data_manh_QTLchr$CHROM)
  
  
  
  
  #fig_med
  
  fig_med <- ggplot() +
    geom_point(data=newrows_size,aes(x=as.numeric(eQTL),y=mediation_estimate),  size=0.1, alpha = 0)+
    geom_point(data=newrows_size,aes(x=as.numeric(endPOS)/1e6,y=mediation_estimate),  size=0.1, alpha = 0) +
    geom_hline(data=data_med, aes(yintercept = q99_med), color = "grey") +
    geom_point(data=subset(data_med,q99_mediator == "Other genes" ),  
               aes(x=eQTL/1e6,
                   y=mediation_estimate),
               color = "gray80",
               size=0.5) +
    geom_point(data=subset(data_med,q99_mediator != "Other genes" ),  
               aes(x=eQTL/1e6,
                   y=mediation_estimate,
                   color =  q99_mediator ),
               size=0.5)  +
    scale_alpha_continuous(range = c(1, 0)) +
    labs(x = "Genomic position (Mb)", 
         y = "Mediation\nestimate",
         color = "Mediator gene") +
    theme_cust + 
    theme(legend.text = ggtext::element_markdown(size=12,  color = "black"),
          strip.text = element_blank(),
          plot.margin = unit(c(0, 0, 2, 0), "mm"),
          panel.spacing = unit(0.01,"line"),
          axis.text.x = element_blank(),
          axis.title.x  = element_blank(),
          legend.position = "none",
          legend.box = "horizontal") +
    scale_color_manual(values=ancestry.colours) +
    guides(color=guide_legend(nrow=2,byrow=TRUE)) +
    facet_grid(.~QTL_Chr,scales = "free" ) +
    scale_x_continuous(breaks=seq(0, 21,5) ) +
    ggrepel::geom_label_repel(data=subset(data_med,q99_mediator != "Other genes" ), aes(label= q99_mediator ,
                                                                                        x=eQTL/1e6,
                                                                                        y=mediation_estimate,
                                                                                        color =  q99_mediator ),  
                              label.size = NA, 
                              size = 3,
                              segment.color = NA,
                              min.segment.length = unit(0, 'lines'), 
                              fill = NA,
                              fontface = "italic" ) 
  
  
  return(fig_med)
  
}

medaplot_aba <- function (trait){
  
  
  
  data_med <- data_fig5_med %>% 
    dplyr::filter(organism_level_trait==trait)
  
  
  q99_med = quantile(data_med$mediation_estimate, probs = 0.99)[[1]]
  
  
  data_manh_QTLchr <- data_fig5_gwa %>% 
    dplyr::filter(organism_level_trait==trait) %>% 
    dplyr::filter(!is.na(startPOS))
  
  
  
  newrows_size <- data.frame(organism_level_trait=NA,
                             QTL_Chr=c("I","II","III","IV","V","X"),
                             eQTL=0,
                             startPOS=0,
                             endPOS = c(14972282, 15173999, 13829314, 17450860, 20914693, 17748731),
                             mediation_estimate=0) %>% 
    dplyr::filter(QTL_Chr %in% data_manh_QTLchr$CHROM)
  
  
  
  
  #fig_med
  
  fig_med <- ggplot() +
    geom_point(data=newrows_size,aes(x=as.numeric(eQTL),y=mediation_estimate),  size=0.1, alpha = 0)+
    geom_point(data=newrows_size,aes(x=as.numeric(endPOS)/1e6,y=mediation_estimate),  size=0.1, alpha = 0) +
    geom_hline(data=data_med, aes(yintercept = q99_med), color = "grey") +
    geom_point(data=subset(data_med,q99_mediator == "Other genes" ),  
               aes(x=eQTL/1e6,
                   y=mediation_estimate),
               color = "gray80",
               size=0.5) +
    geom_point(data=subset(data_med,q99_mediator != "Other genes" ),  
               aes(x=eQTL/1e6,
                   y=mediation_estimate,
                   color =  q99_mediator ),
               size=0.5)  +
    scale_alpha_continuous(range = c(1, 0))   +
    labs(x = "Genomic position (Mb)", 
         y = "Mediation\nestimate",
         color = "Mediator gene") +
    theme_cust + 
    theme( strip.text = element_blank(),
           plot.margin = unit(c(0, 0, 2, 0), "mm"),
           panel.spacing = unit(0.01,"line"),
           axis.text.x = element_blank(),
           legend.title = element_blank(),
           legend.position = "bottom",
           legend.margin = margin(0,0,0,-20),
           legend.text = element_text(size=10,  color = "black",face = "italic" ,
                                      margin = margin(r = 0.1,l = 0.1, t = 0.1,b = 0.1,unit = "pt")),
           legend.spacing = unit(0.001, "cm"), 
           legend.spacing.x = unit(0.001, 'cm'),
           legend.key = element_rect(size = 3, fill = "white", colour = NA), 
           legend.key.size = unit(0.3, "cm"),
           legend.box = "horizontal") +
    scale_color_manual(values=ancestry.colours) +
    guides(color=guide_legend(nrow=4,byrow=TRUE,ncol = 3)) +
    facet_grid(.~QTL_Chr,scales = "free" ) +
    scale_x_continuous(breaks=seq(0, 21,5) )  +
    scale_y_continuous(breaks=seq(0, 0.3,0.1) )
  
  
  
  return(fig_med)
  
}

medaplot_bro <- function (trait){
  
  
  
  data_med <- data_fig5_med %>% 
    dplyr::filter(organism_level_trait==trait)
  
  
  q99_med = quantile(data_med$mediation_estimate, probs = 0.99)[[1]]
  
  
  data_manh_QTLchr <- data_fig5_gwa %>% 
    dplyr::filter(organism_level_trait==trait) %>% 
    dplyr::filter(!is.na(startPOS))
  
  newrows_size <- data.frame(organism_level_trait=NA,
                             QTL_Chr=c("I","II","III","IV","V","X"),
                             eQTL=0,
                             startPOS=0,
                             endPOS = c(14972282, 15173999, 13829314, 17450860, 20914693, 17748731),
                             mediation_estimate=0) %>% 
    dplyr::filter(QTL_Chr %in% data_manh_QTLchr$CHROM)
  
  
  
  
  #fig_med
  
  fig_med <- ggplot() +
    geom_point(data=newrows_size,aes(x=as.numeric(eQTL),y=mediation_estimate),  size=0.1, alpha = 0)+
    geom_point(data=newrows_size,aes(x=as.numeric(endPOS)/1e6,y=mediation_estimate),  size=0.1, alpha = 0) +
    geom_hline(data=data_med, aes(yintercept = q99_med), color = "grey") +
    geom_point(data=subset(data_med,q99_mediator == "Other genes" ),  
               aes(x=eQTL/1e6,
                   y=mediation_estimate),
               color = "gray80",
               size=0.5) +
    geom_point(data=subset(data_med,q99_mediator != "Other genes" ),  
               aes(x=eQTL/1e6,
                   y=mediation_estimate,
                   color =  q99_mediator ),
               size=0.5)  +
    scale_alpha_continuous(range = c(1, 0))  +
    labs(x = "Genomic position (Mb)", 
         y = "Mediation\nestimate",
         color = "Mediator gene") +
    theme_cust + 
    theme( strip.text = element_blank(),
           axis.title.y =  element_blank(),
           plot.margin = unit(c(0, 0, 2, 0), "mm"),
           panel.spacing = unit(0.01,"line"),
           axis.text.x = element_blank(),
           legend.title = element_blank(),
           legend.position = "bottom",
           legend.margin = margin(0,0,0,-27),
           legend.text = element_text(size=10,  color = "black",face = "italic" ,
                                      margin = margin(r = 0.1,l = 0.1, t = 0.1,b = 0.1,unit = "pt")),
           legend.spacing = unit(0.001, "cm"), 
           legend.spacing.x = unit(0.001, 'cm'),
           legend.key = element_rect(size = 3, fill = "white", colour = NA), 
           legend.key.size = unit(0.3, "cm"),
           legend.box = "horizontal") +
    scale_color_manual(values=ancestry.colours) +
    guides(color=guide_legend(nrow=4,byrow=TRUE,ncol = 5)) +
    facet_grid(.~QTL_Chr,scales = "free" ) +
    scale_x_continuous(breaks=seq(0, 21,5) )  
  
  
  
  
  return(fig_med)
  
}


