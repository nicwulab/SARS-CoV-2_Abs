#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(gridExtra)
library(ggforce)
library(ggrepel)
library(ggbeeswarm)
library(sinaplot)
require(cowplot)

plot_CDRH3_cluster_size <- function(CDRH3_cluster, graphname){
  textsize <- 7
  ID_of_interest <- c()
  label_table <- CDRH3_cluster %>%
                   filter(cluster_ID %in% ID_of_interest)
  seed_num <- 4
  p <-  ggplot() +
          geom_point(data=CDRH3_cluster, aes(x=epitope, y=cluster_size, size=num_donors),
                     alpha=0.3, pch=16, fill='black', color='black', position=position_jitter(h=0,w=0.3,seed=seed_num)) +
          geom_text_repel(data=label_table, aes(x=epitope, y=cluster_size,label=cluster_ID),
                          color="black", min.segment.length=0, segment.size=0.2, size=2, force=3, force_pull=2,
                          position=position_jitter(h=0,w=0.3,seed=seed_num), max.overlaps = Inf) +
          scale_size_continuous(name = "# of donors",
                                breaks = c(10,20,30,40,50),
                                range = c(0, 3.5)) +
          theme_cowplot(12) +
          theme(axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title.x=element_blank(),
                axis.title.y=element_text(size=textsize,face="bold"),
                axis.line = element_line(colour = 'black', size = 0),
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                legend.position = "right",
                legend.title    = element_text(size=textsize,face="bold"),
                legend.text=element_text(size=textsize,face="bold"),
                legend.justification='center',
                legend.key.size = unit(0.3,"line")) +
          ylab("cluster size") +
          ylim(0,max(CDRH3_cluster$cluster_size)+6) +
          scale_x_discrete(labels=c("S:RBD" = "RBD",
                                    "S:S2" = "S2",
                                    "S:NTD" = "NTD"))
  ggsave(graphname,p,width=2.5, height=2, bg='white', dpi=1200)
  }

CDRH3_cluster <- read_tsv('result/CDRH3_cluster_summary.tsv') %>%
                   filter(cluster_size >= 2) %>%
                   filter(epitope != "unknown") %>%
                   mutate(epitope=factor(epitope,levels=c('S:RBD','S:NTD','S:S2')))
plot_CDRH3_cluster_size(CDRH3_cluster, 'graph/CDRH3_cluster_size.png')
