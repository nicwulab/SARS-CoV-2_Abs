#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(data.table)
library(qualpalr)
require(cowplot)

plot_Vgenes <- function(Vgene_table, graphname){
  textsize <- 6
  #palette  <- c(brewer.pal(12,"Paired"))
  palette <- qualpal(n = 21, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))
  palette <- palette$hex
  p <-  ggplot() +
          geom_bar(data=Vgene_table, aes(x=Vgene, y=count, fill=Patient), stat='identity', width=0.6) +
          scale_fill_manual(values=palette,drop=FALSE) +
          theme_cowplot(12) +
          theme(plot.title=element_text(size=textsize,face="bold",colour = 'black',hjust=0.5),
                axis.text=element_text(size=textsize-1,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title.x=element_blank(),
                axis.title.y=element_text(size=textsize,face="bold"),
                axis.line = element_line(colour = 'black', size = 0),
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                legend.position = "none",
                legend.title    = element_text(size=textsize,face="bold"),
                legend.text=element_text(size=textsize,face="bold"),
                legend.justification='center',
                legend.key.size = unit(0.3,"line")) +
          ggtitle('Cluster 7') +
          ylab("# of antibodies")
  ggsave(graphname,p,width=1.2, height=1.4, bg='white', dpi=1200)
  }

Ab_table <- read_tsv('result/Ab_info_CDRH3_clustering.tsv')
Ab_table <- Ab_table %>%
              filter(cluster==7) %>%
              mutate(`Heavy V Gene`=mapply(function(x){return (str_split(x, '\\*')[[1]][1])}, `Heavy V Gene`)) %>%
              mutate(`Light V Gene`=mapply(function(x){return (str_split(x, '\\*')[[1]][1])}, `Light V Gene`))
Vgene_table <- Ab_table %>%
                 select(`Patient ID`, `Heavy V Gene`) %>%
                 setnames(c('Patient','Vgene')) %>%
                 mutate(Vgene=mapply(function(x){if (x=='nan'){return ('unknown')}else{return (x)}}, Vgene)) %>%
                 mutate(Vgene=mapply(function(x){if (x=='IGHV3-30-3'){return ('IGHV3-30')}else{return (x)}}, Vgene)) %>%
                 mutate(num=1) %>%
                 group_by(Patient, Vgene) %>%
                 summarise(count=n())
print (Vgene_table)
plot_Vgenes(Vgene_table, 'graph/Vgenes_cluster7.png')
