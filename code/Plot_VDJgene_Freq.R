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
library(ggbeeswarm)
library(sinaplot)
library(readxl)
require(cowplot)

plot_Vgene_usage <- function(data,ref,p_width,graphname){
  palette  <- c(brewer.pal(3,"Accent"))
  textsize <- 7
  p <-  ggplot() +
          geom_bar(data=data, aes(x=gene, y=freq*100, fill=epitope), stat='identity', width=0.6, position=position_dodge()) +
          scale_fill_manual(values=palette,drop=FALSE) +
          geom_errorbar(data = ref,aes(x=V_gene,ymin=(range_low)*100, ymax=(range_high)*100),
                        size=0.25, width=0.5,position=position_dodge(), color='gray30') +
          theme_cowplot(12) +
          theme(axis.text=element_text(size=textsize-1,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5,colour = 'black'),
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
          ylab("frequency (%)")
  ggsave(graphname,p,width=p_width, height=1.3, bg='white', dpi=1200)
}
plot_Dgene_usage <- function(data,ref,graphname){
  palette  <- c(brewer.pal(3,"Accent"))
  textsize <- 7
  p <-  ggplot() +
          geom_bar(data=data, aes(x=gene, y=freq*100, fill=epitope), stat='identity', width=0.6, position=position_dodge()) +
          scale_fill_manual(values=palette,drop=FALSE) +
          geom_errorbar(data = ref,aes(x=D_gene,ymin=(range_low)*100, ymax=(range_high)*100),
                        size=0.25, width=0.5,position=position_dodge(), color='gray30') +
          theme_cowplot(12) +
          theme(axis.text=element_text(size=textsize,face="bold",colour = 'black'),
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
          ylab("frequency (%)")
  ggsave(graphname,p,width=3, height=1.3, bg='white', dpi=1200)
}

set.seed(5)
HVgene_ref  <- read_excel('data/HV_Repertoire_freq.xlsx')
HVgene_level <- sort(HVgene_ref$V_gene)
HVgene_data <- read_tsv('result/HVgene_freq.tsv') %>%
                filter(gene %in% HVgene_level) %>%
                mutate(gene=factor(gene, level=HVgene_level))
plot_Vgene_usage(HVgene_data, HVgene_ref, 6,'graph/HV_gene_usage.png')

LVgene_ref  <- read_excel('data/LV_Repertoire_freq.xlsx')
LVgene_level <- sort(LVgene_ref$V_gene)
LVgene_data <- read_tsv('result/LVgene_freq.tsv') %>%
                filter(gene %in% LVgene_level) %>%
                mutate(gene=factor(gene, level=LVgene_level))
plot_Vgene_usage(LVgene_data, LVgene_ref, 6,'graph/LV_gene_usage.png')

Dgene_ref  <- read_excel('data/D_Repertoire_freq.xlsx')
Dgene_level <- sort(Dgene_ref$D_gene)
Dgene_data <- read_tsv('result/Dgene_freq.tsv') %>%
                filter(gene %in% Dgene_level) %>%
                mutate(gene=factor(gene, level=Dgene_level))
plot_Dgene_usage(Dgene_data, Dgene_ref, 'graph/D_gene_usage.png')
