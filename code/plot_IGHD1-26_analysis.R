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
library(readxl)
require(cowplot)

plot_Vgenes <- function(Vgene_table, graphname, title){
  Vgene_table <- Vgene_table %>%
    mutate(Vgene=mapply(function(x){if (x=='nan'){return ('unknown')}else{return (x)}}, Vgene)) %>%
    mutate(num=1) %>%
    group_by(Patient, Vgene) %>%
    summarise(count=n())
  textsize <- 7
  #palette  <- c(brewer.pal(12,"Paired"))
  palette <- qualpal(n = 21, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex
  p <-  ggplot() +
          geom_bar(data=Vgene_table, aes(x=Vgene, y=count), stat='identity', width=0.6) +
          #scale_fill_manual(values=palette,drop=FALSE) +
          theme_cowplot(12) +
          theme(plot.title=element_text(size=textsize,face="bold",colour = 'black',hjust=0.5),
                axis.text=element_text(size=textsize,face="bold",colour = 'black'),
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
          ggtitle(title) +
          ylab("# of antibodies")
  ggsave(graphname,p,width=3, height=1.5, bg='white', dpi=1200)
  }

plot_Jgenes <- function(Ab_table, graphname){
  textsize <- 7
  IGHJ_table <- Ab_table %>%
    select(`Patient ID`, `Heavy J Gene`) %>%
    setnames(c('Patient','Jgene')) %>%
    mutate(Jgene=mapply(function(x){if (x=='nan'){return ('unknown')}else{return (x)}}, Jgene)) %>%
    mutate(num=1) %>%
    group_by(Jgene) %>%
    summarise(count=n()) %>%
    mutate(Jgene=factor(Jgene, levels=c('IGHJ1','IGHJ2','IGHJ3','IGHJ4','IGHJ5','IGHJ6')))
  print (sum(IGHJ_table$count))
  palette <- qualpal(n = 7, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex
  p <- ggplot(IGHJ_table,aes(x='1',y=count,group=Jgene,fill=Jgene)) +
     #geom_bar(stat="identity", position=position_dodge(), width=1) +
     geom_bar(stat="identity", width=1) +
     scale_fill_manual(values=palette,drop=FALSE) +
     theme_cowplot(12) +
     theme(axis.title=element_text(size=textsize,face="bold"),
           axis.text=element_blank(),
           legend.title=element_blank(),
           legend.text=element_text(size=textsize-1,face="bold", margin = margin(l = -5, unit = "pt")),
           legend.key.size = unit(0.3,"line"),
           legend.position='top') +
     labs(y=expression(bold('count')),x=expression()) +
     coord_polar("y", start=0)
  ggsave(graphname, p, height=1.5, width=2.5,bg='white')
  }

plot_CDRH3_len <- function(Ab_table, graphname){
  print (table(Ab_table$Ab_class))
  Ab_table <- Ab_table %>%
                mutate(CDR_len=mapply(str_length,  CDRH3)) %>%
                group_by(Ab_class, CDR_len) %>%
                summarise(count=n()) %>%
                group_by(Ab_class) %>%
                mutate(freq=100*count/sum(count))
  textsize <- 7
  #palette  <- c(brewer.pal(12,"Dark2"))
  palette <- qualpal(n = 3, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex
  print (filter(Ab_table, CDR_len==14))
  p <-  ggplot() +
          geom_bar(data=Ab_table, aes(x=CDR_len, y=freq, fill=Ab_class), stat='identity', width=0.6,position=position_dodge()) +
          scale_fill_manual(values=palette,drop=FALSE) +
          theme_cowplot(12) +
          theme(plot.title=element_text(size=textsize,face="bold",colour = 'black',hjust=0.5),
                axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title.x=element_text(size=textsize,face="bold"),
                axis.title.y=element_text(size=textsize,face="bold"),
                axis.line = element_line(colour = 'black', size = 0),
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                legend.position = "top",
                legend.title    = element_blank(),
                legend.text=element_text(size=textsize-1,face="bold"),
                legend.justification='right',
                legend.key.size = unit(0.5,"line")) +
          ylab("Frequency (%)") +
          xlab("CDR H3 length (IMGT numbering)")
  ggsave(graphname,p,width=3, height=1.2, bg='white', dpi=1200)
  }

classify_Ab <- function(D, epi){
  if (D=='IGHD1-26' & epi=='S:S2'){
    return ('IGHD1-26 (S2)')
    } 
  else if (D=='IGHD1-26' & epi=='S:S2 Stem Helix'){
    return ('IGHD1-26 (S2)')
    }
  else if (epi=='S:S2' || epi =='S:S2 Stem Helix'){
    return ('non-IGHD1-26 (S2)')
    }
  else if (epi=='S:NTD' || epi=='S:RBD' || epi=='S:S1' || epi=='S:S1 non-RBD'){
    return ('non-S2')
    }
  else{
    return ('discard')
    }
  }

Ab_table <- read_excel('data/SARS-CoV-2-Abs.xlsx')
Ab_table <- Ab_table %>%
              filter(`Heavy D Gene` != '') %>%
              filter(`Protein + Epitope` != '') %>%
              filter(CDRH3 != '') %>%
              mutate(`Heavy V Gene`=mapply(function(x){return (str_split(x, '\\*')[[1]][1])}, `Heavy V Gene`)) %>%
              mutate(`Light V Gene`=mapply(function(x){return (str_split(x, '\\*')[[1]][1])}, `Light V Gene`)) %>%
              mutate(`Heavy J Gene`=mapply(function(x){return (str_split(x, '\\*')[[1]][1])}, `Heavy J Gene`)) %>%
              mutate(`Light J Gene`=mapply(function(x){return (str_split(x, '\\*')[[1]][1])}, `Light J Gene`)) %>%
              mutate(`Heavy D Gene`=mapply(function(x){return (str_split(x, '\\*')[[1]][1])}, `Heavy D Gene`)) %>%
              mutate(Ab_class=mapply(classify_Ab, `Heavy D Gene`, `Protein + Epitope`)) %>%
              filter(Ab_class!='discard')
table_S2_IGHD1_26 <- Ab_table %>%
   filter(Ab_class=='IGHD1-26 (S2)')
table_other <- Ab_table %>%
   filter(Ab_class!='IGHD1-26 (S2)')
HVgene_table_S2_IGHD1_26 <- table_S2_IGHD1_26 %>%
   select(`Patient ID`, `Heavy V Gene`) %>%
   setnames(c('Patient','Vgene')) %>%
   mutate(Vgene=mapply(function(x){if (x=='IGHV3-30-3'){return ('IGHV3-30')}else{return (x)}}, Vgene))
LVgene_table_S2_IGHD1_26 <- table_S2_IGHD1_26 %>%
   select(`Patient ID`, `Light V Gene`) %>%
   setnames(c('Patient','Vgene'))
HJgene_table_S2_IGHD1_26 <- table_S2_IGHD1_26 %>%
   select(`Patient ID`, `Heavy J Gene`) %>%
   setnames(c('Patient','Jgene'))
plot_Vgenes(HVgene_table_S2_IGHD1_26, 'graph/S2_IGHD1-26_HVgenes.png', 'IGHV usage in\nIGHD1-26-encoded S2 antibodies')
plot_Vgenes(LVgene_table_S2_IGHD1_26, 'graph/S2_IGHD1-26_LVgenes.png', 'IGK(L)V usage in\nIGHD1-26-encoded S2 antibodies')
plot_Jgenes(table_S2_IGHD1_26, 'graph/HJgenes_pie_S2_IGHD1-26.png')
plot_Jgenes(table_other, 'graph/HJgenes_pie_other.png')
plot_CDRH3_len(Ab_table, 'graph/S2_IGHD1-26_CDRH3_len.png')

print (length(table_S2_IGHD1_26$CDRH3))
table_S2_IGHD1_26 <- table_S2_IGHD1_26 %>%
  filter(str_length(CDRH3)==14) 
print (paste("# of antibodies:", length(table_S2_IGHD1_26$`Patient ID`), sep=' '))
print (paste("# of donors:", length(unique(table_S2_IGHD1_26$`Patient ID`)), sep=' '))
