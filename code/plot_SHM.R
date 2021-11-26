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

y_nudging <- function(label){
  #if (label == 'Q3R (IGHV3-30/IGKV1-33)'){return (15)}
  if (label == 'F27L (IGHV3-53/IGKV1-9/c4)'){return (-5)}
  if (label == 'F27L (IGHV3-66/IGKV1-9/c1)'){return (-5)}
  if (label == 'Q3R (IGHV3-30/IGKV1-33/c5)'){return (-5)}
  if (label == 'D1A (IGHV3-30/IGKV1-33/c5)'){return (-10)}
  if (label == 'D1A (IGHV3-66/IGKV1-9/c1)'){return (5)}
  if (label == 'M3V (IGHV4-39/IGLV6-57/c7)'){return (-2)}
  if (label == 'F10S (IGHV3-53/IGKV1-9/c1)'){return (-2)}
  if (label == 'S29R (IGHV1-58/IGKV3-20/c3)'){return (0)}
  if (label==F){return (2)}
  else{return (2)}
  }

x_nudging <- function(label){
  if (label == 'F27L (IGHV3-66/IGKV1-9/c1)'){return (-5)}
  if (label == 'F27V (IGHV3-66/IGKV1-9/c1)'){return (-5)}
  if (label == 'F27L (IGHV3-53/IGKV1-9/c4)'){return (-1)}
  if (label == 'M3V (IGHV4-39/IGLV6-57/c7)'){return (15)}
  if (label == 'D1A (IGHV3-30/IGKV1-33/c5)'){return (-2)}
  if (label == 'Q3R (IGHV3-30/IGKV1-33/c5)'){return (15)}
  if (label == 'D1A (IGHV3-30/IGKV1-33/c5)'){return (-2)}
  if (label == 'F10S (IGHV3-53/IGKV1-9/c1)'){return (10)}
  if (label == 'S29R (IGHV1-58/IGKV3-20/c3)'){return (-20)}
  if (label == 'Y58F (IGHV3-53/IGKV3-20/c2)'){return (10)}
  if (label==F){return (2)}
  else{return (2)}
  }

plot_SHM <- function(SHM_table, donor_cutoff, chain_id, graphname, title, x_breaks, x_labels){
  SHM_table   <- filter(SHM_table, chain==chain_id)
  if (chain_id == 'HC'){cutoff_freq <- 0.4}
  if (chain_id == 'LC'){cutoff_freq <- 0.2}
  if (chain_id == 'LC'){
    ##These two mutations are mis-called (IgBLast does not match Kabat numbering) due to a deletion##
    SHM_table <- SHM_table %>%
                   filter(label != 'V28S (IGHV3-53/IGKV3-20)') %>%
                   filter(label != 'S29V (IGHV3-53/IGKV3-20)')
    }
  label_table <- SHM_table %>%
                   filter(total_donor >= donor_cutoff) %>%
                   filter(freq>=cutoff_freq) %>%
                   mutate(nudge_y=mapply(y_nudging, label)) %>%
                   mutate(nudge_x=mapply(x_nudging, label))
  textsize <- 7
  p <-  ggplot() +
          geom_point(data=SHM_table, aes(x=adj_pos, y=freq*100, size=total_donor, color=color),
                     alpha=0.3, pch=16) +
          geom_text_repel(data=label_table, aes(x=adj_pos, y=freq*100,label=label),
                          color="black", min.segment.length=0, segment.size=0.2, size=2, force=10, force_pull=2,
                          seed=3, nudge_x=label_table$nudge_x, nudge_y=label_table$nudge_y, max.overlaps = Inf) +
          scale_color_manual(values = c("black","white")) +
          scale_size_continuous(name = "# of donors in\nthe clonotype",
                                #breaks = c(5,10,15,20,25),
                                breaks = c(10,15,20,25,30),
                                range = c(0.5, 3.5)) +
          theme_cowplot(12) +
          theme(plot.title=element_text(size=textsize+2,hjust=0.5,vjust=0.5,face="bold",colour = 'black'),
                axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title.x=element_text(size=textsize,face="bold"),,
                axis.title.y=element_text(size=textsize,face="bold"),
                axis.line = element_line(colour = 'black', size = 0),
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                legend.position = "right",
                legend.title    = element_text(size=textsize,face="bold"),
                legend.text=element_text(size=textsize,face="bold"),
                legend.justification='center',
                legend.key.size = unit(0.3,"line")) +
          ggtitle(title) +
          ylab("% within clonotype") +
          xlab("position (Kabat numbering)") +
          ylim(0, max(label_table$freq)*100+2) +
          scale_x_continuous(breaks=x_breaks, labels=x_labels)
  ggsave(graphname,p,width=7, height=2.5, bg='white', dpi=1200)
  }

coloring <- function(total_donor, donor_cutoff){
  if (total_donor < donor_cutoff){
    return ('white')
    }
  else{
    return ('black')
    }
  }

set.seed(4)
donor_cutoff <- 9
SHM_table <- read_tsv('result/SHM_frequency.tsv') %>%
                   filter(total_donor >= 7) %>%
                   filter(occurrence >= 2) %>%
                   mutate(color=mapply(coloring, total_donor, rep(donor_cutoff,length(total_donor))))
plot_SHM(SHM_table, donor_cutoff, 'HC', 'graph/SHM_HC_frequency.png', 'Recurring somatic hypermutations in IGHV-encoded region',
         c(10, 21, 31, 48, 60, 74, 86, 96, 109), c(10,20,30,40,50,60,70,80,90))
plot_SHM(SHM_table, donor_cutoff, 'LC', 'graph/SHM_LC_frequency.png', 'Recurring somatic hypermutations in IGK(L)V-encoded region',
         c(9, 21, 37, 49, 59, 76, 88, 100), c(10,20,30,40,50,60,70,80))
