#R code
library(ape)
library(dplyr)
library(readr)
library(scales)
library(ggtree)
library(ggplot2)
library(qualpalr)
library(RColorBrewer)
require(cowplot)

tree  <- read.tree('result/Cluster3_H158K320_LC.tree')
tree  <- root(tree, outgroup="IGKV3-20_germline", resolve.root = TRUE)
class <- read_tsv('result/Cluster3_H158K320_LC_class.tsv') %>%
           mutate(mut=factor(mut, levels=c("Germline",'S29/G92','R29/G92','R29/D92','S29/N92','S29/V92')))
textsize <- 6
palette <- qualpal(n = 6, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex
palette <- c('gray40', palette[5], palette[2], palette[3], palette[4], palette[6])
p  <- ggtree(tree, lwd=0.2)
p  <- p %<+% 
       class + 
       geom_tippoint(aes(color=mut), size=1.2, alpha=0.7, pch=16)
p  <- p +
        scale_color_manual(values=palette,drop=FALSE) +
        ylim(NA, length(class$ID)) +
        theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),
              legend.title    = element_blank(),
              legend.position = "right",
	      legend.text=element_text(size=textsize,face="bold"),
	      legend.justification='top',
              legend.key.size = unit(0.6,"line")) 
ggsave('graph/Cluster3_H158K320_LC_tree.png',p, width=4, height=2, bg='white', dpi=600)
