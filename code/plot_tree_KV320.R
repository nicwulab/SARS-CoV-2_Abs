#R code
library(ape)
library(dplyr)
library(readr)
library(scales)
library(ggtree)
library(ggplot2)
library(RColorBrewer)
require(cowplot)

tree  <- read.tree('result/Cluster3_H158K320_LC.tree')
tree  <- root(tree, outgroup="IGKV3-20_germline", resolve.root = TRUE)
class <- read_tsv('result/Cluster3_H158K320_LC_class.tsv')
textsize <- 6
palette  <- c(brewer.pal(5,"Set1"))
palette  <- c(palette[5], palette[1], palette[4], 'gray40')
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
