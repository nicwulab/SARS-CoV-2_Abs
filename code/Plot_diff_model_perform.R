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

plot_model_perform <- function(data,p_width,graphname){
  palette  <- c(brewer.pal(3,"Set2"))
  textsize <- 7
  p <-  ggplot() +
          geom_bar(data=data, aes(x=Metrics, y=value, fill=Model), stat='identity', width=0.6, position=position_dodge()) +
          scale_fill_manual(values=palette,drop=FALSE,labels = c("RBD/HA","S/HA","RBD/NTD/S2")) +
          theme_cowplot(12) + ylim(0,1)+
          theme(axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title.x=element_blank(),
                axis.title.y=element_text(size=textsize,face="bold"),
                axis.line = element_line(colour = 'black', size = 0),
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                legend.position = "right",
                legend.title    = element_text(size=textsize,face="bold"),
                legend.text=element_text(size=textsize,face="bold"),
                legend.justification='center',
                legend.key.size = unit(0.5,"line")) +
          ylab("Performance")
  ggsave(graphname,p,width=p_width, height=1.5, bg='white', dpi=300)
}

model_level <- c('RBD-HA','Spike-HA','RBD-NTD-S2')
metrics_level <- c('Accuracy','Precision','Recall','ROC AUC')
model_data <- read_excel('result/model_comparison.xlsx') %>%
                mutate(Model=factor(Model, level=model_level)) %>%
                mutate(Metrics=factor(Metrics, level=metrics_level))
plot_model_perform(model_data, 3,'graph/model_comparison.png')
