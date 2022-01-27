# Title     : to plot distribution of clonotype of datasets
# Objective :
# Created by: yiquan
# Created on: 1/18/22
library(ggplot2)
library(readxl)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(data.table)
library(GGally)
library(ggpubr)
library(ggExtra)
require(cowplot)
#define plot function

plot_hist <- function(df,df2,column,x_lab,y1,y2,main_title){
  palette  <- c(brewer.pal(3,"Pastel1"))
  p1 <- ggplot(df) +
         geom_histogram(aes(x=column,fill=state),binwidth=1) +
         scale_fill_manual(values=palette,drop=FALSE) +
         theme_cowplot(12) +
         theme(plot.title=element_text(size=7,face="bold", hjust=0.5),
               axis.title=element_text(size=7,face="bold"),
               axis.text=element_text(size=7,face="bold"),
               legend.key.size=unit(0.1,'in'),
               legend.spacing.x=unit(0.03, 'in'),
               legend.title=element_blank(),
               legend.text=element_text(size=7,face="bold"),
               legend.position='top') +scale_y_continuous(labels = label_number(accuracy = 1))+
     ggtitle(main_title) +
     xlab(x_lab) + ylab("Counts") # coord_cartesian(xlim=c(-2.5,1))

  p2 <- ggplot(df2)+
    geom_bar(aes(x=dataset,y=percent,fill=state),position="stack", stat="identity",width = 0.7)+coord_flip()+
    scale_fill_manual(values=palette,drop=FALSE) +
    theme_cowplot(12) +
         theme(plot.title=element_blank(),
               axis.title=element_blank(),
               axis.text=element_text(size=7,face="bold"),
               # axis.line = element_blank(),
               axis.ticks.y=element_blank(),
               legend.title=element_blank(),
               legend.position = 'none')+scale_y_continuous(labels = scales::percent_format(accuracy = 1))
  p <- p1+annotation_custom(ggplotGrob(p2), xmin = 20, xmax = 5000,
                       ymin = y1, ymax = y2)
  return(p)
  }


#-------------------------------------read dataframe--------------------------------------#
Ab_DF <- read_excel('result/Data S2_v2.xlsx')%>% filter(clonotype != 'NA')
train_df <- Ab_DF %>% filter(dataset=='train set')
train_clono <- train_df$clonotype
train_df$state <- 'overalp with train set'

train_df2 <- train_df %>% group_by(dataset,state)%>%
  dplyr::summarise(count = n()) %>%
  mutate(percent = count/sum(count))
#define not in list operator
`%!in%` <- Negate(`%in%`)

val_df <- Ab_DF %>% filter(dataset=='val set')%>%
  mutate(state=case_when(clonotype %in% train_clono~'overlap with train set',
                         clonotype %!in% train_clono~'unique'))
val_df2 <- val_df %>% group_by(dataset,state)%>%
  dplyr::summarise(count = n()) %>%
  mutate(percent = count/sum(count))
#print the precentage of overlapped state in val set
print(prop.table(table(val_df$state)))
test_df <-Ab_DF %>% filter(dataset=='test set')%>%
  mutate(state=case_when(clonotype %in% train_clono~'overlap with train set',
                         clonotype %!in% train_clono~'unique'))
test_df2 <- test_df %>% group_by(dataset,state)%>%
  dplyr::summarise(count = n()) %>%
  mutate(percent = count/sum(count))
#print the precentage of overlapped state in test set
print(prop.table(table(test_df$state)))
p1<-plot_hist(train_df,train_df2,train_df$clonotype,'Clonotype ID',11,24,'')
p2<-plot_hist(val_df,val_df2,val_df$clonotype,'Clonotype ID',3.2,6.8,'')
p3<-plot_hist(test_df,test_df2,test_df$clonotype,'Clonotype ID',5,11,'Distribution of clonotype ID')

p <- ggarrange(p3,p2,p1,nrow=3,ncol=1,common.legend = TRUE,legend = 'top')
ggsave('graph/Clonotype_distribution.png',p,height=5,width=3.4,bg='white')
