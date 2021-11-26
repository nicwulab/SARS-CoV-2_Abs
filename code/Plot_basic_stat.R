# Title     : to plot basic statistics of Ab table
# Objective :
# Created by: yiquan
# Created on: 10/7/21
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
library(e1071)
library(ggforce)
library(ggbeeswarm)
library(ggExtra)
require(cowplot)
require(ggseqlogo)
#define plot function
##bar plot


#-------------------------------------Function--------------------------------------#
#set factor


reform2heatmap <-function (subclass_df,path,x_lab,y_lab,title){
  subclass_V_df<-data.frame(table(select(subclass_df, HV, LV))) %>%
    filter(LV!='NA')
  subclass_V_df[subclass_V_df==0] <- NA
  #set factor level
  HV_level <- unique(subclass_V_df$HV)
  LV_level <- unique(subclass_V_df$LV)
  #get donor df
  subclass_donor_df <- subclass_df %>%
    group_by(HV,LV,`Patient ID`) %>%
    summarize(donor=n())%>%
    group_by(HV,LV) %>%
    summarize(donor=n())%>%
    filter(donor>=2) %>%
    mutate(HV=factor(as.character(HV),levels=HV_level))%>%
    mutate(LV=factor(as.character(LV),levels=LV_level))
  subclass_donor_df$HV <- as.integer(subclass_donor_df$HV)
  subclass_donor_df$LV <- as.integer(subclass_donor_df$LV)
  plot_heat_map(subclass_V_df,subclass_donor_df,path,x_lab,y_lab,title)
}
plot_bar_counts <- function(DF,column,column2,path,x_lab,y_lab){
  p <-ggplot(DF) +
  theme_classic() +
  geom_bar(aes(x=column,fill=column2))+theme(plot.title=element_blank(),
                axis.text.x=element_text(size=7,angle=45,hjust=0.5,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(size=7,angle=0,hjust=1,vjust=0.5,colour = 'black'),
                axis.title.y=element_text(size=7,face="bold"),
                axis.title.x=element_text(size=7,face="bold"),
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                # legend.title=element_blank(),
                legend.text=element_text(size=7,face="bold"),
                legend.position='right')+
  scale_color_brewer(palette = "PuOr")+ labs(fill = "Neutralization") +
  xlab(x_lab) +
  ylab(y_lab)
  ggsave(path,p,height = 4,width = 6,bg='white')

}

plot_heat_map <- function (heatmap_DF,highlight_df,path,x_lab,y_lab,title){

  V_heatp <-ggplot(heatmap_DF, aes(x = LV, y = HV, fill = Freq)) +
    geom_tile(color = "gray") +
    scale_fill_gradient2(low = "#075AFF", mid = "#FFFFCC", high = "#FF0000", na.value="gray100") +
    geom_text(aes(label = Freq), color = "black", size = 2) +
    coord_fixed()+
    theme(plot.title = element_text(color = "black", size = 9, face = "bold",hjust = 0.5),
          axis.text.x=element_text(size=7,angle=90,hjust=0.5,vjust=0.5,colour = 'black'),
          axis.text.y=element_text(size=7,angle=0,hjust=1,vjust=0.5,colour = 'black'),
          axis.title.y=element_text(size=7,face="bold"),
          axis.title.x=element_text(size=7,face="bold"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.key.size=unit(3,"mm"),
          legend.text=element_text(size=7,face="bold"),
          legend.position='right') +
    geom_rect(data=highlight_df,size=0.2, fill=NA, colour="black", aes(xmin=LV-0.5, xmax=LV+0.5, ymin=HV-0.5, ymax=HV+0.5)) +
    xlab(x_lab) + ylab(y_lab)+ggtitle(title)
  ggsave(path,V_heatp,height = 8,width = 10,bg='white')
}

plot_hist <- function(df,column,path,x_lab,main_title){

  p <- ggplot(df,aes(x=column)) +
         geom_histogram(binwidth=1) +
         theme_cowplot(12) +
         theme(plot.title=element_text(size=7,face="bold", hjust=0.5),
               axis.title=element_text(size=7,face="bold"),
               axis.text=element_text(size=7,face="bold"),
               legend.key.size=unit(0.1,'in'),
               legend.spacing.x=unit(0.03, 'in'),
               legend.title=element_blank(),
               legend.text=element_text(size=7,face="bold"),
               legend.position='right') +
     ggtitle(main_title) +
     scale_fill_manual(values=c('black'),drop=FALSE) + xlab(x_lab) +
          ylab("Counts")# coord_cartesian(xlim=c(-2.5,1))
  ggsave(path,p,height = 2,width = 2,bg='white')
  }


plot_scatter <- function (df,col_x,col_y,Xlab,Ylab,title,path){
  p <- ggplot(df, aes(x=col_x, y=col_y)) +
    geom_point(pch=16, alpha=0.5,color='gray30')+
    theme_cowplot(12) +
    theme(plot.title=element_text(size=9,face="bold", hjust=0.5),
          axis.title=element_text(size=7,face="bold"),
          axis.text=element_text(size=7,face="bold"),
          legend.key.size=unit(0.1,'in'),
          legend.spacing.x=unit(0.03, 'in'),
          legend.title=element_blank(),
          legend.text=element_text(size=9,face="bold"),
          legend.position='right') +
    xlab(Xlab) +
    ylab(Ylab) +
    ggtitle(title)
  ggsave(path,p,height = 3,width = 5,bg='white')
}



#-------------------------------------read dataframe--------------------------------------#
Ab_DF <- read_excel('data/SARS-CoV-2-Abs.xlsx')%>%
  separate(`Heavy V Gene`,into = c("HV", "HV_rest"),sep = '\\*',extra = "merge") %>%
  separate(`Heavy J Gene`,into = c("HJ", "HJ_rest"),sep = '\\*',extra = "merge") %>%
  separate(`Light V Gene`,into = c("LV", "LV_rest"),sep = '\\*',extra = "merge") %>%
  separate(`Light J Gene`,into = c("LJ", "LJ_rest"),sep = '\\*',extra = "merge")
Ab_DF$LV <- str_replace(Ab_DF$LV, "IGLV", "L")
Ab_DF$LV <- str_replace(Ab_DF$LV, "IGKV", "K")
Ab_DF$HV <- str_replace(Ab_DF$HV, "IGHV", "H")


#----------------------------------sub class dataframe--------------------------------------#
RBD_DF <- Ab_DF %>%
  filter(Ab_DF$`Protein + Epitope`=='S:RBD')

#----------------------------------correlation scatter plot--------------------------------------#
Neutral_RBD_DF <- RBD_DF %>%
  filter(grepl("SARS-CoV-2",RBD_DF$`Neutralising Vs`)) %>%
  group_by(HV,LV) %>%
  summarize(Freq=n())
nonNeutral_RBD_DF <- RBD_DF %>%
  filter(grepl("SARS-CoV-2",RBD_DF$`Not Neutralising Vs`)) %>%
  group_by(HV,LV) %>%
  summarize(Freq=n())
RBD_neutral_DF <- merge(Neutral_RBD_DF,nonNeutral_RBD_DF,by=c('HV','LV'),all = TRUE)
RBD_neutral_DF[is.na(RBD_neutral_DF)] <-0

plot_scatter(RBD_neutral_DF,RBD_neutral_DF$Freq.x,RBD_neutral_DF$Freq.y,'# of neutraling antibodies','# of non-neutraling antibodies',
             'RBD antibodies with shared\nIGHV/IGK(L)V pair','graph/RBD_neutralization_sactter.png')
