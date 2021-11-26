# Title     : to plot basic statistics of Ab table
# Created by: yiquan
# Created on: 10/22/21

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



plot_point_heatmap <- function (heatmap_DF,group_col,path,x_lab,y_lab,title){
  p  <- ggplot(heatmap_DF, aes(x=HV,y=LV))+
          geom_point(aes(fill=group_col,size=abs(freq)),alpha = 0.5,group=group_col,color='black',pch=21) +
          theme_classic() + scale_fill_brewer(palette = "Set2")+
          theme(plot.title=element_text(size=7,face="bold",hjust = 0.5),
                text = element_text(size=7,face="bold"),
                legend.position = 'right',
                legend.key.size = unit(1, 'lines'),
                legend.background=element_blank(),
                axis.line = element_line(colour = 'black', size = 0),
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                axis.text.y=element_text(size=6.5,face="bold",color='black'),
                axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=6.5,face="bold",color='black')) + ggtitle(title) +xlab(x_lab)+ylab(y_lab)+ labs(fill='Epitope',size='Freq')
  ggsave(path,p,height = 8,width = 7,bg='white')

}


#-------------------------------------read dataframe--------------------------------------#


Ab_DF <- read_excel('result/SARS-CoV-2-Abs_v29.xlsx')%>%
  separate(`Heavy V Gene`,into = c("HV", "HV_rest"),sep = '\\*',extra = "merge") %>%
  separate(`Heavy J Gene`,into = c("HJ", "HJ_rest"),sep = '\\*',extra = "merge") %>%
  separate(`Heavy D Gene`,into = c("HD", "HD_rest"),sep = '\\*',extra = "merge") %>%
  separate(`Light V Gene`,into = c("LV", "LV_rest"),sep = '\\*',extra = "merge") %>%
  separate(`Light J Gene`,into = c("LJ", "HJ_rest"),sep = '\\*',extra = "merge")
Ab_DF$LV <- str_replace(Ab_DF$LV, "IGLV", "L")
Ab_DF$LV <- str_replace(Ab_DF$LV, "IGKV", "K")
Ab_DF$HV <- str_replace(Ab_DF$HV, "IGHV", "H")
Ab_DF$HD <- str_replace(Ab_DF$HD, "IGHD", "D")

#-------------------------------------V gene neutralization profiling--------------------------------------#

#sort out datafram

neural_V_df<-Ab_DF %>%mutate(neutral=case_when(`Neutralising Vs`=="SARS-CoV-2"~'Yes',
                                               `Not Neutralising Vs`=="SARS-CoV-2"~'No')) %>%
  filter(neutral != 'NA')%>%
  group_by(HV,LV,neutral) %>%
  summarize(counts=n())%>%
  group_by(neutral)%>%
  mutate(freq = counts / sum(counts))
#set factor level
HV_level <- unique(neural_V_df$HV)
LV_level <- unique(neural_V_df$LV)
neural_V_df <- neural_V_df %>%
    mutate(HV=factor(as.character(HV),levels=HV_level))%>%
    mutate(LV=factor(as.character(LV),levels=LV_level))%>%
  filter(HV!='NA')%>%
  filter(LV!='NA')
plot_point_heatmap(neural_V_df,neural_V_df$neutral,'graph/HLV_neutral_heatmap.png','IGHV','IGLV/IGKV','')
#-------------------------------------V gene Epitope profiling--------------------------------------#

#sort out datafram
Ab_DF$`Protein + Epitope` <- str_replace(Ab_DF$`Protein + Epitope`,'S:S Stem Helix','S:S2')
point_HM_V_df<-Ab_DF %>%filter(`Protein + Epitope`=='S:RBD'| `Protein + Epitope`=='S:NTD'|`Protein + Epitope`=='S:S2')%>%
  group_by(HV,LV,`Protein + Epitope`) %>%
  summarize(counts=n())%>%
  group_by(`Protein + Epitope`)%>%
  mutate(freq = counts / sum(counts))
#set factor level
HV_level <- unique(point_HM_V_df$HV)
LV_level <- unique(point_HM_V_df$LV)
point_HM_V_df <- point_HM_V_df %>%
    mutate(HV=factor(as.character(HV),levels=HV_level))%>%
    mutate(LV=factor(as.character(LV),levels=LV_level))%>%
  filter(HV!='NA')%>%
  filter(LV!='NA')
plot_point_heatmap(point_HM_V_df,point_HM_V_df$`Protein + Epitope`,'graph/HLV_epitope_heatmap.png','IGHV','IGLV/IGKV','')
# #find same Vgene + Dgene but different Epitope
# M_df <- point_HM_V_df%>%
#   group_by(HV,LV)%>%
#   summarize(donor=n())%>%
#   filter(donor>=2)
# M_D_df <- merge(M_df,Ab_DF,on=c('HV','LV')) %>%
#   group_by(HV,LV,HD,`Protein + Epitope`)%>%
#   select(HV,LV,HD,`Protein + Epitope`)%>%
#   summarize(counts=n())%>%
#   group_by(HV,LV,HD)%>%
#   summarize(counts=n())%>%
#   filter(counts>=2)
