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


#-------------------------------------Function--------------------------------------#



plot_point_heatmap <- function (heatmap_DF,group_col,path,x_lab,y_lab,title){
  p  <- ggplot(heatmap_DF, aes(x=LV,y=HV))+
          geom_point(aes(fill=group_col,size=abs(freq)),alpha = 0.5,group=group_col,color='black',pch=21) +
          theme_classic() + 
          scale_fill_brewer(palette = "Accent", labels = c("NTD", "RBD", "S2")) +
          scale_size_continuous(name = "Freq (%)",
                                range = c(0.2, 2.7)) +
          theme(plot.title=element_text(size=7,face="bold",hjust = 0.5),
                text = element_text(size=7,face="bold"),
                legend.position = 'right',
                legend.key.size = unit(1, 'lines'),
                legend.background=element_blank(),
                legend.text=element_text(size=7,face="bold",color='black'),
                legend.spacing.y = unit(0.02, 'lines'),
                axis.line = element_line(colour = 'black', size = 0),
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                axis.text.y=element_text(size=5,face="bold",color='black'),
                axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=5,face="bold",color='black')) +
          ggtitle(title) +xlab(x_lab)+ylab(y_lab)+ labs(fill='Epitope',size='Freq')
  ggsave(path,p,height = 4,width = 6,bg='white')

}


#-------------------------------------read dataframe--------------------------------------#


Ab_DF <- read_excel('data/SARS-CoV-2-Abs.xlsx')%>%
  separate(`Heavy V Gene`,into = c("HV", "HV_rest"),sep = '\\*',extra = "merge") %>%
  separate(`Heavy J Gene`,into = c("HJ", "HJ_rest"),sep = '\\*',extra = "merge") %>%
  separate(`Heavy D Gene`,into = c("HD", "HD_rest"),sep = '\\*',extra = "merge") %>%
  separate(`Light V Gene`,into = c("LV", "LV_rest"),sep = '\\*',extra = "merge") %>%
  separate(`Light J Gene`,into = c("LJ", "HJ_rest"),sep = '\\*',extra = "merge")
Ab_DF$LV <- str_replace(Ab_DF$LV, "IGLV", "L")
Ab_DF$LV <- str_replace(Ab_DF$LV, "IGKV", "K")
Ab_DF$HV <- str_replace(Ab_DF$HV, "IGHV", "H")
Ab_DF$HD <- str_replace(Ab_DF$HD, "IGHD", "D")

#-------------------------------------V gene Epitope profiling--------------------------------------#

#sort out datafram
Ab_DF$`Protein + Epitope` <- str_replace(Ab_DF$`Protein + Epitope`,'S:S Stem Helix','S:S2')
Ab_DF$`Protein + Epitope` <- str_replace(Ab_DF$`Protein + Epitope`,'S:S1 non-RBD','S:NTD')
point_HM_V_df<-Ab_DF %>% filter(`Protein + Epitope`=='S:RBD'| `Protein + Epitope`=='S:NTD'|`Protein + Epitope`=='S:S2')%>%
  group_by(HV,LV,`Protein + Epitope`) %>%
  summarize(counts=n())%>%
  group_by(`Protein + Epitope`)%>%
  mutate(freq = counts / sum(counts) * 100)
#set factor level
LVgene_ref  <- read_excel('data/LV_Repertoire_freq.xlsx')
LV_level <- LVgene_ref$V_gene %>%
              str_replace('IGKV','K') %>%
              str_replace('IGLV','L') %>%
              sort()
HVgene_ref  <- read_excel('data/HV_Repertoire_freq.xlsx')
HV_level <- HVgene_ref$V_gene %>%
              str_replace('IGHV','H') %>%
              sort()
point_HM_V_df <- point_HM_V_df %>%
    mutate(HV=factor(as.character(HV),levels=HV_level))%>%
    mutate(LV=factor(as.character(LV),levels=LV_level))%>%
  filter(HV!='NA')%>%
  filter(LV!='NA')
print (table(point_HM_V_df$`Protein + Epitope`))
plot_point_heatmap(point_HM_V_df,point_HM_V_df$`Protein + Epitope`,'graph/HLV_epitope_heatmap.png','IGK(L)V','IGHV','')
