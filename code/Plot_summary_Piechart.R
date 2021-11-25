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
library(ggforce)
library(ggbeeswarm)
library(ggExtra)
require(cowplot)
library(forcats)
library(ggpubr)


plot_pie <- function (data,column,title,set){
  p <- ggplot(data, aes(x="", y=Freq, fill=fct_inorder(column))) +
    geom_col(width=1, color=1,position = position_stack(reverse = TRUE)) +
    coord_polar("y" ) +
    scale_fill_brewer(palette=set)+
    geom_label_repel(aes(y=pos,label = paste0(Freq, "%")), nudge_x = 1,show.legend = FALSE, size=1.5) +
    theme_void() +
    theme(legend.position="bottom",
          legend.key.size = unit(0.25, 'cm'),
          legend.text = element_text(size=6),
          legend.title = element_text(size=7,face="bold"))+
    guides(fill = guide_legend(title = title,nrow = 4,byrow=TRUE))
  return(p)
}
plot_beeswarm <- function (df,column,path){
  p <- ggplot(df,aes(x=column,y=1)) +
    geom_beeswarm(color = "black",size=0.5,
                   groupOnX = FALSE)+
    theme_cowplot(12) +
    theme(plot.title=element_text(size=7,face="bold", hjust=0.5),
          axis.title=element_text(size=7,face="bold"),
          axis.text=element_text(size=7,face="bold"),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_line(colour = 'black', size = 0),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.key.size=unit(0.1,'in'),
          legend.spacing.x=unit(0.03, 'in'),
          legend.title=element_blank(),
          legend.text=element_text(size=6.5,face="bold"),
          legend.position='right') +
     scale_fill_manual(values=c('black'),drop=FALSE) +ylab('Sources')+
          xlab("Counts")# coord_cartesian(xlim=c(-2.5,1))
  p3 <- ggMarginal(p, margins = 'x', size=1,type="histogram",fill = "slateblue", xparams = list(  bins=100))
  ggsave(path,p3,height = 2,width = 6,bg='white')
}


Ab_DF <- read_excel('result/SARS-CoV-2-Abs_v37.xlsx')%>%
  separate(`Heavy V Gene`,into = c("HV", "HV_rest"),sep = '\\*',extra = "merge") %>%
  separate(`Heavy J Gene`,into = c("HJ", "HJ_rest"),sep = '\\*',extra = "merge") %>%
  separate(`Light V Gene`,into = c("LV", "LV_rest"),sep = '\\*',extra = "merge") %>%
  separate(`Light J Gene`,into = c("LJ", "LJ_rest"),sep = '\\*',extra = "merge")
Ab_DF$LV <- str_replace(Ab_DF$LV, "IGLV", "L")
Ab_DF$LV <- str_replace(Ab_DF$LV, "IGKV", "K")
Ab_DF$HV <- str_replace(Ab_DF$HV, "IGHV", "H")

source_count <-data.frame(table(Ab_DF$Sources))
origin_count <-data.frame(table(Ab_DF$Origin))
patient_count <-data.frame(table(Ab_DF$`Patient ID`))
epitope_count <-data.frame(table(Ab_DF$`Protein + Epitope`))
# epitope_count$Var1 <- str_replace(epitope_count$Var1, "S:S2 Stem Helix", "S:S2")
epitope_df<-epitope_count%>%
  mutate(Freq=round(Freq*100 / sum(epitope_count$Freq),1))%>%
  arrange(desc(Freq))
epitope_df$pos = (cumsum(c(0, epitope_df$Freq)) + c(epitope_df$Freq /2 , .01))[1:nrow(epitope_df)]
epitope_p <- plot_pie(epitope_df,epitope_df$Var1,'Epitope','Set3')

origin_df<-origin_count%>%
  mutate(Freq=round(Freq*100 / sum(origin_count$Freq),2))%>%
  arrange(desc(Freq))
origin_df$Var1 <- str_replace(origin_df$Var1, "B-cells;", "")
origin_df$Var1 <- str_replace(origin_df$Var1, "hybridoma;", "")

origin_df$pos = (cumsum(c(0, origin_df$Freq)) + c(origin_df$Freq /2 , .01))[1:nrow(origin_df)]
origin_p <- plot_pie(origin_df,origin_df$Var1,'Origin','Accent')
ggsave('graph/origin_pie.png',origin_p,width = 3.9,height = 2.2,bg='white', dpi=1200)
ggsave('graph/epitope_pie.png',epitope_p,width = 2.7,height = 2.2,bg='white', dpi=1200)

plot_beeswarm(source_count,source_count$Freq,'graph/source_beeswarm.png')