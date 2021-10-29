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

plot_beeswarm <- function (df,column,path){
  p <- ggplot(df,aes(x=group,y=column)) +
    geom_beeswarm(cex = 0.5)+
    theme_cowplot(12) +
    theme(plot.title=element_text(size=7,face="bold", hjust=0.5),
          axis.title=element_text(size=7,face="bold"),
          axis.text=element_text(size=7,face="bold"),
          legend.key.size=unit(0.1,'in'),
          legend.spacing.x=unit(0.03, 'in'),
          legend.title=element_blank(),
          legend.text=element_text(size=7,face="bold"),
          legend.position='right') +
     scale_fill_manual(values=c('black'),drop=FALSE) +xlab('')+
          ylab("Counts")# coord_cartesian(xlim=c(-2.5,1))
  ggsave(path,p,height = 3,width = 2,bg='white')
}

plot_pie <- function (data,column,title,wid,set,path){
  data <- data %>%
  mutate(prop = Freq / sum(data$Freq) *100)
  p <- ggplot(data, aes(x="", y=prop, fill=column)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() +
    theme(legend.position="right")+scale_fill_brewer(palette=set)+ labs(fill  = title)
  ggsave(path,p,height = 3,width = wid,bg='white')
}

plot_scatter <- function (df,col_x,col_y,Xlab,Ylab,path){
  p <- ggplot(df, aes(x=col_x, y=col_y)) +
    geom_point()+
    theme_cowplot(12) +
    theme(plot.title=element_text(size=7,face="bold", hjust=0.5),
          axis.title=element_text(size=7,face="bold"),
          axis.text=element_text(size=7,face="bold"),
          legend.key.size=unit(0.1,'in'),
          legend.spacing.x=unit(0.03, 'in'),
          legend.title=element_blank(),
          legend.text=element_text(size=7,face="bold"),
          legend.position='right') +
    xlab(Xlab)+ylab(Ylab)
  ggsave(path,p,height = 3,width = 3,bg='white')
}



#-------------------------------------read dataframe--------------------------------------#


Ab_DF <- read_excel('result/SARS-CoV-2-Abs_v28-9.xlsx')%>%
  separate(`Heavy V Gene`,into = c("HV", "HV_rest"),sep = '\\*',extra = "merge") %>%
  separate(`Heavy J Gene`,into = c("HJ", "HJ_rest"),sep = '\\*',extra = "merge") %>%
  separate(`Light V Gene`,into = c("LV", "LV_rest"),sep = '\\*',extra = "merge") %>%
  separate(`Light J Gene`,into = c("LJ", "LJ_rest"),sep = '\\*',extra = "merge")
Ab_DF$LV <- str_replace(Ab_DF$LV, "IGLV", "L")
Ab_DF$LV <- str_replace(Ab_DF$LV, "IGKV", "K")
Ab_DF$HV <- str_replace(Ab_DF$HV, "IGHV", "H")


#----------------------------------sub class dataframe--------------------------------------#

neutral_DF <- Ab_DF %>%
  filter(Ab_DF$`Neutralising Vs`=="SARS-CoV-2")
nonneutral_DF <- Ab_DF %>%
  filter(Ab_DF$`Not Neutralising Vs`=="SARS-CoV-2")
crossbinder_DF <- Ab_DF %>%
  filter(grepl(',',`Binds to`))

vaccine_DF <- Ab_DF %>%
  filter(grepl('Vaccinee',Origin))
patient_DF <- Ab_DF %>%
  filter(grepl('Patient',Origin))

nonRBD_DF <- Ab_DF %>%
  filter(Ab_DF$`Protein + Epitope`!='S:RBD')
# nonRBD_DF$CDRH3_Length <- str_length(nonRBD_DF$CDRH3)
# nonRBD_DF$CDRL3_Length <- str_length(nonRBD_DF$CDRL3)

RBD_DF <- Ab_DF %>%
  filter(Ab_DF$`Protein + Epitope`=='S:RBD')
# RBD_DF$CDRH3_Length <- str_length(RBD_DF$CDRH3)
# RBD_DF$CDRL3_Length <- str_length(RBD_DF$CDRL3)

NTD_DF <- Ab_DF %>%
  filter(Ab_DF$`Protein + Epitope`=='S:NTD')

S2_DF <- Ab_DF %>%
  filter(Ab_DF$`Protein + Epitope`=='S:S2'|Ab_DF$`Protein + Epitope`=='S:S2 Stem Helix')

#-------------------------------------Plot--------------------------------------#

reform2heatmap(NTD_DF,'graph/NTD-V_heatmap.png',"IGLV/IGKV", "IGHV",'NTD')
reform2heatmap(S2_DF,'graph/S2-V_heatmap.png',"IGLV/IGKV", "IGHV",'S2')
reform2heatmap(RBD_DF,'graph/RBD-V_heatmap.png',"IGLV/IGKV", "IGHV",'RBD')
reform2heatmap(neutral_DF,'graph/neutral-V_heatmap.png',"IGLV/IGKV", "IGHV",'Neutralizing antibody')
reform2heatmap(nonneutral_DF,'graph/nonneutral-V_heatmap.png',"IGLV/IGKV", "IGHV",'Non-neutralizing antibody')
reform2heatmap(crossbinder_DF,'graph/crossbinder-V_heatmap.png',"IGLV/IGKV", "IGHV",'Cross-binder')
reform2heatmap(vaccine_DF,'graph/vaccine-V_heatmap.png',"IGLV/IGKV", "IGHV",'Vaccine')
reform2heatmap(patient_DF,'graph/patient-V_heatmap.png',"IGLV/IGKV", "IGHV",'Patient')


#----------------------------------count column number--------------------------------------#

source_count <-data.frame(table(Ab_DF$Sources)) %>%
  mutate(group='Sources')
origin_count <-data.frame(table(Ab_DF$Origin)) %>%
  mutate(group='Origin')
patient_count <-data.frame(table(Ab_DF$`Patient ID`)) %>%
  mutate(group='patient')
epitope_count <-data.frame(table(Ab_DF$`Protein + Epitope`))%>%
  mutate(group='epitope')

#rbind multi dataframe
# beeswarm_df <-do.call("rbind",list(source_count,origin_count,patient_count,epitope_count))
plot_beeswarm(source_count,source_count$Freq,'graph/source_beeswarm.png')
plot_pie(epitope_count,epitope_count$Var1,'Epitope',3,'Set3','graph/epitope_piechart.png')
plot_pie(origin_count,origin_count$Var1,'Origin',5,'Set2','graph/origin_piechart.png')

#----------------------------------correlation scatter plot--------------------------------------#
## vaccine vs patient
Freq_patient_DF <- patient_DF %>%
  filter(patient_DF$`Protein + Epitope`=='S:RBD') %>%
  group_by(HV,LV) %>%
  summarize(Freq=n())
Freq_vaccine_DF <- vaccine_DF %>%
  filter(vaccine_DF$`Protein + Epitope`=='S:RBD') %>%
  group_by(HV,LV) %>%
  summarize(Freq=n())
pati_vacc_DF <- merge(Freq_patient_DF,Freq_vaccine_DF,by=c('HV','LV'),all = TRUE)
pati_vacc_DF[is.na(pati_vacc_DF)] <-0

# #statistical for outlier identification
# pati_vacc_DF.model = lm(Freq.y~Freq.x,data = pati_vacc_DF)
# summary(pati_vacc_DF.model)
# yhat = pati_vacc_DF.model$fit
# t = rstudent(pati_vacc_DF.model)
# d=rstandard(pati_vacc_DF.model)
# cooksd = cooks.distance(pati_vacc_DF.model)
# cooksd[abs(cooksd)>1]
# d[abs(d)>6]
# t[abs(t)>6]
# plot(
#           yhat, d, main = "Versus fits (response is time)",
#                             pch="*",
#                             xlab="Fitted Values",
#                             ylab="Standerdized Residuals "
#        )
# abline(h = 0, col="red")
# abline(h = 6, col="red")
# abline(h = -6, col="red")

## Bait by Spike vs Bait by RBD
Bait_RBD_DF <- RBD_DF %>%
  filter(RBD_DF$Bait=='SARS-CoV-2 RBD') %>%
  group_by(HV,LV) %>%
  summarize(Freq=n())
Bait_S_DF <- RBD_DF %>%
  filter(RBD_DF$Bait=='SARS-CoV-2 spike') %>%
  group_by(HV,LV) %>%
  summarize(Freq=n())
RBD_S_DF <- merge(Bait_RBD_DF,Bait_S_DF,by=c('HV','LV'),all = TRUE)
RBD_S_DF[is.na(RBD_S_DF)] <-0

plot_scatter(RBD_S_DF,RBD_S_DF$Freq.x,RBD_S_DF$Freq.y,'Baited by RBD','Baited by Spike','graph/RBD_S_RBD_sactter.png')