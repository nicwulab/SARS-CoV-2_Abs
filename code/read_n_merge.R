# Title     : MERGE Abs
# Objective :
# Created by: yiquan
# Created on: 7/27/21

library(readr)
library(dplyr)
library("tidyr")
library(stats)


files <- list.files(path = "/Users/yiquan/PycharmProjects/antibody project/CoV_Ab_Dataset/result", pattern = "Ab_ls", full.names = T)
merfiles <- sapply(files,read_tsv,simplify = FALSE) %>%
  bind_rows()
df <- merfiles %>%
  separate(sequence_id, c("mAb_name", "GenBank_ID"), "\\+")%>%
  mutate(type=recode(locus,'IGH'="VH",.default ="VL")) %>%
  select(mAb_name,type,GenBank_ID,locus,sequence_alignment,sequence_alignment_aa,v_call,d_call,j_call,cdr3_aa)
df < df[order(df$mAb_name),]
row.names(df)<-NULL
#separate the df by type
df_H <- df %>%
  filter(type=='VH')%>%
  rename(VH_nuc=sequence_alignment,
         VH_AA=sequence_alignment_aa,
         VH_GenBank_ID=GenBank_ID,
         Heavy_V_gene = v_call,
         Heavy_J_gene = j_call,
         Heavy_D_gene = d_call,
         CDRH3_AA=cdr3_aa)
df_L <- df %>%
  filter(type=='VL')%>%
  rename(VL_nuc=sequence_alignment,
         VL_AA=sequence_alignment_aa,
         VL_GenBank_ID=GenBank_ID,
         Light_V_gene = v_call,
         Light_J_gene = j_call,
         CDRL3_AA=cdr3_aa)
DF <- merge(df_H,df_L,by='mAb_name',all = TRUE) %>%
  select(mAb_name,VH_nuc,VH_AA,VL_nuc,VL_AA,Heavy_V_gene,Heavy_J_gene,Heavy_D_gene,Light_V_gene,Light_J_gene,CDRH3_AA,CDRL3_AA,VH_GenBank_ID,VL_GenBank_ID)
write.table(DF, file = "result/Ab_list.tsv", sep = "\t",col.names = NA)