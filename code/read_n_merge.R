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
df_reshape <-reshape(df,
             timevar = 'type',
             idvar = 'mAb_name',
             direction = 'wide')
write.table(df, file = "result/Ab_list.tsv", sep = "\t",col.names = NA)