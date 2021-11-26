# Title     : MERGE Abs
# Objective :
# Created by: yiquan
# Created on: 7/27/21

library(readr)
library(dplyr)
library("tidyr")
library(stats)
library("writexl")

files <- list.files(path = "/Users/yiquan/PycharmProjects/antibody project/CoV_Ab_Dataset/result/igblast_result7", pattern = "tsv.gz", full.names = T)
merfiles <- sapply(files,read_tsv,simplify = FALSE) %>%
  bind_rows()
df <- merfiles %>%
  #separate(sequence_id, c("Name", "GenBank_ID"), "\\+")%>% #only for orignial separate Ab name and Genbank ID
  mutate(type=recode(locus,'IGH'="VH",.default ="VL")) %>%
  select(sequence_id,type,locus,sequence,sequence_alignment_aa,v_call,d_call,j_call,cdr1_aa,cdr2_aa,cdr3_aa) #better use sequence rather than sequence alignment
df < df[order(df$sequence_id),]
row.names(df)<-NULL
#separate the df by type
df_H <- df %>%
  filter(type=='VH')%>%
  distinct(sequence_id, .keep_all = TRUE) %>%
  rename(VH_nuc=sequence,
         VH_AA=sequence_alignment_aa,
         Name=sequence_id,
         #VH_GenBank_ID=GenBank_ID,
         Heavy_V_gene = v_call,
         Heavy_J_gene = j_call,
         Heavy_D_gene = d_call,
         CDRH1_AA=cdr1_aa,
         CDRH2_AA=cdr2_aa,
         CDRH3_AA=cdr3_aa)
df_L <- df %>%
  filter(type=='VL')%>%
  distinct(sequence_id, .keep_all = TRUE)%>%
  rename(VL_nuc=sequence,
         VL_AA=sequence_alignment_aa,
         Name=sequence_id,
         #VL_GenBank_ID=GenBank_ID,
         Light_V_gene = v_call,
         Light_J_gene = j_call,
         CDRL1_AA=cdr1_aa,
         CDRL2_AA=cdr2_aa,
         CDRL3_AA=cdr3_aa)
DF <- merge(df_H,df_L,by='Name',all = TRUE) %>%
  select(Name,VH_nuc,VH_AA,VL_nuc,VL_AA,Heavy_V_gene,Heavy_J_gene,Heavy_D_gene,Light_V_gene,Light_J_gene,CDRL1_AA,CDRL2_AA,CDRL3_AA,CDRH1_AA,CDRH2_AA,CDRH3_AA)
write_xlsx(DF, "result/igblast_result7/HA-S_CDR_table.xlsx")