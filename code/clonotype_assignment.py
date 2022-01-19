#!/usr/bin/python
import os
import sys
import glob
import pandas as pd
from Bio.Seq import Seq
from collections import defaultdict, Counter

def read_ab_info(filename):
  aas = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
  print ('reading: %s' % filename)
  datasheet = pd.read_excel(filename, sheet_name='main')
  Ab_dict = {}
  for index, info_dict in datasheet.iterrows():
    for i in info_dict.keys():
      info_dict[i] = str(info_dict[i])
    Ab_dict[info_dict['Name']] = info_dict
  return Ab_dict
  
def read_cluster_info(filename):
  infile = open(filename, 'r')
  CDRH3_cluster_dict = {}
  for line in infile.readlines():
    if 'cluster' in line: continue
    assert(line.count('|')==1)
    line = line.rstrip().rsplit("\t")
    Ab      = line[0].rsplit('|')[1]
    cluster = line[2]
    CDRH3_cluster_dict[Ab] = cluster 
  infile.close()
  return CDRH3_cluster_dict

def clean_clonotype(clonotype_dict):
  new_dict = {}
  count = 0
  for ID in sorted(clonotype_dict.keys(), key=lambda x:len(clonotype_dict[x]), reverse=True):
     count += 1
     for Ab in clonotype_dict[ID]:
       new_dict[Ab] = count
  return new_dict

def clonotype_assignment(Ab_dict, CDRH3_cluster_dict):
  clonotype_dict = defaultdict(list)
  for Ab in Ab_dict.keys():
    HV = Ab_dict[Ab]['Heavy V Gene'].rsplit('*')[0]
    HJ = Ab_dict[Ab]['Heavy J Gene'].rsplit('*')[0]
    LV = Ab_dict[Ab]['Light V Gene'].rsplit('*')[0]
    LJ = Ab_dict[Ab]['Light J Gene'].rsplit('*')[0]
    CDRH3 = Ab_dict[Ab]['CDRH3']
    if Ab not in CDRH3_cluster_dict.keys(): continue
    CDRH3_cluster = CDRH3_cluster_dict[Ab]
    if 'nan' in [HV, HJ, LV, LJ, CDRH3, CDRH3_cluster]: continue
    clonotype_ID = '|'.join(map(str,[HV, HJ, LV, LJ, len(CDRH3), CDRH3_cluster]))
    clonotype_dict[clonotype_ID].append(Ab)
  clonotype_assign = clean_clonotype(clonotype_dict)
  return clonotype_assign

def write_clonotype(clonotype_assign, filename, outfile_clonotype):
  df = pd.read_excel(filename, sheet_name='main')
  df['clonotype'] = df['Name'].map(clonotype_assign)
  df = df[df.clonotype.notnull()]
  df = df[["Name", "clonotype"]]
  print ('writing: %s' % outfile_clonotype)
  df.to_csv(outfile_clonotype, sep="\t", index=False)

def main():
  filename           = 'data/SARS-CoV-2-Abs.xlsx'
  CDRH3_cluster      = 'result/CDRH3_cluster.tsv'
  outfile_clonotype  = 'result/clonotype_assign.tsv'
  CDRH3_cluster_dict = read_cluster_info(CDRH3_cluster)
  Ab_dict            = read_ab_info(filename)
  clonotype_assign   = clonotype_assignment(Ab_dict, CDRH3_cluster_dict)
  write_clonotype(clonotype_assign, filename, outfile_clonotype)

if __name__ == "__main__":
  main()
