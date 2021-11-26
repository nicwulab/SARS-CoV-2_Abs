#!/usr/bin/python
import os
import sys
import glob
import pandas as pd
from collections import defaultdict, Counter

def extract_Dgene(filename):
  print ('reading: %s' % filename)
  Dgene_dict = defaultdict(list)
  datasheet = pd.read_excel(filename, sheet_name='main')
  for index, info_dict in datasheet.iterrows():
    for i in info_dict.keys():
      info_dict[i] = str(info_dict[i])
    if info_dict['Heavy D Gene'] == 'nan' or info_dict['Heavy D Gene'] == 'NA': continue
    if info_dict['Protein + Epitope'] not in ['S:NTD','S:RBD','S:S2','S:S2 Stem Helix','S:S1 non-RBD']: continue
    epitope = info_dict['Protein + Epitope']
    if epitope == 'S:RBD': epitope = 'RBD'
    if epitope == 'S:NTD': epitope = 'NTD'
    if epitope == 'S:S2': epitope = 'S2'
    if epitope == 'S:S2 Stem Helix': epitope = 'S2'
    if epitope == 'S:S1 non-RBD': epitope = 'NTD'
    Dgenes  = list(set([Dgene.rsplit('*')[0] for Dgene in info_dict['Heavy D Gene'].replace('"','').rsplit(',')]))
    if len(Dgenes) > 1: continue
    Dgene_dict[epitope].append(Dgenes[0])
  for epitope in Dgene_dict.keys():
    Dgene_dict[epitope] = Counter(Dgene_dict[epitope])
  return Dgene_dict

def extract_HVgene(filename):
  print ('reading: %s' % filename)
  HVgene_dict = defaultdict(list)
  datasheet = pd.read_excel(filename, sheet_name='main')
  for index, info_dict in datasheet.iterrows():
    for i in info_dict.keys():
      info_dict[i] = str(info_dict[i])
    if info_dict['Heavy V Gene'] == 'nan' or info_dict['Heavy V Gene'] == 'NA': continue
    if info_dict['Protein + Epitope'] not in ['S:NTD','S:RBD','S:S2','S:S2 Stem Helix','S:S1 non-RBD']: continue
    epitope = info_dict['Protein + Epitope']
    if epitope == 'S:RBD': epitope = 'RBD'
    if epitope == 'S:NTD': epitope = 'NTD'
    if epitope == 'S:S2': epitope = 'S2'
    if epitope == 'S:S2 Stem Helix': epitope = 'S2'
    if epitope == 'S:S1 non-RBD': epitope = 'NTD'
    HVgene  = info_dict['Heavy V Gene'].replace('"','').rsplit('*')[0]
    if HVgene[-1] == 'D':
      HVgene=HVgene[:-2]
    HVgene_dict[epitope].append(HVgene)
  for epitope in HVgene_dict.keys():
    HVgene_dict[epitope] = Counter(HVgene_dict[epitope])
  return HVgene_dict

def extract_LVgene(filename):
  print ('reading: %s' % filename)
  LVgene_dict = defaultdict(list)
  datasheet = pd.read_excel(filename, sheet_name='main')
  for index, info_dict in datasheet.iterrows():
    for i in info_dict.keys():
      info_dict[i] = str(info_dict[i])
    if info_dict['Light V Gene'] == 'nan' or info_dict['Light V Gene'] == 'NA': continue
    if info_dict['Protein + Epitope'] not in ['S:NTD','S:RBD','S:S2','S:S2 Stem Helix','S:S1 non-RBD']: continue
    epitope = info_dict['Protein + Epitope']
    if epitope == 'S:RBD': epitope = 'RBD'
    if epitope == 'S:NTD': epitope = 'NTD'
    if epitope == 'S:S2': epitope = 'S2'
    if epitope == 'S:S2 Stem Helix': epitope = 'S2'
    if epitope == 'S:S1 non-RBD': epitope = 'NTD'
    LVgene  = info_dict['Light V Gene'].replace('"','').rsplit('*')[0]

    LVgene_dict[epitope].append(LVgene)
  for epitope in LVgene_dict.keys():
    LVgene_dict[epitope] = Counter(LVgene_dict[epitope])
  return LVgene_dict

def write_output(Dgene_dict, outfile):
  print ("writing: %s" % outfile)
  outfile = open(outfile, 'w')
  outfile.write("\t".join(['epitope', 'gene', 'freq'])+"\n")
  Dgenes = sorted(list(set([Dgene for epitope in Dgene_dict.keys() for Dgene in Dgene_dict[epitope]])))
  for epitope in Dgene_dict.keys():
    total_count = float(sum(Dgene_dict[epitope].values()))
    for Dgene in Dgenes:
      freq = float(Dgene_dict[epitope][Dgene])/total_count
      outfile.write("\t".join(map(str, [epitope, Dgene, freq]))+"\n")
  outfile.close()
    
def main():
  outfile1 = 'result/Dgene_freq.tsv'
  outfile2 = 'result/HVgene_freq.tsv'
  outfile3 = 'result/LVgene_freq.tsv'
  filename = 'data/SARS-CoV-2-Abs.xlsx'
  HVgene_dict = extract_HVgene(filename)
  Dgene_dict = extract_Dgene(filename)
  LVgene_dict = extract_LVgene(filename)
  write_output(Dgene_dict,outfile1)
  write_output(HVgene_dict, outfile2)
  write_output(LVgene_dict, outfile3)
if __name__ == "__main__":
  main()
