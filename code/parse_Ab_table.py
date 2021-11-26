#!/usr/bin/python
import os
import sys
import glob
import pandas as pd
from Bio.Seq import Seq
from collections import defaultdict, Counter

def extract_basic_info(filename, outfile_CDRH3, outfile_pep, outfile_ref):
  aas = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
  print ('reading: %s' % filename)
  print ('writing: %s' % outfile_CDRH3)
  print ('writing: %s' % outfile_pep)
  print ('writing: %s' % outfile_ref)
  print ('writing: %s' % out_tsv)
  outfile1 = open(outfile_CDRH3, 'w')
  outfile2 = open(outfile_pep, 'w')
  datasheet = pd.read_excel(filename, sheet_name='main')
  datasheet.to_csv(out_tsv, sep="\t", index=False)
  refs = []
  donors = []
  epi    = []
  count_SARS2_Ab = 0
  for index, info_dict in datasheet.iterrows():
    for i in info_dict.keys():
      info_dict[i] = str(info_dict[i])
    ref = info_dict['Sources']
    refs += ref.rsplit(';')
    donors.append(info_dict['Patient ID'])
    if 'SARS-CoV-2' in str(info_dict['Binds to']) or 'SARS-CoV-2' in str(info_dict['Neutralising Vs']):
      count_SARS2_Ab += 1
      epi.append(str(info_dict['Protein + Epitope'])) 
    if str(info_dict['Patient ID']) == 'nan': continue
    if str(info_dict['CDRH3']) == 'nan': continue
    ID     = str(info_dict['Patient ID'])+'|'+str(info_dict['Name'])
    CDRH3  = str(info_dict['CDRH3'])
    CDRL3  = str(info_dict['CDRL3'])
    if 'PDB' in ref: print (ID)
    if '*' in CDRH3 or 'NA' == CDRH3: continue
    outfile1.write(ID+"\t"+CDRH3+"\n")
    HC_pep = str(info_dict['VH'])
    LC_pep = str(info_dict['VL'])
    if '*' in HC_pep: continue
    if '*' in LC_pep: continue
    if HC_pep != 'nan' and HC_pep != '':
      if CDRH3 != '' and CDRH3 != 'nan' and CDRH3 in HC_pep:
        HC_pep = HC_pep.replace('-','').rsplit(CDRH3)[0]
        outfile2.write('>'+ID+'__HC'+"\n"+HC_pep+"\n")
    if LC_pep != 'nan' and LC_pep != '':
      if CDRL3 != '' and CDRL3 != 'nan' and CDRL3 in LC_pep: 
        LC_pep = LC_pep.replace('-','').rsplit(CDRL3)[0]
        outfile2.write('>'+ID+'__LC'+"\n"+LC_pep+"\n")
  outfile1.close()
  outfile2.close()
  outfile3 = open(outfile_ref, 'w')
  outfile3.write("\n".join(sorted(list(set((refs)))))+"\n")
  outfile3.close()
  print ('# of Research Papers: %i' % len(set([ref for ref in refs if 'Patent' not in ref])))
  print ('# of Patents: %i' % len(set([ref for ref in refs if 'Patent' in ref])))
  print ('# of total donors: %i' % len(set(donors)))
  print ('# of SARS2 Abs: %i' % count_SARS2_Ab)
  print ('  ',Counter(epi))

def main():
  filename = 'data/SARS-CoV-2-Abs.xlsx'
  outfile_pep    = 'Fasta/SARS-CoV-2-Ab.pep'
  outfile_CDRH3  = 'result/CDRH3.tsv'
  outfile_ref    = 'result/refs.txt'
  extract_basic_info(filename, outfile_CDRH3, outfile_pep, outfile_ref)
  

if __name__ == "__main__":
  main()
