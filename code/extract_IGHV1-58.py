#!/usr/bin/python
import os
import sys
import glob
import pandas as pd
from Bio.Seq import Seq
from collections import defaultdict

def extract_info(filename, outfile, classfile, KV320_ref):
  print ("writing: %s" % outfile)
  print ("writing: %s" % classfile)
  outfile   = open(outfile, 'w')
  classfile = open(classfile, 'w')
  datasheet = pd.read_csv(filename, sep="\t")
  outfile.write('>IGKV3-20_germline'+"\n"+KV320_ref+"\n")
  classfile.write('ID'+"\t"+"mut"+"\n")
  classfile.write('IGKV3-20_germline'+"\t"+"Germline"+"\n")
  R29_D92 = []
  R29 = []
  for index, info_dict in datasheet.iterrows():
    for i in info_dict.keys():
      info_dict[i] = str(info_dict[i])
    if str(info_dict['Patient ID']) == 'nan': continue
    if str(info_dict['CDRH3']) == 'nan': continue
    ID      = str(info_dict['Patient ID'])+'|'+str(info_dict['Name'])
    cluster = str(info_dict['cluster'])
    LC_pep = str(info_dict['VL'])
    HC_V    = str(info_dict['Heavy V Gene']).rsplit('*')[0]
    LC_V    = str(info_dict['Light V Gene']).rsplit('*')[0]
    CDRH3   = str(info_dict['CDRH3'])
    CDRL3   = str(info_dict['CDRL3'])
    if '*' in CDRH3 or 'NA' == CDRH3: continue
    if '*' in LC_pep: continue
    if cluster=='3' and HC_V == 'IGHV1-58' and LC_V == 'IGKV3-20' and LC_pep != 'nan':
      if ID == 'CV07|HK_CV07-287': 
        seq = str(LC_pep)[0:94]
      else:
        seq = str(LC_pep)[0:95]
      outfile.write('>'+ID+"\n"+seq+"\n")
      resi29 = str(LC_pep)[29]
      resi32 = str(LC_pep)[32]
      resi92 = str(LC_pep)[92]
      if resi29 == 'R' and resi92 == 'D':
        classfile.write(ID+"\t"+"R29/D92"+"\n")
      elif resi29 == 'R':
        classfile.write(ID+"\t"+"R29/G92"+"\n")
      else:
        classfile.write(ID+"\t"+"S92/G92"+"\n")
  outfile.close()
  classfile.close()

def main():
  outfile   = 'Fasta/Cluster3_H158K320_LC.pep'
  classfile = 'result/Cluster3_H158K320_LC_class.tsv'
  filename  = 'result/Ab_info_CDRH3_clustering.tsv'
  KV320_ref = 'EIVLTQSPGTLSLSPGERATLSCRASQSVSSSYLAWYQQKPGQAPRLLIYGASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQYGSS'
  extract_info(filename, outfile, classfile, KV320_ref)

if __name__ == "__main__":
  main()
