#!/usr/bin/python
import os
import sys
from collections import defaultdict

def extract_germline_info(filename):
  Ab_info = defaultdict(dict)
  print ('reading: %s' % filename)
  infile  = open(filename, 'r')
  line_count = 0
  for line in infile.readlines():
    line_count += 1
    if line_count == 1:
      header = line.rstrip().rsplit("\t")
    else:
      info_dict = defaultdict(str)
      for field, entry in zip(header, line.rstrip().rsplit("\t")):
        info_dict[field] = entry
      ID = info_dict['Patient ID']+'|'+info_dict['Name']
      Ab_info[ID] = info_dict 
  infile.close()
  return Ab_info

def reading_kabat_numbering(file_HV_num, Ab_dict):
  print ("reading: %s" % file_HV_num)
  infile = open(file_HV_num,'r')
  line_count = 0
  for line in infile:
    line_count += 1
    if line_count == 1:
      header = line.rstrip().rsplit(",")
    else:
      info_dict = defaultdict(str)
      for field, entry in zip(header, line.rstrip().rsplit(",")):
        info_dict[field] = entry
      Ab_dict[info_dict['Id']] = info_dict
      #if info_dict['52F'] != '-' and info_dict['52F'] != '': print (info_dict['Id'])
  infile.close()
  return (Ab_dict)

def position_to_integer(pos):
  try:
    return (int(pos))
  except:
    return (int(pos[0:-1]))

def compare_ref_and_Ab(SHM_dict, ref_kabat, Ab_kabat, clonotype, donor, chain):
  SHM_list = []
  for pos in ref_kabat.keys():
    if '-' not in pos and '_' not in pos and 'score' != pos and 'Id' != pos:
      if chain == 'HC' and position_to_integer(pos) >= 92: continue
      if chain == 'LC' and position_to_integer(pos) >= 87: continue
      #if position_to_integer(pos) <= 10: continue
      Ab_aa  = Ab_kabat[pos]
      ref_aa = ref_kabat[pos]
      pos = pos.lower()
      if Ab_aa == '':  Ab_aa = '-'
      if ref_aa == '': ref_aa = '-'
      if ref_aa == '-': continue
      if Ab_aa == 'X': continue
      if Ab_aa != ref_aa and Ab_aa != '-':
        mut = ref_aa+pos+Ab_aa
        SHM_dict[clonotype][chain][mut].append(donor)
        SHM_list.append(mut)
        #if clonotype == 'IGHV4-31_IGKV1-33_26':
        #  print (donor, mut)
  return SHM_dict, SHM_list



def call_SHM(Ab_info, Ab_sars2_dict, Ab_ref_dict, shmfile):
  print ('writing: %s' % shmfile)
  outfile = open(shmfile, 'w')
  outfile.write("\t".join(['ID', 'clonotype', 'chain', 'SHM'])+"\n")
  depth_dict = defaultdict(list)
  SHM_dict   = {}
  for Ab in Ab_info.keys():
    HV = Ab_info[Ab]['Heavy V Gene'].rsplit(',')[0]
    LV = Ab_info[Ab]['Light V Gene'].rsplit(',')[0]
    if HV == '' or LV == '': continue
    CDRH3_cluster = Ab_info[Ab]['cluster']
    clonotype = HV.rsplit('*')[0]+'_'+LV.rsplit('*')[0]+'_c'+CDRH3_cluster
    donor = Ab.rsplit('|')[0]
    depth_dict[clonotype].append(donor)
    if clonotype not in SHM_dict.keys():
      SHM_dict[clonotype]   = {'HC':defaultdict(list), 'LC':defaultdict(list)}
    if Ab+'__HC' in Ab_sars2_dict.keys():
      Ab_kabat  = Ab_sars2_dict[Ab+'__HC']
      ref_kabat = Ab_ref_dict[HV]
      if '*' in HV:
        SHM_dict, SHM_list = compare_ref_and_Ab(SHM_dict, ref_kabat, Ab_kabat, clonotype, donor, 'HC')
        outfile.write("\t".join([Ab, clonotype, 'HC', ','.join(SHM_list)])+"\n")
    if Ab+'__LC' in Ab_sars2_dict.keys():
      Ab_kabat  = Ab_sars2_dict[Ab+'__LC']
      ref_kabat = Ab_ref_dict[LV]
      if '*' in LV: 
        SHM_dict, SHM_list = compare_ref_and_Ab(SHM_dict, ref_kabat, Ab_kabat, clonotype, donor, 'LC')
        outfile.write("\t".join([Ab, clonotype, 'LC', ','.join(SHM_list)])+"\n")
  outfile.close()
  return depth_dict, SHM_dict

def writing_SHM(depth_dict, SHM_dict, outfile, pos_dict):
  print ('writing: %s' % outfile)
  outfile = open(outfile, 'w')
  outfile.write("\t".join(['clonotype', 'chain', 'label', 'adj_pos', 'SHM', 'total_donor', 'occurrence', 'freq'])+"\n") 
  for clonotype in depth_dict.keys():
    total_donor = len(set(depth_dict[clonotype]))
    for chain in ['HC', 'LC']:
      for SHM in SHM_dict[clonotype][chain].keys():
        occurrence = len(set(SHM_dict[clonotype][chain][SHM])) 
        freq  = float(occurrence)/float(total_donor)
        pos   = pos_dict[chain].index(SHM[1:-1])
        label = SHM+' ('+'/'.join(clonotype.rsplit('_')[0:3])+')'
        outfile.write("\t".join(map(str, [clonotype, chain, label, pos, SHM, total_donor, occurrence, freq]))+"\n") 
  outfile.close()

def replace_Vgene(V):
  if V == 'IGHV3-30-3':
    return 'IGHV3-30'
  else:
    return V

def main():
  outfile = 'result/SHM_frequency.tsv'
  shmfile = 'result/SHM_antibody.tsv'
  Ab_info = extract_germline_info('result/Ab_info_CDRH3_clustering.tsv')
  Ab_ref_dict = defaultdict(dict)
  Ab_ref_dict = reading_kabat_numbering('result/Human_IGV_gene_kabat_num_H.csv', Ab_ref_dict)
  Ab_ref_dict = reading_kabat_numbering('result/Human_IGV_gene_kabat_num_KL.csv', Ab_ref_dict)
  HC_positions   = [pos.lower() for pos in Ab_ref_dict['IGHV3-53*01'].keys() if '-' not in pos and '_' not in pos
                    and 'score' != pos and 'Id' != pos]
  HC_positions   = list(map(replace_Vgene, HC_positions))
  LC_positions   = [pos.lower() for pos in Ab_ref_dict['IGKV3-20*01'].keys() if '-' not in pos and '_' not in pos
                    and 'score' != pos and 'Id' != pos]
  print ([HC_positions.index(n) for n in map(str,[10,20,30,40,50,60,70,80,90,100])])
  print ([LC_positions.index(n) for n in map(str,[10,20,30,40,50,60,70,80])])
  pos_dict      = {'HC':HC_positions,'LC':LC_positions}
  Ab_sars2_dict = defaultdict(dict)
  Ab_sars2_dict = reading_kabat_numbering('result/SARS-CoV-2-Ab_H.csv', Ab_sars2_dict)
  Ab_sars2_dict = reading_kabat_numbering('result/SARS-CoV-2-Ab_KL.csv', Ab_sars2_dict)
  depth_dict, SHM_dict = call_SHM(Ab_info, Ab_sars2_dict, Ab_ref_dict, shmfile)
  writing_SHM(depth_dict, SHM_dict, outfile, pos_dict)

if __name__ == "__main__":
  main()
