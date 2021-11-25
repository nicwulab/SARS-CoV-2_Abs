from collections import defaultdict
from Bio import SeqIO
import math
""" 
    ########---NOTICE---########
    
    Before running this parser, first run igblat locally.
    igblastn -query result/igblast_result7/test.fasta 
    -germline_db_V imgt_database/human_nuc/IGV.fasta 
    -germline_db_J imgt_database/human_nuc/IGJ.fasta 
    -germline_db_D imgt_database/human_nuc/IGD.fasta 
    -organism human -domain_system kabat 
    -auxiliary_data imgt_database/optional_file/human_gl.aux 
    -out result/igblast_result7/output_file
    
    ########------------########
"""

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    if 'N' in codon:
        aa = 'X'
    else:
        aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep
def CDR_pos_parser(input):
    """ input is the igblast output file """
    infile= open(input,'r')
    CDR1_pos=defaultdict(list)
    CDR2_pos=defaultdict(list)
    name_flag=False
    CDR1_fromto=[]
    CDR2_fromto=[]
    for line in infile.readlines():
        # skip empty rows
        if len(line) < 2: continue
        if line[:5] == "Query":
            name = line.strip()[7:]
            name_flag = True
        elif line[:5] == 'CDR1' + '\t':
            CDR1_fromto=line[5:].split('\t')

        elif line[:5] == 'CDR2' + '\t':
            CDR2_fromto = line[5:].split('\t')

        if (line[:10] == 'Alignments') & name_flag:
            CDR1_pos[name]=CDR1_fromto
            CDR2_pos[name]=CDR2_fromto
            name_flag = False
            CDR1_fromto = []
            CDR2_fromto = []
    return CDR1_pos,CDR2_pos
def fix_aa_pos(start,end):
    diff = end-start
    aa_num = math.ceil(diff/3)
    new_start = end - 3*aa_num
    new_end = end
    return new_start,new_end
def write_output(output,fasta,CDR1_pos,CDR2_pos):
    """write CDR kabat output file as tsv format"""
    # fasta file is the igblast input file
    with open(output, 'w') as f:
        header = "\t".join(['Name', 'CDRL1_kabat_AA', 'CDRL2_kabat_AA'])
        f.write(header + '\n')
        for record in SeqIO.parse(fasta, "fasta"):
            ID = str(record.id)
            seq = str(record.seq)
            CDR1_aa=''
            CDR2_aa = ''
            CDR1_index = CDR1_pos[ID]
            CDR2_index = CDR2_pos[ID]
            if CDR1_index != []:
                CDR1_start, CDR1_end = fix_aa_pos((int(CDR1_index[0]) - 1), int(CDR1_index[1]))
                CDR1_nuc = seq[CDR1_start:CDR1_end]
                CDR1_aa = translation(CDR1_nuc)
            if CDR2_index != []:
                CDR2_start, CDR2_end = fix_aa_pos((int(CDR2_index[0]) - 1), int(CDR2_index[1]))
                CDR2_nuc = seq[CDR2_start:CDR2_end]
                CDR2_aa = translation(CDR2_nuc)
            f.write("\t".join([ID, CDR1_aa, CDR2_aa]) + '\n')

def main():
    input='result/igblast_result7/VL_igblast_output'
    CDR1_pos,CDR2_pos = CDR_pos_parser(input)
    fasta='result/igblast_result7/VL_nuc.fasta'
    output='result/igblast_result7/VL_kabat_CDR12.tsv'
    write_output(output,fasta,CDR1_pos,CDR2_pos)
if __name__ == "__main__":
    main()




