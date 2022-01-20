#!/usr/bin/python3
from Bio import SeqIO

def extract_aa(filename, outfile):
  print ("writing: %s" % outfile)
  outfile = open(outfile, 'w')
  outfile.write('ID'+"\t"+"mut"+"\n")
  for record in SeqIO.parse(filename, "fasta"):
    Ab_ID = str(record.id)
    resi29 = str(record.seq)[29]
    resi32 = str(record.seq)[32]
    resi92 = str(record.seq)[92]
    variant = 'Germline' if 'germline'in Ab_ID else resi29+"29/"+resi92+'92'
    outfile.write(Ab_ID+"\t"+variant+"\n")
  outfile.close()

def main():
  filename  = 'Fasta/Cluster3_H158K320_LC.aln'
  outfile = 'result/Cluster3_H158K320_LC_class.tsv'
  extract_aa(filename, outfile)
  

if __name__ == "__main__":
  main()
