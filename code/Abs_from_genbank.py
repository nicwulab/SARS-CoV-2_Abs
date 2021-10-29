from Bio import Entrez,SeqIO
import pandas as pd

Entrez.email ='yiquan2@illinois.edu'
def Dload_GB(seq_path,output):
    # downolad genbank file from genbank ID
    # seq='data/sequence.seq'
    with open(seq_path, 'r') as f:
        IDs = f.read().replace('\n', ',').split(',')
        IDs = ', '.join(set(IDs))
        with Entrez.efetch(
                db="nucleotide", rettype="gb", retmode="text", id=IDs
        ) as handle:
            with open(output, 'w') as o:
                o.write('Genbank ID' + "\t" + 'sequence' + "\t" + "PMID" + "\t" + "Description" "\t" + "References" + "\n")
                for record in SeqIO.parse(handle, "gb"):
                    #if record.annotations['molecule_type'] == 'DNA': continue
                    ID = record.id
                    print('processing %s' % ID)
                    Des = record.description
                    PMID = record.annotations['references'][0].pubmed_id
                    Ref = record.annotations['references']
                    Seq = record.seq
                    o.write(str(ID) + "\t" + str(Seq) + "\t" + str(PMID) + "\t" + Des + "\t" + str(Ref) + "\n")
def extract_gb(gb_path,output):
    records = SeqIO.parse(gb_path, "genbank")
    with open(output,'w') as o:
        o.write('Genbank ID' + "\t" + 'sequence' + "\t" + "PMID" + "\t" + "Description" "\t" + "References" + "\n")
        for record in records:
            if record.annotations['molecule_type'] == 'DNA': continue
            ID = record.id
            print('processing %s' % ID)
            Des = record.description
            PMID = record.annotations['references'][0].pubmed_id
            Ref = record.annotations['references']
            Seq = record.seq
            o.write(str(ID) + "\t" + str(Seq) + "\t" + str(PMID) + "\t" + Des + "\t" + str(Ref) + "\n")

def main():
    seq_file = 'data/ab_GBIDv20.seq'
    output = 'result/Abs_from_GBID.tsv'
    # gb_file = 'data/sequence_v2.gb'
    Dload_GB(seq_file,output)


if __name__ == "__main__":
    main()