from Bio import Entrez
import pandas as pd
Entrez.email ='yiquan2@illinois.edu'

def get_seq(id_ls, name_ls,output):
    ids = name_ls
    with open (output,'w') as f:
        for idx, ID in enumerate(id_ls):
            handle = Entrez.efetch(db="nucleotide", id=ID, rettype="fasta", retmode="text")
            record = handle.read()
            record=record[0]+ids[idx]+'+'+record[1:]
            print('writing %s' % ids[idx])
            f.write(record)
def write_seq(name_ls, seq_ls,output):
    with open(output, 'w') as f:
        for idx, seq in enumerate(seq_ls):
            name = '>'+name_ls[idx]+'+'
            f.write(name+'\n'+seq+'\n')
def main():
    output = 'result/Genbank_ID.fasta'
    df=pd.read_excel('data/CoVmAb_info_to_extract_v1.xlsx')
    df_id = df[df['nuc seq / GenBank ID'].str.len()<10]
    df_seq1 = df[df['nuc seq / GenBank ID'].str.len()>10]
    df_seq2 = df[df['additional nuc seq / GenBank ID'].str.len() > 10]
    # writing sequence1 to fasta
    name_ls = df_seq1['Ab name'].tolist()
    seq_ls = df_seq1['nuc seq / GenBank ID'].tolist()
    write_seq(name_ls,seq_ls,'result/sequence1.fasta')
    # writing sequence2 to fasta
    name_ls2 = df_seq2['Ab name'].tolist()
    seq_ls2 = df_seq2['additional nuc seq / GenBank ID'].tolist()
    write_seq(name_ls2, seq_ls2, 'result/sequence2.fasta')
    # writing genbank id sequence
    name_ls=df_id['Ab name'].tolist()
    id_ls=df_id['nuc seq / GenBank ID'].tolist()
    get_seq(id_ls,name_ls,output)

if __name__ == "__main__":
    main()