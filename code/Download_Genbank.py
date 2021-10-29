from Bio import Entrez
import pandas as pd
Entrez.email ='yiquan2@illinois.edu'

def id_get_seq(id_ls, name_ls, output):
    ids = name_ls
    with open (output,'w') as f:
        for idx, ID in enumerate(id_ls):
            ID=ID.strip()
            handle = Entrez.efetch(db="nucleotide", id=ID, rettype="fasta", retmode="text")
            record = handle.read()
            record=record[0]+ids[idx]+'+'+record[1:]
            print('writing %s' % ids[idx])
            f.write(record)
def excel_write_seq(name_ls, seq_ls, output):
    with open(output, 'w') as f:
        for idx, seq in enumerate(seq_ls):
            seq=str(seq)
            seq=seq.strip()
            name = '>'+name_ls[idx]
            f.write(name+'\n'+seq+'\n')
def main():
    df=pd.read_excel('data/CoVmAb_info_to_extract_v2.xlsx')
    df_id1 = df[df['nuc seq / GenBank ID'].str.len()<11]
    df_id2 = df[df['additional nuc seq / GenBank ID'].str.len() < 11]
    df_seq1 = df[df['nuc seq / GenBank ID'].str.len()>11]
    df_seq2 = df[df['additional nuc seq / GenBank ID'].str.len() > 11]

    # writing VH sequence1 from excel to fasta
    name_ls1 = df_seq1['Ab name'].tolist()
    seq_ls1 = df_seq1['nuc seq / GenBank ID'].tolist()
    excel_write_seq(name_ls1, seq_ls1, 'result/sequence1.fasta')
    # writing VL sequence2 from excel to fasta
    name_ls2 = df_seq2['Ab name'].tolist()
    seq_ls2 = df_seq2['additional nuc seq / GenBank ID'].tolist()
    excel_write_seq(name_ls2, seq_ls2, 'result/sequence2.fasta')
    # writing VH genbank id sequence
    name_ls3=df_id1['Ab name'].tolist()
    id_ls3=df_id1['nuc seq / GenBank ID'].tolist()
    id_get_seq(id_ls3, name_ls3, 'result/sequence3.fasta')
    # writing VL genbank id sequence
    name_ls4 = df_id2['Ab name'].tolist()
    id_ls4 = df_id2['additional nuc seq / GenBank ID'].tolist()
    id_get_seq(id_ls4, name_ls4, 'result/sequence4.fasta')

if __name__ == "__main__":
    main()