from Download_Genbank import id_get_seq,excel_write_seq
from Bio import Entrez
import pandas as pd

Entrez.email ='yiquan2@illinois.edu'
def main():
    df=pd.read_excel('result/SARS-CoV-2-Abs_v35.xlsx')

    df_VH_seq = df[df['VH_nuc'].notnull()]
    print(len(df_VH_seq))
    df_VL_seq = df[df['VL_nuc'].notnull()]
    print(len(df_VL_seq))


    # writing VH sequence from excel to fasta
    name_ls1 = df_VH_seq['Name'].str.replace(' ', '_').tolist()
    seq_ls1 = df_VH_seq['VH_nuc'].tolist()
    excel_write_seq(name_ls1, seq_ls1, 'result/igblast_result5/VH_nuc.fasta')
    # writing VL sequence from excel to fasta
    name_ls2 = df_VL_seq['Name'].str.replace(' ', '_').tolist()
    seq_ls2 = df_VL_seq['VL_nuc'].tolist()
    excel_write_seq(name_ls2, seq_ls2, 'result/igblast_result5/VL_nuc.fasta')


if __name__ == "__main__":
    main()