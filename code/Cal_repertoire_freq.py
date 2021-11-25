import pandas as pd
from collections import defaultdict,Counter

donor_df=pd.read_excel('data/baseline_database/cAb_donor_info.xlsx')
donor_dic=pd.Series(donor_df.donor.values,index=donor_df.SRA_id).to_dict()

def cal_repertoire_freq(input,outfile,donor_dict,gene):
    infile=open(input,'r')
    Ab_database_dict=defaultdict(list)
    df = pd.DataFrame(columns=[gene])
    report_freq=100000
    for iter, line in enumerate(infile.readlines()):
        Gene=''
        if gene == 'V_gene':
            Gene= line.rsplit(' ')[1].rsplit('=')[1].rsplit('*')[0]
        elif gene == 'D_gene':
            Gene= line.rsplit(' ')[3].rsplit('=')[1].rsplit('*')[0]
        if Gene == 'not_found':continue
        if Gene == 'N/A':continue
        if Gene == 'NA':continue
        SRA_ID= line.rsplit(' ')[0].rsplit('_')[0].rsplit('>')[1]
        donor=donor_dict[SRA_ID]
        Ab_database_dict[donor].append(Gene)
        if iter % report_freq == 0:
            print('Finish processing: {}'.format(float(iter)))
    print('writing file:{}'.format(outfile))
    for donor,gene_ls in Ab_database_dict.items():
        gene_dict=Counter(gene_ls)
        count_col=str(donor)+'_Counts'
        V_df=pd.DataFrame.from_dict(gene_dict, orient='index',columns=[count_col])
        V_df[gene] = V_df.index
        freq_col=str(donor)+'_Freq'
        V_df[freq_col]=V_df[count_col]/V_df[count_col].sum()
        df=df.merge(V_df,on=gene,how='outer')
    df.to_excel(outfile)


cal_repertoire_freq('data/baseline_database/Abs_info/Nature_Briney_name.txt','data/baseline_database/Briney_HV_Repertoire_freq.xlsx',donor_dic,'V_gene')
cal_repertoire_freq('data/baseline_database/Abs_info/Nature_Briney_name.txt','data/baseline_database/Briney_D_Repertoire_freq.xlsx',donor_dic,'D_gene')
cal_repertoire_freq('data/baseline_database/Abs_info/Cinque_heavy_name.txt','data/baseline_database/Cinque_HV_Repertoire_freq.xlsx',donor_dic,'V_gene')
cal_repertoire_freq('data/baseline_database/Abs_info/Cinque_heavy_name.txt','data/baseline_database/Cinque_D_Repertoire_freq.xlsx',donor_dic,'D_gene')

