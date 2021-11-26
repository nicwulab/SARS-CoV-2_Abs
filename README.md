This README describes the analyzes in "A large-scale systematic survey of SARS-CoV-2 antibodies
reveals recurring molecular features"


## SARS-CoV-2 Antibodies


## Dependencies ##
* python=3.9
* [Igblast](https://github.com/ncbi/igblast)
* [PyIR](https://github.com/crowelab/PyIR)
* [BioPython](https://github.com/biopython/biopython)
* [Pandas](https://pandas.pydata.org/)
* [Openpyxl](https://openpyxl.readthedocs.io/en/stable/)
* [Distance](https://pypi.org/project/Distance/)
* [ANARCI](https://github.com/oxpig/ANARCI)
* [Logomaker](https://logomaker.readthedocs.io/en/latest/)

## Installation ##
Install everything dependencies by conda:

```conda create -n Abs -c bioconda -c anaconda -c conda-forge python=3.9 biopython pandas openpyxl distance logomaker igblast anarci```


## Local igblast set up

```conda activate Abs```

PyIR: An IgBLAST wrapper and parser[https://github.com/crowelab/PyIR]

```pip3 install crowelab_pyir```

Database set up in pyir library directory

```pyir setup```

- Manuualy install amino acid database from imgt

1. Sequence download from  http://www.imgt.org/vquest/refseqh.html#VQUEST
2. Copy and paste, save as fasta(save all V gene in one file)
3. Clean data (raw edit_imgt_file.pl can be found on igblast-1.17.1xxx/bin)

```edit_imgt_file.pl imgt_database/human_prot/imgt_raw/IGV.fasta > imgt_database/human_prot/IGV.fasta```

4. Create database (use "-dbtype prot" for protein sequence, use "-dbtype nucl" for DNA sequence)

```makeblastdb -parse_seqids -dbtype prot -in imgt_database/human_prot/IGV.fasta```

- Run PyIR for igBlast

  see [PyIr.py](./code/_PyIR_.py)

- Run local igblast and CDR parser

## Baseline VDJ setup

## CDR H3 clustering analysis

1. Extract information from the antibody dataset for downstream analyses   
```python3 code/parse_Ab_table.py```   
    - Input file:
      - [./data/SARS-CoV-2-Abs.xlsx](./data/SARS-CoV-2-Abs.xlsx)
    - Output files:
      - [./Fasta/SARS-CoV-2-Ab.pep](./Fasta/SARS-CoV-2-Ab.pep)
      - [./result/CDRH3.tsv](./result/CDRH3.tsv)
      - [./result/refs.txt](./result/refs.txt)

2. Clustering CDR H3 sequences   
```python3 code/CDRH3_clustering_optimal.py```   
    - Input file:
      - [./result/CDRH3.tsv](./result/CDRH3.tsv)
    - Output file:
      - [./result/CDRH3_cluster.tsv](./result/CDRH3_cluster.tsv)

3. Analyzing CDR H3 clustering results   
```python3 code/analyze_CDRH3_cluster.py```   
    - Input files:
      - [./data/SARS-CoV-2-Abs.xlsx](./data/SARS-CoV-2-Abs.xlsx)
      - [./result/CDRH3_cluster.tsv](./result/CDRH3_cluster.tsv)
    - Output file:
      - [./result/Ab_info_CDRH3_clustering.tsv](./result/Ab_info_CDRH3_clustering.tsv)
      - [./result/CDRH3_cluster_summary.tsv](./result/CDRH3_cluster_summary.tsv)

## Identification of recurring somatic hypermutation (SHM)
1. Numbering SARS2 antibody sequences according to Kabat numbering   
```ANARCI --scheme kabat --csv -i Fasta/SARS-CoV-2-Ab.pep -o result/SARS-CoV-2-Ab```   
    - Input file:
      - [./Fasta/SARS-CoV-2-Ab.pep](./Fasta/SARS-CoV-2-Ab.pep)
    - Output files:
      - [./result/SARS-CoV-2-Ab_H.csv](./result/SARS-CoV-2-Ab_H.csv)
      - [./result/SARS-CoV-2-Ab_KL.csv](./result/SARS-CoV-2-Ab_KL.csv)

2. Numbering germline sequences according to Kabat numbering   
```ANARCI --scheme kabat --csv -i imgt_database/human_prot/IGV.fasta -o result/Human_IGV_gene_kabat_num```
    - Input files:
      - [./imgt_database/human_prot/IGV.fasta](./imgt_database/human_prot/IGV.fasta)
    - Output files:
      - [./result/Human_IGV_gene_kabat_num_H.csv](./result/Human_IGV_gene_kabat_num_H.csv)
      - [./result/Human_IGV_gene_kabat_num_KL.csv](./result/Human_IGV_gene_kabat_num_KL.csv)

3. Calling SHMs   
```python3 code/SHM_analysis.py```   
    - Input files:
      - [./result/Ab_info_CDRH3_clustering.tsv](./result/Ab_info_CDRH3_clustering.tsv)
      - [./result/SARS-CoV-2-Ab_H.csv](./result/SARS-CoV-2-Ab_H.csv)
      - [./result/SARS-CoV-2-Ab_KL.csv](./result/SARS-CoV-2-Ab_KL.csv)
      - [./result/Human_IGV_gene_kabat_num_H.csv](./result/Human_IGV_gene_kabat_num_H.csv)
      - [./result/Human_IGV_gene_kabat_num_KL.csv](./result/Human_IGV_gene_kabat_num_KL.csv)
    - Output files:
      - [./result/SHM_antibody.tsv](./result/SHM_antibody.tsv)
      - [./result/SHM_frequency.tsv](./result/SHM_frequency.tsv)

## Deep learning model for antigen identification

## Plotting
1. Plot CDR H3 cluster size   
```Rscript code/plot_CDRH3_cluster_summary.R```   
    - Input file:
      - [./result/CDRH3_cluster_summary.tsv](./result/CDRH3_cluster_summary.tsv)
    - Output file:
      - [./graph/CDRH3_cluster_size.png](./graph/CDRH3_cluster_size.png)

2. Generate sequence logo for different CDR H3 clusters   
```python3 code/CDRH3_seqlogo.py```   
    - Input file:
      - [./result/Ab_info_CDRH3_clustering.tsv](./result/Ab_info_CDRH3_clustering.tsv)
    - Output files:
      - ./CDRH3_seqlogo/*.png

3. Plot SHM   
```Rscript code/plot_SHM.R```   
    - Input file:
      - [./result/SHM_frequency.tsv](./result/SHM_frequency.tsv)
    - Output files:
      - [./graph/SHM_HC_frequency.png](./graph/SHM_HC_frequency.png)
      - [./graph/SHM_LC_frequency.png](./graph/SHM_LC_frequency.png)
