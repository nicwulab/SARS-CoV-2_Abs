# Sequence analysis of SARS-CoV-2 antibodies
This README describes the analyses in:   
[A large-scale systematic survey reveals recurring molecular features of public antibody responses to SARS-CoV-2](https://www.cell.com/immunity/fulltext/S1074-7613(22)00142-X)

## Contents

* [Local igblast setup](#local-igblast-setup)   
* [Analysis of VDJ gene usage and V gene pairing](#Analysis-of-VDJ-gene-usage-and-V-gene-pairing)   
* [Baseline VDJ setup](#baseline-vdj-setup)   
* [CDR H3 clustering analysis](#cdr-h3-clustering-analysis)    
* [Identification of recurring somatic hypermutation (SHM)](#identification-of-recurring-somatic-hypermutation-shm)   
* [Analysis of recurring SHM in IGHV1-58/IGKV3-20 antibodies](#analysis-of-recurring-shm-in-ighv1-58igkv3-20-antibodies)
* [Clonotype assignment](#Clonotype-assignment)
* [Deep learning model for antigen identification](#deep-learning-model-for-antigen-identification)   
* [Plotting](#plotting)  

## Input files 

* [./data/SARS-CoV-2-Abs.xlsx](./data/SARS-CoV-2-Abs.xlsx): List of antibodies collected from publications and patents
* [./data/HV_Repertoire_freq.xlsx](./data/HV_Repertoire_freq.xlsx): Baseline IGHV usage (see [Baseline VDJ setup](#Baseline-VDJ-setup))
* [./data/LV_Repertoire_freq.xlsx](./data/LV_Repertoire_freq.xlsx): Baseline IGK(L)V usage (see [Baseline VDJ setup](#Baseline-VDJ-setup))
* [./data/D_Repertoire_freq.xlsx](./data/D_Repertoire_freq.xlsx): Baseline IGHD usage (see [Baseline VDJ setup](#Baseline-VDJ-setup))

## Dependencies ##
* [python3](https://www.python.org/downloads/)
* [Igblast](https://github.com/ncbi/igblast)
* [MAFFT](https://mafft.cbrc.jp/alignment/software/)
* [FastTree](http://www.microbesonline.org/fasttree/)
* [PyIR](https://github.com/crowelab/PyIR)
* [BioPython](https://github.com/biopython/biopython)
* [Pandas](https://pandas.pydata.org/)
* [Openpyxl](https://openpyxl.readthedocs.io/en/stable/)
* [Distance](https://pypi.org/project/Distance/)
* [ANARCI](https://github.com/oxpig/ANARCI)
* [Logomaker](https://logomaker.readthedocs.io/en/latest/)
* [R](https://www.r-project.org/)

## Dependencies Installation ##
Install dependencies by conda:

```
conda create -n Abs -c bioconda -c anaconda -c conda-forge \
  python=3.9 \
  biopython \
  pandas \
  openpyxl \
  distance \
  logomaker \
  igblast \
  anarci \
  mafft \
  fasttree
```

## Local igblast setup

Before analysis, do:

``conda activate Abs``

### PyIR: An IgBLAST wrapper and parser

``pip3 install crowelab_pyir``

Database set up in pyir library directory

``pyir setup``

### Manually install IMGT REF database

1. Sequence download from  http://www.imgt.org/vquest/refseqh.html#VQUEST

2. Copy and paste, save as fasta (save all V gene in one file; all D gene in one file; all J gene in one file)

3. Clean data (raw edit_imgt_file.pl can be found on igblast-1.17.1xxx/bin)

``edit_imgt_file.pl imgt_database/human_prot/imgt_raw/IGV.fasta > imgt_database/human_prot/IGV.fasta``

4. Create database (use "-dbtype prot" for protein sequence, use "-dbtype nucl" for DNA sequence). For example:

``makeblastdb -parse_seqids -dbtype prot -in imgt_database/human_prot/IGV.fasta``   

``makeblastdb -parse_seqids -dbtype nucl -in imgt_database/human_nuc/IGV.fasta``

5. Run PyIR for igBlast

- see [PyIR.py](./code/_PyIR_.py)

### Run local igblast and CDR parser

1. Run local igblast on kabat numbering system

```
igblastn -query result/test.fasta \
  -germline_db_V imgt_database/human_nuc/IGV.fasta \
  -germline_db_J imgt_database/human_nuc/IGJ.fasta \
  -germline_db_D imgt_database/human_nuc/IGD.fasta \
  -organism human -domain_system kabat \
  -auxiliary_data imgt_database/optional_file/human_gl.aux \
  -out result/igblast_output
```

2. Parse igblast output

- Using [CDR_Parser.py](./code/CDR_Parser.py) for igblast_output

## Baseline VDJ setup

1. Download antibody repertoire data for healthy donors from [cAb-Rep](https://www.frontiersin.org/articles/10.3389/fimmu.2019.02365/full)

2. [Cal_repertoire_freq.py](./code/Cal_repertoire_freq.py) is used to establish the baseline germline usage frequency

## Analysis of VDJ gene usage and V gene pairing
1. Extract VDJ gene usage   
``python3 code/VDJgene_freq_analysis.py``   
    - Input file:
      - [./data/SARS-CoV-2-Abs.xlsx](./data/SARS-CoV-2-Abs.xlsx)
    - Output files:
      - [./result/Dgene_freq.tsv](./result/Dgene_freq.tsv)
      - [./result/HVgene_freq.tsv](./result/HVgene_freq.tsv)
      - [./result/LVgene_freq.tsv](./result/LVgene_freq.tsv)

## CDR H3 clustering analysis

1. Extract information from the antibody dataset for downstream analyses   
``python3 code/parse_Ab_table.py``   
    - Input file:
      - [./data/SARS-CoV-2-Abs.xlsx](./data/SARS-CoV-2-Abs.xlsx)
    - Output files:
      - [./Fasta/SARS-CoV-2-Ab.pep](./Fasta/SARS-CoV-2-Ab.pep)
      - [./result/CDRH3.tsv](./result/CDRH3.tsv)
      - [./result/refs.txt](./result/refs.txt)

2. Clustering CDR H3 sequences   
``python3 code/CDRH3_clustering_optimal.py``   
    - Input file:
      - [./result/CDRH3.tsv](./result/CDRH3.tsv)
    - Output file:
      - [./result/CDRH3_cluster.tsv](./result/CDRH3_cluster.tsv)

3. Analyzing CDR H3 clustering results   
``python3 code/analyze_CDRH3_cluster.py``   
    - Input files:
      - [./data/SARS-CoV-2-Abs.xlsx](./data/SARS-CoV-2-Abs.xlsx)
      - [./result/CDRH3_cluster.tsv](./result/CDRH3_cluster.tsv)
    - Output file:
      - [./result/Ab_info_CDRH3_clustering.tsv](./result/Ab_info_CDRH3_clustering.tsv)
      - [./result/CDRH3_cluster_summary.tsv](./result/CDRH3_cluster_summary.tsv)

## Identification of recurring somatic hypermutation (SHM)
1. Numbering SARS2 antibody sequences according to Kabat numbering   
``ANARCI --scheme kabat --csv -i Fasta/SARS-CoV-2-Ab.pep -o result/SARS-CoV-2-Ab``   
    - Input file:
      - [./Fasta/SARS-CoV-2-Ab.pep](./Fasta/SARS-CoV-2-Ab.pep)
    - Output files:
      - [./result/SARS-CoV-2-Ab_H.csv](./result/SARS-CoV-2-Ab_H.csv)
      - [./result/SARS-CoV-2-Ab_KL.csv](./result/SARS-CoV-2-Ab_KL.csv)

2. Numbering germline sequences according to Kabat numbering   
``ANARCI --scheme kabat --csv -i imgt_database/human_prot/IGV.fasta -o result/Human_IGV_gene_kabat_num``   
    - Input files:
      - [./imgt_database/human_prot/IGV.fasta](./imgt_database/human_prot/IGV.fasta)
    - Output files:
      - [./result/Human_IGV_gene_kabat_num_H.csv](./result/Human_IGV_gene_kabat_num_H.csv)
      - [./result/Human_IGV_gene_kabat_num_KL.csv](./result/Human_IGV_gene_kabat_num_KL.csv)

3. Calling SHMs   
``python3 code/SHM_analysis.py``   
    - Input files:
      - [./result/Ab_info_CDRH3_clustering.tsv](./result/Ab_info_CDRH3_clustering.tsv)
      - [./result/SARS-CoV-2-Ab_H.csv](./result/SARS-CoV-2-Ab_H.csv)
      - [./result/SARS-CoV-2-Ab_KL.csv](./result/SARS-CoV-2-Ab_KL.csv)
      - [./result/Human_IGV_gene_kabat_num_H.csv](./result/Human_IGV_gene_kabat_num_H.csv)
      - [./result/Human_IGV_gene_kabat_num_KL.csv](./result/Human_IGV_gene_kabat_num_KL.csv)
    - Output files:
      - [./result/SHM_antibody.tsv](./result/SHM_antibody.tsv)
      - [./result/SHM_frequency.tsv](./result/SHM_frequency.tsv)

## Analysis of recurring SHM in IGHV1-58/IGKV3-20 antibodies
1. Extracting light chain sequences from IGHV1-58/IGKV3-20 antibodies
``python3 code/extract_IGHV1-58_VL.py``
    - Input file:
      - [./result/Ab_info_CDRH3_clustering.tsv](./result/Ab_info_CDRH3_clustering.tsv)
    - Output file:
      - [./Fasta/Cluster3_H158K320_LC.pep](./Fasta/Cluster3_H158K320_LC.pep)

2. Multiple sequence alignment
``mafft Fasta/Cluster3_H158K320_LC.pep > Fasta/Cluster3_H158K320_LC.aln``

3. Identifying amino acid variants at VL residues 29 and 92
``code/extract_IGHV1-58_aa.py``
    - Input file:
      - [./Fasta/Cluster3_H158K320_LC.aln](./Fasta/Cluster3_H158K320_LC.aln)
    - Output file:
      - [./result/Cluster3_H158K320_LC_class.tsv](./result/Cluster3_H158K320_LC_class.tsv)
4. Constructing phylogenetic tree
``FastTree Fasta/Cluster3_H158K320_LC.aln > result/Cluster3_H158K320_LC.tree``

## Clonotype assignment
1. Antibodies with the same IGHV, IGK(L)V, IGHJ, IGK(L)J, and belong to the same CDR H3 cluster will be assigned to the same clonotype
    - Input files:
      - [./data/SARS-CoV-2-Abs.xlsx](./data/SARS-CoV-2-Abs.xlsx)
      - [./result/CDRH3_cluster.tsv](./result/CDRH3_cluster.tsv)
    - Output file:
      - [./result/clonotype_assign.tsv](./result/clonotype_assign.tsv)

## Deep learning model for antigen identification

Deep learning model is under [CoV_Encoder](./code/CoV_Encoder)

* Classifier for S vs HA: [./code/CoV_Encoder/S_HA_classifier.ipynb](code/CoV_Encoder/S_HA_classifier.ipynb)
* Classifier for RBD vs HA: [./code/CoV_Encoder/RBD_HA_classifier.ipynb](./code/CoV_Encoder/RBD_HA_classifier.ipynb)
* Classifier for RBD vs NTD vs S2: [./code/CoV_Encoder/RBD_S2_NTD_classifier_VH.ipynb](./code/CoV_Encoder/RBD_S2_NTD_classifier_VH.ipynb)

- Input files:
  - [./code/CoV_Encoder/data](./code/CoV_Encoder/data)  
- Output files: 
  - [./code/CoV_Encoder/result](./code/CoV_Encoder/result)

## Plotting
1. Plot summary statistics   
``Rscript code/Plot_summary_Piechart.R``   
    - Input file:
      - [./data/SARS-CoV-2-Abs.xlsx](./data/SARS-CoV-2-Abs.xlsx)
    - Output file:
      - [./graph/origin_pie.png](./graph/origin_pie.png)
      - [./graph/epitope_pie.png](./graph/epitope_pie.png)
      - [./graph/source_beeswarm.png](./graph/source_beeswarm.png)

2. Plot VDJ gene usage   
``Rscript code/Plot_VDJgene_Freq.R``   
    - Input files:
      - [./result/HVgene_freq.tsv](./result/HVgene_freq.tsv)
      - [./result/LVgene_freq.tsv](./result/LVgene_freq.tsv)
      - [./result/Dgene_freq.tsv](./result/Dgene_freq.tsv)
      - [./data/HV_Repertoire_freq.xlsx](./data/HV_Repertoire_freq.xlsx)
      - [./data/LV_Repertoire_freq.xlsx](./data/LV_Repertoire_freq.xlsx)
      - [./data/D_Repertoire_freq.xlsx](./data/D_Repertoire_freq.xlsx)
    - Output files:
      - [./graph/HV_gene_usage.png](./graph/HV_gene_usage.png)
      - [./graph/LV_gene_usage.png](./graph/LV_gene_usage.png)
      - [./graph/D_gene_usage.png](./graph/D_gene_usage.png)

3. plot IGHV/IGK(L)V pairing frequency   
``Rscript code/Plot_point_heatmap.R``   
    - Input files:
      - [./data/SARS-CoV-2-Abs.xlsx](./data/SARS-CoV-2-Abs.xlsx)
      - [./data/HV_Repertoire_freq.xlsx](./data/HV_Repertoire_freq.xlsx)
      - [./data/LV_Repertoire_freq.xlsx](./data/LV_Repertoire_freq.xlsx)
    - Output file:
      - [./graph/HLV_epitope_heatmap.png](./graph/HLV_epitope_heatmap.png)

4. Plot CDR H3 cluster size   
``Rscript code/plot_CDRH3_cluster_summary.R``   
    - Input file:
      - [./result/CDRH3_cluster_summary.tsv](./result/CDRH3_cluster_summary.tsv)
    - Output file:
      - [./graph/CDRH3_cluster_size.png](./graph/CDRH3_cluster_size.png)

5. Generate sequence logo for different CDR H3 clusters   
``python3 code/CDRH3_seqlogo.py``   
    - Input file:
      - [./result/Ab_info_CDRH3_clustering.tsv](./result/Ab_info_CDRH3_clustering.tsv)
    - Output files:
      - ./CDRH3_seqlogo/*.png

6. Plot IGHV gene usage for CDR H3 cluster 7   
``code/plot_cluster7_Vgenes.R``   
    - Input file:
      - [./result/Ab_info_CDRH3_clustering.tsv](./result/Ab_info_CDRH3_clustering.tsv)
    - Output file:
      - [./graph/Vgenes_cluster7.png](./graph/Vgenes_cluster7.png)

7. Plot analysis results for IGHD1-26-encoded S2 antibodies   
``Rscript code/plot_IGHD1-26_analysis.R``   
    - Input file:
      - [./data/SARS-CoV-2-Abs.xlsx](./data/SARS-CoV-2-Abs.xlsx)
    - Output file:
      - [./graph/S2_IGHD1-26_HVgenes.png](./graph/S2_IGHD1-26_HVgenes.png)
      - [./graph/S2_IGHD1-26_LVgenes.png](./graph/S2_IGHD1-26_LVgenes.png)
      - [./graph/HJgenes_pie_S2_IGHD1-26.png](./graph/HJgenes_pie_S2_IGHD1-26.png)
      - [./graph/HJgenes_pie_other.png](./graph/HJgenes_pie_other.png)
      - [./graph/S2_IGHD1-26_CDRH3_len.png](./graph/S2_IGHD1-26_CDRH3_len.png)
 
8. Plot SHM   
``Rscript code/plot_SHM.R``   
    - Input file:
      - [./result/SHM_frequency.tsv](./result/SHM_frequency.tsv)
    - Output files:
      - [./graph/SHM_HC_frequency.png](./graph/SHM_HC_frequency.png)
      - [./graph/SHM_LC_frequency.png](./graph/SHM_LC_frequency.png)

9. Plot the number of neutralizing vs non-neutralizing antibodies in each IGHV/IGK(L)V pair   
``Rscript code/Plot_basic_stat.R``   
    - Input file:
      - [./data/SARS-CoV-2-Abs.xlsx](./data/SARS-CoV-2-Abs.xlsx)
    - Output file:
      - [./RBD_neutralization_sactter](./RBD_neutralization_sactter)

10. Plot phylogenetic tree for the light chains of IGHV1-58/IGKV3-20 antibodies
``Rscript code/plot_tree_KV320.R``
    - Input files:
      - [./result/Cluster3_H158K320_LC.tree](./result/Cluster3_H158K320_LC.tree)
      - [./result/Cluster3_H158K320_LC_class.tsv](./result/Cluster3_H158K320_LC_class.tsv)
    - Output file:
      - [./graph/Cluster3_H158K320_LC_tree.png](./graph/Cluster3_H158K320_LC_tree.png)
