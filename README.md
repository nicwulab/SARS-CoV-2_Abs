This README describes the scripts used for the sequence analysis in "A large-scale systematic survey of SARS-CoV-2 antibodies
reveals recurring molecular features"


## SARS-CoV-2 Antibodies


## Dependencies ##
* python=3.9
* [igblast](https://github.com/ncbi/igblast)
* [PyIR](https://github.com/crowelab/PyIR)
* [BioPython](https://github.com/biopython/biopython)
## Installation ##
Install everything dependencies by conda:

```conda create -n Abs -c bioconda -c anaconda -c conda-forge python=3.9 biopython igblast```


### Local igblast set up

```conda activate Abs(for Mac)```

PyIR:An IgBLAST wrapper and parser[https://github.com/crowelab/PyIR]

```pip3 install crowelab_pyir```

database set up in pyir library directory

```pyir setup```

- manuualy intall amino acid database from imgt

1. sequence download from  http://www.imgt.org/vquest/refseqh.html#VQUEST
2. copy and paste, save as fasta(save all V gene in one file)
3. clean data (raw edit_imgt_file.pl can be found on igblast-1.17.1xxx/bin)

```edit_imgt_file.pl imgt_database/human_prot/imgt_raw/IGV.fasta > imgt_database/human_prot/IGV.fasta```

4. make database( For protein sequence,DNA using -dbtype nucl )

```makeblastdb -parse_seqids -dbtype prot -in imgt_database/human_prot/IGV.fasta```

- run PyIR for igBlast

  see [PyIr.py](./code/_PyIR_.py)

- run local igblast and CDR parser

### Baseline VDJ setup

### CDR H3 clustering analysis

### Identification of recurring somatic hypermutation (SHM)![image](https://user-images.githubusercontent.com/53042063/143530364-0e7c69e1-f61d-40a7-97b1-464df17dedbf.png)

### Deep learning model for antigen identification

### Plotting

