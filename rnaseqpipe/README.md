# rnaseqpipe
### _A pipeline for analyzing RNAseq data_
Developed in Institut Germans Trias i Pujol, high throughput genomics and bioinformatics department. 

### Description 
rnasepipe is a bioinformatic pipeline that maps and counts genomic variants from pair-end and single-end reads. 
It contains several modules. 
- Prep: it creates three differnt folders:
  - Info: useful information about the files used, the parameters for the softwares have used. 
  - Data: for each samples it will generate a folder
  - Report: what each module has done to samples
- QC : check input quality 
- trimm : deletes the adapters. It will re-calculate the quality after this module
- Map : map the input reads to the reference genome
- Count: create a matix with the variants found in each sample

### Softwares 
This pipeline includes many open source software. In each module we can use:
- Map module: HISAT2, Salmon, Kallisto, STAR

### Installation 