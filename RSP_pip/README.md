# RSP
### _A pipeline for analyzing RNAseq data_
Developed in Institut Germans Trias i Pujol, high throughput genomics and bioinformatics department. 

### Description 
RSP is a bioinformatic pipeline that maps and counts genomic variants from pair-end and single-end reads. 
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

RSP will require python v3.7 and java (we tested in openjdk 14 2020-03-17).

The RSP python pipeline will be available in `pip` and also available using `conda`.

RSP depends on multiple third party software that we have listed here [TODO: Add link].

We encourage you to install RSP and all dependencies using the `conda` environment we created and following these instructions. 

To create a new conda environment, install third party software, install RSP and missing dependencies, do as follows:

1) Get requirements file from RSP git repo  [TODO: Add link].

```sh
wget https://raw.githubusercontent.com/HCGB-IGTP/XICRA/master/XICRA_pip/devel/conda/environment.yml
```

or available under the git repository in:

```sh
ls RNA_seq_pipe/RSP_pip/devel/conda/environment.yml
```

2) Create environment named RSP and install required packages using conda: 

```sh
conda env create -f environment.yml
```

3) Activate environment and install RSP
```sh
## activate
conda activate RSP

## install latest python code
pip install RSP
```

Or install developer python code by doing:

```sh
## get git repo
git clone https://github.com/HCGB-IGTP/RNA_seq_pipe.git
cd RNA_seq_pipe/RSP_pip
sh devel/pypi/test_module.sh
```

To check everything is fine, try executing the `config` module:
```sh
RSP config
```

To additionally check RSP is working try to run a test example from https://github.com/griffithlab/rnaseq_tutorial

There are instructions to download this dataset within RSP git repository in: 

```
cd toy_dataset/
sh get_toyDataset.sh
sh execute_test.sh
```

or available to download by executing the following command:

```
RSP test
sh execute_test.sh
```

