#!/usr/bin/env python3
############################################################
## Author: Jose F. Sanchez & Mireia Marin                 ##
## Copyright (C) 2022                                     ##
## High Content Genomics and Bioinformatics IGPT Unit     ## 
## Lauro Sumoy Lab, IGTP, Spain                           ##
############################################################

import os
import argparse
import sys
from subprocess import Popen, PIPE

## import my modules
from HCGB import functions
from RSP.config import set_config
from RSP.scripts import samtools


### INDEXING FUNCTION ###

def salmon_index(path_reference, reference_genome, index):

    #check if the index exist in the path_reference (even inside folders in the path_reference)
    
    command = "find " + path_reference + " -name duplicate_clusters.tsv | wc -l" #the indexing function creates a folder that contains this file (duplicate_clusters.tsv)
    
    stdout = Popen(command, shell=True, stdout=PIPE).stdout 
    output = stdout.read() #keeps the output in the terminal as a python variable
    
    reference_abs_path = os.path.join(path_reference, reference_genome) #path reference genome
    index_abs_path = os.path.join(path_reference, index) #path index
    
    selection = "0".encode()
   
    if selection in output:
        indexing = "salmon index -t "+ reference_abs_path + " -i "+ index_abs_path
        print("Indexing the transcriptome")
        print("Path to the index:", index_abs_path)
        print(indexing)
        
        os.system(indexing)
    
    else:
        print("Transcriptome is already indexed")
        print("Path to the index:", index_abs_path)
    
def salmon_quant(path_reference, index, reads_list, output, threads):
    
    path_index = os.path.join(path_reference, index) #path to the genome index
    
    #loop to check number of reads
         
    for i in reads_list: #for every list inside reads_list
        if len(i) == 3: #if there are three elements it's a pair-end analysis
            print("Pair-end analysis")
            read1 = i[1] #reads forward 
            read2 = i[2] #read reverse
            outputs_name = i[0] #variable that will go through the functions, sample name
            folder_results = os.path.join(output, outputs_name) #path to results folder + sample
            os.mkdir(folder_results)
            
            path_results = os.path.join(folder_results, outputs_name)          
            print ("Sample name", outputs_name, "\nRead forward:", read1, "\nRead reverse:", read2)

            quant = "salmon quant -i" + path_index + " -p "+ threads + " -l A -1 " + read1 + " -2 " + read2 + " -o " + path_results
            
            print(quant)
            os.system(quant)
            
        else:
            print("Single-end analysis")
            single_read = i[1] #single end analysis
            outputs_name = i[0] #variable that will go through the functions 
            folder_results = os.path.join(output, outputs_name) #path to results folder + sample
            os.mkdir(folder_results)
            
            path_results = os.path.join(folder_results, outputs_name)
            print ("Sample name", outputs_name, "\nSingle read:", single_read)
            
            quant = "salmon quant -i" + path_index + " -p "+ threads + " -l A -r " + single_read + " -o " + path_results
            
            print(quant)
            os.system(quamnt)
           

### MAIN FUNCTION ###  
    
def main():
    
    parser = argparse.ArgumentParser (description = 'Indexation and mapping using HISAT2') #all arguments are mandatory 
    parser.add_argument('-p', '--path_reference', help = 'Directory where the reference genome is, it will be the index directory too ', required = "TRUE")
    parser.add_argument('-g', '--reference_genome', help = 'Name of the file containing the reference genome', required = "TRUE")
    parser.add_argument('-i', '--index', help= 'It will create a folder called as your argument', required = "TRUE")
    parser.add_argument('-r', '--path_reads', help='Path to your txt with the format: $PATH/sample_name;$PATH/read1;$PATH/read or $PATH/sample_name;$PATH/read1. IT IS VERY IMPORTANT TO MANTAIN THE ORDER OF THE ELEMENTS AND THAT THEY ARE SEPARATED BY ;', required = "TRUE")
    parser.add_argument('-o', '--output', help= 'Path where the outputs will be created', required = "TRUE")
    parser.add_argument('-t', '--threads', help = 'Choose how many threads you want to use to execute HISAT2')

        
    args = parser.parse_args() #interprets the arguments 
    
    # python3 salmon.py -p /home/mireia/IGTP/salmon -g cds.fa.gz -i index -r /home/mireia/IGTP/reads/reads2 -o results
    
    ## inputs for other functions (can't be called as args.)
    
    path_reference = args.path_reference
    reference_genome = args.reference_genome
    index = args.index
    path_reads = args.path_reads
    output = args.output 
    threads = args.threads

  
    lines = open(path_reads).readlines() #read the txt with the information of the reads 
    reads_list = [] #empty list
    for i in lines: #for every sample
        reads_list.append(i.strip().split(';')) #split the different camps
        
    # lines.close() #close document
        
    print("Reads in main", reads_list)
    
    


    
if __name__ == '__main__':
    main()
    
