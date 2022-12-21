#!/usr/bin/env python3
############################################################
## Author: Jose F. Sanchez & Mireia Marin                 ##
## Copyright (C) 2022                                     ##
## High Content Genomics and Bioinformatics IGPT Unit     ## 
## Lauro Sumoy Lab, IGTP, Spain                           ##
############################################################
"""
Calls multiQC to generate HTML statistics reports.
"""
## useful imports
import os
import io
import sys
from io import open
from sys import argv
from termcolor import colored

import os 
import argparse
import glob

## import my modules
from HCGB import functions
import HCGB.functions.system_call_functions as HCGB_sys

from RSP.config import set_config
from RSP.scripts import samtools

####INDEXING FUNCTION########################################################################################################################
def hisat2_index(path_reference, reference_genome, index, threads, Debug):

    check_index = glob.glob(os.path.join(path_reference, index+"*.ht2")) #save all the files with extension .ht2 into a list 
    
    #check sixe of the list 
    if len(check_index) > 0: #if there is any .ht2 file, we assume the indexation is present 
        print("The genome is already indexed. Index files:", check_index)

    else: #if there is no .ht2 in the path_reference --> build index 
        
        hisat2_build = set_config.get_exe('hisat2-build')

        reference_abs_path = reference_genome #path reference genome
        index_abs_path = os.path.join(path_reference, index) #path index
        indexing = hisat2_build + ' -p ' + str(threads) + ' ' + reference_abs_path + ' ' + index_abs_path #bash command 
        print("Path to the index:", index_abs_path) 

        ## system call & return
        code = HCGB_sys.system_call(indexing)
        if (code):    
            return(index_abs_path)

####MAPPING FUNCTION#########################################################################################################################
def hisat2_mapping(sample_name, index_path_reference, reads_list, output, threads, extra_params, Debug):
        
    outputs_name = sample_name #variable that will go through the functions, sample name
    path_results = os.path.join(output, outputs_name)
    path_sam = path_results + ".sam" # safe sam in results folder
    errLog = path_results + ".err" 
    outLog = path_results + ".log" 
    outSummary = path_results + ".summary"
        
    #loop to check number of reads
    hisat2 = set_config.get_exe('hisat2')
       
    if len(reads_list) == 2: #if there are two elements it's a pair-end analysis
        read1 = reads_list[0] #reads forward 
        read2 = reads_list[1] #read reverse
        
        if (Debug):
            print (colored("**DEBUG: Sample name **", 'yellow'))
            print (outputs_name)

            print (colored("**DEBUG: read1 **", 'yellow'))
            print (read1)
            
            print (colored("**DEBUG: read2 **", 'yellow'))
            print (read2)

        ## hisat call
        #hisat2 paired end mapping command
        mapping = hisat2 + " -x " + index_path_reference + " -p " + str(threads) 
        mapping = mapping + " -1 " + read1 + " -2 " + read2 
        mapping = mapping + " -S " + path_sam
        mapping = mapping + " --rg-id " + sample_name + ' --rg ' + sample_name
        mapping = mapping + " --new-summary --summary-file " + outSummary
        mapping = mapping + extra_params
        mapping = mapping + " > " +  outLog + " 2> " + errLog
        
        ## system call & return
        code = HCGB_sys.system_call(mapping)
        if (code):
            return(samtools.sam_to_bam(path_results, path_sam, threads, Debug))
        else:
            print("Some error occurred during mapping PE reads...") 
        
    else:
        print("Single-end analysis")
        single_read = reads_list[0] #single end analysis

        if (Debug):
            print (colored("**DEBUG: Sample name **", 'yellow'))
            print (outputs_name)

            print (colored("**DEBUG: single_read **", 'yellow'))
            print (single_read)
            
            
        ## mapping call
        mapping = hisat2 + " -x " + index_path_reference + " -p " + str(threads) 
        mapping = mapping + " -U " + single_read
        mapping = mapping + " -S " + path_sam 
        mapping = mapping + " --rg-id " + sample_name + ' --rg ' + sample_name
        mapping = mapping + " --new-summary --summary-file " + outSummary
        mapping = mapping + extra_params
        mapping = mapping + " > " +  outLog + " 2> " + errLog
        
        ## system call & return
        code = HCGB_sys.system_call(mapping)
        if (code):
            return(samtools.sam_to_bam(path_results, path_sam, threads, Debug))
        else:
            print("Some error occurred during mapping SE reads...") 
            

####MAIN FUNCTION##############################################################################################################################
def main():
    
    parser = argparse.ArgumentParser (description = 'Indexation and mapping using HISAT2') #all arguments are mandatory 
    parser.add_argument('-p', '--path_reference', help = 'Directory where the reference genome is, it will be the index directory too ', required = "TRUE")
    parser.add_argument('-g', '--reference_genome', help = 'Name of the file containing the reference genome', required = "TRUE")
    parser.add_argument('-i', '--index', help= 'Name you want your indexes to have, by default the program check if the index already exists in the path_reference', required = "TRUE")
    parser.add_argument('-r', '--path_reads', help='Path to your txt with the format: $PATH/sample_name;$PATH/read1;$PATH/read or $PATH/sample_name;$PATH/read1. IT IS VERY IMPORTANT TO MANTAIN THE ORDER OF THE ELEMENTS AND THAT THEY ARE SEPARATED BY ;', required = "TRUE")
    parser.add_argument('-o', '--output', help= 'Path where the sam will be created', required = "TRUE")
    parser.add_argument('-t', '--threads', help = 'Choose how many threads you want to use to execute HISAT2')
    parser.add_argument('-c', '--gtf', help = 'Path to your gtf')

    
    args = parser.parse_args() #interprets the arguments 
    
    # python3 hisat2.py -p reference -g chr22.fa -i index -r /home/mireia/IGTP/reads/reads2 -o results
    
    ## inputs for other functions (can't be called as args.)
    
    path_reference = args.path_reference
    reference_genome = args.reference_genome
    index = args.index
    path_reads = args.path_reads
    output = args.output
    threads = args.threads
    gtf = args.gtf
  
    #function calling
  
    if function == "Indexing": #if in --function they choose indexing just the index function will run
        hisat2_index(path_reference, reference_genome, index)
    elif function == "Mapping":#if in --function they choose mapping just the mapping function will run
        hisat2_mapping(path_reference, index, reads_list, output, threads,gtf)
    else: #if they don't specify the function, both will run
        hisat2_index(path_reference, reference_genome, index)
        hisat2_mapping(path_reference, index, reads_list, output, threads,gtf)



    
if __name__ == '__main__':
    main()
    
    
