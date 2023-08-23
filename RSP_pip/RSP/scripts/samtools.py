#!/usr/bin/env python3
############################################################
## Author: Jose F. Sanchez & Mireia Marin                 ##
## Copyright (C) 2022                                     ##
## High Content Genomics and Bioinformatics IGPT Unit     ## 
## Lauro Sumoy Lab, IGTP, Spain                           ##
############################################################

import os
import sys
## import my modules
from HCGB import functions
from RSP.config import set_config
import HCGB.functions.system_call_functions as HCGB_sys
from builtins import str
from termcolor import colored

#### SAM TO BAM FUNCTION ####################################
def sam_to_bam(path_results, path_sam, threads, Debug):
    
    samtools_path = set_config.get_exe('samtools', Debug=Debug)
    
    if (Debug):
        print (colored("**DEBUG: samtools_path **", 'yellow'))
        print (samtools_path)

        print (colored("**DEBUG: path_results **", 'yellow'))
        print (path_results)

        print (colored("**DEBUG: path_sam **", 'yellow'))
        print (path_sam)
        
        print (colored("**DEBUG: threads **", 'yellow'))
        print (str(threads))

    print("+ Converting SAM to BAM")
    path_bam = path_results + ".bam" #safe bam into results folder
    sam_to_bam = samtools_path + " view -bS "+ path_sam + " --threads " + str(threads) + " > " + path_bam  #samtools sam to bam command

    ## system call & return
    code = HCGB_sys.system_call(sam_to_bam)
    if code!="OK":
        print("An error occurred when converting sam file: " + path_sam)
        return(False)

    return(path_bam)


#### BAM TO SORTED_BAM FUNCTION ###############
def bam_to_sorted_bam(path_results, path_bam, threads, Debug):
    
    samtools_path = set_config.get_exe('samtools')
    
    print("+ Converting BAM to Sorted BAM")
    sorted_bam_name = path_results +".sorted.bam" #safe sorted bam into results folder
    bam_to_sorted_bam = samtools_path + " sort "+ path_bam + " --threads " + str(threads) + " -o " + sorted_bam_name #samtools bam to sorted bam command

    ## system call & return
    code = HCGB_sys.system_call(bam_to_sorted_bam)
    if code!="OK":
        print("An error occurred when converting bam file: " + sorted_bam_name)
        return(False)
    
    return(sorted_bam_name)

#### SAM TO SORTED_BAM FUNCTION ###############
def sam_to_sorted_bam(path_results, path_sam, threads, Debug):
    print()

    tmp_bam = sam_to_bam(path_results, path_sam, threads, Debug)
    sorted_bam = bam_to_sorted_bam(path_results, tmp_bam, threads, Debug)

    return(sorted_bam)
   

def main():

    ## control if options provided or help
    if len(sys.argv) > 3:
        print ("+ Convert SAM-BAM files")
    else:
        print ("*** ATTENTION: Provide the following arguments:***\npython samtools.py path_results path_sam num_threads [sam2bam|sam2sortedbam|bam2sortedbam]")
        exit()


    if sys.argv[4]=="sam2sortedbam":
        sorted_bam = sam_to_sorted_bam(sys.argv[1], sys.argv[2], sys.argv[3], True)

    if sys.argv[4]=="bam2sortedbam":
        bam_to_sorted_bam(sys.argv[1], sys.argv[2], sys.argv[3], True)



if __name__ == '__main__':
    main()
    
