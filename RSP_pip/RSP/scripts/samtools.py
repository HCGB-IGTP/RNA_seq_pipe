#!/usr/bin/env python3
############################################################
## Author: Jose F. Sanchez & Mireia Marin                 ##
## Copyright (C) 2022                                     ##
## High Content Genomics and Bioinformatics IGPT Unit     ## 
## Lauro Sumoy Lab, IGTP, Spain                           ##
############################################################

import os
## import my modules
from HCGB import functions
from RSP.config import set_config
import HCGB.functions.system_call_functions as HCGB_sys


####SAM TO BAM FUNCTION########################################################################################################################
def sam_to_bam(path_results, path_sam, threads, gtf):
    
    samtools_path = set_config.get_exe('samtools')
    print("Converting SAM to BAM")
    path_bam = path_results + ".bam" #safe bam into results folder
    sam_to_bam = samtools_path + " view -bS "+ path_sam + " --threads " + threads + " > " + path_bam  #samtools sam to bam command

    ## system call & return
    code = HCGB_sys.system_call(sam_to_bam)
    return(code)


####BAM TO SORTED_BAM FUNCTION########################################################################################################################
def bam_to_sorted_bam(path_results, path_bam, threads, gtf):
    
    samtools_path = set_config.get_exe('samtools')
    
    print("Converting BAM to Sorted BAM")
    sorted_bam_name = path_results +".sorted.bam" #safe sorted bam into results folder
    bam_to_sorted_bam = samtools_path + " sort "+ path_bam + " --threads " + threads + " -o " + sorted_bam_name #samtools bam to sorted bam command

    ## system call & return
    code = HCGB_sys.system_call(bam_to_sorted_bam)
    return(code)


if __name__ == '__main__':
    main()
    
