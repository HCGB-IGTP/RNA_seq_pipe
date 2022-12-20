#!/usr/bin/env python3
############################################################
## Author: Jose F. Sanchez                                ##
## Copyright (C) 2022                                     ##
## High Content Genomics and Bioinformatics IGPT Unit     ## 
## Lauro Sumoy Lab, IGTP, Spain                           ##
############################################################
'''
Calls cutadapt to trim raw reads
'''
## useful imports
import time
import io
import os
import re
import sys
from sys import argv
from io import open
from termcolor import colored

## import my modules
from HCGB import functions
from RSP.config import set_config

#############################################
def cutadapt(cutadapt_exe, reads, path, sample_name, num_threads, min_len_given, Debug, adapters, extra):
    """
    Executes cutadapt sofware for each sample cutting the adapters of each read
    
    :param cutadapt_exe: to call cutadapt software
    :param reads: name of the fastqc files of the sample to be trimmed
    :param path: path to the sample folder to store the results
    :param sample_name: name of the sample to be trimmed
    :param num_threads: number of CPUs to use.
    :param Debug: show additional message for debugging purposes.
    :param adapters: dictionary with the introduced adapters
    :param extra: provided extra options for cutadapt trimming process
    
    :type cutadapt_exe: string
    :type reads: string
    :type path: string
    :type sample_name: string
    :type num_threads: string
    :type Debug: boolean
    :type adapters: dictionary
    :type extra: string
    
    :returns: the trimmed files
    """
    logfile = os.path.join(path, sample_name + '.cutadapt.log')
    
    if (len(reads) == 2):
        if not adapters['adapter_a'] or not adapters['adapter_A']:
             print ("** ERROR: Missing adapter information")
             exit()
        
        if extra:
            o_param = os.path.join(path, sample_name + '_temp1_trim_R1.fastq')
            p_param = os.path.join(path, sample_name + '_temp1_trim_R2.fastq')
        else:
            p_param = os.path.join(path, sample_name + '_trim_R2.fastq')
            o_param = os.path.join(path, sample_name + '_trim_R1.fastq')
        
        ## paired-end mode, 15 bps as the min length cutoff
        cmd = '%s -j %s -m %s -a %s -A %s -o %s -p %s %s %s > %s' %(cutadapt_exe,  
                                                                       num_threads, min_len_given, 
                                                                       adapters['adapter_a'], 
                                                                       adapters['adapter_A'], o_param, 
                                                                       p_param, reads[0], reads[1], logfile)
    elif (len(reads) == 1):
        if not adapters['adapter_a']:
             print ("** ERROR: Missing adapter information")
             exit()

        if extra:
            o_param = os.path.join(path, sample_name + '_temp1_trim.fastq')
        else:
            o_param = os.path.join(path, sample_name + '_trim.fastq')
        
        ## single-end mode:
        cmd = '%s -j %s -m %s -a %s -o %s %s > %s' %(cutadapt_exe, num_threads, 
                                                     min_len_given,
                                                     adapters['adapter_a'], 
                                                     o_param, reads[0], logfile)    
    else:
        print ('** Wrong number of files provided for sample: %s...' %sample_name)
        return(False)

    ##
    code = functions.system_call_functions.system_call(cmd)

    ## if additional options, run a second cutadapt command
    ## to ensure this options take effect.
    if (extra):
        if (len(reads) == 2):
            o_param2 = os.path.join(path, sample_name + '_trim_R1.fastq')
            p_param2 = os.path.join(path, sample_name + '_trim_R2.fastq')
        
            ## paired-end mode
            extra_cmd = '%s %s -j %s -m %s -a %s -A %s -o %s -p %s %s %s >> %s' %(cutadapt_exe, 
                                                                                  extra, 
                                                                                  num_threads, 
                                                                                  min_len_given, 
                                                                                  adapters['adapter_a'], 
                                                                                  adapters['adapter_A'],
                                                                                  o_param2, p_param2, 
                                                                                  o_param, p_param, logfile)
        
        elif (len(reads) == 1):
            o_param2 = os.path.join(path, sample_name + '_trim.fastq')
            ## single-end mode:
            extra_cmd = '%s %s -j %s -m %s -a %s -o %s %s >> %s' %(cutadapt_exe, 
                                                                   extra, 
                                                                   num_threads, 
                                                                   min_len_given,
                                                                   adapters['adapter_a'], 
                                                                   o_param2, o_param, logfile)    
        
        code2 = functions.system_call_functions.system_call(extra_cmd)
        
        ## remove: o_param p_param
        if (len(reads) == 2):        
            os.remove(p_param)
       
        os.remove(o_param)
        return (code2)
    
    else:
        return (code)
