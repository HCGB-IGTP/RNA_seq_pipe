#!/usr/bin/env python3
############################################################
## Author: Jose F. Sanchez & Mireia Marin                 ##
## Copyright (C) 2022                                     ##
## High Content Genomics and Bioinformatics IGPT Unit     ## 
## Lauro Sumoy Lab, IGTP, Spain                           ##
############################################################

import os 

def featurecounts (path_results, sorted_bam_name, threads, gtf):

    output = os.path.join(path_results, "Counts_Matrix.txt") 
       
    featurecounts_command = "featureCounts -p -O -T " + threads + " -a " + gtf + " -o " + output + " " + sorted_bam_name 
    
    print(featurecounts_command)
    
    os.system(featurecounts_command) 

