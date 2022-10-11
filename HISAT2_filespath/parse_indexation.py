'''
Created on Oct 10, 2022

@author: mireia
'''

import os 
import argparse


def hisat2_index():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path_reference')
    parser.add_argument('-r','--reference')
    parser.add_argument('-i', '--index')

    args = parser.parse_args()
    
    os.chdir(args.path_reference) #change to the working directory 
    
    indexing = "hisat2-build "+args.reference+" "+args.index #bash command 
    
    os.system(indexing)

hisat2_index() 
    