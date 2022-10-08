'''
Created on Oct 8, 2022

@author: mireia
'''

import os 
import sys


def hisat2_index():
    
    path_reference = sys.argv[1] #first argument --> working directory 
    reference = sys.argv[2]  #second argument--> reference genome name 
    index = sys.argv[3] #third argument --> ouput name
   
    os.chdir(path_reference) #change to the working directory 
    
    indexing = "hisat2-build "+reference+" "+index #bash command 
    
    os.system(indexing)

hisat2_index()