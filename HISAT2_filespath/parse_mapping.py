'''
Created on Oct 10, 2022

@author: mireia
'''

import argparse
import os 

def hisat2_mapping():
    
    parser = argparse.ArgumentParser (description = 'Indexation and mapping using HISAT2')
    parser.add_argument('-p', '--path_index', help = 'Directory where the reference genome is, it will be the index directory too ')
    parser.add_argument('-i', '--index', help= 'Name you want your indexes to have')
    parser.add_argument('-r', '--path_reads', help='Path to your reads')
    parser.add_argument('-o', '--output', help= 'Path where the sam stored')

    args = parser.parse_args() #interprets the arguments 
        
    index_join = os.path.join(args.path_index, args.index)
    
    os.chdir(args.path_reads)
    
    reads = os.listdir() #saves files in current directory as a variable 

    reads_1 = [] #list with reads forward
    reads_2 = [] #list with reads reverse  
    
    """
    This loop must be changed when using real named reads --> NOT FINAL VERSION
    
    """
    
    for r in reads: 
        if r.endswith('.fastq'): #extract the reads from all the files in the directory
            if "1" in r: #reads -1
                reads_1.append(r)
            else: # reads -2
                reads_2.append(r) 
        
    """
    """
    
    for i,j in zip (reads_1, reads_2): #for every couple of reads we run the hisat2 command
        forward = i #-1 read
        reverse = j  # -2 read
        name=i+j+".sam" #name of output sam 
        name=name.replace(".fastq", "") #aesthetic name
        path_sam = os.path.join(args.output,name)
        
        mapping = "hisat2 -x"+index_join+" -1"+forward+" -2"+reverse+" -S"+path_sam #hisat2 mapping command 
        
        os.system(mapping)
        
hisat2_mapping()
    