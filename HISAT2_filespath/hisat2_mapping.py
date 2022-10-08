'''
Created on Oct 8, 2022

@author: mireia
'''

import sys
import os 

def hisat2_mapping():
    

    path_index = sys.argv[1] #index path 
    index = sys.argv[2] #index name
    index_join = os.path.join(path_index,index) #creates correct path to index

    path_reads = sys.argv[3] #reads path
    path_output = sys.argv[4] #output path 
    
    os.chdir(path_reads)
    
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
        path_sam = os.path.join(path_output,name)
        
        mapping = "hisat2 -x"+index_join+" -1"+forward+" -2"+reverse+" -S"+path_sam #hisat2 mapping command 
        
        os.system(mapping)
        
hisat2_mapping()
    


    
    
    
