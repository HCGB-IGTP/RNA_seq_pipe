'''
Created on Oct 8, 2022

@author: mireia
'''

import os 
import sys

#########################################################
def hisat2_index(path_reference, reference, index):

    ## OLD
    #os.chdir(path_reference) #change to the working directory 
    #indexing = "hisat2-build "+reference+" "+index #bash command 
    
    ## new variable
    reference_abs_path = os.path.join(path_reference, reference)
    index_abs_path = os.path.join(path_reference, index)
    indexing = "hisat2-build " + reference_abs_path + " " + index_abs_path #bash command 

    os.system(indexing)
    
    return(index_abs_path)    
    

 #########################################################
 def hisat2_mapping(index_join, path_reads, path_output):   
    
    os.chdir(path_reads)

    reads = os.listdir() #saves files in current directory as a variable 

    reads_1 = [] #list with reads forward
    reads_2 = [] #list with reads reverse  

    """
    !!!!!!!!!!!!! This loop must be changed when using real named reads !!!!!!!!!!!

    """
    
    ## get reads
    for r in reads: 
        if r.endswith('.fastq'): #extract the reads from all the files in the directory
            if "1" in r: #reads -1
                reads_1.append(r)
            else: # reads -2
                reads_2.append(r) 

    """
    """

    ## create mapping
    for i,j in zip (reads_1, reads_2): #for every couple of reads we run the hisat2 command
        forward = i #-1 read
        reverse = j  # -2 read
        name=i+j+".sam" #name of output sam 
        name=name.replace(".fastq", "") #aesthetic name
        path_sam = os.path.join(path_output,name)

        mapping = "hisat2 -x"+index_join+" -1"+forward+" -2"+reverse+" -S"+path_sam #hisat2 mapping command 

        os.system(mapping)
   
#########################################################
def main():  
    ## sys,argv[0] script
    ## sys.argv[1] #first argument --> working directory 
    ## sys.argv[2]  #second argument--> reference genome name 
    ## sys.argv[3] #third argument --> index name
    ## sys.argv[4] # reads folder
    ## sys.argv[5] #output
    
    ## check length(sys.argv) == 6
    ## si error: 
        ## help: 
        ## python hisat2_indexation_mapping.py /dir/folder/with/reference_folder file_name_reference index_name reads_folder output_folder
    
    ## index
    path_reference = sys.argv[1] #first argument --> working directory 
    reference = sys.argv[2]  #second argument--> reference genome name 
    index = sys.argv[3] #third argument --> index name
    
    ### call index
    index_abs_path = hisat2_index(path_reference, reference, index)

    ## map reads
    path_reads = sys.argv[4] #reads path
    path_output = sys.argv[5] #output path 

    hisat2_mapping(index_abs_path, path_reads, path_output)
    
    
#########################################################
if __name__ == '__main__':
    main()
       
