'''
Created on Oct 10, 2022

@author: mireia
'''

import os 
import argparse

def hisat2_index():
               
    os.chdir(args.path_index) #change to the working directory 

    indexing = 'hisat2-build '+args.reference_genome+' '+args.index #bash command 

    os.system(indexing)


def hisat2_mapping():

    ## TODO:
    ## abs_path_index, output_folder, name_sample, list_reads
    ## list_reads = [read_R1, read_R2]
    ## name_sample = sample1, sample2, sample3...

    path_sam = os.path.join(args.output_folder, name_sample )

    ## control length(list_reads)
    ## if==1
        # single-end

    ## ATTENTION
    mapping = 'hisat2 -x ' + abs_path_index + ' -1'+ list_reads[0] + ' -S'+path_sam #hisat2 mapping command 

    ## if==2
        # paired-end 
    mapping = 'hisat2 -x ' + abs_path_index + ' -1'+ list_reads[0] +' -2' + list_reads[1] + ' -S'+path_sam #hisat2 mapping command 


    os.system(mapping)
    return(path_sam)

def get_reads(folder_reads):
    ## get reads
    for r in reads: 
        if r.endswith('.fastq'): #extract the reads from all the files in the directory
            if "1" in r: #reads -1
                reads_1.append(r)
            else: # reads -2
                reads_2.append(r) 
                
     return([reads1, reads2])



def main():
    
    parser = argparse.ArgumentParser (description = 'Indexation and mapping using HISAT2')
    parser.add_argument('-p', '--path_index', help = 'Directory where the reference genome is, it will be the index directory too ')
    parser.add_argument('-g', '--reference_genome', help = 'Name of the file containing the reference genome')
    parser.add_argument('-i', '--index', help= 'Name you want your indexes to have')
    
    #parser.add_argument('-r', '--path_reads', help='Path to your reads')
    
    ## option:
    parser.add_argument('-r1', '--path_reads1', help='Path to your reads')
    parser.add_argument('-r2', '--path_reads2', help='Path to your reads')
    
    parser.add_argument('-o', '--output', help= 'Path where the sam stored')

    args = parser.parse_args() #interprets the arguments 
    
    ## 
    list_of_Reads = get_reads(args.path_reads)
    
    
    ## TODO
    hisat2_index()
    
    ## for each reads
        hisat2_mapping()
        ## samtools...
    

if __name__ == '__main__':
    main()

    
