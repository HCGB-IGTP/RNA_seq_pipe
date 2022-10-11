'''
Created on Oct 10, 2022

@author: mireia
'''

import os 
import argparse


def main():
    
    parser = argparse.ArgumentParser (description = 'Indexation and mapping using HISAT2')
    parser.add_argument('-p', '--path_index', help = 'Directory where the reference genome is, it will be the index directory too ')
    parser.add_argument('-g', '--reference_genome', help = 'Name of the file containing the reference genome')
    parser.add_argument('-i', '--index', help= 'Name you want your indexes to have')
    parser.add_argument('-r', '--path_reads', help='Path to your reads')
    parser.add_argument('-o', '--output', help= 'Path where the sam stored')

    args = parser.parse_args() #interprets the arguments 
    
    def hisat2_index():
               
        os.chdir(args.path_index) #change to the working directory 
        
        indexing = 'hisat2-build '+args.reference_genome+' '+args.index #bash command 
        
        os.system(indexing)
       
    
    def hisat2_mapping():
        
        index_join = os.path.join(args.path_index,args.index) #creates correct path to index
    
        os.chdir(args.path_reads)
        
        reads = os.listdir() #saves files in current directory as a variable 
    
        reads_1 = [] #list with reads forward
        reads_2 = [] #list with reads reverse  
        
        """
        !!!!!!!!!!!!! This loop must be changed when using real named reads !!!!!!!!!!!
        
        """
        
        for r in reads: 
            if r.endswith('.fastq'): #extract the reads from all the files in the directory
                if '1' in r: #reads -1
                    reads_1.append(r)
                else: # reads -2
                    reads_2.append(r) 
            
        """
        """
        
        for i,j in zip (reads_1, reads_2): #for every couple of reads we run the hisat2 command
            forward = i #-1 read
            reverse = j  # -2 read
            name=i+j+'.sam' #name of output sam 
            name=name.replace('.fastq', '') #aesthetic name
            path_sam = os.path.join(args.output,name)
            
            mapping = 'hisat2 -x'+index_join+' -1'+forward+' -2'+reverse+' -S'+path_sam #hisat2 mapping command 
            
            os.system(mapping)
   
    hisat2_index()
    hisat2_mapping()    

main()
    
if __name__ == '__main__':
    pass

    
