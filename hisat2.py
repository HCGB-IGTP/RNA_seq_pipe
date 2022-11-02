import os 
import argparse
import glob
from samtools import *

####INDEXING FUNCTION########################################################################################################################

def hisat2_index(path_reference, reference_genome, index):

    check_index = glob.glob(os.path.join(path_reference, "*.ht2")) #save all the files with extension .ht2 into a list 
    
    #check sixe of the list 
    if len(check_index) > 0: #if there is any .ht2 file, we assume the indexation is present 
        print("The genome is already indexed. Index files:", check_index)
    
    else: #if there is no .ht2 in the path_reference --> build index 

        reference_abs_path = os.path.join(path_reference, reference_genome) #path reference genome
        index_abs_path = os.path.join(path_reference, index) #path index
        indexing = 'hisat2-build ' + reference_abs_path + ' ' + index_abs_path #bash command 
        print("Path to the index:", index_abs_path) 
        print(indexing)
        os.system(indexing)

####MAPPING FUNCTION#########################################################################################################################

def hisat2_mapping(path_reference, index, reads_list, output, threads):
        
    path_index = os.path.join(path_reference, index) #path to the genome index
    
    #loop to check number of reads
       
    for i in reads_list: #for every list inside reads_list
        if len(i) == 3: #if there are three elements it's a pair-end analysis
            print("Pair-end analysis")
            read1 = i[1] #reads forward 
            read2 = i[2] #read reverse
            outputs_name = i[0] #variable that will go through the functions, sample name
            print ("Sample name", outputs_name, "Read forward:", read1, "Read reverse:", read2)


            path_results = os.path.join(output, outputs_name) #path to results folder + sample name 
            path_sam = path_results + ".sam" # safe sam in results folder
                    
            mapping = "hisat2 -x " + path_index + " -p " + threads + " -1 " + read1 + " -2 " + read2 + " -S " + path_sam #hisat2 paired end mapping command  
            print(mapping)
            os.system(mapping)
            
            SamToBam(path_results, path_sam) 
            
        else:
            print("Single-end analysis")
            single_read = i[1] #single end analysis
            outputs_name = i[0] #variable that will go through the functions 
            print ("Sample name", outputs_name, "Single read:", single_read)

            name_sam = outputs_name +".sam" #sample name

            
            path_results = os.path.join(output, outputs_name) #path to results folder + sample name 
            path_sam = path_results + ".sam" # safe sam in results folder
            
            mapping = "hisat2 -x " + path_index + " -p "+ threads + " -U " + single_read + " -S " + path_sam #hisat2 single end mapping command  
            
            print(mapping)
            os.system(mapping)
            SamToBam(path_results, path_sam)


####MAIN FUNCTION##############################################################################################################################

def main():
    
    parser = argparse.ArgumentParser (description = 'Indexation and mapping using HISAT2') #all arguments are mandatory 
    parser.add_argument('-p', '--path_reference', help = 'Directory where the reference genome is, it will be the index directory too ', required = "TRUE")
    parser.add_argument('-g', '--reference_genome', help = 'Name of the file containing the reference genome', required = "TRUE")
    parser.add_argument('-i', '--index', help= 'Name you want your indexes to have', required = "TRUE")
    parser.add_argument('-r', '--path_reads', help='Path to your txt with the format: $PATH/sample_name;$PATH/read1;$PATH/read or $PATH/sample_name;$PATH/read1. IT IS VERY IMPORTANT TO MANTAIN THE ORDER OF THE ELEMENTS AND THAT THEY ARE SEPARATED BY ;', required = "TRUE")
    parser.add_argument('-o', '--output', help= 'Path where the sam will be created', required = "TRUE")
    parser.add_argument('-t', '--threads', help = 'Choose how many threads you want to use to execute HISAT2')
    
    args = parser.parse_args() #interprets the arguments 
    
    # python3 hisat2.py -p reference -g chr22.fa -i index -r reads/reads2 -o results
    
    ## inputs for other functions (can't be called as args.)
    
    path_reference = args.path_reference
    reference_genome = args.reference_genome
    index = args.index
    path_reads = args.path_reads
    output = args.output
    threads = args.threads
  
    lines = open(path_reads).readlines() #read the txt with the information of the reads 
    reads_list = [] #empty list
    for i in lines: #for every sample
        reads_list.append(i.strip().split(';')) #split the different camps
        
    # lines.close() #close document
        
    print("Reads in main", reads_list)
    
    
    #function calling
    
    hisat2_index(path_reference, reference_genome, index)
    hisat2_mapping(path_reference, index, reads_list, output, threads)


    
if __name__ == '__main__':
    main()
    
    
