import os
import argparse


### INDEXING FUNCTION ###

def kallisto_index(path_reference, reference_genome, index, kmers):

    
    reference_abs_path = os.path.join(path_reference, reference_genome) #path reference genome
    index_abs_path = os.path.join(path_reference, index) #path index
    
    #check if the index has already been created
    
    if os.path.exists(index_abs_path) == False:
        indexing = "kallisto index -i " + index_abs_path + path_reference #kmers option can be added

        print("Indexing the transcriptome")
        print("Path to the index:", index_abs_path)
        print(indexing)
        #os.system(indexing)
    
    else:
        print("Transcriptome is already indexed")
        print("Path to the index:", index_abs_path)
    
def kallisto_quant(path_reference, index, reads_list, output, threads):
    
    path_index = os.path.join(path_reference, index) #path to the genome index
    
    #loop to check number of reads
         
    for i in reads_list: #for every list inside reads_list
        if len(i) == 3: #if there are three elements it's a pair-end analysis
            print("Pair-end analysis")
            read1 = i[1] #reads forward 
            read2 = i[2] #read reverse
            outputs_name = i[0] #variable that will go through the functions, sample name
            path_results = os.path.join(output, outputs_name) #path to results folder + sample            
            print ("Sample name", outputs_name, "Read forward:", read1, "Read reverse:", read2)

            quant =  "kallisto quant -i " + path_index + " -t " + threads + " -o " + ouput + read1 + read2 

            
            print(quant)
            #os.system(quant)
            
        else:
            print("Single-end analysis")
            single_read = i[1] #single end analysis
            outputs_name = i[0] #variable that will go through the functions 
            path_results = os.path.join(output, outputs_name) #path to results folder + sample
            print ("Sample name", outputs_name, "Single read:", single_read)
            
            quant = "kallisto quant -i " + path_index + " -t " + threads + " -o " + ouput + read1             
            print(quant)
            #os.system(quamnt)
           

### MAIN FUNCTION ###  
    
def main():
    
    parser = argparse.ArgumentParser (description = 'Indexation and mapping using HISAT2') #all arguments are mandatory 
    parser.add_argument('-p', '--path_reference', help = 'Directory where the reference genome is, it will be the index directory too ', required = "TRUE")
    parser.add_argument('-g', '--reference_genome', help = 'Name of the file containing the reference genome', required = "TRUE")
    parser.add_argument('-i', '--index', help= 'Name of the folder containing the index (it must be in the same path of the reference) or name you want to give to the folder that will be created to store the index', required = "TRUE")
    parser.add_argument('-r', '--path_reads', help='Path to your txt with the format: $PATH/sample_name;$PATH/read1;$PATH/read or $PATH/sample_name;$PATH/read1. IT IS VERY IMPORTANT TO MANTAIN THE ORDER OF THE ELEMENTS AND THAT THEY ARE SEPARATED BY ;', required = "TRUE")
    parser.add_argument('-o', '--output', help= 'Path where the outputs will be created', required = "TRUE")
    parser.add_argument('-t', '--threads', type = int, help = 'Choose how many threads you want to use to execute HISAT2')
    parser.add_argument('-f', '--function', choices = ['Indexing', 'Quantifying'], help = 'You can choose to execute the indexing function, the quantifying function or, by default, both. Possible choices for this argument: Indexing / Quantifying')
    parser.add_argument ('-k', '--kmers', type = int, help = 'Desired k-mers length when indexing. Default: 31')
        
    args = parser.parse_args() #interprets the arguments 
    
    
    
    ## inputs for other functions (can't be called as args.)
    
    path_reference = args.path_reference
    reference_genome = args.reference_genome
    index = args.index
    path_reads = args.path_reads
    output = args.output 
    threads = args.threads
    function = args.function
    kmers =args.kmers
    
  
    lines = open(path_reads).readlines() #read the txt with the information of the reads 
    reads_list = [] #empty list
    for i in lines: #for every sample
        reads_list.append(i.strip().split(';')) #split the different camps
        
    # lines.close() #close document
        
    print("Reads in main", reads_list)
    
    
    #function calling
    
    if function == "Indexing": #if in --function they choose indexing just the index function will run
        kallisto_index(path_reference, reference_genome, index, kmers)
    elif function == "Quantifying":#if in --function they choose mapping just the mapping function will run
        kallisto_quant(path_reference, index, reads_list, output, threads)
    else: #if they don't specify the function, both will run
        kallisto_index(path_reference, reference_genome, index)
        kallisto_quant(path_reference, index, reads_list, output, threads)

    
if __name__ == '__main__':
    main()
    
