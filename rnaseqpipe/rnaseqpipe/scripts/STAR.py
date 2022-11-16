import os
import argparse



### INDEXING FUNCTION ###

def star_index(path_reference, reference_genome, index, threads, gtf):
    
    reference_abs_path = os.path.join(path_reference, reference_genome) #path reference genome
    
    """
    #Check if reference genome and gtf are compressed
    
    fasta = (".fa", ".fasta")
    if reference_abs_path.endswith(fasta) == False:  
        print("Your genome must be descompressed")
        gzip_genome = "gzip -k -d " + reference_abs_path
        os.system(gunzip_genome)
    
    if gtf.endswith(".gtf") == False:
        print("Your GTF file must be descompressed")
        gunzip_gtf = "gunzip -k -d " + gtf
        os.system(gunzip_gtf)
     
     """
             
    index_abs_path = os.path.join(path_reference, index) #path index
    
    #check if the index exist by searching the folder called as index
        
    if os.path.isdir(index_abs_path) == True: 
    
        print("The genome is already indexed")
       
        print("Path to the index:", index_abs_path)


        else:
            indexing = "STAR --runThreadN " + threads + " --runMode genomeGenerate --genomeDir " + index_abs_path + " --genomeFastaFiles " + reference_abs_path + " --sjdbGTFfile " + gtf  #index the
            print(indexing)
            print("Path to the index:", index_abs_path)
    
    
    
    #if the folder didn't exist we create it and index the genome
    
    else: 
        os.mkdir(index_abs_path) #create directory where the indexes will be stored.
        
    indexing = "STAR --runThreadN " + threads + " --runMode genomeGenerate --genomeDir " + index_abs_path + " --genomeFastaFiles " + reference_abs_path + " --sjdbGTFfile " + gtf
    print(indexing)
    
    #os.system(star) 
    
def star_mapping(path_reference, index, reads_list, output, threads):
    
    path_index = os.path.join(path_reference, index) #path to the genome index
    
    #loop to check number of reads
         
    for i in reads_list: #for every list inside reads_list
        if len(i) == 3: #if there are three elements it's a pair-end analysis
            print("Pair-end analysis")
            read1 = i[1] #reads forward 
            read2 = i[2] #read reverse
            outputs_name = i[0] #variable that will go through the functions, sample name
            folder_results = os.path.join(output, outputs_name) #path to results folder + sample
            os.mkdir(folder_results)
            
            path_results = os.path.join(folder_results, outputs_name)        
            
            print ("Sample name", outputs_name, "Read forward:", read1, "Read reverse:", read2)

            mapping = "STAR --genomeDir " + path_index + " --runThreadN " + threads + " --readFilesIn " + read1 + " " + read2 + " --outFileNamePrefix " + path_results 
            #two fq separated by space --> paired end analysis
                        
            print(mapping)
            #os.system(mapping)
            #SamToBam(path_results, path_results, threads, gtf) 
        else:
            print("Single-end analysis")
            single_read = i[1] #single end analysis
            outputs_name = i[0] #variable that will go through the functions 
            folder_results = os.path.join(output, outputs_name) #path to results folder + sample
            os.mkdir(folder_results)
            
            path_results = os.path.join(folder_results, outputs_name)
            print ("Sample name", outputs_name, "Single read:", single_read)
            
            mapping = "STAR --genomeDir " + path_index + " --runThreadN " + threads + " --readFilesIn " + read1 + " --outFileNamePrefix " + path_results
            #one fq --> single end analysis 
            
            print(mapping)
            #os.system(mapping)
            #SamToBam(path_results, path_results, threads, gtf) 
           

### MAIN FUNCTION ###  
    
def main():
    
    parser = argparse.ArgumentParser (description = 'Indexation and mapping using HISAT2') #all arguments are mandatory 
    parser.add_argument('-p', '--path_reference', help = 'Directory where the reference genome is, it will be the index directory too ', required = "TRUE")
    parser.add_argument('-g', '--reference_genome', help = 'Name of the file containing the reference genome', required = "TRUE")
    parser.add_argument('-i', '--index', help= 'It will create a folder called as your argument', required = "TRUE")
    parser.add_argument('-r', '--path_reads', help='Path to your txt with the format: $PATH/sample_name;$PATH/read1;$PATH/read or $PATH/sample_name;$PATH/read1. IT IS VERY IMRTANT TO MANTAIN THE ORDER OF THE ELEMENTS AND THAT THEY ARE SEPARATED BY ;', required = "TRUE")
    parser.add_argument('-o', '--output', help= 'Path where the outputs will be created', required = "TRUE")
    parser.add_argument('-t', '--threads', help = 'Choose how many threads you want to use to execute HISAT2')
    parser.add_argument('-f', '--function', choices = ['Indexing', 'Mapping'], help = 'You can choose to execute the indexing function, the quantifying function or, by default, both. Possible choices for this argument: Indexing / Quantifying')
    parser.add_argument ('-b', '--gtf', help = 'Path to the file with annotated transcripts in the standard GTF format.')

        
    args = parser.parse_args() #interprets the arguments 
    
    # python3 /home/mireia/IGTP/STAR/star.py -p /home/mireia/IGTP/STAR -g /home/mireia/IGTP/salmon/cds.fa.gz -i index -r /home/mireia/IGTP/reads/reads2 -o results -f Indexing -b /home/mireia/IGTP/STAR/reference.gtf.gz
    
    ## inputs for other functions (can't be called as args.)
    
    path_reference = args.path_reference
    reference_genome = args.reference_genome
    index = args.index
    path_reads = args.path_reads
    output = args.output 
    threads = args.threads
    function = args.function
    gtf = args.gtf
    
  
    lines = open(path_reads).readlines() #read the txt with the information of the reads 
    reads_list = [] #empty list
    for i in lines: #for every sample
        reads_list.append(i.strip().split(';')) #split the different camps
        
    # lines.close() #close document
        
    print("Reads in main", reads_list)
    
    
    #function calling
    
    if function == "Indexing": #if in --function they choose indexing just the index function will run
        star_index(path_reference, reference_genome, index, threads, gtf)
    elif function == "Quantifying":#if in --function they choose mapping just the mapping function will run
        star_mapping(path_reference, index, reads_list, output, threads, gtf)
    else: #if they don't specify the function, both will run
        star_index(path_reference, reference_genome, index, threads, gtf)
        star_mapping(path_reference, index, reads_list, output, threads, gtf)

    
if __name__ == '__main__':
    main()
