import os 
import argparse

####INDEXING FUNCTION########################################################################################################################

def hisat2_index(path_reference, reference_genome, index):
               
    reference_abs_path = os.path.join(path_reference, reference_genome) #path reference genome
    index_abs_path = os.path.join(path_reference, index) #path index
    indexing = 'hisat2-build ' + reference_abs_path + ' ' + index_abs_path #bash command 
    
    print("Path to the index:", index_abs_path)
    
    #os.system(indexing)

####MAPPING FUNCTION#########################################################################################################################

def hisat2_mapping(path_reference, index, reads_list, output):
        
    #loop to check number of reads
    
    for i in reads_list:
        if len(i) == 2: #if there is two position in the line there is a pair-end analysis
            print("Pair-end analysis")
            read1 = i[0] #reads forward 
            read2 = i[1] #read reverse

            """ No se com posar els noms dels outputs pq vagen cambiant en cada read
            sam = read1 + read2 + ".sam" #output name
            sam = sam.replace(".fastq", "")
            """
            
            path_sam = os.path.join(output, "hola.sam") #sSOLES FUNCIONA QUAN LI DONES NOM A NOM PER A L'OUTPUT

            path_index = os.path.join(path_reference, index) #index path
            
            mapping = "hisat2 -x " + path_index + " -1 " + read1 + " -2 " + read2 + " -S " + path_sam #hisat2 paired end mapping command  
            #print(mapping)
            os.system(mapping)
            SamToBam(path_sam) #IT SHOULD TAKE SAM NAME AND OUTPUT (PATH TO RESULTS)
            
        else:
            print("Single-end analysis")
            single_read = i[0] #single end analysis
            """
            sam = single_read + ".sam"
            sam = sam.replace(".fastq", "")
            path_sam = os.path.join(output, sam)
            """
            path_index = os.path.join(path_reference, "hola")

            mapping = "hisat2 -x " + path_index + " -U " + single_read + " -S " + path_sam #hisat2 paired end mapping command  
            print(mapping)
            SamToBam(output, path_sam)


####SAM TO BAM FUNCTION########################################################################################################################

def SamToBam(path_sam):
    print("Converting SAM to BAM")
    bam_name = "eg2.bam"
    sam_to_bam = "samtools view -bS "+ path_sam + " > " + bam_name
    os.system(sam_to_bam)
    BamToSortedBam(bam_name)


####BAM TO SORTED BAM FUNCTION#################################################################################################################

def BamToSortedBam(bam_name):
    print("Converting BAM to Sorted BAM")
    bam_to_sorted_bam = "samtools sort "+ bam_name + " -o  eg2.sorted.bam"
    os.system(bam_to_sorted_bam)

####MAIN FUNCTION##############################################################################################################################

def main():
    
    parser = argparse.ArgumentParser (description = 'Indexation and mapping using HISAT2') #all arguments are mandatory 
    parser.add_argument('-p', '--path_reference', help = 'Directory where the reference genome is, it will be the index directory too ', required = "TRUE")
    parser.add_argument('-g', '--reference_genome', help = 'Name of the file containing the reference genome', required = "TRUE")
    parser.add_argument('-i', '--index', help= 'Name you want your indexes to have', required = "TRUE")
    parser.add_argument('-r', '--path_reads', help='Path to your txt with the reads names', required = "TRUE")
    parser.add_argument('-o', '--output', help= 'Path where the sam will be created', required = "TRUE")
    
    args = parser.parse_args() #interprets the arguments 
    
    # python3 hisat2.py -p /home/mireia/HISAT2/reference -g chr22.fa -i index -r /home/mireia/HISAT2/reads/reads -o /home/mireia/HISAT2/results
    
    ## inputs for other functions (can't be called as args.)
    
    path_reference = args.path_reference
    reference_genome = args.reference_genome
    index = args.index
    path_reads = args.path_reads
    output = args.output 
    
    with open(args.path_reads, 'r') as f: #from files creates list
        reads_list = [line.strip().split(';') for line in f] #arguments in line are separated by a tab
 
    #function calling
    
    hisat2_index(path_reference, reference_genome, index)
    hisat2_mapping(path_reference, index, reads_list, output)


    
if __name__ == '__main__':
    main()
    
    

