import os
import sys
import glob
sys.path.append("../scripts") #añade en el path lo que hay en la carpeta scripts
from salmon import *
from hisat2 import *
from STAR import *
from kallisto import *



def module_map (path_reference, reference_genome, index, path_reads, output, threads, gtf, kmers, software): #two inputs path reads and software  ??
 
    lines = open(path_reads).readlines() #read the txt with the information of the reads 
    reads_list = [] #empty list
    for i in lines: #for every sample
        reads_list.append(i.strip().split(';')) #split the different camps

    print("Reads in main", reads_list)
      
    if software == "HISAT2":
        hisat2_index(path_reference, reference_genome, index)
        hisat2_mapping(path_reference, index, reads_list, output, threads, gtf)
    elif software == "STAR":
        star_index(path_reference, reference_genome, index, threads, gtf)
        star_mapping(path_reference, index, reads_list, output, threads)
    elif software == "Kallisto":
        kallisto_index(path_reference, reference_genome, index, kmers)
        kallisto_quant(path_reference, index, reads_list, output, threads)
    elif software =="Salmon":
        salmon_index(path_reference, reference_genome, index)
        salmon_quant(path_reference, index, reads_list, output, threads)



def main():

    parser = argparse.ArgumentParser (description = 'Map module with different softwares')

    parser.add_argument('-p', '--path_reference', help = 'Directory where the reference genome is, it will be the index directory too ', required = "TRUE")
    parser.add_argument('-g', '--reference_genome', help = 'Name of the file containing the reference genome', required = "TRUE")
    parser.add_argument('-i', '--index', help= 'Name of the index files. They will be created in the same folder of the reference genome', required = "TRUE")
    parser.add_argument('-r', '--path_reads', help='Path to your txt with the format: $PATH/sample_name;$PATH/read1;$PATH/read or $PATH/sample_name;$PATH/read1. IT IS VERY IMRTANT TO MANTAIN THE ORDER OF THE ELEMENTS AND THAT THEY ARE SEPARATED BY ;', required = "TRUE")
    parser.add_argument('-o', '--output', help= 'Path where the outputs will be created', required = "TRUE")
    parser.add_argument('-t', '--threads', help = 'Choose how many threads you want to use to execute HISAT2')
    parser.add_argument ('-b', '--gtf', help = 'Path to the file with annotated transcripts in the standard GTF format.')
    parser.add_argument ('-k', '--kmers', type = int, help = 'Desired k-mers length when indexing. Default: 31', default=31)
    parser.add_argument ('-s', '--software', choices = ["HISAT2","STAR","Salmon","Kallisto"], required = "TRUE")

    args = parser.parse_args()  # interprets the arguments

    path_reference = args.path_reference
    reference_genome = args.reference_genome
    index = args.index
    path_reads = args.path_reads
    output = args.output
    threads = args.threads
    kmers = args.kmers
    gtf = args.gtf
    software = args.software

    #module_map(path_reference, reference_genome, index, path_reads, output, threads, gtf, kmers, software)

    parser.set_defaults(func=rnaseqpip.modules.map.module_map)

    
if __name__ == '__main__':
    main()