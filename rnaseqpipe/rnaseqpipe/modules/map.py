import os
import sys
import glob
from rnaseqpipe.scripts import * 

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


