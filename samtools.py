import os

####SAM TO BAM FUNCTION########################################################################################################################

def SamToBam(path_results, path_sam):

    print("Converting SAM to BAM")

    path_bam = path_results + ".bam" #safe bam into results folder
    sam_to_bam = "samtools view -bS "+ path_sam + " > " + path_bam #samtools sam to bam command
    print(sam_to_bam)
    os.system(sam_to_bam)
    
    
    BamToSortedBam(path_results, path_bam)


####BAM TO SORTED BAM FUNCTION#################################################################################################################

def BamToSortedBam(path_results, path_bam):
    
    print("Converting BAM to Sorted BAM")
    sorted_bam_name = path_results +".sorted.bam" #safe sorted bam into results folder
    bam_to_sorted_bam = "samtools sort "+ path_bam + " -o " + sorted_bam_name #samtools bam to sorted bam command
    print(bam_to_sorted_bam)
    os.system(bam_to_sorted_bam)



if __name__ == '__main__':
    main()
    
