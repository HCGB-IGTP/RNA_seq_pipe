import os


####SAM TO BAM FUNCTION########################################################################################################################

def sam_to_bam(path_results, path_sam, threads, gtf):

    print("Converting SAM to BAM")

    path_bam = path_results + ".bam" #safe bam into results folder
    sam_to_bam = "samtools view -bS "+ path_sam + " > " + path_bam #samtools sam to bam command
    print(sam_to_bam)
    
    os.system(sam_to_bam)


####BAM TO SORTED_BAM FUNCTION########################################################################################################################


def bam_to_sorted_bam(path_results, path_bam, threads, gtf):
    
    print("Converting BAM to Sorted BAM")
    sorted_bam_name = path_results +".sorted.bam" #safe sorted bam into results folder
    bam_to_sorted_bam = "samtools sort "+ path_bam + " -o " + sorted_bam_name #samtools bam to sorted bam command
    print(bam_to_sorted_bam)
    
    os.system(bam_to_sorted_bam)
    



if __name__ == '__main__':
    main()
    
