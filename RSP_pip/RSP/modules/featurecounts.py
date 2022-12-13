import os
from RSP.scripts import * 

def run_featurecounts(path_results, path_sam, threads, gtf):

    ## checks before running FeatureCounts
    
    check_sam = path_sam.endswith("sorted.bam")

    if check_sam.endswith(".sam"): #check if input is a sam 
        
        sam_to_bam(path_results, path_sam, threads, gtf) #converting sam to bam 
        bam_to_sorted_bam(path_results, path_bam, threads, gtf) #converting bam to sorted bam 

        print("Converting SAM to BAM and BAM to sorted.BAM")
        
    elif check_sam == False:
        
        bam_to_sorted_bam(path_results, path_bam, threads, gtf) #converting bam to sorted bam 
        print("Converting BAM to sorted.BAM")
        
    else:
        print("The input was already a sorted BAM")
        
        
    ## Create counts matrix
    
    featurecounts(path_results, sorted_bam_name, threads, gtf) #call feature counts
    
   
   
if __name__ == '__main__':
    main()
    
