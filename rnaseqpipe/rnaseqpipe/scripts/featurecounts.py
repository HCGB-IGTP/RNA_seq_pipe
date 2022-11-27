import os 

def featurecounts (path_results, sorted_bam_name, threads, gtf):

    output = os.path.join(path_results, "Counts_Matrix_Output.txt") 
       
    featurecounts = "featureCounts -p -O -T " + threads + " -a " + gtf + " -o " + output + " " + sorted_bam_name 
    
    print(featurecounts)
    
    os.system(featurecounts) 

