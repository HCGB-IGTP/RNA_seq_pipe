import os 

def featurecounts (path_results, sorted_bam_name, threads, gtf):

    featurecounts = "featureCounts -p -O -T " + threads + " -a " + gtf + " -o " + path_results + " " + sorted_bam_name 
    print(featurecounts)
    os.system(featurecounts) 

