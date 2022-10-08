'''
Created on Sep 28, 2022

@author: mireia
'''

import os
import sys



# os.mkdir("NEW FOLDER") #use os. before the function 

# path = "/home/mireia/IGTP/22"

# os.chdir(path) #change directory 

def mapper_histat():
    print()


def main():

    ## take arguments
    # folder with reads
    # output results
    # reference folder
    # reference name

    print(sys.argv) #LISTA
    
    
    
    
    
    
    exit()
    
    reads = os.listdir() #saves files in current directory as a variable 
    
    reads_1 = [] #list with reads forward
    reads_2 = [] #list with reads reverse
    
    for r in reads: 
        if r.endswith('.fastq'): #extract the reads from all the files in the directory
            if "1" in r: #reads -1
                reads_1.append(r)
            else: # reads -2
                reads_2.append(r) 
        
    
    
    for i,j in zip (reads_1, reads_2): #for every couple of reads we run the hisat2 command
        forward = i #-1 read
        reverse = j  # -2 read
        name=i+j+".sam" #name of output sam 
        name=name.replace(".fastq", "") #aesthetic name 
        
        # call mapper_histat
        #sam_output = mapper_histat(ref, read1, read2, output_path_name)
        mapping = "hisat2 -x chr22_index -1"+forward+" -2"+reverse+" -S"+name #hisat2 mapping command 
        print(mapping)
        #os.system(mapping)


######
if __name__== "__main__":
    main()





## os.path()
# os.path.join(folder, name_index)
# os.path.join(results, name_sam)





















