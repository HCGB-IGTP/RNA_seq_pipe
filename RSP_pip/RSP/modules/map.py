#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                      ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
"""
Maps sequence reads to reference genome
"""
## import useful modules
import os
import sys
import re
import time
from io import open
import shutil
import concurrent.futures
from termcolor import colored

## import my modules
from RSP.scripts import multiQC_report
from RSP.scripts import hisat2
from RSP.scripts import kallisto
from RSP.scripts import STAR
from RSP.scripts import salmon

from RSP.config import set_config
from RSP.data import data_files
from RSP import __version__ as pipeline_version
from RSP.modules import help_RSP
from RSP.modules import qc

from HCGB import sampleParser
import HCGB.functions.info_functions as HCGB_info
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.main_functions as HCGB_main

###############################################3
def run_map(options):
    """Main function of the module, organizes the mapping process.
    
    There are several alternative softwares to map the reads.
    
    Some of them might have special configuration details. 
    
    :param options: input parameters introduced by the user. See RSP map -h.

    :returns: None
    """
    ## init time
    start_time_total = time.time()

    ##################################
    ### show help messages if desired    
    ##################################
    if (options.help_format):
        ## help_format option
        help_RSP.help_fastq_format()
    elif (options.help_project):
        ## information for project
        help_RSP.project_help()
        exit()
    elif (options.help_multiqc):
        ## information for Multiqc
        help_RSP.multiqc_help()
        exit()
        
    ################################################
    ## Set defaults
    ################################################
    global Debug
    if (options.debug):
        Debug = True
    else:
        Debug = False
        
    ### set as default paired_end mode
    if (options.single_end):
        options.pair = False
    else:
        options.pair = True
    
    ## absolute path for in & out
    input_dir = os.path.abspath(options.input)
    outdir=""

    ## set mode: project/detached
    if (options.detached):
        outdir = os.path.abspath(options.output_folder)
        options.project = False
    else:
        options.project = True
        outdir = input_dir        

    ################################################
    ## Let's go
    ################################################
    HCGB_aes.pipeline_header('RSP')
    HCGB_aes.boxymcboxface("Mapping samples")
    print ("--------- Starting Process ---------")
    HCGB_time.print_time()
    
    ################################################
    ## Check reference genome, annotations, CDS
    ################################################
    print ('+ Getting reference genome information provided... ')

    if HCGB_files.is_non_zero_file(options.ref_genome):
        options.ref_genome = os.path.abspath(options.ref_genome)
        print ('+ Reference genome provided:... ')
        print('\t' + options.ref_genome)
    else:
        print("Fasta reference genome provided does not exist")
        exit()
        
    if HCGB_files.is_non_zero_file(options.ref_annot):
        options.ref_annot = os.path.abspath(options.ref_annot)
        print ('+ Reference genome annotation provided:... ')
        print('\t' + options.ref_annot)
    else:
        print("Reference genome annotation provided does not exist")
        exit()

    ## if desired, use a different folder to store index files
    if (options.index_folder):
        path_reference = os.path.abspath(options.index_folder)
        HCGB_files.create_folder(path_reference)
    else:
        path_reference = os.path.dirname(options.ref_genome)

    print ('+ Reference index_folder to store index files:... ')
    print('\t' + path_reference) 
     
    reference_genome = options.ref_genome ## reference_genome
    index_ref_name = options.ref_index ## index_name
    
    print ('+ Reference index name provided: ' + index_ref_name)
    print('') 
    
                                                 
     ## debug message
    if (Debug):
        print (colored("**DEBUG: path_reference **", 'yellow'))
        print(path_reference)
        
        print (colored("**DEBUG: index_ref_name **", 'yellow'))
        print(index_ref_name)
    
        print (colored("**DEBUG: reference_genome **", 'yellow'))
        print(reference_genome)

    ################################################
    ## Check parameters for each software
    ################################################
    map_params = {}
    
    print ("--------- Check paramaters ---------\n")

    
    ##----------------------------
    ## set hisat2 params
    ##----------------------------
    map_params_hisat2 = {}
    if "hisat2" in options.soft_name:
        
        print ('+ Setting HISAT2 parameters: ')
    
        print ('+ Checking index for HISAT2: ')
        ## Prepare or check index
        abs_path_index =check_index("hisat2", path_reference, reference_genome, 
                    index_ref_name + '_hisat2', threads=options.threads, Debug=Debug)
        
        # use default if not provided
        map_params_hisat2 = {
            'index':abs_path_index
        }
        map_params["hisat2"] = map_params_hisat2
        
        
    ##----------------------------
    ## set salmon params
    ##----------------------------
    map_params_salmon = {}
    if "salmon" in options.soft_name:
        print ('+ Setting salmon parameters: ')
    
        print ('+ Checking index for salmon: ')
        ## Prepare or check index
        abs_path_index =check_index("salmon", path_reference, reference_genome, 
                    index_ref_name + '_salmon', threads=options.threads, Debug=Debug)
        
        # use default if not provided
        map_params_salmon = {
            'index':abs_path_index
        }
        map_params["salmon"] = map_params_salmon
        
    ##----------------------------
    ## set kallisto params
    ##----------------------------
    map_params_kallisto = {}
    if "kallisto" in options.soft_name:
        print ('+ Setting kallisto parameters: ')
    
        print ('+ Checking index for kallisto: ')
        ## Prepare or check index
        abs_path_index =check_index("kallisto", path_reference, reference_genome, 
                    index_ref_name + '_kallisto', threads=options.threads, Debug=Debug)
        
        # use default if not provided
        map_params_kallisto = {
            'index':abs_path_index
        }
        map_params["kallisto"] = map_params_kallisto
        
    ##----------------------------
    ## set star params
    ##----------------------------
    map_params_star = {}
    if "star" in options.soft_name:
        print ('+ Setting star parameters: ')
    
        print ('+ Checking index for star: ')
        ## Prepare or check index
        abs_path_index =check_index("star", path_reference, reference_genome, 
                    index_ref_name + '_star', threads=options.threads, Debug=Debug)
        
        # use default if not provided
        map_params_star = {
            'index':abs_path_index
        }
        map_params["star"] = map_params_star


    ################################################
    ## get files
    ################################################
    print ("\n--------- Get reads to map ---------\n")

    print ('+ Getting files from input folder... ')
    
    ## get files: use trimmed reads or not trimmed if specified
    if options.pair:
        options.pair = False ## set paired-end to false for further prepocessing
        if options.noTrim:
            print ('+ Mode: fastq.\n+ Extension: ')
            print ("[ fastq, fq, fastq.gz, fq.gz ]\n")
            pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "fastq", ["fastq", "fq", "fastq.gz", "fq.gz"], options.debug)
        else:
            print ('+ Mode: join.\n+ Extension: ')
            print ("[_trim.*fastq]\n")
            pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "trim", ['_trim'], options.debug)
    else:
        if options.noTrim:
            print ('+ Mode: fastq.\n+ Extension: ')
            print ("[ fastq, fq, fastq.gz, fq.gz ]\n")
            pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "fastq", ["fastq", "fq", "fastq.gz", "fq.gz"], options.debug)
        else:
            print ('+ Mode: trim.\n+ Extension: ')
            print ("[_trim.fastq]\n")
            pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "trim", ['_trim'], options.debug)
    
    ## debug message
    if (Debug):
        print (colored("**DEBUG: pd_samples_retrieve **", 'yellow'))
        print (pd_samples_retrieved)

    ## generate output folder, if necessary
    print ("\n+ Create output folder(s):")
    if not options.project:
        HCGB_files.create_folder(outdir)
    ## for samples
    outdir_dict = HCGB_files.outdir_project(outdir, options.project, pd_samples_retrieved, "map", options.debug)
    
    ## optimize threads
    name_list = set(pd_samples_retrieved["new_name"].tolist())
    threads_job = HCGB_main.optimize_threads(options.threads, len(name_list)) ## threads optimization
    max_workers_int = int(options.threads/threads_job)

    ## debug message
    if (Debug):
        print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
        print (colored("**DEBUG: max_workers " +  str(max_workers_int) + " **", 'yellow'))
        print (colored("**DEBUG: cpu_here " +  str(threads_job) + " **", 'yellow'))
    ##########################

    ##########################
    ## Let's start processing
    ##########################
    print ("+ Mapping reads for each sample and software specified...")    
    
    # Group dataframe by sample name
    sample_frame = pd_samples_retrieved.groupby(["new_name"])
    
    ## send for each sample
    ## use software specified

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
            ## sample_name, path_reference, reference_genome, index_name, reads_list, main_output, threads, parameters, software_list, Debug
            commandsSent = { executor.submit(module_map, 
                                             name, ## sample_name,
                                             path_reference, ## path_reference, 
                                             reference_genome, ## reference_genome,
                                             index_ref_name, ## index_name,
                                             sorted(cluster["sample"].tolist()), ## reads_list 
                                             outdir_dict[name], ## main_output
                                             threads_job, ## threads 
                                             map_params, ## dictionary parameters
                                             options.soft_name, ## software list
                                             Debug): name for name, cluster in sample_frame }
    
            for cmd2 in concurrent.futures.as_completed(commandsSent):
                details = commandsSent[cmd2]
                try:
                    data = cmd2.result()
                except Exception as exc:
                    print ('***ERROR:')
                    print (cmd2)
                    print('%r generated an exception: %s' % (details, exc))


###############################################3
def check_index(soft_name, path_reference, reference_genome, index_ref_name, threads, Debug):
    
    ## check if previously mapped and succeeded
    filename_stamp = path_reference + '/.success_index_' + soft_name
    if os.path.isfile(filename_stamp):
        stamp = HCGB_time.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s]" %(stamp, 'index hisat2'), 'yellow'))
        
        return (os.path.join(path_reference, index_ref_name))
        
    else:
        if soft_name=="hisat2":
            code_returned= hisat2.hisat2_index(path_reference, 
                                          reference_genome, index_ref_name, threads, Debug)
        
        if soft_name=="salmon":
            #salmon_index(path_reference, reference_genome, index)
            code_returned= hisat2.hisat2_index(path_reference, 
                                          reference_genome, index_ref_name, threads, Debug)
        
        if soft_name=="kallisto":
            #kallisto_index(path_reference, reference_genome, index, kmers)
            code_returned= hisat2.hisat2_index(path_reference, 
                                          reference_genome, index_ref_name, threads, Debug)
        
        if soft_name=="star":
            #star_index(path_reference, reference_genome, index, threads, gtf)
            code_returned= hisat2.hisat2_index(path_reference, 
                                          reference_genome, index_ref_name, threads, Debug)
        
        if code_returned:
            HCGB_time.print_time_stamp(filename_stamp)
            return (os.path.join(path_reference, index_ref_name))

        else:
            print ('** Sample %s failed...' %sample_name)

###############################################
def module_map (sample_name, path_reference, reference_genome, index_ref_name, reads_list, main_output, 
                threads, map_params, software_list, Debug): 
    
    ##-------------------------------------
    ## send commands for each software
    ##-------------------------------------
    
    ##-------------------------------------
    ## HISAT2
    ##-------------------------------------
    if "hisat2" in software_list:
        output = HCGB_files.create_subfolder("hisat2", main_output)
        
        ## TODO: Add additional parameters saved in dictionary map_params['hisat2']
        extra_params=""
        
        ## check if previously mapped and succeeded
        filename_stamp = output + '/.success_hisat2'
        if os.path.isfile(filename_stamp):
            stamp = HCGB_time.read_time_stamp(filename_stamp)
            print (colored("\tA previous command generated results on: %s [%s -- %s: %s]" %(stamp, sample_name, 'map', 'hisat2'), 'yellow'))
        else:
            ## create call to mapping
            code_returned= hisat2.hisat2_mapping(sample_name, map_params['hisat2']['index'], 
                                                 reads_list, output, threads, 
                                                 extra_params, Debug)
            if code_returned:
                HCGB_time.print_time_stamp(filename_stamp)
            else:
                print ('** Sample %s failed...' %sample_name)

    ##-------------------------------------
    ## star
    ##-------------------------------------
    if "star" in software_list:
        output = HCGB_files.create_subfolder("STAR", main_output)
        #star_index(path_reference, reference_genome, index, threads, gtf)
        #star_mapping(path_reference, index, reads_list, output, threads)
    
    ##-------------------------------------
    ## kallisto
    ##-------------------------------------
    if "kallisto" in software_list:
        output = HCGB_files.create_subfolder("kallisto", main_output)
        #kallisto_index(path_reference, reference_genome, index, kmers)
        #kallisto_quant(path_reference, index, reads_list, output, threads)
    
    ##-------------------------------------
    ## salmon
    ##-------------------------------------
    if "salmon" in software_list:
        output = HCGB_files.create_subfolder("salmon", main_output)
        #salmon_index(path_reference, reference_genome, index)
        #salmon_quant(path_reference, index, reads_list, output, threads)


