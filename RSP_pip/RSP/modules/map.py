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
from RSP.scripts import STAR_caller    
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

    ## multimapping:
    if options.no_multiMapping:
        multimapping = False
    else:
        multimapping = True
    


    ################################################
    ## Let's go
    ################################################
    HCGB_aes.pipeline_header('RSP')
    HCGB_aes.boxymcboxface("Mapping samples")
    print ("--------- Starting Process ---------")
    HCGB_time.print_time()
    

    ## HCGB_time.timestamp
    start_time_partial = HCGB_time.timestamp(start_time_total)
    
    ################################################
    ## Check reference genome, annotations, CDS
    ################################################
    print ('+ Getting reference genome information provided... ')

    # Check genome fasta file provided exists
    if HCGB_files.is_non_zero_file(options.ref_genome):
        options.ref_genome = os.path.abspath(options.ref_genome)
        print ('+ Reference genome provided:... ')
        print('\t' + options.ref_genome)
    else:
        print("Fasta reference genome provided does not exist")
        exit()
    # Check genome annotation file provided exists
    if HCGB_files.is_non_zero_file(options.ref_annot):
        options.ref_annot = os.path.abspath(options.ref_annot)
        print ('+ Reference genome annotation provided:... ')
        print('\t' + options.ref_annot)
    else:
        print("Reference genome annotation provided does not exist")
        exit()

    ## if desired, provide a shared folder with pre-computed indexed files
    if (options.index_folder):
        path_reference = os.path.abspath(options.index_folder)
    else:
        path_reference = os.path.abspath(options.ref_folder)

    print ('+ Reference folder to store index files:... ')
    HCGB_files.create_folder(path_reference)
    print('\t' + path_reference) 
     
    reference_genome = options.ref_genome ## reference_genome
    index_ref_name = options.ref_name ## index_name
    
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
                    index_ref_name + '_hisat2', threads=options.threads, 
                    extra_index=options.extra_index, index_folder=options.index_folder, Debug=Debug)
        
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
                    index_ref_name + '_salmon', threads=options.threads, 
                    extra_index=options.extra_index, index_folder=options.index_folder, Debug=Debug)
        
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
                    index_ref_name + '_kallisto', threads=options.threads, 
                    extra_index=options.extra_index, index_folder=options.index_folder, Debug=Debug)

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
        abs_path_index = check_index("star", path_reference, reference_genome, 
                    index_ref_name + '_star', threads=options.threads, 
                    extra_index=options.extra_index, index_folder=options.index_folder, limitGenomeGenerateRAM=options.limitGenomeGenerateRAM,
                    Debug=Debug)
        
        # use default if not provided
        map_params_star = {
            'index':abs_path_index
        }
        map_params["star"] = map_params_star


   ## debug message
    if (Debug):
        print (colored("**DEBUG: map_params **", 'yellow'))
        print(map_params)


    ## HCGB_time.timestamp
    start_time_partial = HCGB_time.timestamp(start_time_partial)
   

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
    
    ## send for each sample and each software specified unless STAR

    ## if star, it is better to load genome, map all samples and then remove genome from memory
    if "star" in options.soft_name:
        print(options.soft_name)    
        options.soft_name = options.soft_name.remove("star")
        print(options.soft_name)

        ## Create call for STAR only
        mapReads_module_STAR(options, pd_samples_retrieved, outdir_dict, Debug, 
                    max_workers_int, threads_job, start_time_partial, outdir, multimapping, map_params["star"]["index"])
    




    else:

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
def check_index(soft_name, path_reference, reference_genome, index_ref_name, threads, extra_index, index_folder, limitGenomeGenerateRAM, Debug):
    
    """Checks the index genome folder for each software 
    
    :param soft_name: Software name to use in the mapping
    :param path_reference: Main folder to store the results
    :param index_ref_name: Name to use to create index
    :param threads: number of threads to do the computation
    :param reference_genome: path to the genome directory
    :param extra_index: Addtional string options to include in the index call
    :param index_folder: Folder provided with indexed files.
    :param limitGenomeGenerateRAM: limit RAM bytes to be used in the computation
    :param Debug: Print debugging messages or not

    :type soft_name: string
    :type path_reference: string
    :type index_ref_name: string
    :type threads: int 
    :type reference_genome: string
    :type extra_index: string
    :type index_folder: string
    :type limitGenomeGenerateRAM: int
    :type Debug: boolean

    :returns: genomeDir
    """

    if (index_folder): # folder is provided containing pre-computed results for a single software: eg. STAR
        print("+ Check index folder provided: ")
        print(index_folder)
        genomeDir=index_folder
    else:
        path_reference = HCGB_files.create_subfolder(soft_name, path_reference)
        genomeDir=path_reference
        ## /path/to/any/folder/
        ## -------------------- star
        ## -------------------- hisat2


    ## check if previously mapped and succeeded
    filename_stamp = path_reference + '/.success_index_' + soft_name
    
    if os.path.isfile(filename_stamp):
        stamp = HCGB_time.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s]" %(stamp, 'index ', soft_name), 'yellow')
        print("+ Let's check everything is alright...'")
    
    if soft_name=="hisat2":
        code_returned= hisat2.hisat2_index(path_reference, 
                                        reference_genome, index_ref_name, threads, Debug)
        
    if soft_name=="salmon":
        # Fix
        #salmon_index(path_reference, reference_genome, index)
        code_returned= hisat2.hisat2_index(path_reference, 
                                        reference_genome, index_ref_name, threads, Debug)
        
    if soft_name=="kallisto":
        # Fix
        #kallisto_index(path_reference, reference_genome, index, kmers)
        code_returned= hisat2.hisat2_index(path_reference, 
                                        reference_genome, index_ref_name, threads, Debug)
        
    if soft_name=="star":
        code_returned= STAR_caller.check_index(path_reference, 
                                        reference_genome, index_ref_name, threads, extra_index, limitGenomeGenerateRAM, Debug)
        
    if code_returned:
        HCGB_time.print_time_stamp(filename_stamp)
        
    else:
        print ('** Software %s failed to index genome provided...' %soft_name)

    return genomeDir


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



#########################################
def mapReads_module_STAR(options, pd_samples_retrieved, outdir_dict, Debug, 
                    max_workers_int, threads_job, start_time_partial, outdir, multimapping, genomeDir):
    
    """Organizes the mapping of the samples, executed in parallel.

    First, checks if the files needed to do the mapping (and posterior
    classification of the reads) have been provided by the user:
    fasta sequence + annotation or STAR index directory. 

    Then, sends the mapping in parallel for each sample calling mapReads_caller().    
    
    Finally, generate the MultiQC report of the mapping for each sample.

    Original code from XICRA code, modified and updated for RSP.
    
    :param options: input parameters introduced by the user. See XICRA biotype -h.
    :param pd_samples_retrieved: data frame with the information of the samples
    :param outdir_dict: dictionary with the names of the samples and their files
    :param Debug: show extra information of the process
    :param max_workers_int: number of workers for each thread
    :param threads_job: number of threads to do the computation
    :param start_time_partial: time of the beggining of the process
    :param outdir: directory to store the results

    :type Debug: boolean
    :type max_workers_int: int
    :type threads_job: int 
    :type start_time_partial: int
    :type outdir: string

    :returns: None
    """

    # Group dataframe by sample name
    sample_frame = pd_samples_retrieved.groupby(["new_name"])
    
    ## options
    STAR_exe = set_config.get_exe("STAR", Debug=Debug)
    cwd_folder = os.path.abspath("./")
    folder=files_functions.create_subfolder('STAR_files', cwd_folder)

    ## For many samples it will have to load genome index in memory every time.
    ## For a unique sample it will not matter. Take care genome might stay in memory.
    ## Use before loop option LoadAndExit and then:
        ## in loop
        ## Use option LoadAndKeep, set shared memory > 30 Gb
    ## when finished loop Remove memory        
    
    ## remove previous reference genome from memory
    print ("+ Remove genome in memory from previous call... (if any)")
    #STAR_caller.remove_Genome(STAR_exe, options.genomeDir, folder, options.threads)
    
    ## load reference genome
    STAR_caller.load_Genome(folder, STAR_exe, genomeDir, options.threads)

    ## functions.time_functions.timestamp
    start_time_partial = HCGB_time.timestamp(start_time_partial)
    
    print ("+ Mapping sequencing reads for each sample retrieved...")

    ## send for each sample
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
        commandsSent = { executor.submit(mapReads_caller_STAR, sorted(cluster["sample"].tolist()), 
                                         outdir_dict[name], name, threads_job, STAR_exe, 
                                         options.genomeDir, options.limitRAM, Debug, multimapping): name for name, cluster in sample_frame }

        for cmd2 in concurrent.futures.as_completed(commandsSent):
            details = commandsSent[cmd2]
            try:
                data = cmd2.result()
            except Exception as exc:
                print ('***ERROR:')
                print (cmd2)
                print('%r generated an exception: %s' % (details, exc))

    print ("\n\n+ Mapping reads has finished...")
    
    ## functions.time_functions.timestamp
    start_time_partial = HCGB_time.timestamp(start_time_partial)

    ## remove reference genome from memory
    STAR_caller.remove_Genome(STAR_exe, genomeDir, folder, options.threads)
    
    ## functions.time_functions.timestamp
    start_time_partial = HCGB_time.timestamp(start_time_partial)

    ## retrieve mapping files
    if options.detached:
        input_dir = outdir
    else:
        input_dir = os.path.abspath(options.input)
    
    results_SampleParser = sampleParser.files.get_files(options, input_dir, "map", ["Aligned.sortedByCoord.out.bam"], options.debug, bam=True)
    del results_SampleParser['dirname']
    del results_SampleParser['ext']
    del results_SampleParser['tag']
    #del results_SampleParser['new_name']

    results_SampleParser = results_SampleParser.set_index('name')
    mapping_results = results_SampleParser.to_dict()['sample']

    ## Create mapping report    
    if (options.skip_report):
        print ("+ No report generation...")
    else:
        print ("\n+ Generating a report using MultiQC module.")
        outdir_report = HCGB_files.create_subfolder("report", outdir)
        map_outdir_report = HCGB_files.create_subfolder("map", outdir)

        ## get subdirs generated and call multiQC report module
        givenList = []
        print ("+ Detail information for each sample could be identified in separate folders:")
        
        ## call multiQC report module
        givenList = [ v for v in outdir_dict.values() ]
        my_outdir_list = set(givenList)

        ## debug message
        if (Debug):
            print (colored("\n**DEBUG: my_outdir_list for multiqc report **", 'yellow'))
            print (my_outdir_list)
            print ("\n")
        
        map_report = HCGB_files.create_subfolder("STAR", map_outdir_report)
        multiQC_report.multiQC_module_call(my_outdir_list, "STAR", map_report,"-dd 2")
        print ('\n+ A summary HTML report of each sample is generated in folder: %s' %map_report)

    return(start_time_partial, mapping_results)

#################################
def mapReads_caller_STAR(files, folder, name, threads, STAR_exe, genomeDir, limitRAM_option, Debug, multimapping):
    """Mapping of a given sample with STAR

    First, checks if the trimmed unjoined files exist for the sample and also
    if the calculation has not been done previously. 

    Executes STAR and generate the BAM file of the sample with mapReads script. 
    
    :param files: unjoined trimmed files of the sample
    :param folder: sample folder to store the results
    :param name: sample name
    :param threads: number of threads to do the computation
    :param STAR_exe: check the STAR software is available
    :param genomeDir: path to the genome directory to do the mappig
    :param limitRAM_option: limit RAM bytes to be used in the computation
    :param Debug: show extra information of the process
    :param multimapping: Flag to say whether to use multimapping reads or not
    
    :type folder: string
    :type name: string
    :type threads: int 
    :type start_exe: boolean
    :type genomeDir: string
    :type limitRAM_option: int
    :type Debug: boolean
    :type multimapping: boolean

    :returns: None
    """

    ## check if previously joined and succeeded
    filename_stamp = folder + '/.success'
    if os.path.isfile(filename_stamp):
        stamp = time_functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'STAR'), 'yellow'))
    else:
        ##
        if Debug:
            print ("\n** DEBUG: mapReads_caller options **\n")
            print ("folder: " + folder) 
            print ("name: " + name)
            print ("threads: " + str(threads))
            print ("STAR_exe: " + STAR_exe) 
            print ("genomeDir: " + genomeDir) 
            print ("limitRAM_option: " + str(limitRAM_option))
            print ("files: ")
            print (files)
            
        # Call STAR
        code_returned = STAR_caller.mapReads("LoadAndKeep", files, folder, name, STAR_exe, genomeDir, limitRAM_option, threads, Debug, multimapping)
        
        if (code_returned):
            time_functions.print_time_stamp(filename_stamp)
        else:
            print ("+ Mapping sample %s failed..." %name)
    
    ## return results
    return()
    