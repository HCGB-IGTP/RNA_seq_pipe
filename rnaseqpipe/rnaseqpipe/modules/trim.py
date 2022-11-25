#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                      ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
"""
Trims sequence adapters within fastq reads.
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
from rnaseqpipe.scripts import multiQC_report
from rnaseqpipe.scripts import cutadapt_caller
from rnaseqpipe.scripts import trimmomatic_call

from rnaseqpipe.config import set_config
from rnaseqpipe.data import data_files
from rnaseqpipe import __version__ as pipeline_version
from rnaseqpipe.modules import help_RSP
from rnaseqpipe.modules import qc

from HCGB import sampleParser
import HCGB.functions.info_functions as HCGB_info
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.main_functions as HCGB_main


##############################################
def run_trim(options):
    """Main function of the module, organizes the trimming process.

    First, checks if the adapter(s) sequence(s) have been provided by the user:
    - adapters_a or adapters_A
    If there is no adapter sequence provided the process will be stopped.  

    If the adapters have been introduced, it calls cutadapt_caller() for each sample in parallel.
    Finally, generates a report using MultiQC module if desired.
    
    :param options: input parameters introduced by the user. See rnaseqpipe trim -h.

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
    elif (options.help_trimm_adapters):
        ## help on trimm adapters
        help_RSP.print_help_adapters()
        exit()
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
    HCGB_aes.pipeline_header('rnaseqpipe')
    HCGB_aes.boxymcboxface("Trimming samples")
    print ("--------- Starting Process ---------")
    HCGB_time.print_time()

    ##########################
    ## Software:
    ## Trimmomatic or Cutadapt    
    ##########################
    if options.software == "trimmomatic":
        
        #-----------------------------
        ## Trimmomatic parameters:
        #-----------------------------
        
        # Trimming adapters
        if (options.adapters):
            # Adapter file provided
            options.adapters = os.path.abspath(options.adapters)
            print("\t- Adapters file provided...")
        else:
            # Get default adpaters file
            print("\t- Default Trimmomatic adapters (v0.39) will be used...")
            options.adapters = data_files.data_list("available_Trimmomatic_adapters")
    
        # use default if not provided
        trim_params = {
            "ILLUMINACLIP": options.ILLUMINACLIP,
            "LEADING": str(options.LEADING),
            "TRAILING": str(options.TRAILING),
            "SLIDINGWINDOW": options.SLIDINGWINDOW,
            "MINLEN": str(options.MINLEN),
            "adapters": options.adapters
        }
        adapters_dict = {
            "adapters": options.adapters
        }
        
    ##########################
    elif options.software == "cutadapt":
        
        # Trimming adapters
    
        ## check adapters provided
            ## options.adapters_a
            ## options.adapters_A
            ## options.extra
            
        ## no adapters provided
        if (not options.adapters_a and not options.adapters_A and not options.extra):
            print (colored("** ERROR: No adapter trimming options provided...", 'red'))
            print ("Please provide any option")
            exit()
        
        ## create dictionary with 
        adapters_dict = {}
        if (options.adapters_a):
            adapters_dict['adapter_a'] = options.adapters_a
        
        if (options.adapters_A):
            adapters_dict['adapter_A'] = options.adapters_A
        else:
            options.adapters_A = ""
        
        ## set default
        if not options.min_len_read:
            options.min_len_read=15
    
    
        # use default if not provided
        trim_params = {
            "adapters_a": options.adapters_a,
            "adapters_A": options.adapters_A,
            "min_len_read": options.min_len_read
        }
    ##########################
    
    ################################################
    ## get files
    ################################################
    print ('+ Getting files from input folder... ')
    print ('+ Mode: fastq.\n+ Extension: ')
    print ("[ fastq, fq, fastq.gz, fq.gz ]\n")
    pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "fastq", ["fastq", "fq", "fastq.gz", "fq.gz"], options.debug)
    
    ## debug message
    if (Debug):
        print (colored("**DEBUG: pd_samples_retrieve **", 'yellow'))
        print (pd_samples_retrieved)

        print (colored("**DEBUG: adapters_dict **",'yellow'))
        print (adapters_dict)

    ## generate output folder, if necessary
    print ("\n+ Create output folder(s):")
    if not options.project:
        HCGB_files.create_folder(outdir)
    ## for samples
    outdir_dict = HCGB_files.outdir_project(outdir, options.project, pd_samples_retrieved, "trim", options.debug)
    
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
    print ("+ Trimming adapters for each sample retrieved...")    
    
    # Group dataframe by sample name
    sample_frame = pd_samples_retrieved.groupby(["new_name"])
    
    ## send for each sample
    ## use software trimmomatic or cutadapt
    
    ###############################
    if options.software == "trimmomatic":
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
            commandsSent = { executor.submit(trimmo_module, sorted(cluster["sample"].tolist()), 
                                             outdir_dict[name], name, threads_job, Debug, 
                                             trim_params, options.adapters): name for name, cluster in sample_frame }
    
            for cmd2 in concurrent.futures.as_completed(commandsSent):
                details = commandsSent[cmd2]
                try:
                    data = cmd2.result()
                except Exception as exc:
                    print ('***ERROR:')
                    print (cmd2)
                    print('%r generated an exception: %s' % (details, exc))

    ###############################
    elif options.software == "cutadapt":
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
            commandsSent = { executor.submit(cutadapt_module. sorted(cluster["sample"].tolist()), 
                                             outdir_dict[name], name, threads_job, 
                                             options.min_read_len, Debug, adapters_dict, options.extra): name for name, cluster in sample_frame }
    
            for cmd2 in concurrent.futures.as_completed(commandsSent):
                details = commandsSent[cmd2]
                try:
                    data = cmd2.result()
                except Exception as exc:
                    print ('***ERROR:')
                    print (cmd2)
                    print('%r generated an exception: %s' % (details, exc))

    ###############################
    print ("\n\n+ Trimming samples has finished...")
    ###############################
    
    ## functions.time_functions.timestamp
    start_time_partial = HCGB_time.timestamp(start_time_total)

    ## get files generated and generate symbolic link
    if not options.project:
        dir_symlinks = HCGB_files.create_subfolder('link_files', outdir)
        files2symbolic = []
        folders = os.listdir(outdir)

        ## debug message
        if (Debug):
            print (colored("**DEBUG: generate symbolic links for each file in " + dir_symlinks + "**", 'yellow'))
        
        for fold in folders:
            if fold.endswith(".log"):
                continue
            else:
                this_folder = outdir + '/' + fold
                subfiles = os.listdir(this_folder)
                for files in subfiles:
                    files_search = re.search(r".*trim_R\d{1}.*", files) ## only paired-end. Todo: single end
                    if files_search:
                        files2symbolic.append(this_folder + '/' + files)
    
        HCGB_files.get_symbolic_link(files2symbolic, dir_symlinks)

    ###############################
    ## Report & Summary
    ###############################
    if (options.skip_report):
        print ("+ No report generation...")
    else:
        print ("\n+ Generating a report using MultiQC module.")
        outdir_report = HCGB_files.create_subfolder("report", outdir)
    
        ## call multiQC report module
        givenList = [ v for v in outdir_dict.values() ]
        my_outdir_list = set(givenList)
        
        ## debug message
        if (Debug):
            print (colored("\n**DEBUG: my_outdir_list for multiqc report **", 'yellow'))
            print (my_outdir_list)
            print ("\n")

        trimm_report = HCGB_files.create_subfolder("trim", outdir_report)
        multiQC_report.multiQC_module_call(my_outdir_list, "Cutadapt", trimm_report,"")
        print ('\n+ A summary HTML report of each sample is generated in folder: %s' %trimm_report)
        
        ## QC analysis for trimmed reads
        if (Debug):
            print (colored("** Beginning FAStQC analysis **", 'red'))

        ## functions.time_functions.timestamp
        start_time_partial = HCGB_time.timestamp(start_time_partial)

    ## create FASTQC calling for trimmed reads
    pd_samples_retrieved_trimmed = sampleParser.files.get_files(options, input_dir, "trim", ['_trim'], options.debug)
    outdir_dict_fq = qc.fastqc(pd_samples_retrieved_trimmed, outdir, options, "trimmed", start_time_partial, Debug)
    qc.multiQC_rep(options, outdir, outdir_dict_fq, "trimmed")

    ################################################
    ## dump information and parameters
    ################################################
    ## samples information dictionary
    samples_info = {}
    samples_frame = pd_samples_retrieved_trimmed.groupby('new_name')
    for name, grouped in samples_frame:
        samples_info[name] = grouped['sample'].to_list()
    
    info_dir = HCGB_files.create_subfolder("info", outdir)
    print("+ Dumping information and parameters")
    runInfo = { "module":"trim", "time":time.time(),
                "RSP version":pipeline_version,
                'sample_info': samples_info,
                "outdir_dict": outdir_dict,
                 "outdir_dict_fq": outdir_dict_fq,
               'trim_params': trim_params }
    
    HCGB_info.dump_info_run(info_dir, 'trim', options, runInfo, options.debug)
    ################################################
    
    ################################################
    print ("\n*************** Finish *******************")
    start_time_partial = HCGB_time.timestamp(start_time_total)

    print ("\n+ Exiting trim module.")
    exit()

#############################################
def cutadapt_module(list_reads, sample_folder, name, threads, min_read_len, Debug, adapters, extra):
    """ Checks if the trimming process have been done previously. If not, it executes it
    calling cutadapt()    
    
    :param list_reads: name of the fastqc files of the sample to be trimmed
    :param sample_folder: path to the sample folder to store the results
    :param name: Name of the sample to be analyzed
    :param threads: number of CPUs to use.
    :param Debug: show additional message for debugging purposes.
    :param adapters: dictionary with the introduced adapters
    :param extra: provided extra options for cutadapt trimming process

    :type list_reads: string
    :type sample_folder: string
    :type name: string
    :type threads: string
    :type Debug: boolean
    :type adapters: dictionary
    :type extra: string
    

    :returns: None
    """
    
    ## check if previously trimmed and succeeded
    filename_stamp = sample_folder + '/.success_cut'
    if os.path.isfile(filename_stamp):
        stamp = HCGB_time.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'cutadapt'), 'yellow'))
    else:
        # Call cutadapt
        cutadapt_exe = set_config.get_exe('cutadapt')
        code_returned = cutadapt_caller.cutadapt(cutadapt_exe, list_reads, sample_folder, name, threads, min_read_len, Debug, adapters, extra)
        if code_returned:
            HCGB_time.print_time_stamp(filename_stamp)
        else:
            print ('** Sample %s failed...' %name)

    
#############################################
def trimmo_module(list_reads, sample_folder, name, threads, Debug, trimmomatic_params, adapters):
    
    """ This functions generates a trimmomatic call using java and trimmomatic 
    It checks if the trimming process have been done previously. If not, it executes it
    calling trimmomatic_caller.trimmo_call
    
    Checks if adapter file exists
    Returns code from trimmomatic_caller.trimmo_call: OK/FAIL
    
    :param list_reads: name of the fastqc files of the sample to be trimmed
    :param sample_folder: path to the sample folder to store the results
    :param name: Name of the sample to be analyzed
    :param threads: number of CPUs to use.
    :param Debug: show additional message for debugging purposes.
    :param trimmomatic_params: dictionary with the Trimmomatic parameters
    
    :type list_reads: string
    :type sample_folder: string
    :type name: string
    :type threads: string
    :type Debug: boolean
    :type trimmomatic_params: dictionary
    
    :returns: None
    """
    ## check if it exists
    if os.path.isfile(trimmomatic_params['adapters']):
        ## debug message
        if (Debug):
            print (colored("**DEBUG: trimmomatic_adapters file exists **", 'yellow'))
            print (trimmomatic_params['adapters'])
    else:
        ## rise error & exit
        print (colored("***ERROR: Trimmomatic adapters file does not exist: " + trimmomatic_params['adapters'],'red'))
        exit()
    
    ## check if previously trimmed and succeeded
    filename_stamp = sample_folder + '/.success_trimmo'
    if os.path.isfile(filename_stamp):
        stamp = HCGB_time.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s]" %(stamp, name), 'yellow'))
    
    else:

        ## get exe
        trimmomatic_jar = set_config.get_exe('trimmomatic')
        java_path = set_config.get_exe('java')
        if (Debug):
            print (colored("**DEBUG: trimmomatic executable **", 'yellow'))
            print (trimmomatic_jar)
            print(java_path)
    
        ## call: prints success if it works
        code_trim = trimmomatic_call.trimmo_call(java_path, sample_folder, name, list_reads, 
                           trimmomatic_jar, threads, trimmomatic_params, Debug)
        if code_trim:
            HCGB_time.print_time_stamp(filename_stamp)
        else:
            print ('** Sample %s failed...' %name)
