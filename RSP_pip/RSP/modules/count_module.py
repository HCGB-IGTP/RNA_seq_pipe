#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                      ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
"""
Counts sequence reads to reference genome
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
import pandas as pd
import csv

## import my modules
from RSP.scripts import multiQC_report, featurecounts, generate_matrix

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
def run_count(options):
	"""Main function of the module, organizes the counting process.
	
	There are several alternative softwares to count the reads.
	
	Some of them might have special configuration details. 
	
	:param options: input parameters introduced by the user. See RSP count -h.

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
	HCGB_aes.boxymcboxface("Counting reads for samples")
	print ("--------- Starting Process ---------")
	HCGB_time.print_time()
	
	## HCGB_time.timestamp
	start_time_partial = HCGB_time.timestamp(start_time_total)
	
	################################################
	## Check reference genome, annotations, CDS
	################################################
	print ('+ Getting reference genome information provided... ')

	# Check genome annotation file provided exists
	if HCGB_files.is_non_zero_file(options.ref_annot):
		options.ref_annot = os.path.abspath(options.ref_annot)
		print ('+ Reference genome annotation provided:... ')
		print('\t' + options.ref_annot)
	else:
		print("Reference genome annotation provided does not exist")
		exit()

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


	## Get bam files	
	print ('+ Mode: Sorted Bam files.\n+ Extension: ')
	print ("[ sortedByCoord.out.bam ]\n")
	
	## reorder sample bam files available
	bam_file_dict = {}
	files_dict = {
		"star": "Aligned.sortedByCoord.out.bam",
		"hisat2": "xx",
		"salmon": "x",
		"kallisto": "x"
	}


	## 
	sample_frame = pd_samples_retrieved.groupby(["new_name"])
	map_outdir_dict = HCGB_files.outdir_project(outdir, options.project, pd_samples_retrieved, "map", options.debug)


	## for each software create dictionary with files
	for soft_name2check in options.soft_name:
		bam_file_dict[soft_name2check] = {}
		for name, cluster in sample_frame:
			bam_file_dict[soft_name2check][name] = os.path.join(map_outdir_dict[name], soft_name2check, files_dict[soft_name2check])

	## debug message
	if (Debug):
		print (colored("**DEBUG: bam_file_dict **", 'yellow'))
		print (bam_file_dict)

	## generate output folder, if necessary
	print("\n+ Create output folder(s):")
	if not options.project:
		files_functions.create_folder(outdir)

	## for samples
	counts_outdir_dict = HCGB_files.outdir_project(outdir, options.project, pd_samples_retrieved, "counts", options.debug)
	
	## debug message
	if (Debug):
		print (colored("**DEBUG: counts_outdir_dict **", 'yellow'))
		print (counts_outdir_dict)

	# time stamp
	start_time_partial = HCGB_time.timestamp(start_time_total)

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
	## Let's start processing
	##########################
	print ("+ Counting reads for each sample and software specified...")    
	
	# Group dataframe by sample name
	#sample_frame = pd_samples_retrieved.groupby(["new_name"])
	
	for soft_name2check in options.soft_name:

		with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
				## sample_name, path_reference, reference_genome, index_name, reads_list, main_output, threads, parameters, software_list, Debug
				commandsSent = { executor.submit(gene_count_caller, 
													counts_outdir_dict[name], ## main_output,
													options.ref_annot, ## reference_genome annotation,
													bam_file_dict[soft_name2check][name], ## bam file to use,
													name,
													threads_job, ## threads 
													multimapping,
													options.stranded,
													soft_name2check, ## software list
													Debug): name for name, cluster in sample_frame }
				## (path, gtf_file, bam_file, name, threads, allow_multimap, stranded, , Debug)

				for cmd2 in concurrent.futures.as_completed(commandsSent):
					details = commandsSent[cmd2]
					try:
						data = cmd2.result()
					except Exception as exc:
						print ('***ERROR:')
						print (cmd2)
						print('%r generated an exception: %s' % (details, exc))

	## create report/summary for each software
	results_dict_soft = {}
	results_count={}
	outdir_report = HCGB_files.create_subfolder("report", outdir)
	module_outdir_report = HCGB_files.create_subfolder("counts", outdir_report)
	results_fold_dict_soft = {}

	for soft in options.soft_name:
		results_dict_soft[soft]={} 
		results_fold_dict_soft[soft]={}
		for name, cluster in sample_frame:
			results_dict_soft[soft][name] = os.path.join(counts_outdir_dict[name], soft, "featureCount.out")
			results_fold_dict_soft[soft][name] = os.path.join(counts_outdir_dict[name], soft)


		# Create count report    
		if (options.skip_report):
			print ("+ No report generation...")
		else:
			multiQC_report.create_module_report(main_outdir=outdir, soft_name=soft, 
				outdir_dict_given=results_fold_dict_soft[soft], module_given="counts", 
				options2multiqc="-dd 3")

		## for each software create count matrix
		all_counts_matrix_soft = generate_matrix.generate_matrix(results_dict_soft[soft], "Geneid")
		
		## dump data in folder provided
		csv_outfile = os.path.join(module_outdir_report, 'counts_RNAseq_' + soft + '.csv')
		all_counts_matrix_soft.to_csv(csv_outfile, quoting=csv.QUOTE_NONNUMERIC)
		
		## Dump count files into report folder
		print("Save counts in file: " + csv_outfile)
	

	## debug message
	if (Debug):
		print (colored("**DEBUG: results_dict_soft **", 'yellow'))
		print (results_dict_soft)
		
	
	## Finish count

	################################################
	## dump information and parameters
	################################################
	## samples information dictionary
	samples_info = {}
	samples_frame = pd_samples_retrieved.groupby('new_name')
	for name, grouped in samples_frame:
		samples_info[name] = grouped['sample'].to_list()
	
	info_dir = HCGB_files.create_subfolder("info", outdir)
	print("+ Dumping information and parameters")
	runInfo = { "module":"count", "time":time.time(),
				"RSP version":pipeline_version,
				'sample_info': samples_info,
				"outdir_dict": counts_outdir_dict,
				"results_count": results_count}
	
	HCGB_info.dump_info_run(info_dir, 'count', options, runInfo, options.debug)
	################################################
	
	################################################
	print ("\n*************** Finish *******************")
	start_time_partial = HCGB_time.timestamp(start_time_total)

	print ("\n+ Exiting count module.")
	
	return(runInfo)

################################################################################## 
def gene_count_caller(output_folder, gtf_file, bam_file, name2use, threads, multimapping, stranded, soft_name2use, Debug):

	## create subfolder for this mapping software
	output_folder = HCGB_files.create_subfolder(soft_name2use, output_folder)

	## check if previously counted and succeeded
	filename_stamp = output_folder + '/.success_featureCounts'
	if os.path.isfile(filename_stamp):
		stamp = HCGB_time.read_time_stamp(filename_stamp)
		print (colored("\tA previous command generated results on: %s [%s -- %s: %s]" %(stamp, name2use, 'count', 'featureCounts'), 'yellow'))
	else:
		## create call to mapping
		code_returned = featurecounts.get_counts_gene(output_folder, gtf_file, bam_file, name2use, str(threads), multimapping, str(stranded), Debug)

		if not os.path.isfile(code_returned):
			print ('** Sample %s failed...' %name2use)
