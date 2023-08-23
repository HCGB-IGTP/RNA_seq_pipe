#!/usr/bin/env python3
############################################################
## Author: Jose F. Sanchez                                ##
## Copyright (C) 2022                                     ##
## High Content Genomics and Bioinformatics IGPT Unit     ## 
## Lauro Sumoy Lab, IGTP, Spain                           ##
############################################################

"""
Calls multiQC to generate HTML statistics reports.
"""
## useful imports
import os
import io
import sys
from io import open
from sys import argv
from termcolor import colored

## import my modules
from HCGB import functions
from RSP.config import set_config

############
def multiQC_module_call(givenList, name, path, option):
    """
    Prepares files for multiQC report generation.
    
    :param givenList: List of folder to search for multiQC report.
    :param name: Name to include in the html report.
    :param path: Absolute path for the output folder.
    :param option: Some options to provide to multiQC_call.
    
    :type givenList: list
    :type name: string
    :type path: string
    :type option: string
    
    .. seealso:: This function depends on other RSP functions called:
    
        - :func:`RSP.scripts.functions.printList2file`
        
        - :func:`RSP.scripts.multiQC_report.multiQC_call`
    
    """
    pathFile = path + '/' + 'samples.txt'
    functions.main_functions.printList2file(pathFile, givenList)
    multiQC_call(pathFile, name, path, option)    
    
############
def multiQC_call(pathFile, name, folder, option):
    """
    multiQC_ report generation call.
    
    :param pathFile: File containing list of files to include in report.
    :param name: Name to include in the html report.
    :param folder: Absolute path for the output folder.
    :param option: Options to provide to multiQC call.
    
    :type pathFile: string
    :type name: string 
    :type folder: string 
    :type option: string
    
    :returns: :func:`RSP.scripts.functions.system_call_functions.system_call` output (OK/FALSE)
        
    .. seealso:: This function depends on other RSP functions called:
    
        - :func:`RSP.scripts.functions.system_call_functions.system_call`
    
    """
    multiqc_bin = set_config.get_exe("multiqc")
    ## set options for call
    cmd = "%s --force -o %s -n %s -l %s -p -i 'MultiQC report' -b 'HTML report generated for multiple samples and steps' %s" %(multiqc_bin, folder, name, pathFile, option)
    
    ## if a report was previously generated in the folder 
    ## force to delete and generate a new one
    return(functions.system_call_functions.system_call(cmd))

#################################
def create_module_report(main_outdir, soft_name, outdir_dict_given, module_given, options2multiqc):
	
	print ("\n+ Generating a report using MultiQC module for results obtained with software: " + soft_name)
	outdir_report = HCGB_files.create_subfolder("report", main_outdir)
	module_outdir_report = HCGB_files.create_subfolder(module_given, outdir_report)

	## get subdirs generated and call multiQC report module
	givenList = []
	print ("+ Detailed information for each sample could be identified in separate folders.")
		
	## call multiQC report module
	givenList = [ v for v in outdir_dict_given.values() ]
	my_outdir_list = set(givenList)

	## debug message
	if (Debug):
		print (colored("\n**DEBUG: my_outdir_list for multiqc report **", 'yellow'))
		print (my_outdir_list)
		print ("\n")
		
	module_given_report = HCGB_files.create_subfolder(soft_name, module_outdir_report)
	multiQC_module_call(my_outdir_list, soft_name, module_given_report, options2multiqc)
	print ('\n+ A summary HTML report of each sample for software %s is generated in folder: %s' %(soft_name, module_given_report))

