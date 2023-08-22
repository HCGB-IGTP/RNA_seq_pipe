#!/usr/bin/env python3
############################################################
## Author: Jose F. Sanchez                                ##
## Copyright (C) 2022                                     ##
## High Content Genomics and Bioinformatics IGPT Unit     ## 
## Lauro Sumoy Lab, IGTP, Spain                           ##
############################################################

## useful imports
import time
import io
import os
import re
import sys
from sys import argv
import subprocess
from termcolor import colored
import concurrent.futures

## import my modules
from RSP.config import set_config

## import HCGB
from HCGB.functions import system_call_functions, main_functions, time_functions
from HCGB.functions import files_functions, math_functions
from HCGB.functions.aesthetics_functions import debug_message

## plots
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from pandas.plotting import table

#####################
def help_info():
	'''Provide information on RNA biotype analysis'''
	print ("** TODO: Add description on RNAbiotype analysis and details")
	exit()

#######################################################################
def pie_plot_results(RNAbiotypes_stats_file, name, folder, Debug):
	
	##
	filename_stamp_plot = folder + '/.success_plot'
	if os.path.isfile(filename_stamp_plot):
		stamp = time_functions.read_time_stamp(filename_stamp_plot)
		print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'plot results'), 'yellow'))
	else:
		
		# PLOT and SHOW results
		RNAbiotypes_stats = main_functions.get_data(RNAbiotypes_stats_file, '\t', 'header=None')
	
		# create plot
		plt.figure(figsize=(16,8))
		df_genetype_2 = pd.DataFrame({'Type':RNAbiotypes_stats[0], 
									'Count':RNAbiotypes_stats[1]}).sort_values(by=['Count'])
	
		## get total count
		df_genetype_ReadCount_sum = df_genetype_2['Count'].sum()
	
		## filter 1% values
		minimun = df_genetype_ReadCount_sum * 0.01
		df_genetype_filter_greater = df_genetype_2[ df_genetype_2['Count'] >= minimun ]
		df_genetype_filter_smaller = df_genetype_2[ df_genetype_2['Count'] < minimun ]
	
		## create %values
		df_genetype_2['Percentage'] = (df_genetype_2['Count']/df_genetype_ReadCount_sum*100).round(3)
		
		## merge and generate Other class
		df_genetype_filter_smaller_sum = df_genetype_filter_smaller['Count'].sum() ## total filter smaller
		df_genetype_filter_greater2 = df_genetype_filter_greater.append({
			'Count':df_genetype_filter_smaller_sum, 
			'Type':'Other'}, ignore_index=True)
	
		## Create Pie Plot
		ax1 = plt.subplot(121, aspect='equal')
		df_genetype_filter_greater2.plot.pie(
			y = 'Count', 
			ax=ax1, 
			autopct='%1.2f%%', 
			shadow=False, 
			labels=df_genetype_filter_greater2['Type'], 
			legend = False)
	
		# plot table
		ax2 = plt.subplot(122)
		plt.axis('off')
		tbl = ax2.table(
			cellText=df_genetype_2.values, 
			colLabels=df_genetype_2.columns,
			loc='center', rowLoc='left', cellLoc='center', 
			)
		tbl.auto_set_font_size(True)
		#tbl.set_fontsize(12)
		tbl.scale(1.1,1.1)
	
		## set PDF name
		name_figure = os.path.join(folder, name + '_RNAbiotypes.pdf')
	
		## generate image
		plt.savefig(name_figure)		
		plt.close(name_figure)
		plt.close()
		
		## print time stamps
		time_functions.print_time_stamp(filename_stamp_plot)
		filename_stamp_all = folder + '/.success_all'
		time_functions.print_time_stamp(filename_stamp_all)
		
#######################################################################
def main():
	
	## ARGV
	if len (sys.argv) < 6:
		print ("\nUsage:")
		print ("python3 %s bam_file folder gtf_file threads name featureCount_bin multimapping[True/False]\n" %os.path.realpath(__file__))
		exit()
	
	bam_file = os.path.abspath(argv[1])
	folder = os.path.abspath(argv[2])
	gtf_file = os.path.abspath(argv[3])
	threads = argv[4]
	name = argv[5]
	featureCount_exe = argv[6]
	multimapping= argv[7]

	## Debug
	Debug=True
	
	## variables
	biotype_all(featureCount_exe, folder, gtf_file, bam_file, name, threads, Debug, multimapping)
	## plot results
	RNAbiotypes_stats_file = os.path.join(folder, name + '_RNAbiotype.tsv')
	if files_functions.is_non_zero_file(RNAbiotypes_stats_file):
		pie_plot_results(RNAbiotypes_stats_file, name, folder, Debug)
			
	
######
if __name__== "__main__":
	main()



#===============================================================================
# Example STAR Log Log.final.out
#===============================================================================
#
# 								 Started job on |	Oct 20 10:44:42
# 							 Started mapping on |	Oct 20 10:44:42
# 									Finished on |	Oct 20 11:13:25
# 	   Mapping speed, Million of reads per hour |	63.38
# 
# 						  Number of input reads |	30332431
# 					  Average input read length |	65
# 									UNIQUE READS:
# 				   Uniquely mapped reads number |	11259635
# 						Uniquely mapped reads % |	37.12%
# 						  Average mapped length |	47.60
# 					   Number of splices: Total |	0
# 			Number of splices: Annotated (sjdb) |	0
# 					   Number of splices: GT/AG |	0
# 					   Number of splices: GC/AG |	0
# 					   Number of splices: AT/AC |	0
# 			   Number of splices: Non-canonical |	0
# 					  Mismatch rate per base, % |	0.06%
# 						 Deletion rate per base |	0.00%
# 						Deletion average length |	1.00
# 						Insertion rate per base |	0.00%
# 					   Insertion average length |	1.36
# 							 MULTI-MAPPING READS:
# 		Number of reads mapped to multiple loci |	0
# 			 % of reads mapped to multiple loci |	0.00%
# 		Number of reads mapped to too many loci |	18685152
# 			 % of reads mapped to too many loci |	61.60%
# 								  UNMAPPED READS:
#   Number of reads unmapped: too many mismatches |	268064
# 	   % of reads unmapped: too many mismatches |	0.88%
# 			Number of reads unmapped: too short |	20122
# 				 % of reads unmapped: too short |	0.07%
# 				Number of reads unmapped: other |	99458
# 					 % of reads unmapped: other |	0.33%
# 								  CHIMERIC READS:
# 					   Number of chimeric reads |	0
# 							% of chimeric reads |	0.00%
#===============================================================================
