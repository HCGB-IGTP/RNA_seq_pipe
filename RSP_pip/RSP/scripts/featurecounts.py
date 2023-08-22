#!/usr/bin/env python3
############################################################
## Author: Jose F. Sanchez & Mireia Marin                 ##
## Copyright (C) 2022                                     ##
## High Content Genomics and Bioinformatics IGPT Unit     ## 
## Lauro Sumoy Lab, IGTP, Spain                           ##
############################################################

#####################
def featurecounts_call(path, gtf_file, bam_file, name, threads, allow_multimap, stranded, option_featureCount, Debug):
	
	## option_featureCount: RNAbiotype, Gene count


	featureCount_exe = set_config.get_exe('featureCounts')

	## folder for results
	if not os.path.isdir(path):
		files_functions.create_folder(path)

	out_file = os.path.join(path, 'featureCount.out')
	logfile = os.path.join(path, name + '_RNAbiotype.log')

	filename_stamp_all = path + '/.success_all'
	if os.path.isfile(filename_stamp_all):
		stamp = time_functions.read_time_stamp(filename_stamp_all)
		print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, option_featureCount), 'yellow'))
		return()

	else:
		filename_stamp_featureCounts = path + '/.success_featureCounts'
		if os.path.isfile(filename_stamp_featureCounts):
			stamp = time_functions.read_time_stamp(filename_stamp_featureCounts)
			print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'featureCounts'), 'yellow'))
		else:
            
			## debugging messages
			if Debug:
				print ("** DEBUG:")
				print ("featureCounts system call for sample: " + name)
				print ("out_file: " + out_file)
				print ("logfile: " + logfile)
				print("option_featureCount: " + option_featureCount)
		
			## Mode
			if (option_featureCount=="RNAbiotype"):
			
				## Allow multimapping
				if allow_multimap:
					cmd_featureCount = ('%s -s %s -M -O -T %s -p -t exon -g transcript_biotype -a %s -o %s %s 2> %s' %(
						featureCount_exe, stranded, threads, gtf_file, out_file, bam_file, logfile)
					)
				else:
					cmd_featureCount = ('%s -s %s --largestOverlap -T %s -p -t exon -g transcript_biotype -a %s -o %s %s 2> %s' %(
						featureCount_exe, stranded, threads, gtf_file, out_file, bam_file, logfile)
					)
			
			elif (option_featureCount=="Gene Count"):
				## "-t exon":   Specify feature type in GTF annotation.
				##                              `exon' by default. Features used for read
				##                              counting will be extracted from annotation using the provided value.

				## "-g gene_name":      Specify attribute type in GTF annotation. `gene_id' by
				##                      default. Meta-features used for read counting will be
				##                      extracted from annotation using the provided value.

				## inicialmente se hizo con gene_name pero tras ver resultados entendemos que es mejor a nivel de gene_id ya que será IDs unicos
				
				## Allow multimapping
				if allow_multimap:
					cmd_featureCount = ('%s -p -t exon -g gene_id -s %s -M -O -T %s -p -a %s -o %s %s 2> %s' %(
						featureCount_exe, stranded, threads, gtf_file, out_file, bam_file, logfile)
						## -T threads
						## -s stranded
						## -M multimapping
					)
				else:
					cmd_featureCount = ('%s -p -t exon -g gene_id -s %s --largestOverlap -T %s -a %s -o %s %s 2> %s' %(
						featureCount_exe, stranded, threads, gtf_file, out_file, bam_file, logfile)
					)
		
			## send command for feature count
			## system call
			cmd_featureCount_code = system_call_functions.system_call(cmd_featureCount, False, True)
			if not cmd_featureCount_code:
				print("** ERROR: featureCount failed for sample " + name)
				exit()
				
			## print time stamp
			time_functions.print_time_stamp(filename_stamp_featureCounts)
		
	return (out_file)

#####################
def biotype_count(path, gtf_file, bam_file, name, threads, Debug, allow_multimap, stranded):
	
	out_file = featurecounts_call(path, gtf_file, bam_file, name, threads, allow_multimap, stranded, 'RNAbiotype', Debug)

	## parse results
	(extended_Stats_file, RNAbiotypes_stats_file) = parse_featureCount(out_file, path, name, bam_file, Debug)
		
	## debugging messages
	if Debug:
		print ("** DEBUG:")
		print ("extended_Stats: " + extended_Stats_file)
		print (main_functions.get_data(extended_Stats_file, '\t', 'header=None'))
		print ("RNAbiotypes_stats: " + RNAbiotypes_stats_file)
		print (main_functions.get_data(RNAbiotypes_stats_file, '\t', 'header=None'))

	return (out_file, extended_Stats_file, RNAbiotypes_stats_file)

#######################################################################
def parse_featureCount(out_file, path, name, bam_file, Debug):
	"""
	Parses featureCount results from STAR for RNAbiotype analysis.
	
	:param out_file: Name provided to featureCount for output results.
	:param path:
	:param name:
	
	"""

	## file names
	out_tsv_file_name = out_file + '.tsv'
	RNA_biotypes_file_name = os.path.join(path, name + '_RNAbiotype.tsv')

	##
	filename_stamp_parse = path + '/.success_parse'
	if os.path.isfile(filename_stamp_parse):
		stamp = time_functions.read_time_stamp(filename_stamp_parse)
		print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'parse results'), 'yellow'))
	else:
	
		## debugging messages
		if Debug:
			print ("** DEBUG:")
			print ("Parse results for sample: " + name)
			
		## parse results
		out_tsv_file = open(out_tsv_file_name, 'w')
		RNA_biotypes_file = open(RNA_biotypes_file_name, 'w')
		tRNA_count = 0
		
		##########################################
		### read count file
		##########################################
		count_file = open(out_file)
		count_file_text = count_file.read()
		count_file_lines = count_file_text.splitlines()	
	
		for line in count_file_lines:
			if line.startswith('#'):
				continue
			elif line.startswith('Geneid'):
				continue
			else:
				ID = line.split('\t')[0]
				count = int(line.split('\t')[-1])
				string2write_raw = "%s\t%s\n" %(ID, count)
				out_tsv_file.write(string2write_raw)
	
				tRNA_search = re.search(r".*tRNA", ID)
				if tRNA_search:
					tRNA_count = int(tRNA_count) + int(count)				
				elif (count > 0):
					RNA_biotypes_file.write(string2write_raw)
		
		## count and summary tRNA
		string2write = "tRNA\t%s\n" %tRNA_count
		RNA_biotypes_file.write(string2write)
		RNA_biotypes_file.close()
		
		##########################################
		### read summary count file
		##########################################
		summary_count_file = open(out_file + '.summary')
		summary_count_file_text = summary_count_file.read()
		summary_count_file_lines = summary_count_file_text.splitlines()	
	
		for line in summary_count_file_lines:
			if line.startswith('Status'):
				continue
			elif line.startswith('Assigned'):
				continue
			else:
				## adds Unassigned_Ambiguity
				## adds Unassigned_NoFeatures
				ID = line.split('\t')[0]
				count = int(line.split('\t')[-1])
	
				## skip empty entries
				if count == 0:
					continue
				string2write_raw = "%s\t%s\n" %(ID, count)
				out_tsv_file.write(string2write_raw)
	
		##########################################
		## get mapping statistics according to mapping software
		##########################################
		count_multi = 0
		count_unmap = 0
		mapping_folder = os.path.dirname(bam_file)
		mapping_stats = mapping_folder + '/Log.final.out'
		
		## ATTENTION: See example at the end of this file
		
		## -------------------------------- ##
		### STAR mapping		
		## -------------------------------- ##
		if files_functions.is_non_zero_file(mapping_stats):
			## debugging messages
			if Debug:
				print ("** DEBUG:")
				print ("STAR mapping available for sample: " + name)
				print ("mapping_folder: " + mapping_folder)
	
			mapping_stats_file = open(mapping_stats)
			mapping_stats_file_text = mapping_stats_file.read()
			mapping_stats_file_lines = mapping_stats_file_text.splitlines()	
	
			for line in mapping_stats_file_lines:
				multi_search = re.search(r".*Number of reads mapped to", line)
				unmap_search = re.search(r".*Number of reads unmapped", line)
				input_search = re.search(r".*input reads.*", line)
			
				if Debug:
					debug_message("mapping_stats_file Line")
					debug_message(line)
			
				if input_search:
					total_input_reads = int(line.split('\t')[-1])
					if Debug:
						print("total_input_reads")
						print(type(total_input_reads))
						print(total_input_reads)
							
				if multi_search:
					count_tmp = int(line.split('\t')[-1])
					count_multi = count_multi + count_tmp
					if Debug:
						print("count_tmp")
						print(type(count_tmp))
						print(count_tmp)
						print("count_multi")
						print(type(count_multi))
						print(count_multi)
	
				if unmap_search:
					count_reads_tmp = int(line.split('\t')[-1])
					#count_reads = math_functions.percentage(perc_tmp, total_input_reads)
					count_unmap = count_unmap + count_reads_tmp
					if Debug:
						print("perc_tmp")
						print(type(count_reads_tmp))
						print(count_reads_tmp)
						
						print("count_unmap")
						print(type(count_unmap))
						print(count_unmap)
		
		else:
	
			## -------------------------------- ##
			## tophat
			## -------------------------------- ##
	
			mapping_stats = mapping_folder + '/align_summary.txt' 
			count_map = 0
			total_input_reads = 0
			
			if files_functions.is_non_zero_file(mapping_stats):
				## debugging messages
				if Debug:
					print ("** DEBUG:")
					print ("tophat mapping available for sample: " + name)
					print ("mapping_folder: " + mapping_folder)
				
				mapping_stats_file = open(mapping_stats)
				mapping_stats_file_text = mapping_stats_file.read()
				mapping_stats_file_lines = mapping_stats_file_text.splitlines()	
	
				for line in mapping_stats_file_lines:
					map_search2 = re.search(r"Aligned.*\:\s+(\d+).*", line)
					input_search2 = re.search(r".*Input.*\:\s+(\d+).*", line)
					if input_search2:
						total_input_reads = input_search2.group(1)
					if map_search2:
						count_map = map_search2.group(1)
		
				####
				count_unmap = int(total_input_reads) - int(count_map)
	
			else:
				## other
				print ("Neither tophat or STAR..., no mapping statistics")
	
		### print mapping stats
		string2write_unmap = "unmapped\t%s\n" %count_unmap
		out_tsv_file.write(string2write_unmap)
		
		## close files
		out_tsv_file.close()
		summary_count_file.close()
		mapping_stats_file.close()
		count_file.close()
		## print timestamp
		time_functions.print_time_stamp(filename_stamp_parse)

	return(out_tsv_file_name, RNA_biotypes_file_name)

###########################
def get_counts_gene():

	out_file = featurecounts_call(path, gtf_file, bam_file, name, threads, allow_multimap, stranded, 'Gene Count', Debug)
	return(out_file)

	