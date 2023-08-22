


############################################################
def generate_matrix(dict_files, index_name):
	######################################################
	### For a dictionary containing names as keys and files as values, 
	###	generates a count matrix with the classification from featurecounts
	###		
	### It returns a dataframe containing for each index generated
	### count values for each sample (sample_name) in columns
	###
	#########################################################
	all_data = pd.DataFrame()
	for key,values in dict_files.items():
		print ('+ Reading information from sample: ', key)	
		data = pd.read_csv(values, sep='\t', header=None, names=[index_name, key])
		
		## skip if file is empty
		if data.empty:
			continue

		## get info, generate unique name and merge for samples
		data = data.set_index(index_name)
		all_data = pd.concat([all_data, data], axis=1, sort=True)
	
	return (all_data)
	##
