
import pandas as pd


############################################################
def generate_matrix(dict_files, index_name, Debug=False):
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
		#print ('+ File: ' + values)

		## check if file values exists
		data = pd.read_csv(values, sep='\t', skiprows=1)
		if Debug:
			print ("\n**DEBUG: values content file **")
			print (data)
		
		## skip if file is empty
		if data.empty:
			continue

		## get info, generate unique name and merge for samples
		data = data.set_index(index_name)
		data = data.iloc[:,-1:] ## subset
		data.columns = [*data.columns[:-1], key] ## rename
		all_data = pd.concat([all_data, data], axis=1, sort=True)

	return (all_data)
	##
