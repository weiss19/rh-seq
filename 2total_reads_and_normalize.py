from sys import argv
import pandas as pd
import numpy as np
import re

### Usage ###
# python total_reads_and_normalize.py output_folder file_location/*.fastq_pooled_reads_clean 
# calculates total number of mapped reads for each file. calculates average total mapped reads across all files. 
# normalizes read counts "n" to average total mapped reads / file-specific total mapped reads


def parse_file(filename, sep='\t'): # reads in file, returning a pandas dataframe
	
	#finding the filename identifier from the filename
	file_string = str(filename)
	p = re.compile('[^/]+[^\.](?=\.)') # matches everything before the period and after the last slash, ie. the identifier of the file
	m = p.search(file_string) # finds regex in file name
	file_identifier = m.group() # prints match in string format
	
	#reading in file as a pandas dataframe and adding file identifier to new column
	df = pd.read_csv(filename, sep=sep)
	df['tech_rep_name'] = file_identifier
	return df
	

def total_mapped_reads(df):
	#sums up n for a given dataframe and adds a new column where total_n is in every row
	total_mapped_reads = float(df['n'].sum())
	df['total_mapped_reads'] = total_mapped_reads
	return df, total_mapped_reads

def average_mapped_reads(list):
	#averages all the total mapped reads calculated across all the files 
	average = np.mean(list)
	return float(average)

def normalize_reads(df, average):
	#divides each 'n' by file-specific total mapped reads, multiplies by average-total-mapped-reads
	df['all_files_average'] = average
	df['normalized_read'] = ( df['n'] / df['total_mapped_reads'] ) * df['all_files_average']
	return df

if __name__ == '__main__':

	output_folder = argv[1].strip('/')
	files = argv[2:]

	list_of_dfs = []
	list_of_total_mapped_reads = []

	for each_file in files:
		df = parse_file(each_file)
		df2, file_average_mapped_reads = total_mapped_reads(df)
		list_of_total_mapped_reads.append(file_average_mapped_reads)
		list_of_dfs.append(df2)

	all_files_average = average_mapped_reads(list_of_total_mapped_reads)

	for each_df in list_of_dfs:
		df3 = normalize_reads(each_df, all_files_average)
		file_id = df3['tech_rep_name'][1]
		df3.to_csv(str(output_folder)+str('/')+str(file_id)+'.normalized_pooled_reads', sep='\t', index=False)
		#saves output df to "file_identifier.normalized_pooled_reads" in output directory specified






