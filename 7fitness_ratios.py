from sys import argv
import sys
import pandas as pd
import numpy as np
import math
import re

### USAGE ###
# python fitness_ratios.py output_folder file_location/filtered_file.filtered_inserts
# for each insert, calculates log2(39/28) using final normalized read counts


def parse_file(filename, sep='\t'): 
	df = pd.read_csv(filename, sep=sep)
	file_string = str(filename)
	p = re.compile('[^/]+[^\.](?=\.)') # matches everything before the period and after the last slash, ie. the identifier of the file
	m = p.search(file_string) # finds regex in file name
	file_identifier = m.group() # prints match in string format
	return df, file_identifier


if __name__ == '__main__':

	output_folder = argv[1].strip('/')
	files = argv[2:]

	file_id_dict = {}
	df_list = []

	for each_file in files:
		df, file_identifier = parse_file(each_file)
		df.set_index('ID', inplace=True, drop=False)
		df_list.append(df)
		file_id_dict[file_identifier] = df

	for file_identifier, each_df in file_id_dict.iteritems():
		columns = each_df.columns.values
		read_columns = [col for col in columns if '_averaged_reads' in col]
		for each_col in read_columns:
			if each_col.startswith('28'):
				reads_28 = each_col
			if each_col.startswith('39'):
				reads_39 = each_col
			if each_col.startswith('T0'):
				reads_T0 = each_col
		div39_28 = each_df[reads_39] / each_df[reads_28]
		each_df['39_28_log2'] = np.log2(div39_28)
		each_df.to_csv(str(output_folder)+str('/')+str(file_identifier)+'.insert_ratios', sep='\t', index=False)



		


