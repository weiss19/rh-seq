import pandas as pd
from sys import argv
import sys

### USAGE ###
# python filter_inserts.py output_folder file_location/*.normalized_averaged_bioreps
# combines data from all bioreps into a single dataframe, filters inserts based on variables below
# can decide to filter based on coefficient of variation between bioreps, on the final average read count, or on both
# equal to cutoff value passes
# will generate a new file with filtering info in the title, to keep track of which file had which filtering scheme (coeffvar#_read#_...)
# keeps the insert if passes filtering metric in either 39C or 28C, doesn't need to pass in both
# when insert is found in one temperature but not the other, substitute a 1
# also saves a file before filtering with all the coding inserts "all_data.unfiltered_inserts_coding"
# if you want to try multiple filtering schemes, this is the point where you will have to run through the rest of the analysis separately for every separate filtered file

coeff_var_cutoff = 1.5 # cutoff for coefficient of variation
norm_reads_cutoff = 10.0 # cutoff for normalized read count
filter_strategy = 'both' # can be 'coeff' or 'reads' or 'both'


def parse_file(filename, sep='\t'): 
	df = pd.read_csv(filename, sep=sep)
	return df

global test
test = 0.0

def fillna_and_combine(grp):
	# keeps track of how many inserts have gone through, prints out so you can see how long the script is taking to run
	global test
	test += 1
	if test % 10000 == 0:
		print test

	grp['ID'] = grp.index
	filled = grp.fillna(method='ffill').fillna(method='bfill')
	return filled.iloc[0]

def make_and_populate_new_columns(grp):
	
	condition = grp.ix[0, 'condition']

	if condition != 'T0':
		column_name = condition+'_br_averaged_reads'
		grp[column_name] = grp['br_averaged_reads']
		column_name2 = condition+'_br_coeffvar'
		grp[column_name2] = grp['br_coeffvar']
	return grp

if __name__ == '__main__':

	output_folder = argv[1].strip('/')
	files = argv[2:]

	file_dict = {}
	df_list = []

	for each_file in files:
		df = parse_file(each_file)
		df.set_index('ID', inplace=True, drop=False)
		condition_ID = df['condition'][1]
		if condition_ID == 'T0':
			df.drop('bio_rep_ID', axis=1, inplace=True)
		file_dict[condition_ID] = df
		df_list.append(df)
	
	concatd_df = pd.concat(df_list, axis=0)
	all_columns = list(concatd_df.columns.values)
	drop_columns = [col for col in all_columns if '_n_' in col]
	drop_columns.extend([col for col in all_columns if '_tr_' in col])
	concatd_df.drop(drop_columns, axis=1, inplace=True)
	concatd_df.condition = concatd_df.condition.astype(str)
	concatd_df.set_index('condition', inplace=True, drop=False)
	grouped2 = concatd_df.groupby('condition', sort=False).apply(make_and_populate_new_columns)
	grouped2.set_index('ID', inplace=True, drop=False)
	grouped2.drop(['br_averaged_reads', 'br_coeffvar', 'condition'], axis=1, inplace=True)
	grouped2['NEWID'] = grouped2['ID']
	grouped_tmp = grouped2.groupby('NEWID', sort=False)
	group3 = grouped_tmp.apply(fillna_and_combine)
	columns = list(group3.columns.values)
	to_fill_with_one = [col for col in columns if 'reads' in col]
	group3[to_fill_with_one] = group3[to_fill_with_one].fillna(value=1.0) # leave empty "coeffvar" cells as NaN
	group3.to_csv(str(output_folder)+str('/')+'all_data.unfiltered_inserts_coding', sep='\t', index=False) 
	
	if filter_strategy == 'coeff':
		columns_to_filter_by = [col for col in columns if 'br_coeffvar' in col] # filter by either of the "conditions" coeffvar
		# if br coeffvar of 39 or 28 is lower than or equal to cutoff, keep
		filtered = group3[(group3[columns_to_filter_by[0]] <= coeff_var_cutoff) | (group3[columns_to_filter_by[1]] <= coeff_var_cutoff)]
		#will ignore NaN, it won't count either way and will never evaluate true
		print len(group3), 'before filtering'
		print len(filtered), 'after filtering'
		filtered.to_csv(str(output_folder)+str('/')+str(coeff_var_cutoff)+'_coeffvarcutoff.filtered_inserts', sep='\t', index=False)

	if filter_strategy == 'reads':
		columns_to_filter_by = [col for col in columns if 'br_averaged_reads' in col]
		# if norm averaged reads in 39 or 28 is higher than or equal to read cutoff, keep
		filtered = group3[(group3[columns_to_filter_by[0]] >= norm_reads_cutoff) | (group3[columns_to_filter_by[1]] >= norm_reads_cutoff)]
		print len(group3), 'before filtering'
		print len(filtered), 'after filtering'
		filtered.to_csv(str(output_folder)+str('/')+str(norm_reads_cutoff)+'_readscutoff.filtered_inserts', sep='\t', index=False)

	if filter_strategy == 'both':
		columns_to_filter_by = [col for col in columns if 'br_coeffvar' in col]
		columns_to_filter_by2 = [col for col in columns if 'br_averaged_reads' in col]
		filtered = group3[(group3[columns_to_filter_by[0]] <= coeff_var_cutoff) | (group3[columns_to_filter_by[1]] <= coeff_var_cutoff)]
		filtered2 = filtered[(filtered[columns_to_filter_by2[0]] >= norm_reads_cutoff) | (filtered[columns_to_filter_by2[1]] >= norm_reads_cutoff)]
		print len(group3), 'before filtering'
		print len(filtered), 'after filtering for coeffvar'
		print len(filtered2), 'after filtering for coeffvar and read count'
		filtered2.to_csv(str(output_folder)+str('/')+str(coeff_var_cutoff)+'_'+str(norm_reads_cutoff)+'_coeffvarANDreadcutoff.filtered_inserts', sep='\t', index=False)
