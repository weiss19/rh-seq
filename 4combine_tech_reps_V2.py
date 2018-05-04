import sys
from sys import argv
import pandas as pd

### USAGE ###
# python combine_tech_reps_V2.py output_folder file_location/*.normalized_pooled_reads_coding
# combines tech reps into a single file, one for each biorep
# e.g. 28A1, 28A2, 28A3 become 28A
# finds all inserts from across all files, takes average read count from total number of tech reps (replaces 0 with 1)
# calculates coefficient of variation across each insert for the average it calculated. takes average of all of those and prints out for each bio rep


def parse_file(filename, sep='\t'): 
	df = pd.read_csv(filename, sep=sep)
	file_identifier = df['tech_rep_name'][1]
	if file_identifier.startswith('T0'):
		bio_rep_id = 'T0'
	else:
		bio_rep_id = file_identifier[:-1].strip('_')  # chops off the last character, removing info for 
	return df, bio_rep_id

def make_and_populate_new_columns(grp):
	column2_name = grp.ix[0, 'tech_rep_name']+'_norm_n'
	grp[column2_name] = grp['normalized_read']
	return grp

def combine_rows_by_summing(grp):
	# will spit out into the command line every 10,000 inserts, just to keep track of how much longer it needs to finish running
	global test
	test += 1
	if test % 10000 == 0:
		print test

	column_names = list(grp.columns.values)
	column_names_to_sum = [name for name in column_names if '_norm_n' in name]
	if len(grp) > 1: # when more than one row is found, ie. when the insert is found in more than one techrep
		for each_column in column_names_to_sum:
			grp[each_column+'_summed'] = grp[each_column].sum()
	
	if len(grp) == 1: # when the insert was found only in one techrep
		for each_column in column_names_to_sum:
			grp[each_column+'_summed'] = grp[each_column]
	grp_1 = grp.iloc[0]	
	return grp_1

if __name__ == '__main__':

	output_folder = argv[1].strip('/')
	files = argv[2:]

	bio_rep_df_dict = {}
	for each_file in files:
		df, bio_rep_id = parse_file(each_file)

		if bio_rep_id in bio_rep_df_dict: # sort files into a dict based on which condition/biorep they are a part of
			bio_rep_df_dict[bio_rep_id].append(df)
		if bio_rep_id not in bio_rep_df_dict: 
			bio_rep_df_dict[bio_rep_id] = [df]

	for each_biorep, list_of_dfs in bio_rep_df_dict.iteritems():
		global test
		test = 0.0

		# pandas magic to combine the files and make it look pretty
		concatd_df = pd.concat(list_of_dfs) # merges dataframes from same biorep
		idx_concat_df = concatd_df.set_index('tech_rep_name', drop=False)
		grouped = idx_concat_df.groupby('tech_rep_name', sort=False).apply(make_and_populate_new_columns).fillna(value=0.0)
		grouped.drop(['n', 'normalized_read', 'tech_rep_name', 'total_mapped_reads'], axis=1, inplace=True)
		id_idx_grouped = grouped.set_index('ID', drop=False)
		column_names = list(id_idx_grouped.columns.values)
		column_names_to_sum = [name for name in column_names if '_n' in name]
		id_grouped = id_idx_grouped.groupby('ID', sort=False).apply(combine_rows_by_summing)
		id_grouped.replace(to_replace=0.0, value=1.0, inplace=True)
		for each_column in column_names_to_sum:
			id_grouped.drop(each_column, 1, inplace=True)
		columns_to_average = [column for column in list(id_grouped.columns.values) if 'norm_n' in column]
		id_grouped['averaged_tech_rep_reads'] = id_grouped[columns_to_average].mean(axis=1)
		id_grouped['tr_coeffvar'] = (id_grouped[columns_to_average].std(axis=1)) / (id_grouped[columns_to_average].mean(axis=1))
		mean_coeffvar = id_grouped['tr_coeffvar'].mean()

		print 'The mean coeffvar of ', each_biorep, ' is ', mean_coeffvar

		# save the combined data to a new file
		id_grouped.to_csv(str(output_folder)+str('/')+str(each_biorep)+'.normalized_averaged_techreps', sep='\t', index=False)




