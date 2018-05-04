from sys import argv
import pandas as pd
import sys
import re

### USAGE ###
# python combine_bio_reps.py output_folder file_location/*.normalized_average_techreps
# combines biological replicates into one file, e.g. 28A, 28B, 28C into 28
# takes average of read counts across reps, if insert isn't found in all reps, sub in value of 1. average is across all reps.
# calculates coefficient of variation across the bio rep values. this value is what is used to filter inserts later.


def parse_file(filename, sep='\t'): 
	df = pd.read_csv(filename, sep=sep)
	
	file_string = str(filename)
	p = re.compile('[^/]+[^\.](?=\.)') # matches everything before the period and after the last slash, ie. the identifier of the file
	m = p.search(file_string) # finds regex in file name
	file_identifier = m.group() # prints match in string format
	if file_identifier == 'T0':
		condition_ID = 'T0'
	else:
		condition_ID = file_identifier[:-1]  # chops off the last character, e.g. 39C becomes 39
	df['bio_rep_ID'] = file_identifier
	df['condition'] = condition_ID
	return df, condition_ID

def make_and_populate_new_columns(grp):
	column_name = grp.ix[0, 'bio_rep_ID']+'_n_av'
	grp[column_name] = grp['averaged_tech_rep_reads']
	column_name2 = grp.ix[0, 'bio_rep_ID']+'_tr_coeff'
	grp[column_name2] = grp['tr_coeffvar']
	return grp

def combine_rows_by_summing(grp):
	global test
	test += 1
	if test % 100 == 0:
		print test
	column_names = list(grp.columns.values)
	column_names_to_sum = column_names_to_sum = [name for name in column_names if '_n_av' in name or 'n_summed' in name or 'tr_coeff' in name]
	for each_column in column_names_to_sum:
		grp[each_column+'_summed'] = grp[each_column].sum()  # a sum of number plus 0
	grp_1 = grp.iloc[0] # all rows should be the same now, just return the first one
	return grp_1

def return_first(grp):
	grp_1 = grp.iloc[0]
	return grp_1	

if __name__ == '__main__':

	output_folder = argv[1].strip('/')
	files = argv[2:]

	condition_df_dict = {}

	for each_file in files:
		df, condition_id = parse_file(each_file)

		if condition_id in condition_df_dict: # sort files into a dict based on which condition/biorep they are a part of
			condition_df_dict[condition_id].append(df)
		if condition_id not in condition_df_dict: 
			condition_df_dict[condition_id] = [df]

	for each_condition, list_of_df in condition_df_dict.iteritems():
		global test
		test = 0.0

		# our first set of experiments didn't have "bio reps" of T0, only tech reps, so it is handled slightly differently here.
		if each_condition == 'T0':
			for each_df in list_of_df:
				each_df['condition'] = 'T0'
				each_df.rename(columns={'tr_coeffvar': 'T0_tr_coeffvar', 'averaged_tech_rep_reads': 'T0_averaged_reads'}, inplace=True)
				each_df.to_csv(str(output_folder)+str('/')+str(each_condition)+'.normalized_averaged_bioreps', sep='\t', index=False)
				print 'done T0'
		else:
			# pandas magic to get data into correct format
			concatd_df = pd.concat(list_of_df, axis=0).fillna(value=0.0) # merges dataframes from same biorep
			idx_df = concatd_df.set_index('bio_rep_ID', drop=False, verify_integrity=False)
			#group by bio rep and then move all biorep specific row info to new row
			grouped = idx_df.groupby('bio_rep_ID', sort=False).apply(make_and_populate_new_columns).fillna(value=0.0)
			columns_to_drop = []
			columns_to_drop.extend(['averaged_tech_rep_reads', 'bio_rep_ID', 'tr_coeffvar'])
			grouped.drop(columns_to_drop, axis=1, inplace=True)
			idx_grouped = grouped.set_index('ID', drop=False)
			column_names = list(idx_grouped.columns.values)
			column_names_to_sum = [name for name in column_names if '_n_av' in name or 'n_summed' in name or 'tr_coeff' in name]
			id_grouped = idx_grouped.groupby('ID', sort=False)[column_names_to_sum].sum()
			columns_to_add = list((set(column_names) - set(column_names_to_sum)))
			to_add = idx_grouped.groupby('ID', sort=False)[columns_to_add].first()	
			merged = pd.concat([to_add, id_grouped], axis=1)
			merged.replace(to_replace=0.0, value=1.0, inplace=True)
			columns_to_average = [column for column in list(merged.columns.values) if 'n_av' in column]
			merged['br_averaged_reads'] = merged[columns_to_average].mean(axis=1)
			merged['br_coeffvar'] = (merged[columns_to_average].std(axis=1)) / (merged[columns_to_average].mean(axis=1))
			mean_coeffvar = merged['br_coeffvar'].mean()

			print 'The mean coeffvar of ', each_condition, ' is ', mean_coeffvar

			merged.to_csv(str(output_folder)+str('/')+str(each_condition)+'.normalized_averaged_bioreps', sep='\t', index=False)




