import sys
from sys import argv
import pandas as pd
import numpy as np
import re

### USAGE ###
# python organize_and_filter_genes.py output_folder file_location/filtered_insert_ratio_file.insert_ratios
# organizes inserts into which gene they are in, counts how many inserts per allele of each gene
# filters inserts according to number of inserts required in both alleles before we want to test that gene (change below)
# for mann whitney U test, want at least 5 inserts (I think!) in each allele
# prints out the number of inserts before and after filtering
# generates a file specific to this filtering number, added to the file name at the end: "_#"
# for all inserts within a given allele, calculate the mean of all the log2(39/28) for each insert, "mean39_28"

number_inserts_per_allele_needed_to_test = 5.0

def parse_file(filename, sep='\t'): 
	df = pd.read_csv(filename, sep=sep)
	file_string = str(filename)
	p = re.compile('[^/]+[^\.](?=\.)') # matches everything before the period and after the last slash, ie. the identifier of the file
	m = p.search(file_string) # finds regex in file name
	file_identifier = m.group() # prints match in string format
	return df, file_identifier

def calc_mean_insert_ratios(grp):
	grp['mean39_28'] = np.mean(grp['39_28_log2']) # calculate the geometric mean of each allele's ratios
	return grp

def check_for_singles(grp):
	number_alleles_found = len(grp['alleles'].unique())
	if number_alleles_found == 1:
		print grp
		print 'found one'
	return grp

class GeneObject(object):

	def __init__(self, gene_name):
		"""
		Each GeneObject is a gene (not specific to sp or sc). Keeps track of how many inserts per allele there are.
		"""
		self.gene = gene_name
		self.number_sp_inserts = 0.0
		self.number_sc_inserts = 0.0

if __name__ == '__main__':

	output_folder = argv[1].strip('/')
	files = argv[2:]

	file_id_dict = {}
	df_list = []

	for each_file in files:
		df, file_identifier = parse_file(each_file)
		df_list.append(df)
		file_id_dict[file_identifier] = df

	for file_identifier, each_df in file_id_dict.iteritems():
		each_df['gene'] = each_df.loc[:, 'annotation'].copy()
		each_df['gene'] = each_df['gene'].map(lambda x: x[2:])
		each_df['allele'] = each_df.loc[:, 'annotation'].copy()
		each_df['allele'] = each_df['allele'].map(lambda x: x[0:2])
		grouped = each_df.groupby('annotation', sort=False).apply(calc_mean_insert_ratios)
		grouped.set_index('gene', inplace=True, drop=False)
		unique_genes = list(grouped.gene.unique())
		gene_dict = {key: None for key in unique_genes}
		for each_gene in unique_genes:
			geneobj = GeneObject(each_gene)
			gene_dict[each_gene] = geneobj
		for index, row in grouped.iterrows():
			row_gene_name = row['gene']
			row_allele = row['allele']
			gene_object = gene_dict[row_gene_name]
			if row_allele == 'sc':
				gene_object.number_sc_inserts += 1.0
			if row_allele == 'sp':
				gene_object.number_sp_inserts += 1.0
			gene_dict[row_gene_name] = gene_object

		filtered_genes = []

		for each_gene, gene_info in gene_dict.iteritems():
			number_sc = float(gene_info.number_sc_inserts)
			number_sp = float(gene_info.number_sp_inserts)
			if number_sc >= number_inserts_per_allele_needed_to_test and number_sp >= number_inserts_per_allele_needed_to_test:
				filtered_genes.append(each_gene)

		filtered_gene_df = grouped[grouped['gene'].isin(filtered_genes)]

		print len(grouped)
		print len(filtered_gene_df)

		filtered_gene_df.to_csv(str(output_folder)+str('/')+str(file_identifier)+'_'+str(number_inserts_per_allele_needed_to_test)+'.filtered_gene_inserts', sep='\t', index=False)

		