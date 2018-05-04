import pandas as pd
import numpy as np
import sys
from sys import argv
import re
from scipy import stats
import statsmodels.stats.multitest as smm
import math
import matplotlib.pyplot as plt


### USAGE ###
# finally, computes the Mann-Whitney U test, compares the log2(39/28) values of inserts in the Scer allele vs the Spar allele 
# test is two-tailed (want to see if scer > spar or spar > scer), so multiply pvals by 2
# applies multiple testing correction via Benjamini Hochberg
# saves a file with each gene that was tested, pval and adjusted pval
# prints out number of genes tested
# makes old bar plots of any genes tested below a certain adjusted pval cutoff, if you want (default is turned off)


def parse_file(filename, sep='\t'): 
	df = pd.read_csv(filename, sep=sep)
	file_string = str(filename)
	p = re.compile('[^/]+[^\.](?=\.)') # matches everything before the period and after the last slash, ie. the identifier of the file
	m = p.search(file_string) # finds regex in file name
	file_identifier = m.group() # prints match in string format
	return df, file_identifier

class GeneObject(object):

	def __init__(self, df):
		self.df = df
		self.mwu_test = None
		self.two_tailed_pval = None

def make_plots(gene_df):
	
	gene_name = gene_df['gene'].iloc[0]
	sc_df = gene_df.loc[gene_df['allele'] == 'sc']
	sp_df = gene_df.loc[gene_df['allele'] == 'sp']
	sc_ratios = sc_df['39_28_log2'].tolist()
	sp_ratios = sp_df['39_28_log2'].tolist()
	geo_mean_sc_3928 = np.mean(sc_ratios)
	geo_mean_sp_3928 = np.mean(sp_ratios)
	number_sc_inserts = len(sc_df)
	number_sp_inserts = len(sp_df)
	total_inserts = number_sc_inserts + number_sp_inserts
	ticklocats = np.arange(total_inserts)
	fitness_data = sc_ratios + sp_ratios

	abs_fitness = []
	for each_fitness in fitness_data:
		abs_fitness.append(math.fabs(each_fitness))
	highest_fitness = max(abs_fitness)
	axis_limit = highest_fitness + 1
	negative_limit = axis_limit * (-1)
	colours = []
	for x in range(number_sc_inserts):
		colours.append('#00A5CD')  # blue
	for x in range(number_sp_inserts):
		colours.append('#FAB900')  # yellow

	
	# plotting magic

	fig, ax = plt.subplots(figsize = (total_inserts + 8, 10))
	ax.bar(ticklocats, fitness_data, color = colours, width = 0.75)
	fig.suptitle(str(gene_name), size = 40)
	ax.set_xticks([])
	ax.yaxis.set_ticks_position('left')
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.set_ylabel('Log2 ratio 39/28', size = 40)
	ax.set_xlabel('Inserts: Sc = yellow, Sp = blue', size = 40)
	ax.set_ylim([negative_limit, axis_limit])
	ax.tick_params(labelsize=25)
	fig.savefig(str(output_folder)+str('/')+str(file_identifier)+str(gene_name)+'_insert_ratios_plot.eps', format='eps', dpi=1000)
	plt.close()

	mean_dict = {
	'Sc_3928': geo_mean_sc_3928, 
	'Sp_3928': geo_mean_sp_3928
	}

	plt.bar(range(len(mean_dict)), mean_dict.values(), align='center')
	plt.xticks(range(len(mean_dict)), mean_dict.keys())
	plt.savefig(str(output_folder)+str('/')+str(file_identifier)+str(gene_name)+'_geometric_means_plot_6.eps', format='eps', dpi=1000)
	plt.close()

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
		unique_genes = each_df.gene.unique()
		unique_genes_dict = {key: None for key in unique_genes}
		another_gene_dict = {key: None for key in unique_genes}
		pval_dict = {key: None for key in unique_genes}
	
		for each_gene in unique_genes:
			gene_df = each_df.loc[each_df['gene'] == each_gene]  # generate new df for each gene that is just its inserts from both alleles
			sc_df = gene_df.loc[gene_df['allele'] == 'sc']
			sp_df = gene_df.loc[gene_df['allele'] == 'sp']
			appended_df = pd.concat([sc_df, sp_df])
			another_gene_dict[each_gene] = appended_df
			gene_object = GeneObject(gene_df)
			sc_inserts_ratios = [float(i) for i in sc_df['39_28_log2'].tolist()]
			sp_inserts_ratios = [float(i) for i in sp_df['39_28_log2'].tolist()]
			mwu_test = stats.mannwhitneyu(sc_inserts_ratios, sp_inserts_ratios)   # mwu test
			two_tailed_pval = mwu_test[1] * 2.0
			gene_object.mwu_test = mwu_test
			gene_object.two_tailed_pval = two_tailed_pval
			unique_genes_dict[each_gene] = gene_object
			pval_dict[each_gene] = two_tailed_pval

		gene_list = []
		pvalue_list = []
		for each_gene, each_pval in pval_dict.iteritems():
			pvalue_list.append(each_pval)
			gene_list.append(each_gene)

		#multiple testing correction
		fdrbh_output = smm.multipletests(pvalue_list, method='fdr_bh') # benjamini hochberg method
		adjusted_pvalues = fdrbh_output[1].tolist()
		zipped_pvals = zip(pvalue_list, adjusted_pvalues)
		zipped_gene_info = zip(gene_list, zipped_pvals)
		unpacked_gene_info = [(gene, pval, adjusted_pval) for (gene, (pval, adjusted_pval)) in zipped_gene_info]
		pval_df = pd.DataFrame(unpacked_gene_info, columns=['gene', 'pval', 'adjusted_pval'])
		print len(pval_df), 'number of genes tested'

		pval_df.sort_values(by='pval', inplace=True)
		pval_df.to_csv(str(output_folder)+str('/')+str(file_identifier)+'.mwu_test_results', sep='\t', index=False)
		genes_to_graph_df = pval_df.loc[pval_df['adjusted_pval'] <= 0.05]  # set the pval cutoff for genes you want to make plots of
		genes_to_graph_list = genes_to_graph_df['gene'].tolist()

		# you can uncomment below if you want to make bargraphs of gene data where adjusted pval is beneath a certain cutoff
		# we swtiched to other ways of visualizing the data so I don't really use this function anymore

		# for each_gene in genes_to_graph_list:
		# 	gene_df = each_df.loc[each_df['gene'] == each_gene] # get all the inserts in that gene
		# 	make_plots(gene_df)


