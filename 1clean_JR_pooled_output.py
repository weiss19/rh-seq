
from sys import argv

### Usage ###
# python clean_JR_pooled_output.py file_location/annotation_file file_location/pooled_read_file_to_clean 
# saves new clean file to same place as the un-clean file, with _clean at the end 
# makes dict of annotation file gene : chromosome mapping. chucks any inserts that contradict this (insert maps to part a chromosome, but is annotated with a gene that isn't there)


def parse_files(filename, remove_header, sep='\t'):
	"""
	Read in each file and return a list of rows. JR pooled reads input or the clean annotation input.
	"""	
	rows = []
	with open(filename) as f:
		if remove_header == True:
			f.readline() # skip header line if necessary
		for each_row in f:
			row_data = each_row.rstrip().split(sep) # remove trailing whitespace and split by tab delimitor \t
			rows.append(row_data) 
	return rows

def map_genes_to_sc_chr(annotation_data):
	"""
	For every gene, find the chromosome it maps to on sc.
	"""

	roman_numerals_as_ints = {
	'I': '1', 
	'II': '2', 
	'III': '3', 
	'IV': '4',
	'V': '5',
	'VI': '6',
	'VII': '7',
	'VIII': '8',
	'IX': '9',
	'X': '10',
	'XI': '11',
	'XII': '12',
	'XIII': '13',
	'XIV': '14',
	'XV': '15',
	'XVI': '16'
	}

	gene_chr_dict = {}

	for each_gene in annotation_data:
		species = each_gene[0][0:2]
		if species == 'sp':  
			if each_gene[0][2:6] == 'Scer':  #account for weird sp genes that are labeled this way
				gene = each_gene[0][2:]
				if each_gene[0][7:9].isdigit() == True:
					gene_chr_dict[gene] = each_gene[0][7:9]
				else:
					gene_chr_dict[gene] = each_gene[0][7:8]
			else:
				continue # don't process paradoxus genes
		if species == 'sc':
			if each_gene[1] == 'sc2-micron' or each_gene[1] == 'scchrMito':  # don't process genes in mito DNA or 2-micro plasmid
				continue
			if each_gene[0] != 'sc':	# gets rid of "empty" genes that just say sc ...
				gene = each_gene[0][2:]
				scaffold_roman = each_gene[1][5:]
				scaffold_int = roman_numerals_as_ints[scaffold_roman]  # convert roman to int
				gene_chr_dict[gene] = scaffold_int

	return gene_chr_dict

def clean_file(pooled_reads_file_data, genes_to_chromosomes):
	"""
	Go through all inserts in pooled_reads file, populate dict of genes and corresponding chr mapped to.
	Check to make sure it matches chr from genes_to_chromosomes. If not, delete insert (don't add it to new clean file)
	"""
	cleaned_filename = file_to_clean+'_clean'
	wf = open(cleaned_filename, 'w')
	wf.writelines("ID\tscaffold\tstrand\tlocation\tannotation\tn\trel_loc\trel_prop\tgene_length\n")

	number_chucked = 0.0
	bad_genes_dict = {}
	bad_inserts_list = {}
	keep_inserts_list = {}

	for each_insert in pooled_reads_file_data:
		ID = each_insert[0]
		if ID[0:2] != 'sp' and ID[0:2] != 'sc':  # inserts in plasmid
			keep_inserts_list[ID] = []
		if ID[0:2] == 'sc':  # keep all sc inserts
			keep_inserts_list[ID] = []
		if ID[0:2] == 'sp': # all sp inserts
			if ID[0:2] != 'sp':
				print ID
				raise Exception('this isnt sp!')
			if each_insert[4] == 'NC':
				keep_inserts_list[ID] = []  # keep NC 
			else:
				if each_insert[4] != 'sp':  # not odd case where gene is just "sp", keep down below though
					gene = each_insert[4][2:]
					scaffold_num = each_insert[1][2:]
					if gene in genes_to_chromosomes:
						correct_scaffold = genes_to_chromosomes[gene]
						#print scaffold_num
						#print correct_scaffold
						if scaffold_num != correct_scaffold:
							bad_inserts_list[ID] = []
						if scaffold_num == correct_scaffold:
							keep_inserts_list[ID] = []
				else:
					keep_inserts_list[ID] = []
		if ID not in keep_inserts_list and ID not in bad_inserts_list:
			print ID

	print len(pooled_reads_file_data), 'total'
	print len(bad_inserts_list), 'bad'
	print len(keep_inserts_list), 'good'

	for each_insert in pooled_reads_file_data:
		ID = each_insert[0]
		if ID in keep_inserts_list:
			wf.writelines(each_insert[0]+'\t'+each_insert[1]+'\t'+each_insert[2]+'\t'+each_insert[3]+'\t'+each_insert[4]+'\t'+each_insert[5]+'\t'+each_insert[6]+'\t'+each_insert[7]+'\t'+each_insert[8]+'\n')

	wf.close()


annotation_file = argv[1]
file_to_clean = argv[2]

annotations = parse_files(annotation_file, remove_header = False)
file_data = parse_files(file_to_clean, remove_header = True)

genes_to_chromosomes = map_genes_to_sc_chr(annotations)
clean_file(file_data, genes_to_chromosomes)




