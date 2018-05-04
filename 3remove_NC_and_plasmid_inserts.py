from sys import argv
import pandas as pd

### USAGE ###
# python remove_NC_and_plasmid_inserts.py output_directory/ file_location/*.normalized_pooled_reads
# keeps inserts only in genic regions (removes any annotated as non-coding or plasmid)
# also saves another .txt file with summary of how many inserts were in each category

def parse_file(filename, sep='\t'): # reads in file, returning a pandas dataframe
	
	#reading in file as a pandas dataframe and pulls out file ID
	df = pd.read_csv(filename, sep=sep)
	file_id = df['tech_rep_name'][1]

	return df, file_id

def find_coding_inserts(df):  
# removes all inserts with annotations 'NC' and 'plasmid', saves info about how many were removed to a new text file
	not_plasmid = df['annotation'] != 'plasmid'
	not_NC = df['annotation'] != 'NC'
	coding_rows = df[not_plasmid & not_NC]
	total_inserts = len(df) - 1.0
	total_coding = len(coding_rows) - 1.0
	total_noncoding = total_inserts - total_coding

#write info on inserts removed to a new txt file stored in the output directory
	wf.writelines('In '+str(file_id)+', the total number of inserts is: '+str(total_inserts)+', the number of coding inserts is: '+str(total_coding)+' and the number of non-coding inserts removed was: '+str(total_noncoding)+str('\n'))
	
	return coding_rows

output_folder = argv[1].strip('/')
files = argv[2:]
wf = open(str(output_folder)+'/coding_vs_noncoding_stats.txt', 'w')
for each_file in files:
	data, file_id = parse_file(each_file) #parse
	coding_inserts = find_coding_inserts(data) #filter
	coding_inserts.to_csv(str(output_folder)+str('/')+str(file_id)+'.normalized_pooled_reads_coding', sep='\t', index=False) #save coding inserts to new file
wf.close()


