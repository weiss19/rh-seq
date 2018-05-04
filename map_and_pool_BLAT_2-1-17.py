import regex
import numpy as np
import sys
import subprocess as sp

# PARAMETERS # 


mapping_genome = "/usr2/people/carlyweiss/Amended_Genomes/Amended_Genomes/Concatenated/D1373_Z1.fa"  # location of the .udb database
min_ID = str(95) # minimum identity cutoff for mapping reads.  Must be a decimel proportion (0-1)
min_id_pooling = .95 # minimum identity for reads at the pooling stage.  Can be different then at mapping.  Must be (0-1)
gff = '/usr2/people/carlyweiss/SparScerinfo/YS2+CBS432+plasmid_clean' ## GFF to use for annotations.  
min_read_cutoff = 0  ## this is only for reporting summary statistics 

min_id_pooling*=100

# HELP #

if len(sys.argv) == 1:
    print "USAGE: python map_PB_reads.py fastq_file out_directory"
    exit()
# INPUT # 

read_file = sys.argv[1]  
fastq_filename = read_file.split("/")[-1]

out_dir = sys.argv[2]

# BEGIN FUNCTIONS # 

def Filter_Reads(fastq_file):  ## filter out reads that do not have a TN sequence 

    wf = open(out_dir+fastq_filename+'_parsed_reads','w')  # outfile for the trucated reads that will be mapped

    tn_pattern = regex.compile('(?e)(CAGACTATCTTTCTAGGGTTAA){e<=2}')  # last 2 bps of the left arm of the transposon.  searches for this pattern allowing 2 mismatchs for sequencing errors

    # counts for summary stats #

    tot_reads = 0.0  
    reads_with_tn = 0
    line_count = 0
    total_Ns = 0
    total_bases = 0
    reads_too_short = 0
    too_high_Ns = 0
    f = open(fastq_file) # the file with the reads
    
    for line in f:

        # parses fastq file # 

        line_count+=1
        if line_count % 4 == 1:
            header = line
        elif line_count % 4 == 2:
            read = line.strip()
        elif line_count % 4 == 0:
            qual = line
            tot_reads+=1

            tn_match_data = tn_pattern.search(read)  # searches for the tn pattern
            if tn_match_data == None:
                continue  #skips the read if it doesn't have a tn pattern

            total_bases+=len(read)
            total_Ns+=read.count('N') # counts the numer of Ns per read
            
            reads_with_tn+=1  ## counts reads that have the tn sequence

            end_match = tn_match_data.end()  # finds the location in the read that the last tn base matches to
    
            genome_read =  read[tn_match_data.end()-4:]  # gets the sequence  of the read outside the transposon (still contains TTAA)

            if len(genome_read) < 50: # if the length of the parsed read is less than 50, throw it out.  This is appropriate for 150bp reads and should be changed otherwise
                reads_too_short+=1
                continue 
            
            if genome_read.count('N') / float(len(genome_read)) > .2:
                
                too_high_Ns+=1
                continue

        
            # writes a new fastq file that has the filtered read and only the regions TTAA--genome.  The TN sequence is removed # 

            wf.writelines(">"+header)
            wf.writelines(genome_read+"\n")
           # wf.writelines("+\n")
           # wf.writelines(qual[tn_match_data.end()-4:])
            
    wf.close()    

    prop_Ns = float(total_Ns) / total_bases # gets the proportion of Ns in the fastq file as an indicator for overall sequence quality


    wf = open(out_dir+fastq_filename+"_mapping_stats", 'w')  # the file that will contain summary stats of the mapping and additional analysis

    # write mapping statistics" 

    print "total_reads: ", str(tot_reads)    
    wf.writelines("total_reads: "+str(tot_reads)+"\n")
    print "reads with tn: ", str(reads_with_tn), " ("+str(100*reads_with_tn/tot_reads)+"%)"
    wf.writelines("reads with TN: "+str(reads_with_tn)+" ("+str(100*reads_with_tn/tot_reads)+"%)\n")
    print "proportion bases in TN containing reads that are 'N': ", str(prop_Ns)
    wf.writelines("proportion bases in TN containing reads that are 'N': "+ str(prop_Ns)+"\n")
    print "TN contaning reads too short that were discarded: ", str(reads_too_short), " ("+str(100*reads_too_short/float(reads_with_tn))+"%)"
    wf.writelines("TN contaning reads too short that were discarded: "+str(reads_too_short)+" ("+str(100*reads_too_short/float(reads_with_tn))+"%)\n")
    print "TN containing reads with too many Ns: ", str(too_high_Ns), " ("+str(100*too_high_Ns/float(reads_with_tn))+"%)"
    wf.writelines("TN containing reads with too many Ns: "+str(too_high_Ns)+" ("+str(100*too_high_Ns/float(reads_with_tn))+"%)\n")
 


    wf.close()

    tot_parsed_reads = reads_with_tn - reads_too_short - too_high_Ns

    return tot_parsed_reads


def Map_Reads(): # uses usearch to map the reads

#    cmd = ["/opt/bin/usearch8.1.1756_i86linux64", "-threads", nthreads, "-usearch_local",  out_dir+fastq_filename+"_parsed_reads", "-db", db, "-strand", "both", "-id", min_id, "-alnout", out_dir+fastq_filename+"_aligned_reads", "-maxhits", "2", "-notmatched", out_dir+fastq_filename+"_unmapped_reads", "-userout", out_dir+fastq_filename+"_mapped_reads", "-userfields", "query+target+id+alnlen+qstrand+tstrand+tlot+thit", "-mincols", "50", "-query_cov", ".95"]  # the command used to call usearch

    cmd = ["blat", mapping_genome, out_dir+fastq_filename+"_parsed_reads", out_dir+fastq_filename+"_mapped_reads", "-minIdentity=95", "-tileSize=12", "-out=blast8" ]
    
    sp.call(cmd)

def Filter_for_ID_and_multimapping_and_pool(num_parsed_reads):  # filters our mapped reads that map below the ID threshold, those that map to multiple locations, then pools reads into insertion sites

## min_id_pooling

    mapped_file = out_dir+fastq_filename+"_mapped_reads" # where the maped reads are

    first = 'y'

    f = open(mapped_file)


    read_dict = {}
    total_mapped_reads = 0 ###  CHANGE THIS BACK TO 0
    reads_above_identity_threshold = 0
    mulitmapped_reads = 0
    
    all_read_list = set()
    for line in f:

        
        line = line.strip().split("\t")

        
        read = line[0]
        
        if read not in all_read_list:  ## counts all reads that map at least once
            all_read_list.add(read)
            total_mapped_reads+=1


        ID = float(line[2])
        if ID < min_id_pooling:  ## filters out mapped reads below the minimum identity for pooling                                                                                                
            continue

        len_align = float(line[3])
        start_align = float(line[6]) 

        if len_align < 50 or start_align > 3:
            continue

#        print line


        reads_above_identity_threshold+=1  # this counts reads that are above the identity threshold and length of aligment threshold and in which the alignment starts within 3pbs of the read TTAA

        start = line[8]
        end = line[9]

        strand = '+'
        if start < end:
            strand = '-'
        

        insertion_loc = start
        
        loc = line[1]+"__"+strand+"__"+insertion_loc  ## This is the identifier for the insertion. scaffold+strand+position

        # searches for reads that map to more than 1 locatiion # 

        if read in read_dict:  
            read_dict[read].append([loc, ID])  ## appends the read location and %ID of mapping in case of multiple mappings of the same read
        #    mulitmapped_reads+=2  # as long as the usearch multimapping paramter is set to a maximim of 2 reads, this will work.  
        else:
            read_dict[read] = [[loc, ID ]]



    # filters our reads that map to more than one location # 


    print "done making the read_dict"

    loc_dict = {}  # 
    for read in read_dict.keys():
        if len(read_dict[read]) > 2:
            del read_dict[read]
            continue
        if len(read_dict[read]) == 2:

 #           print read_dict[read]


            if read_dict[read][0][1] == read_dict[read][1][1]:
                del read_dict[read]
                continue
            if read_dict[read][0][1] > read_dict[read][1][1]:
                read_dict[read] = [read_dict[read][0]]
            elif read_dict[read][0][1] < read_dict[read][1][1]:
                read_dict[read] = [read_dict[read][1]]

  #          print read_dict[read]
                
            
#        print read_dict

        read_dict[read] = read_dict[read][0][0]  ## appends just the insertion location

 #       print read_dict[read]



        if read_dict[read] in loc_dict:
            loc_dict[read_dict[read]]+=1
        else:
            loc_dict[read_dict[read]] = 1


  #  for loc in loc_dict:
  #      print loc, loc_dict[loc]

            

    print "total_mapped reads: ", str(total_mapped_reads), " ("+str(float(100*total_mapped_reads/num_parsed_reads))+"% of parsed reads)"
    print "mapped reads passing identity cutoff: "+str(reads_above_identity_threshold), "("+str(100*reads_above_identity_threshold/float(total_mapped_reads))+"%)"
    print "remaining reads mapping to one location: "+str(len(read_dict)), "("+str(100*len(read_dict)/float(reads_above_identity_threshold))+"%)"

    wf = open(out_dir+fastq_filename+"_mapping_stats", 'a')

    wf.writelines("total mapped reads: "+str(total_mapped_reads)+" ("+str(float(100*total_mapped_reads/num_parsed_reads))+"% of parsed reads)\n")
    wf.writelines("mapped reads passing identity cutoff: "+str(reads_above_identity_threshold)+" "+str(100*reads_above_identity_threshold/float(total_mapped_reads))+"% of mapped reads\n")
    wf.writelines("remaining reads mapping to one location: "+str(len(read_dict))+" "+str(100*len(read_dict)/float(reads_above_identity_threshold))+"%\n")

    wf.close()

    return loc_dict, len(read_dict)


def Combine_near_mappings(loc_dict, reads_remaining):  # # combines insertion sites that are within 3 bases of each other.  reads are assigned to the site with the initial max number of reads

    split_loc_dict = {}  # will hold hierarchical data on each insertion site.  keys for nested dictionaries are scaffold, strand, position and value is # reads mapping there. 

    for full_location in loc_dict:  ## loc dict holds the identifier for an inserion site as key and reads mapping to that site as a value
        chrom = full_location.split("__")[0]
        strand = full_location.split("__")[1]
        pos = int(full_location.split("__")[2])

        # initialize the dictionary #

        if chrom not in split_loc_dict:
            split_loc_dict[chrom] = {'+' : {}, '-' : {}}

        if pos not in split_loc_dict[chrom][strand]:
            split_loc_dict[chrom][strand][pos] = loc_dict[full_location]

    reads_moved = 0
    
    # sorts the insertion positions, and combines reads forward, then reverses the sorting and combines forward again.  #

    for chrom in split_loc_dict:
        for strand in split_loc_dict[chrom]:

            sorted_positions = sorted(split_loc_dict[chrom][strand])
            first ='y'
            for pos in sorted_positions:
                if first == 'y':
                    first = 'n'
                    last = pos
                    continue

                if int(pos) - int(last) < 4:

                    if split_loc_dict[chrom][strand][pos] >= split_loc_dict[chrom][strand][last]:
                        split_loc_dict[chrom][strand][pos]+=split_loc_dict[chrom][strand][last]
                        reads_moved+=split_loc_dict[chrom][strand][last]
                        del split_loc_dict[chrom][strand][last]

                last = pos
            sorted_positions = sorted(split_loc_dict[chrom][strand])
            sorted_positions.reverse()

            first ='y'
            for pos in sorted_positions:
                if first == 'y':
                    first = 'n'
                    last = pos
                    continue

                if abs(int(pos) - int(last)) < 4:
                    if split_loc_dict[chrom][strand][pos] >= split_loc_dict[chrom][strand][last]:
                        split_loc_dict[chrom][strand][pos]+=split_loc_dict[chrom][strand][last]
                        reads_moved+=split_loc_dict[chrom][strand][last]
                        del split_loc_dict[chrom][strand][last]

                last = pos

    print "remaining reads moved to higher peak: ", str(reads_moved), "("+str(100*float(reads_moved)/reads_remaining)+"%)"
    

    wf = open(out_dir+fastq_filename+"_mapping_stats", 'a')
    wf.writelines("remaining reads moved to higher peak: "+str(reads_moved)+" "+str(100*float(reads_moved)/reads_remaining)+"%\n")
    wf.close()

    return split_loc_dict

def Annotate_insetions(split_loc_dict, mapped_reads):

    out_filename = out_dir+fastq_filename+"_pooled_reads"  # the final output, this will hold the pooled insertion table                                                                      
    wf = open(out_filename,'w')

    wf.writelines("ID\tscaffold\tstrand\tlocation\tannotation\tn\trel_loc\trel_prop\tgene_length\n")


    sc_insertions = 0
    sp_insertions = 0

    sc_genic_insertions = 0
    sp_genic_insertions = 0

    tot_insertions = 0
    tot_min_insertions = 0.0
    plasmid_reads = 0 # reads mapping to the plasmid
    Rtn_reads = 0 # reads mapping to the tn right border
    Ltn_reads = 0 # reads mapping to the tn left border

    gff_dict = {}
    f = open(gff)

    for line in f:

        line = line.split("\t")
        chrom = line[1]
        gene = line[0]
        start = int(line[4])
        end = int(line[5])
        strand = line[3]
        type = line[2]

        # put the annotation information in a dictionary #                                                                                                                                    

        if chrom not in gff_dict:
            gff_dict[chrom] = {}

        gff_dict[chrom][gene] = [start, end, strand]


# Search through the dictionary for insertions that fall within genes                                                                                                                     
    for chrom in split_loc_dict:
        for strand in split_loc_dict[chrom]:
            for pos in split_loc_dict[chrom][strand]:

                if split_loc_dict[chrom][strand][pos] > min_read_cutoff:
                    tot_min_insertions+=1
                    if chrom[:2] == 'sc':
                        sc_insertions+=1
                    elif chrom[:2] == 'sp':
                        sp_insertions+=1
                tot_insertions+=1

                # set defaults for noncoding                                                                                                                                                 

                insertion_type = 'NC'
                gene_length = -1
                relative_insertion_site = -1
                prop_gene_insertion_in = -1

                for gff_chrom in gff_dict:
                    if gff_chrom != chrom:
                        continue
                    for gene in gff_dict[gff_chrom]:
                        if pos >= gff_dict[gff_chrom][gene][0] and pos <= gff_dict[gff_chrom][gene][1]:  # if the insertion falls within a gene                                              
                            insertion_type = gene
                            gene_length = gff_dict[gff_chrom][gene][1] - gff_dict[gff_chrom][gene][0]+1
                            if gff_dict[gff_chrom][gene][2] == '+':
                                relative_insertion_site = pos - gff_dict[gff_chrom][gene][0]

                            else:
                                relative_insertion_site = gff_dict[gff_chrom][gene][1] - pos

                            if gene[:2] == 'sp' and split_loc_dict[chrom][strand][pos] > min_read_cutoff:
                                sp_genic_insertions+=1
                            elif gene[:2] == 'sc'  and split_loc_dict[chrom][strand][pos] > min_read_cutoff:
                                sc_genic_insertions+=1

                            prop_gene_insertion_in = relative_insertion_site / float(gene_length)

                wf.writelines(chrom+"_"+strand+"_"+str(pos)+"\t"+chrom+"\t"+strand+"\t"+str(pos)+"\t"+insertion_type+"\t"+str(split_loc_dict[chrom][strand][pos])+"\t"+str(relative_insertion_site)+"\t"+str(prop_gene_insertion_in)+"\t"+str(gene_length)+"\n")

    wf.close()



    
    f = open(out_filename)
    for line in f:
        if line[:2] != 'pl':
            continue

        line = line.split("\t")
        if line[4] == 'TN_right_arm':
            Rtn_reads+=int(line[5]) # reads mapping to the tn right border        
        elif line[4] == 'TN_left_arm':

            Ltn_reads+=int(line[5]) # reads mapping to the tn right border        
        elif line[4] == 'NC':
            plasmid_reads+=int(line[5])

    f.close()

    
    tot_genic_insertions = float(sc_genic_insertions+sp_genic_insertions)

    wf = open(out_dir+fastq_filename+"_mapping_stats", 'a')
    
    wf.writelines("reads mapping to NC plasmid backbone: "+str(plasmid_reads)+" ("+str(100*plasmid_reads/mapped_reads)+"%) of mapped reads\n")
    wf.writelines("reads mapping to TN right border: "+str(Rtn_reads)+" ("+str(100*Rtn_reads/mapped_reads)+"%) of mapped reads\n")
    wf.writelines("reads mapping to TN left border: "+str(Ltn_reads)+" ("+str(100*Ltn_reads/mapped_reads)+"%) of mapped reads\n")

    wf.writelines("total insertions: "+str(tot_insertions)+"\n")
    wf.writelines("total insertions with >"+str(min_read_cutoff)+" reads: "+str(tot_min_insertions)+"\tscer: "+str(sc_insertions)+" ("+str(100*sc_insertions/tot_min_insertions)+"%)"+" spar: "+str(sp_insertions)+"("+str(100*sp_insertions/tot_min_insertions)+"%)\n")


    if tot_genic_insertions > 0:
        


        wf.writelines("OF  THESE:\n")
        wf.writelines("total genic insertions: "+str(sp_genic_insertions + sc_genic_insertions)+" ("+str(100*(sp_genic_insertions + sc_genic_insertions)/tot_min_insertions)+"% of insertions)\n")
        wf.writelines("Scer genic insertions: "+str(sc_genic_insertions)+" ("+str(100*sc_genic_insertions/tot_genic_insertions)+"% of genic insertions)\n")
        wf.writelines("Spar genic insertions: "+str(sp_genic_insertions)+" ("+str(100*sp_genic_insertions/tot_genic_insertions)+"% of genic insertions)\n")

    else:
        wf.writelines("no genic insertions\n")


    wf.close()



                  

#### START PROGRAM ####

num_parsed_reads = Filter_Reads(read_file)  ## filters out reads that don't have tn sequence, writes the genomic portion of the remaining reads to a new file
Map_Reads()  ## maps reads

print "done mapping"

loc_dict, reads_remaining = Filter_for_ID_and_multimapping_and_pool(num_parsed_reads)  # filters out reads below the identity threshold, that map to multiple locations and then pools insertions 


print "done filtering"
loc_dict = Combine_near_mappings(loc_dict, reads_remaining) # combine insertions within 3bp of each other

print "done combine near mappings"

Annotate_insetions(loc_dict, reads_remaining) # identify insertions in genes, and write the final pooled outfile








