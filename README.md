# rh-seq
Weiss et al RH-seq code repository
This repository houses the final code used in analyzing (link to paper)

Input: single end .fastq files from an RH-seq experiment, filtered according to JGI standards.
First, map all fastqs with map_and_pool_BLAT-1-27-16.py.
Then, run each script in order from 1-9, using the output of the previous script as input into the next script.


