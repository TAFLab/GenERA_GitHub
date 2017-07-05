#!/usr/bin/env python
#__author__ = 'quentinf'

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||| GenERA PIPELINE 

# NAME: GenERA_1.py

# OPERATIONS:
# - reads gDNA and cDNA SAM files for a list of amplicons provided in "./META/input.csv" file
# - maps read deletions using the CIGAR
# - counts the number of reads having a deletion at nucleotide: nucleotide deletion profile
# - saves in "./META/nucleotide_deletion_pattern.csv"

# INPUT:
# - "./META/input.csv": contains a line for each amplicons:
# (format of each line) gDNA SAM file, cDNA SAM file, genome sequence, seed start, seed end,MRE start, MRE end

# OUTPUT:
# - "./META/nucleotide_deletion_pattern.csv": contains a table with the following columns
#   - Sample type (gDNA, cDNA)
#   - position (nucleotide along the reference genome)
#   - deletion stack (one column per amplicon): counts the number of reads with deletion at nucleotide i (specified by the position column).

# - "./META/gDNA_cDNA_pairs.csv": contains meta information for each amplicon
#   - Name
#   - sam_gDNA
#   - sam_cDNA
#   - ref_genome
#   - read_count_gDNA
#   - read_count_cDNA
#   - seed_start
#   - seed_end

# DEPENDENCIES: GenERA_functions.py

# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

from GenERA_functions import *

INPUT_FILE = "../META/input.csv"
OUTPUT_FILE = "../META/nucleotide_deletion_pattern.csv"
OUTPUT_FILE_2 = "../META/gDNA_cDNA_pairs.csv"
PATH_DATA = "../DATA/SAM/" #needs to be changed depending on the system architecture

NAMES = list()
NAMES_gDNA = list()
NAMES_cDNA = list()
L_REF_GENOME = list()
SEED = list()

DELETION_DENT_gDNA = list()
READ_COUNT_gDNA = list()

DELETION_DENT_cDNA = list()
READ_COUNT_cDNA = list()

# reads input file and extract information:
with open(INPUT_FILE, "r") as fo_input:
        for line in fo_input: # read line by line
            print(line)
            # processes line and save information:
            line_split = line.split("\n")
            line_split = line_split[0].split(",")
            FILENAME_INPUT_gDNA = PATH_DATA + line_split[0]
            FILENAME_INPUT_cDNA = PATH_DATA + line_split[1]
            REF_GENOME = line_split[2]
            SEED.append([line_split[3],line_split[4]])

            NAMES_gDNA.append(line_split[0])
            NAMES_cDNA.append(line_split[1])
            L_REF_GENOME.append(REF_GENOME)

            line_split = line_split[0].split(".")
            NAMES.append(line_split[0])

            # computes the nucleotide deletion profile:
            deletion_dent_gDNA = mreCRISPR_deletion_dent(FILENAME_INPUT_gDNA,REF_GENOME)
            deletion_dent_cDNA = mreCRISPR_deletion_dent(FILENAME_INPUT_cDNA,REF_GENOME)


            READ_COUNT_gDNA.append(deletion_dent_gDNA['count'])
            DELETION_DENT_gDNA.append(deletion_dent_gDNA['dent'])

            READ_COUNT_cDNA.append(deletion_dent_cDNA['count'])
            DELETION_DENT_cDNA.append(deletion_dent_cDNA['dent'])

fo_input.close()

genome_length = len(DELETION_DENT_gDNA[0])
print(genome_length)


#///////////////////////////// writes nucleotide deletion profile information in output file
print('writing normalised dent count: count/nb_read*mean(nb_read)')

fo_ouput = open(OUTPUT_FILE,'w')

# writes header
header = 'SampleType,Position'
for i in range(0,len(NAMES)):
    header += ',' + NAMES[i]

header += '\n'
fo_ouput.write(header)

# writes nucleotide deletion profile information
for j in range(0,genome_length):

    line = 'gDNA,' + str(j+1) #/!\ to be compatible with R indexing
    for i in range(0,len(DELETION_DENT_gDNA)):
        # use if the nucleotide deletion profile is to be normalised
        # dent = float(DELETION_DENT_gDNA[i][j])/float(READ_COUNT_gDNA[i])*(float(READ_COUNT_gDNA[i])+float(READ_COUNT_cDNA[i]))/2
        dent = DELETION_DENT_gDNA[i][j]
        line += ',' + str(dent)
    line += '\n'
    fo_ouput.write(line)

for j in range(0,genome_length):

    line = 'cDNA,' + str(j+1) #/!\ to be compatible with R indexing
    for i in range(0,len(DELETION_DENT_cDNA)):
        # use if the nucleotide deletion profile is to be normalised
        # dent = float(DELETION_DENT_cDNA[i][j])/float(READ_COUNT_cDNA[i])*(float(READ_COUNT_gDNA[i])+float(READ_COUNT_cDNA[i]))/2
        dent = DELETION_DENT_cDNA[i][j]
        line += ',' + str(dent)
    line += '\n'
    fo_ouput.write(line)
fo_ouput.close()
#/////////////////////////////

#///////////////////////////// writes meta data
fo_ouput = open(OUTPUT_FILE_2,'w')

header = 'Name,sam_gDNA,sam_cDNA,ref_genome,read_count_gDNA,read_count_cDNA,seed_start,seed_end\n'
fo_ouput.write(header)

for i in range(0,len(NAMES)):
    line  = str(NAMES[i]) + ',' +\
            str(NAMES_gDNA[i]) + ',' + \
            str(NAMES_cDNA[i]) + ',' + \
            L_REF_GENOME[i] + ',' + \
            str(READ_COUNT_gDNA[i]) + ',' +\
            str(READ_COUNT_cDNA[i]) + ',' +\
            str(SEED[i][0]) + ',' +\
            str(SEED[i][1]) + '\n'
    fo_ouput.write(line)

fo_ouput.close()
#/////////////////////////////