#!/usr/bin/env python
#__author__ = 'quentinf'

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||| GenERA PIPELINE

# NAME: GenERA_3.py

# FUNCTION: For each amplicon extract from SAM files wild-type and all UDP (unique deletion patterns) in the zone
# - Using zone coordinates, extract for each amplicon the list of unique cigars and their corresponding read counts
# for both gDNA and cDNA
# - cDNA and gDNA cigars are then matched if possible and compiled into a common file

# INPUT:
# - "./META/gDNA_cDNA_pairs.csv": see GenERA_1.py output
# - gDNA SAM files
# - cDNA SAM files
# - "./META/ROI.csv": ROI delineation for each amplicon see GenERA_2.R output

# OUTPUT:
# - './DATA/CSV/A_RAW_CSV_FILES/'+ L_file_pairs[i][0] + '_GenERA.csv': csv file cataloguing the cDNA and gDNA cigars and their respective read counts (matched when possible)

# DEPENDENCIES: mreCRISPR_functions.py

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||| MAIN

from GenERA_functions import *

# reads gDNA,cDNA sam files
INPUT_FILE = "../META/gDNA_cDNA_pairs.csv"

L_file_pairs = list() # contains all pair of gDNA, cDNA sam files

with open(INPUT_FILE, "r") as fo_input:
        for line in fo_input: # reads line by line
            line_split = line.split('\n')
            V_name_gDNA_cDNA = line_split[0].split(',')

            if V_name_gDNA_cDNA[0] != 'Name': # avoid .SAM header
                L_file_pairs.append(V_name_gDNA_cDNA)# line from gDNA_cDNA pairs file
                print(V_name_gDNA_cDNA)


# reads zone for each pairs
INPUT_FILE = "../META/ROI.csv"
L_zone_name = list()
L_zone = list() # contains all pair of gDNA, cDNA sam files

with open(INPUT_FILE, "r") as fo_input:
        for line in fo_input: # reads line by line
            line_split = line.split('\n')
            line_split = line_split[0].split(',')

            if line_split[1] != 'name': # avoid .SAM header
                L_zone_name.append(line_split[1])
                L_zone.append([int(line_split[2]),int(line_split[3])])

print(len(L_file_pairs),len(L_zone),len(L_zone_name))

# for each file pairs extract all Unique Deletion Patterns
PATH_DATA = "../DATA/SAM/"

for i in range(0,len(L_file_pairs)):
    print('index:',i)
    if L_file_pairs[i][0] == L_zone_name[i]: # check file compatibility

        filename_gDNA = PATH_DATA + L_file_pairs[i][1]
        filename_cDNA = PATH_DATA + L_file_pairs[i][2]
        REF_GENOME = L_file_pairs[i][3]
        ZONE = L_zone[i]
        #/!\ index start at zero
        #/!\ a:b stop at b-1
        ZONE[0] += -1


        print('gDNA file:',filename_gDNA)
        print('cDNA file:',filename_cDNA)
        print('Genome:', REF_GENOME)
        print('Genome length:', len(REF_GENOME))


        SEED = [int(L_file_pairs[i][6])-1,int(L_file_pairs[i][7])] #/!\ specifies the seed
        print(SEED)

        # computes list of unique cigars in zone with their associated read counts for gDNA and then cDNA
        DATA_gDNA = mreCRISPR_cigars_and_paired_count_full(filename_gDNA,REF_GENOME,ZONE,SEED)
        print(len(DATA_gDNA['CIGARs']),len(DATA_gDNA['CIGARs_seed']),len(DATA_gDNA['CIGARs_count']))
        DATA_cDNA = mreCRISPR_cigars_and_paired_count_full(filename_cDNA,REF_GENOME,ZONE,SEED)
        print(len(DATA_cDNA['CIGARs']),len(DATA_cDNA['CIGARs_seed']),len(DATA_cDNA['CIGARs_count']))


        OUTPUT_FILE = '../DATA/CSV/A_RAW_CSV_FILES/'+ L_file_pairs[i][0] + '_GenERA.csv'
        mreCRISPR_pair_and_mirScore_new(DATA_gDNA,DATA_cDNA, OUTPUT_FILE)