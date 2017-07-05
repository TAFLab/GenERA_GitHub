#!/usr/bin/env python
#__author__ = 'quentinf'

# FUNCTIONS ////////////////////////////////////////////////////////////////////////////////////////////////////////////
# test if a character is a number
def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def vector_sum(x):
    s = 0
    for i in range(0,len(x)):
        s += x[i]
    return s

# this function takes a string cigar "10M20D7M" and converts it to arrays of letters ['M','D','M'] and number [10,20,7]
def mreCRISPR_parse_cigar(input_cigar):
    X = list(input_cigar) # splits the string into array of char

    list_num = []
    list_letter = []

    idx_last_letter = -1;
    idx_current = 0;

    for x in X:
        if not is_number(x):
            list_letter.append(x)
            idx_1 = idx_last_letter + 1
            idx_2 = idx_current
            list_num.append(int(''.join(X[idx_1:idx_2])))

            idx_last_letter = idx_current
        idx_current += 1

    return {'list_letter': list_letter, 'list_num': list_num, 'canonical': set(list_letter) == set(['M', 'D', 'M'])}


# this function takes a .sam read and registers to the reference genome
# '-' not sequenced nucleotides
# '*' deleted nucleotides
# if the cigars contains H or I the variable 'read_ok' will be False
# the returned read should be the same length as the reference genome
def mreCRISPR_register_read(sam_position,sam_cigar,sam_read, genome_length):

    sam_cigar_parsed = mreCRISPR_parse_cigar(sam_cigar)
    registered_read = '-' * sam_position
    read_ok = True
    read_contains_deletion = False
    read_weird = False

    cursor_read = 0
    cursor_registered_red = sam_position;

    for i in range(0,len(sam_cigar_parsed['list_letter'])):

        type = sam_cigar_parsed['list_letter'][i]

        if type == 'M':
            registered_read += sam_read[cursor_read: cursor_read + sam_cigar_parsed['list_num'][i]]
            cursor_read += sam_cigar_parsed['list_num'][i]
        elif type == 'D':
            registered_read += '*' * sam_cigar_parsed['list_num'][i]
            read_contains_deletion = True
        elif type == 'S':
            cursor_read += sam_cigar_parsed['list_num'][i]
            read_weird = True
        else:
            read_ok = False

    registered_read += '-' * (genome_length - len(registered_read))
    return{'read':registered_read,'read_ok':read_ok, 'deletion':read_contains_deletion}


# This function computes the cigar in the zone of interest
def mreCRISPR_cigar_in_zone(ref_zone_seq, registered_read_zone):

    cigar = list() # cigar letters
    cigar_count = list() # cigar counts

    # example:
    # cigar ['+','D','-']
    # cigar_count [10,7,23]

    N = len(ref_zone_seq)
    type = ''
    count = 0
    for i in range(0,N):
        if registered_read_zone[i] == ref_zone_seq[i]:
            if type == '':
                type = '+'
                count = 1
            elif type == '+':
                count += 1
            else:
                cigar.append(type)
                cigar_count.append(count)
                type = '+'
                count = 1
        elif registered_read_zone[i] == '*':
            if type == '':
                type = 'D'
                count = 1
            elif type == 'D':
                count += 1
            else:
                cigar.append(type)
                cigar_count.append(count)
                type = 'D'
                count = 1
        elif registered_read_zone[i] == '-':
            if type == '':
                type = '#'
                count = 1
            elif type == '#':
                count += 1
            else:
                cigar.append(type)
                cigar_count.append(count)
                type = '#'
                count = 1
        else:
            if type == '':
                type = '-'
                count = 1
            elif type == '-':
                count += 1
            else:
                cigar.append(type)
                cigar_count.append(count)
                type = '-'
                count = 1
    cigar.append(type)
    cigar_count.append(count)

    cigar_txt = ''
    D_max = 0
    D_max_txt = ''
    D_max_num = [0,0]
    D_count = 0
    full_read = True

    for i in range(0,len(cigar)):
        cigar_txt += str(cigar_count[i]) + cigar[i]

        if cigar[i] == '#':
            full_read = False

        # if cigar[i] == 'D':
        #     D_count += 1
        #     if (i>0) & (i<(len(cigar)-1)):
        #         if (cigar[i-1] == '+') & (cigar[i+1] == '+'):
        #             if cigar_count[i] > D_max:
        #                 D_max_txt = str(sum(cigar_count[0:i])) + 'D' + str(sum(cigar_count[0:i])+cigar_count[i])
        #                 D_max_num = [sum(cigar_count[0:i]),sum(cigar_count[0:i])+cigar_count[i]]
        #                 D_max = cigar_count[i]



    # print(cigar,cigar_count,cigar_txt,D_max_txt)
    return{'cigar':cigar,
           'cigar_count':cigar_count,
           'cigar_txt':cigar_txt,
           'full_read':full_read,
           'D_max':D_max_txt,
           'D_count':D_count,
           'length_D': D_max,
           'D_max_num': D_max_num}


# This function sets the condition to select the reads
def mreCRISPR_condition(registered_cigar_zone):

    # reminder, the input contains:
    # {'cigar':cigar,'cigar_count':cigar_count,'cigar_txt':cigar_txt,'D_max':D_max_txt,'D_count':D_count,'length_D': D_max}


    # should return True or False
    is_the_read_good = False

    # assess if the read is to be used
    # is_the_read_good = registered_cigar_zone['cigar'] == ['+','D','+']
    # is_the_read_good = registered_cigar_zone['cigar_txt'] == '31+26D13+'
    # is_the_read_good = registered_cigar_zone['D_count'] == 1
    # is_the_read_good = (registered_cigar_zone['D_max_num'][1]-registered_cigar_zone['D_max_num'][0]) == 1

    #grab all clean deletions in zone 1 with strictly inside zone 2 (clean deletions that only delete inside the seed)
    if registered_cigar_zone['cigar'] == ['+','D','+']:
        #print(registered_cigar_zone['cigar']); print(registered_cigar_zone['cigar_count'])
        X = registered_cigar_zone['D_max_num']
        if (20<=X[1]<=27) or (20<= X[0]<=27) or (X[0]<=20 and X[1]>=27):
            is_the_read_good = True
    #         print("<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>")

    # grab all deletions that touch the seed
    # if registered_cigar_zone['cigar'] == ['+','D','+']:
    #     print(registered_cigar_zone['cigar']); print(registered_cigar_zone['cigar_count'])
    #     X = registered_cigar_zone['D_max_num']
    #     if (20 <= X[1] <= 51) or (20 <= X[0] <= 51) or ((X[0] <= 20) and (X[1] >= 51)):
    #         is_the_read_good = True
    #        print("<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>")

    # gral all the deletions in zone_1 without overlap the zone 2 (deletions that doesn't touch the seed)
    # if registered_cigar_zone['cigar'] == ['+','D','+']:
    #     X = registered_cigar_zone['D_max_num']
    #     if (X[1]<20) or (X[0]>27):
    #         is_the_read_good = True

    # gral all the wild type reads (pure WT)
    #if registered_cigar_zone['cigar'] == ['+']:
    #     print(registered_cigar_zone['cigar']); print(registered_cigar_zone['cigar_count'])
    #     is_the_read_good = True

    # return the verdict
    return is_the_read_good



# MAIN FUNCTIONS ///////////////////////////////////////////////////////////////////////////////////////////////////////

# This function computes the nucleotide deletion profile
def mreCRISPR_deletion_dent(FILENAME_INPUT,REF_GENOME):

    # FILENAME_INPUT: sam file
    # REF_GENOME: sequence of the amplicon

    ref_genome_length = len(REF_GENOME) # size of the ref genome
    DELETION_dent = [0]*len(REF_GENOME)
    READ_count = 0

    # reads a line from the .sam file
    # use scigar, position, read to build a registered read
    with open(FILENAME_INPUT, "r") as fo_input:
        for line in fo_input: # reads line by line

            line_split = line.split("\t") # converts the string "line" to an array of string

            if len(line_split) > 10: # avoid the .SAM header

                READ_count += 1

                # grad information from line:
                sam_position = int(line_split[3]) - 1  # python indexes start at 0
                sam_cigar = line_split[5]
                sam_read = line_split[9]

                # computes the registered read
                OUTPUT = mreCRISPR_register_read(sam_position,sam_cigar,sam_read,ref_genome_length)
		# print(OUTPUT['read'])

                # updates the nucleotide deletion profile
                for i in range(0,len(REF_GENOME)): #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< NEW
                    if OUTPUT['read'][i] == '*': #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< NEW
                        DELETION_dent[i] += 1 #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< NEW

    fo_input.close()
    return {'count':READ_count, 'dent': DELETION_dent}




# this function runs the program: extracts reads and counts the cigars
def mreCRISPR_cigars_and_count(FILENAME_INPUT,REF_GENOME,ZONE):

    # FILENAME_INPUT: sam file
    # REF_GENOME: sequence of the amplicon
    # ZONE: [a,b] boundary of the region of interest

    ref_genome_length = len(REF_GENOME) # size of the ref genome
    ref_zone_seq = REF_GENOME[ZONE[0]:ZONE[1]] # sequence of the ROI
    print(ref_zone_seq)

    CIGARs = list() # list of unique CIGARs
    CIGARs_count = list() # corresponding count for each CIGAR
    DELETION_dent = [0]*len(REF_GENOME) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< NEW

    fo_output_name = FILENAME_INPUT + '.mreCRISPR.sam' # [file][output] with the reads that pass the condition
    fo_output = open(fo_output_name, "w")

    # reads a line from the .sam file
    # uses cigar, position, read to build a registered read
    # crops using ZONE
    # creates a custom cigar for that region
    # decides if the read is to be used
    # if it is, adds it to the list of unique cigars

    with open(FILENAME_INPUT, "r") as fo_input:
        for line in fo_input: # reads line by line

            line_split = line.split("\t") # converts the string "line" to an array of string

            if len(line_split) > 10: # avoid the .SAM header

                # grads information from line:
                sam_position = int(line_split[3]) - 1  # python indexes start at 0
                sam_cigar = line_split[5]
                sam_read = line_split[9]

                # computes the registered read
                OUTPUT = mreCRISPR_register_read(sam_position,sam_cigar,sam_read,ref_genome_length)

                # updates the nucleotide deletion profile
                for i in range(0,len(REF_GENOME)): #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< NEW
                    if OUTPUT['read'][i] == '*': #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< NEW
                        DELETION_dent[i] += 1 #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< NEW

                # crops the registered read with the ZONE
                registered_read_zone = OUTPUT['read'][ZONE[0]:ZONE[1]]
                print(registered_read_zone)

                if OUTPUT['read_ok']: # read_ok = False: when I or H in the CIGAR
                    # computes custom cigar for the cropped registered read:
                    registered_cigar_zone = mreCRISPR_cigar_in_zone(ref_zone_seq, registered_read_zone)


                    if mreCRISPR_condition(registered_cigar_zone): # chooses if read should be considered

                        fo_output.write(line)

                        if registered_cigar_zone['cigar_txt'] in CIGARs:
                            CIGARs_count[CIGARs.index(registered_cigar_zone['cigar_txt'])] += 1
                        else:
                            CIGARs.append(registered_cigar_zone['cigar_txt'])
                            CIGARs_count.append(1)

    print('DELETION profile for the amplicon:')#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< NEW
    print(DELETION_dent)#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< NEW

    #R code #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< NEW
    # library('ggplot2')
    # X = c()
    # X = data.frame(x=c(1:length(X)),y=X)
    # g = ggplot(X,aes(x=x,y=y))
    # g = g + geom_point() + geom_line
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< NEW

    fo_input.close()
    fo_output.close()
    return{'CIGARs':CIGARs,'CIGARs_count':CIGARs_count}

    return{}

# this function runs the program: extracts reads and counts the cigars
def mreCRISPR_cigars_and_paired_count_full(FILENAME_INPUT,REF_GENOME,ZONE,SEED):

    # FILENAME_INPUT: sam file
    # REF_GENOME: sequence of the amplicon
    # ZONE: [a,b] boundary of the region of interest

    ref_genome_length = len(REF_GENOME) # size of the ref genome
    ref_zone_seq = REF_GENOME[ZONE[0]:ZONE[1]] # sequence of the ROI
    ref_seed = REF_GENOME[SEED[0]:SEED[1]]

    seed_coord_in_zone = [SEED[0]-ZONE[0],0]
    seed_coord_in_zone[1] = seed_coord_in_zone[0] + (SEED[1] - SEED[0]) - 1

    print(ref_zone_seq)

    CIGARs = list() # list of unique CIGARs
    CIGARs_seed = list()
    CIGARs_count = list() # corresponding count for each CIGAR
    SEQUENCE_seed = list()
    SEQUENCE_zone = list()

    # reads a line from the .sam file
    # uses cigar, position, read to build a registered read
    # crops using ZONE
    # creates a custom cigar for that region
    # decides if the read is to be used
    # if it is, adds it to the list of unique cigars

    with open(FILENAME_INPUT, "r") as fo_input:
        for line in fo_input: # reads line by line

            line_split = line.split("\t") # converts the string "line" to an array of string

            if len(line_split) > 10: # avoid the .SAM header

                # grabs information from line:
                sam_position = int(line_split[3]) - 1  # python indexes start at 0
                sam_cigar = line_split[5]
                sam_read = line_split[9]

                # computes the registered read
                OUTPUT = mreCRISPR_register_read(sam_position,sam_cigar,sam_read,ref_genome_length)

                #print(OUTPUT['read'])

                # crops the registered read with the ZONE
                registered_read_zone = OUTPUT['read'][ZONE[0]:ZONE[1]]
                #print('ZONE:',len(registered_read_zone))
                registered_read_seed = OUTPUT['read'][SEED[0]:SEED[1]]
                #print('ZONE:',SEED,len(registered_read_seed))

                #print(registered_read_zone)

                if OUTPUT['read_ok']: # read_ok = False: when I or H in the CIGAR
                    # computes custom cigar for the cropped registered read:
                    registered_cigar_zone = mreCRISPR_cigar_in_zone(ref_zone_seq, registered_read_zone)
                    registered_cigar_seed = mreCRISPR_cigar_in_zone(ref_seed, registered_read_seed)

                    #print(registered_cigar_zone['cigar_txt'])

                    if registered_cigar_zone['cigar_txt'] in CIGARs:
                        CIGARs_count[CIGARs.index(registered_cigar_zone['cigar_txt'])] += 1
                    else:
                        CIGARs.append(registered_cigar_zone['cigar_txt'])
                        CIGARs_seed.append(registered_cigar_seed['cigar_txt'])
                        SEQUENCE_seed.append(registered_read_seed)
                        SEQUENCE_zone.append(registered_read_zone)
                        CIGARs_count.append(1)

    fo_input.close()

    return{'CIGARs':CIGARs,
           'CIGARs_seed':CIGARs_seed,
           'CIGARs_count':CIGARs_count,
           'SEQ_seed':SEQUENCE_seed,
           'SEQ_zone':SEQUENCE_zone,
           'seed_coord_in_zone':seed_coord_in_zone}

    return{}


# this function pairs cDNA and gDNA and compute mirScore
def mreCRISPR_pair_and_mirScore(DATA_gDNA,DATA_cDNA,sort_type, output_plot_file):
    # computes overlap and sort by either mirScore or SeqDepth
    CIGARs = list()
    CIGARs_seed = list()
    SEQUENCE_seed = list()
    SEQUENCE_zone = list()
    CIGARs_COUNT_gDNA = list()
    CIGARs_COUNT_cDNA = list()
    MIRSCORE = list()

    for i in range(0,len(DATA_gDNA['CIGARs'])):
        if DATA_gDNA['CIGARs'][i] in DATA_cDNA['CIGARs']:
            CIGARs.append(DATA_gDNA['CIGARs'][i])
            CIGARs_seed.append(DATA_gDNA['CIGARs_seed'][i])
            SEQUENCE_seed.append(DATA_gDNA['SEQ_seed'][i])
            SEQUENCE_zone.append(DATA_gDNA['SEQ_zone'][i])
            CIGARs_COUNT_gDNA.append(DATA_gDNA['CIGARs_count'][i])

            idx_cDNA = DATA_cDNA['CIGARs'].index(DATA_gDNA['CIGARs'][i])
            CIGARs_COUNT_cDNA.append(DATA_cDNA['CIGARs_count'][idx_cDNA])

            MIRSCORE.append(float(DATA_cDNA['CIGARs_count'][idx_cDNA])/float(DATA_gDNA['CIGARs_count'][i]))

    print('overlap', len(CIGARs))
    # sort
    #mirScore or SeqDepth

    idx_sort = []

    if sort_type == "SeqDepth":
        print('sorting on depth')
        idx_sort = sorted(range(len(CIGARs_COUNT_cDNA)), key=lambda k: CIGARs_COUNT_cDNA[k])
        idx_sort = list(reversed(idx_sort))

    elif sort_type == "mirScore":
        print('sorting on mirScore')
        idx_sort = sorted(range(len(MIRSCORE)), key=lambda k: MIRSCORE[k])
        idx_sort = list(reversed(idx_sort))

    CIGARs = [CIGARs[i] for i in idx_sort]
    CIGARs_seed = [CIGARs_seed[i] for i in idx_sort]
    SEQUENCE_seed = [SEQUENCE_seed[i] for i in idx_sort]
    SEQUENCE_zone = [SEQUENCE_zone[i] for i in idx_sort]
    CIGARs_COUNT_gDNA = [CIGARs_COUNT_gDNA[i] for i in idx_sort]
    CIGARs_COUNT_cDNA = [CIGARs_COUNT_cDNA[i] for i in idx_sort]
    MIRSCORE = [MIRSCORE[i] for i in idx_sort]

    fo_output = open(output_plot_file, "w")
    for i in range(0,len(CIGARs)):
        line = CIGARs[i] + ","
        line += CIGARs_seed[i] + ","
        line += SEQUENCE_seed[i] + ","
        line += SEQUENCE_zone[i] + ","
        line += str(CIGARs_COUNT_gDNA[i]) + ","
        line += str(CIGARs_COUNT_cDNA[i]) + ","
        line += str(MIRSCORE[i]) + ","
        line += str(DATA_gDNA['seed_coord_in_zone'][0]+1) + ","
        line += str(DATA_gDNA['seed_coord_in_zone'][1]+1) + "\n"
        fo_output.write(line)
    fo_output.close()

# this function pairs cDNA and gDNA and compute mirScore
def mreCRISPR_pair_and_mirScore_new(DATA_gDNA, DATA_cDNA, output_plot_file):
    
	# computes overlap and sort by either mirScore or SeqDepth
	CIGARs = list() # contains the list of all unique cigars over the zone
	CIGARs_seed = list() #
	SEQUENCE_seed = list() #
	SEQUENCE_zone = list() #
	CIGARs_COUNT_gDNA = list() 
	CIGARs_COUNT_cDNA = list()

	# finds the complete set:
	TOTAL_CIGARs = []
	for i in range(0,len(DATA_gDNA['CIGARs'])):
		if DATA_gDNA['CIGARs'][i] in CIGARs:
			print('match')
		else:
			CIGARs.append(DATA_gDNA['CIGARs'][i])
			CIGARs_seed.append(DATA_gDNA['CIGARs_seed'][i])
			SEQUENCE_seed.append(DATA_gDNA['SEQ_seed'][i])
			SEQUENCE_zone.append(DATA_gDNA['SEQ_zone'][i])

	for i in range(0,len(DATA_cDNA['CIGARs'])):
		if DATA_cDNA['CIGARs'][i] in CIGARs:
			print('match')
		else:
			CIGARs.append(DATA_cDNA['CIGARs'][i])
			CIGARs_seed.append(DATA_cDNA['CIGARs_seed'][i])
			SEQUENCE_seed.append(DATA_cDNA['SEQ_seed'][i])
			SEQUENCE_zone.append(DATA_cDNA['SEQ_zone'][i])

	print('stats:',len(CIGARs),len(DATA_gDNA['CIGARs'])+ len(DATA_cDNA['CIGARs']))


	# fills out counts
	for i in range(0,len(CIGARs)):
		if CIGARs[i] in DATA_gDNA['CIGARs']:
			idx_gDNA = DATA_gDNA['CIGARs'].index(CIGARs[i])
			CIGARs_COUNT_gDNA.append(DATA_gDNA['CIGARs_count'][idx_gDNA])
		else:
			CIGARs_COUNT_gDNA.append(0)
		
		if CIGARs[i] in DATA_cDNA['CIGARs']:
			idx_cDNA = DATA_cDNA['CIGARs'].index(CIGARs[i])
			CIGARs_COUNT_cDNA.append(DATA_cDNA['CIGARs_count'][idx_cDNA])
		else:
			CIGARs_COUNT_cDNA.append(0)
		
							
	fo_output = open(output_plot_file, "w")
	line = 'cigar,cigar_seed,seq_seed,seq_zone,count_gDNA,count_cDNA,seed_start,seed_end\n'
	fo_output.write(line)

	for i in range(0,len(CIGARs)):
		line = CIGARs[i] + ","
		line += CIGARs_seed[i] + ","
		line += SEQUENCE_seed[i] + ","
		line += SEQUENCE_zone[i] + ","
		line += str(CIGARs_COUNT_gDNA[i]) + ","
		line += str(CIGARs_COUNT_cDNA[i]) + ","
		line += str(DATA_gDNA['seed_coord_in_zone'][0]+1) + ","
		line += str(DATA_gDNA['seed_coord_in_zone'][1]+1) + "\n"
		fo_output.write(line)
	fo_output.close()
