# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||| GenERA PIPELINE

# NAME: GenERA_5.R

# FUNCTION:
# - load reformated file
# - STATS: total_gDNA_count
# - STATS: total_cDNA_count
# > remove all cigars in zone with ‘-‘ (searched in sequence in the zone)
# - STATS: total_gDNA_count
# - STATS: total_cDNA_count
# 
# - STATS: seed_start_in_zone
# - STATS: seed_end_in_zone
# - STATS: zone_length
# 
# - STATS: perfect_wt_gDNA_count (=0 if no wild-type cigar)
# - STATS: perfect_wt_cDNA_count (=0 if no wild-type cigar)
# 
# > remove cigars with too many mutations
# - compute mutation frequency of all cigars
# - keep only mutation frequency < 0.1 cigars
# - STATS: post_mutated_valid_gDNA_count
# - STATS: post_mutated_valid_cDNA_count
# - STATS: nb_UDP_pre_merge
# 
# > merge similar cigars after ignoring mutations
# - convert seq_zone and seq_seed to X***X profiles
# - pair reads based on seq_zone_XDX
# - compute DATA_merge (seq_zone_XDX, 
#                       seq_seed_XDX, 
#                       count_gDNA,
#                       count_cDNA,
#                       seed_start,
#                       seed_end,
#                       zone_length)
# - remove all unpaired UDPs
# - STATS: nb_UDP_post_merge
# - STATS: wt_gDNA_count
# - STATS: wt_cDNA_count
# 
# > UDP stats
# - STATS: nb_deletion_UDP_post_merge
# - STATS: deletion_UDP_post_merge_gDNA_count
# - STATS: deletion_UDP_post_merge_cDNA_count
# 
# - STATS: nb_deletion_out_seed_UDP_post_merge (outside seed)
# - STATS: deletion_out_seed_UDP_post_merge_gDNA_count
# - STATS: deletion_out_seed_UDP_post_merge_cDNA_count
# 
# - STATS: nb_deletion_seed_UDP_post_merge (touch seed)
# - STATS: deletion_seed_UDP_post_merge_gDNA_count
# - STATS: deletion_seed_UDP_post_merge_cDNA_count
# 
# - STATS: nb_full_deletion_seed_UDP_post_merge (complete seed deletion)
# - STATS: full_deletion_seed_UDP_post_merge_gDNA_count
# - STATS: full_deletion_seed_UDP_post_merge_cDNA_count
# 
# - STATS: nb_deletion_mre_only_UDP_post_merge ***
# - STATS: deletion_mre_only_UDP_post_merge_gDNA_count
# - STATS: deletion_mre_only_UDP_post_merge_cDNA_count

# INPUT:
# - output from GenERA_3.py (reformatted with GenERA_4.R)

# OUTPUT:
# - stats.csv file

# DEPENDENCIES: 
# - GenERA_functions.R

# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

setwd("~/Desktop/GenERA_GitHub/CODE") # sets path to working directory; /!\ changes based on system architecture

rm(list = ls()) # clear all variables
source('GenERA_functions.R') # load functions

# fetch files 
path_to_files = "../DATA/CSV/A_REFORMATED_CSV_FILES/"
LIST_FILES = list.files(path = path_to_files, pattern = '*csv')

GLOBAL_STATS = c() # output vector
STATS_FILES = c()

for (f in c(1:length(LIST_FILES))){
  
  file = LIST_FILES[f]
  print(file)
  DATA = read.csv(paste(path_to_files,file, sep=''),header = FALSE)
  colnames(DATA) <- c('cigar_zone','cigar_seed','seq_seed','seq_zone','count_gDNA','count_cDNA','score','seed_start','seed_end')
  
  if(nrow(DATA)>=2){
  
  
  DATA_origin = DATA
  L_zone = length_string(as.character(DATA$seq_zone))
  L_seed = length_string(as.character(DATA$seq_seed))
  
  STATS = c()
  STATS_NAME = c()
  #-------------------------------------------------------------- get total read counts
  STATS = c(STATS, sum(DATA$count_gDNA))
  STATS_NAME = c(STATS_NAME,'total_gDNA_count')
  STATS = c(STATS, sum(DATA$count_cDNA))
  STATS_NAME = c(STATS_NAME,'total_cDNA_count')
  #-------------------------------------------------------------- get total post removal of incomplete reads
  DATA = crop(DATA,'-','seq_zone')
  STATS = c(STATS, sum(DATA$count_gDNA))
  STATS_NAME = c(STATS_NAME,'valid_gDNA_count')
  STATS = c(STATS, sum(DATA$count_cDNA))
  STATS_NAME = c(STATS_NAME,'valid_cDNA_count')
  #-------------------------------------------------------------- get stats before fusion with mutated
  info_zone_length = length_string(as.character(DATA[1,'seq_zone']))
  info_is_there_perfect_wt = nrow(DATA[DATA$cigar_zone == paste(as.character(info_zone_length),'+',sep=''),]) > 0
  
  STATS = c(STATS, DATA$seed_start[1], DATA$seed_end[1], info_zone_length)
  STATS_NAME = c(STATS_NAME,'seed_start_in_zone','seed_end_in_zone','zone_length')
  
  if(info_is_there_perfect_wt){
    data = DATA[DATA$cigar_zone == paste(as.character(info_zone_length),'+',sep=''),]
    STATS = c(STATS, data$count_gDNA)
    STATS_NAME = c(STATS_NAME,'perfect_wt_gDNA_count')
    STATS = c(STATS, data$count_cDNA)
    STATS_NAME = c(STATS_NAME,'perfect_wt_cDNA_count')
  }else{
    STATS = c(STATS, 0)
    STATS_NAME = c(STATS_NAME,'perfect_wt_gDNA_count')
    STATS = c(STATS, 0)
    STATS_NAME = c(STATS_NAME,'perfect_wt_cDNA_count')
  }
  #-------------------------------------------------------------- merge data and remove mutations
  # compute the mutation frequency
  mutation_frequency = rep(0,nrow(DATA))
  if(nrow(DATA)>0){
    for(i in c(1:nrow(DATA))){
      mutation_frequency[i] = compute_mutation_percent_from_cigar(as.character(DATA[i,'cigar_zone']))
    }
  }
  DATA$mutation_freq = mutation_frequency
  
  # how many too mutated reads are removed
  DATA = DATA[DATA$mutation_freq < 0.1,]
  STATS = c(STATS, sum(DATA$count_gDNA))
  STATS_NAME = c(STATS_NAME,'post_mutated_valid_gDNA_count')
  STATS = c(STATS, sum(DATA$count_cDNA))
  STATS_NAME = c(STATS_NAME,'post_mutated_valid_cDNA_count')
  STATS = c(STATS, nrow(DATA))
  STATS_NAME = c(STATS_NAME,'nb_UDP_pre_merge')
  
  # convert seq_zone and seq_seed to X***X profiles
  DATA$seq_zone_XDX = convert_to_XDX_seq(as.character(DATA$seq_zone))
  DATA$seq_seed_XDX = convert_to_XDX_seq(as.character(DATA$seq_seed))
  
  # pair reads based on seq_zone_XDX
  XDX_unique = unique(DATA$seq_zone_XDX)
  XDX_unique_length = length(XDX_unique)
  
  merge_seq_zone_XDX = rep('',XDX_unique_length)
  merge_seq_seed_XDX = rep('',XDX_unique_length)
  merge_count_gDNA = rep(0,XDX_unique_length)
  merge_count_cDNA = rep(0,XDX_unique_length)
  
  for(i in c(1:XDX_unique_length)){
    
    DATA_merge = DATA[DATA$seq_zone_XDX == XDX_unique[i],]
    
    merge_seq_zone_XDX[i] = XDX_unique[i]
    merge_seq_seed_XDX[i] = DATA_merge[1,'seq_seed_XDX']
    merge_count_gDNA[i] = sum(DATA_merge[,'count_gDNA'])
    merge_count_cDNA[i] = sum(DATA_merge[,'count_cDNA'])
  }
  
  if(identical(merge_seq_zone_XDX, character(0))){
    DATA_merge = data.frame()
  }else{
    DATA_merge = data.frame(seq_zone_XDX = merge_seq_zone_XDX, 
                            seq_seed_XDX = merge_seq_seed_XDX, 
                            count_gDNA = merge_count_gDNA,
                            count_cDNA = merge_count_cDNA,
                            seed_start = DATA[1,'seed_start'],
                            seed_end = DATA[1,'seed_end'],
                            zone_length = length_string(as.character(DATA[1,'seq_zone'])))
    
    # remove all unmatched
    DATA_merge = DATA_merge[(DATA_merge$count_gDNA > 0) * (DATA_merge$count_cDNA >0) > 0,]
  }
  
  STATS = c(STATS, nrow(DATA_merge))
  STATS_NAME = c(STATS_NAME,'nb_UDP_post_merge')
  
  
  data_merge_wt = DATA_merge[DATA_merge$seq_zone_XDX == paste(rep('X', info_zone_length),collapse = ''),]
  
  if(nrow(data_merge_wt)>0){
    STATS = c(STATS, data_merge_wt$count_gDNA)
    STATS_NAME = c(STATS_NAME,'wt_gDNA_count')
    STATS = c(STATS, data_merge_wt$count_cDNA)
    STATS_NAME = c(STATS_NAME,'wt_cDNA_count')
  }else{
    STATS = c(STATS, 0)
    STATS_NAME = c(STATS_NAME,'wt_gDNA_count')
    STATS = c(STATS, 0)
    STATS_NAME = c(STATS_NAME,'wt_cDNA_count')
  }
  
  #-------------------------------------------------------------- give stats on UDPs
  data_merge_deletion = keep(DATA_merge,'*','seq_zone_XDX') # deletion only
  STATS = c(STATS, nrow(data_merge_deletion))
  STATS_NAME = c(STATS_NAME,'nb_deletion_UDP_post_merge')
  STATS = c(STATS, sum(data_merge_deletion$count_gDNA))
  STATS_NAME = c(STATS_NAME,'deletion_UDP_post_merge_gDNA_count')
  STATS = c(STATS, sum(data_merge_deletion$count_cDNA))
  STATS_NAME = c(STATS_NAME,'deletion_UDP_post_merge_cDNA_count')
  
  data_merge_deletion = crop(data_merge_deletion,'*','seq_seed_XDX') # deletion out only
  STATS = c(STATS, nrow(data_merge_deletion))
  STATS_NAME = c(STATS_NAME,'nb_deletion_out_seed_UDP_post_merge')
  STATS = c(STATS, sum(data_merge_deletion$count_gDNA))
  STATS_NAME = c(STATS_NAME,'deletion_out_seed_UDP_post_merge_gDNA_count')
  STATS = c(STATS, sum(data_merge_deletion$count_cDNA))
  STATS_NAME = c(STATS_NAME,'deletion_out_seed_UDP_post_merge_cDNA_count')
  
  data_merge_deletion = keep(DATA_merge,'*','seq_seed_XDX') # deletion in only
  STATS = c(STATS, nrow(data_merge_deletion))
  STATS_NAME = c(STATS_NAME,'nb_deletion_seed_UDP_post_merge')
  STATS = c(STATS, sum(data_merge_deletion$count_gDNA))
  STATS_NAME = c(STATS_NAME,'deletion_seed_UDP_post_merge_gDNA_count')
  STATS = c(STATS, sum(data_merge_deletion$count_cDNA))
  STATS_NAME = c(STATS_NAME,'deletion_seed_UDP_post_merge_cDNA_count')
  
  # \/
  data_merge_deletion_full_seed = data_merge_deletion[data_merge_deletion$seq_seed_XDX == paste(rep('*',L_seed),collapse=''),]
  STATS = c(STATS, nrow(data_merge_deletion_full_seed))
  STATS_NAME = c(STATS_NAME,'nb_full_deletion_seed_UDP_post_merge')
  STATS = c(STATS, sum(data_merge_deletion_full_seed$count_gDNA))
  STATS_NAME = c(STATS_NAME,'full_deletion_seed_UDP_post_merge_gDNA_count')
  STATS = c(STATS, sum(data_merge_deletion_full_seed$count_cDNA))
  STATS_NAME = c(STATS_NAME,'full_deletion_seed_UDP_post_merge_cDNA_count')
  # /\
  
  #-------------------------------------------------------------- retrain data to MRE only
  #\/
  info_mre_start_m1 = DATA_merge$seed_start[1] - 14 - 1 
  info_mre_end = DATA_merge$seed_end[1] 
  #/\
  
  if((info_mre_start_m1>0) && (nrow(data_merge_deletion)>0)){
    
    is_this_type_of_UDP = rep(FALSE,nrow(data_merge_deletion))
    
    for(i in c(1:nrow(data_merge_deletion))){ # note the use of data_merge_deletion
      current_UDP = as.character(data_merge_deletion[i,'seq_zone_XDX'])
      current_UDP = string_to_char_array(current_UDP)
      
      is_this_type_of_UDP[i] = (current_UDP[info_mre_start_m1] == 'X') && (current_UDP[info_mre_end] == 'X')
    }
    
    data_merge_deletion$is_this_type_of_UDP = is_this_type_of_UDP
    data_merge_deletion = data_merge_deletion[data_merge_deletion$is_this_type_of_UDP,]
    
    STATS = c(STATS, nrow(data_merge_deletion))
    STATS_NAME = c(STATS_NAME,'nb_deletion_mre_only_UDP_post_merge')
    STATS = c(STATS, sum(data_merge_deletion$count_gDNA))
    STATS_NAME = c(STATS_NAME,'deletion_mre_only_UDP_post_merge_gDNA_count')
    STATS = c(STATS, sum(data_merge_deletion$count_cDNA))
    STATS_NAME = c(STATS_NAME,'deletion_mre_only_UDP_post_merge_cDNA_count')
    
  }else{
    
    STATS = c(STATS, 0)
    STATS_NAME = c(STATS_NAME,'nb_deletion_mre_only_UDP_post_merge')
    STATS = c(STATS, 0)
    STATS_NAME = c(STATS_NAME,'deletion_mre_only_UDP_post_merge_gDNA_count')
    STATS = c(STATS, 0)
    STATS_NAME = c(STATS_NAME,'deletion_mre_only_UDP_post_merge_cDNA_count')
    
  }
  
  GLOBAL_STATS = rbind(GLOBAL_STATS,STATS)
  STATS_FILES = c(STATS_FILES, file)
  
  }else{
    print(paste('>>> REMOVE: ', file, sep = ''))
  }
  
}

GLOBAL_STATS = as.data.frame(GLOBAL_STATS)
colnames(GLOBAL_STATS) <- STATS_NAME
rownames(GLOBAL_STATS) <- STATS_FILES

write.csv(GLOBAL_STATS, '../META/global_stats.csv')



