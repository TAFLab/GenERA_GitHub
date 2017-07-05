# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||| GenERA PIPELINE

# NAME: GenERA_6.R

# FUNCTION:
# - remove the incomplete reads (reads with dilution on 5' or 3' of the ROI)
# - remove all cigars with too many mutations
# - compute DATA_merge = data.frame(seq_zone_XDX, 
#                                   seq_seed_XDX, 
#                                   count_gDNA,
#                                   count_cDNA,
#                                   seed_start,
#                                   seed_end,
#                                   zone_length)
# - remove unmatched
# - add count_min, deletion (boolean), deletion_seed (boolean) columns
# - compute UDP score
# - compute normalised UDP score (UNS)
# - boolean check if the deletion is contained or not within the MRE
# - check for de-novo MRE: (now done in mreCRISPR_7.R)
#   - map back the UDP on the reference genome
# - get genome sequence post deletion
# - search for and flag creation or removal of miRNA seed from top10 miRNA in S2+ cells

# INPUT:
# - output from GenERA_3.py (created with GenERA_4.R) in ./DATA/CSV/A_REFORMATED_CSV_FILES/

# OUTPUT:
# - csv of DATA_merge in ./DATA/CSV/A_FINAL_CSV_FILES/

# DEPENDENCIES: GenERA_functions.R

# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

setwd("~/Desktop/GenERA_GitHub/CODE") # sets path to working directory; /!\ changes based on system architecture

rm(list = ls()) # clear all variables
source('GenERA_functions.R') # load functions
GENOMES = read.csv('../META/GENOMES.csv',header = TRUE)

# list files in directory of interest
path_to_files = "../DATA/CSV/A_REFORMATED_CSV_FILES/"
LIST_FILES = list.files(path = path_to_files, pattern = '*csv')

INFO_NB_ROW = c() # contains stats pre and post merging
info_wt = rep(FALSE,length(LIST_FILES))
LOG = c() # store log informations

path_output_files = '../DATA/CSV/A_FINAL_CSV_FILES/'

for (f in c(1:length(LIST_FILES))){
  
  # load file and extract data 
  file = LIST_FILES[f]
  print(file)
  DATA = read.csv(paste(path_to_files,file, sep=''),header = FALSE)
  colnames(DATA) <- c('cigar_zone','cigar_seed','seq_seed','seq_zone','count_gDNA','count_cDNA','score','seed_start','seed_end')
  
  info_nb_row = rep(0,4)
  
  if(nrow(DATA) >= 2){
  
    info_nb_row[1] = nrow(DATA) # number of unique cigars before cleaning and merging
    
    #----------------------------------- remove the incomplete reads
    
    if(nrow(DATA)>0){
      nb_row_before = nrow(DATA)
      DATA = crop(DATA,'-','seq_zone') #<<< remove the incomplete reads
      nb_row_after = nrow(DATA)
      print(paste(c('removing \'-\'', nb_row_before, nb_row_after), collapse = ' '))
    }else{
      print('nothing before removing \'-\'')
    } 
    
    info_nb_row[2] = nrow(DATA) # number after cleaning
    
    #----------------------------------- remove hyper mutated reads
    
    if(nrow(DATA)>0){
      
      DATA$mutation_percent = compute_mutation_percent(DATA)
      
      nb_row_before = nrow(DATA)
      # remove all cigars with too many mutations
      DATA = DATA[DATA$mutation_percent < 0.1,] #<<< remove too mutated
      nb_row_after = nrow(DATA)
      print(paste(c('removing >10% mutated', nb_row_before, nb_row_after), collapse = ' '))
      
    }else{print('cannot compute mutation frequency')}
    
    #----------------------------------- generate UDP and compile read count for each UDP
    
    info_nb_row[3] = nrow(DATA) # number after removal of too mutated
    
    if(nrow(DATA)>0){
      
      # convert seq_zone and seq_seed to X***X signatures
      DATA$seq_zone_XDX = convert_to_XDX_seq(as.character(DATA$seq_zone))
      DATA$seq_seed_XDX = convert_to_XDX_seq(as.character(DATA$seq_seed))
      
      # further pair reads based on seq_zone_XDX
      XDX_unique = unique(DATA$seq_zone_XDX)
      XDX_unique_length = length(XDX_unique)
      
      merge_seq_zone_XDX = rep('',XDX_unique_length)
      merge_seq_seed_XDX = rep('',XDX_unique_length)
      merge_ref = rep('',XDX_unique_length)
      merge_count_gDNA = rep(0,XDX_unique_length)
      merge_count_cDNA = rep(0,XDX_unique_length)
      
      for(i in c(1:XDX_unique_length)){
        
        DATA_merge_temp = DATA[DATA$seq_zone_XDX == XDX_unique[i],]
        
        merge_seq_zone_XDX[i] = XDX_unique[i]
        merge_seq_seed_XDX[i] = DATA_merge_temp[1,'seq_seed_XDX']
        merge_count_gDNA[i] = sum(DATA_merge_temp[,'count_gDNA'])
        merge_count_cDNA[i] = sum(DATA_merge_temp[,'count_cDNA'])   
        
      }
      
      DATA_merge = data.frame(seq_zone_XDX = merge_seq_zone_XDX, 
                              seq_seed_XDX = merge_seq_seed_XDX, 
                              count_gDNA = merge_count_gDNA,
                              count_cDNA = merge_count_cDNA,
                              seed_start = DATA[1,'seed_start'],
                              seed_end = DATA[1,'seed_end'],
                              zone_length = length_string(as.character(DATA[1,'seq_zone'])))
      
      DATA_merge = DATA_merge[(DATA_merge$count_gDNA>0)*(DATA_merge$count_cDNA>0)>0,] #<<< remove the unmatched
      
      # add further information
      DATA_merge$count_min = compute_min(DATA_merge,c('count_gDNA','count_cDNA'))
      DATA_merge$deletion = is_mutated(as.character(DATA_merge$seq_zone_XDX))
      DATA_merge$deletion_seed = is_mutated(as.character(DATA_merge$seq_seed_XDX))
      
      nb_row_before = nrow(DATA)
      nb_row_after = nrow(DATA_merge)
      print(paste(c('merging UDPs', nb_row_before, nb_row_after), collapse = ' '))
      
      info_nb_row[4] = nrow(DATA_merge)
      
      # add scores
      DATA_merge$score = DATA_merge$count_cDNA / DATA_merge$count_gDNA
      wt_mirScore = DATA_merge[!DATA_merge$deletion,'score']
      
      if(length(wt_mirScore) > 0){
        
        DATA_merge$score_norm_divide = DATA_merge$score / wt_mirScore
        DATA_merge$score_norm_minus = DATA_merge$score - wt_mirScore
        
        # boolean check if the deletion is contained or not within the MRE
        MRE_start = GENOMES$mir_start[match(file,as.character(GENOMES$filename))] - GENOMES$zone_start[match(file,as.character(GENOMES$filename))]+1
        MRE_start = max(MRE_start,1)
        
        MRE_end = GENOMES$mir_end[match(file,as.character(GENOMES$filename))] - GENOMES$zone_start[match(file,as.character(GENOMES$filename))]+1
        MRE_end = min(MRE_end,length_string(as.character(DATA_merge$seq_zone_XDX[1])))
        
        DATA_merge$MRE_start = MRE_start
        DATA_merge$MRE_end = MRE_end
        DATA_merge$mre_only = is_mre_only(DATA_merge,c(MRE_start:MRE_end))

        # cycle through the deletion types and search for de-novo motif -- is now done in a later script
        DATA_merge$DNM = rep('missing',nrow(DATA_merge))
        
        # write the merge file in folder
        write.csv(DATA_merge,paste(path_output_files,file,sep=''))

        info_wt[f] = TRUE
      }else{
        print('No wt!')
      }
      
      
    }else{print('cannot generate XDX')}  
    
    INFO_NB_ROW = rbind(INFO_NB_ROW,info_nb_row)
  
  }
  else{
    print(paste('>> remove: ', file, sep = ''))
    INFO_NB_ROW = rbind(INFO_NB_ROW,info_nb_row)
  }
  
  
}

INFO_NB_ROW = as.data.frame(INFO_NB_ROW)
LOG = data.frame(file = LIST_FILES, info_wt)
LOG = cbind(LOG, INFO_NB_ROW)
colnames(LOG) <- c('file','contains_WT','nb_origin','post_incomplete', 'post_mutant', 'post_merge')

write.csv(LOG, '../META/merging_log.csv')