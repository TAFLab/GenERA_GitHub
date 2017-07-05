# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||| GenERA PIPELINE

# NAME: GenERA_7.R

# FUNCTION:
# - check all UDP footprints for overlap with predicted regulatory elements.
# - check all UDP footprints for de-novo deletion/creation of MREs for top 10 S2 miRNAs (including miR-184) 

# INPUT:
# - merge UDP files from GenERA_6.R foud in ./DATA/CSV/A_FINAL_CSV_FILES/
# - ./META/GENOMES.csv contains the reference genome sequence for each amplicon
# - ./META/predicted_RREs.csv contains the name and coordinates of predicted regulatory elements for each amplicon 
# - ./META/TOP10_miRNA.csv contains the name, seed sequence, and expression levels of the top 10 most expressed miRNA in S2 cells


# OUTPUT:
# - enhanced UDP files (stored in './DATA/CSV/A_FINAL_enhanced_CSV_FILES/')

# DEPENDENCIES: mreCRISPR_functions.R

# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

setwd("~/Desktop/GenERA_GitHub/CODE") # sets path to working directory; /!\ changes based on system architecture

rm(list = ls()) # clear all variables
source('GenERA_functions.R') # load functions

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD METADATA
path_meta_data = "../META/"
# metadata -- genomes
GENOMES = read.csv(paste(path_meta_data,'GENOMES.csv', sep = ''), header = TRUE)

# metadata -- features
FEATURES = read.csv(paste(path_meta_data,"predicted_RREs.csv",sep = '')) # features to be mapped
FEATURES$name = paste(as.character(FEATURES$name),'_GenERA.csv', sep = '') # change name
FEATURES$name = paste('data_',as.character(FEATURES$name), sep = '')
FEATURES = FEATURES[FEATURES$type != 'AUB',] # filter out particular features

# metadata -- miRNA
MIRNA = read.csv(paste(path_meta_data, 'TOP10_miRNA.csv', sep = ''))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD UDP FILES
path_mreCRISPR_data = "../DATA/CSV/A_FINAL_CSV_FILES/"
path_mreCRISPR_data_enhanced = '../DATA/CSV/A_FINAL_enhanced_CSV_FILES/'
LIST_FILES = list.files(path = path_mreCRISPR_data, pattern = '*.csv')


COUNT.UDP = 0
# contains the creation or deletion of motif, if is chinese plot UDP and which amplicon
DNM.data.frame = data.frame() 

for(file in LIST_FILES){
  
  current_DATA = read.csv(paste(c(path_mreCRISPR_data,file),collapse = ''),header = TRUE)
  COUNT.UDP = COUNT.UDP + nrow(current_DATA)
  current_GENOME = GENOMES[GENOMES$filename == file, ]
  current_FEATURES = FEATURES[FEATURES$name == file, ]
  
  # look for feature deletions ---------------------------------------------
  amplicon_genome = as.character(current_GENOME$genome)
  amplicon_genome_array = string_to_char_array(amplicon_genome)
  amplicon_genome_length = length_string(amplicon_genome)
  
  amplicon_zone_start = current_GENOME$zone_start
  amplicon_zone_end = current_GENOME$zone_end
  
  amplicon_zone_footprint = rep(0,amplicon_genome_length)
  amplicon_zone_footprint[c(amplicon_zone_start:amplicon_zone_end)] = 1
  
  FEATURES_DELETED = rep('none',nrow(current_DATA))
  
  if(nrow(current_FEATURES)>0){
    for(k in c(1:nrow(current_FEATURES))){
      
      feature_footprint = rep(0,amplicon_genome_length)
      feature_footprint[c(current_FEATURES$feature_start[k]:current_FEATURES$feature_end[k])] = 1
      
      if(sum(amplicon_zone_footprint * feature_footprint) > 0){ # the feature is within the zone
        
        for(i in c(1:nrow(current_DATA))){
          
          # get the udp
          udp = string_to_char_array(as.character(current_DATA$seq_zone_XDX[i]))
          udp_bool = udp == '*'
          udp_binary = rep(0,current_DATA$zone_length[1])
          udp_binary[udp_bool] = 1
          
          # map to genome
          udp_footprint = rep(0,amplicon_genome_length)
          udp_footprint[c(c(amplicon_zone_start:amplicon_zone_end))] = udp_binary
          
          # check if interferes with feature
          if (sum(udp_footprint * feature_footprint) > 0){ # the udp overlap with the feature
            
            if(FEATURES_DELETED[i] == 'none'){
              FEATURES_DELETED[i] = as.character(current_FEATURES$type[k])
            }else{
              FEATURES_DELETED[i] = paste(FEATURES_DELETED[i], as.character(current_FEATURES$type[k]), sep = '|')
            }
          }
        }
      }
    }
  }
  
  current_DATA$FeatureDeletion = FEATURES_DELETED
  
  # miRNA prior ----------------------------------------------------------------
  amplicon_seed_start = amplicon_zone_start + current_DATA$seed_start[1] - 1
  amplicon_seed_end = amplicon_zone_start + current_DATA$seed_end[1] - 1
  
  amplicon_seed_footprint = rep(0,amplicon_genome_length)
  amplicon_seed_footprint[c(amplicon_seed_start:amplicon_seed_end)] = 1
  
  # map all MRE ----------------------------------------------------------------
  info.name = c()
  info.feature_start = c()
  info.feature_end = c()
  
  for(i in c(1:nrow(MIRNA))){
    x = gregexpr(pattern = as.character(MIRNA$seed_seq[i]), amplicon_genome)
    if(x[[1]][1] != -1){
      
      for(j in c(1:length(x[[1]]))){
        info.name = c(info.name, as.character(MIRNA$name[i]))
        info.feature_start = c(info.feature_start, x[[1]][j])
        info.feature_end = c(info.feature_end, (x[[1]][j] + attr(x[[1]], 'match.length')[j] - 1))
      }
    }
  }
  
  current_MIRNA_PRIOR = data.frame(name = info.name, feature_start = info.feature_start, feature_end = info.feature_end)
  
  current_DATA$has_seed = sum(as.character(current_MIRNA_PRIOR$name) == 'miR-184') 
  current_DATA$has_seed_foreign = sum(as.character(current_MIRNA_PRIOR$name) != 'miR-184')
  
  # search for differencial MRE motifs ***
  
  DNM = rep('none',nrow(current_DATA))
  DNM_exclu = rep('none',nrow(current_DATA))

  for(i in c(1:nrow(current_DATA))){
    
    # get the udp
    udp = string_to_char_array(as.character(current_DATA$seq_zone_XDX[i]))
    udp_binary = rep(0,current_DATA$zone_length[1])
    udp_binary[udp == '*'] = 1
    
    # map to genome
    udp_footprint = rep(0,amplicon_genome_length)
    udp_footprint[c(amplicon_zone_start:amplicon_zone_end)] = udp_binary
    udp_seq_post_deletion = paste(amplicon_genome_array[udp_footprint == 0], collapse = '')
    
    info.before = rep(0, nrow(MIRNA))
    info.after = rep(0, nrow(MIRNA))

    for(k in c(1:nrow(MIRNA))){
      
      x_before = gregexpr(pattern = as.character(MIRNA$seed_seq[k]), amplicon_genome)
      x_after = gregexpr(pattern = as.character(MIRNA$seed_seq[k]), udp_seq_post_deletion)
      
      if(x_before[[1]][1] != -1){
        info.before[k] = length(x_before[[1]])
      }
      if(x_after[[1]][1] != -1){
        info.after[k] = length(x_after[[1]])
      }
    }
    
    current_MIRNA = MIRNA
    current_MIRNA$count_before = info.before
    current_MIRNA$count_after = info.after
    current_MIRNA$count_discrepency = current_MIRNA$count_after - current_MIRNA$count_before
    
    current_MIRNA_interest = current_MIRNA[current_MIRNA$count_discrepency != 0,]
    
    if(nrow(current_MIRNA_interest)>0){
      
      for(k in c(1:nrow(current_MIRNA_interest))){
        to_write = paste(c(as.character(current_MIRNA_interest$name[k]),
                           as.character(current_MIRNA_interest$count_discrepency[k])), collapse = ',')
        
        if(DNM[i] == 'none'){
          DNM[i] = to_write
        }else{
          DNM[i] = paste(DNM[i], to_write, sep = '|')
        }
         
        
        if(to_write != "miR-184,-1"){
          DNM_exclu[i] = paste(DNM_exclu[i], to_write, sep = '|')
          print(paste(c(to_write, as.character(current_DATA$count_min[i]), as.character(current_DATA$deletion_seed[i])), collapse = '___'))
        
          
          dnm = data.frame(file = file, 
                           udp.count.min = current_DATA$count_min[i], 
                           udp.seed.del = current_DATA$deletion_seed[i],
                           udp.normalised.score = current_DATA$score_norm_divide[i],
                           udp.full.seq = udp_seq_post_deletion,
                           event = to_write,
                           miR184.prior = current_MIRNA$count_before[current_MIRNA$name == 'miR-184'])
          DNM.data.frame = rbind(DNM.data.frame, dnm)  
        }
        
      }
    }
  }
  
  current_DATA$DNM = DNM
  current_DATA$DNM_exclu = DNM_exclu
  current_DATA$seed_seq = paste(amplicon_genome_array[amplicon_seed_footprint == 1], collapse = '')
  write.csv(current_DATA, file = paste(path_mreCRISPR_data_enhanced, file, sep = ''))
}


# # characterise those cases of miR184 motif creation
# load("~/Desktop/CORR/INTERCEPT/DATA_intercept.RData")
# I = INTERCEPT[INTERCEPT$has.miRscore,]
# 
# path_data = "./DATA_20160114/"
# path_data_new = './DATA_20160120_chinese_plot/'
# LIST_FILES_to_move = list.files(path = path_data, pattern = '*csv')
# 
# 
# for(file in LIST_FILES_to_move){
#   if(file %in% I$name.miRscore.file){
#     DATA = read.csv(paste(c(path_data,file),collapse = ''),header = TRUE)
#     write.csv(DATA, file = paste(path_data_new, file, sep = ''))
#   }
# }

# I = INTERCEPT[INTERCEPT$has.miRscore,]
# DNM.data.frame$in.chinese.plot = as.character(DNM.data.frame$file) %in% as.character(I$name.miRscore.file)
# DNM.data.frame = DNM.data.frame[DNM.data.frame$event == 'miR-184,1',]
# DNM.data.frame = DNM.data.frame[DNM.data.frame$in.chinese.plot,]
# DNM.data.frame = DNM.data.frame[DNM.data.frame$udp.count.min > 10,]