# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||| GenERA PIPELINE

# NAME: GenERA_4.R

# FUNCTION:
# - reformat output from GenERA_3.py (add one extra column) to be compatible with downstream R scripts 

# INPUT:
# - output files from GenERA_3.py (raw CSV files)

# OUTPUT:
# - reformated output from GenERA_3.py (reformated CSV files)

# DEPENDENCIES: none

# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

setwd("~/Desktop/GenERA_GitHub/CODE") # sets path to working directory; /!\ changes based on system architecture 

rm(list = ls()); # clears variables
library('ggplot2') # loads graphic library

path_to_files = "../DATA/CSV/A_RAW_CSV_FILES/" 
LIST_FILES = list.files(path = path_to_files, pattern = '*csv')


path_to_reformated_files = "../DATA/CSV/A_REFORMATED_CSV_FILES/" 

# goes through all files and modify format
for (f in c(1:length(LIST_FILES))){
  
  INFO_NB_ROW = c()
  INFO_gDNA_count = c()
  INFO_cDNA_count = c()
  INFO_names = c()
  
  # loads file
  file = LIST_FILES[f]
  print(file)
  DATA = read.csv(paste(path_to_files,file, sep=''),header = TRUE)
  
  # renames columns 
  colnames(DATA) <- c('cigar_zone','cigar_seed','seq_seed','seq_zone','count_gDNA','count_cDNA','seed_start','seed_end')
  DATA_2 = data.frame(a = DATA[,1], b = DATA[,2], c = DATA[,3], d = DATA[,4], e = DATA[,5], f = DATA[,6], g = rep(0,nrow(DATA)), h = DATA[,7], i = DATA[,8])
  
  # writes reformated data 
  write.table(DATA_2, paste(path_to_reformated_files,file, sep=''), col.names = FALSE, row.names = FALSE, sep = ',')
  
}




