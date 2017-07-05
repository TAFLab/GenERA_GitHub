# contains all the R functions for mreCRISPR pipeline

read_cigars_letter <- function(cigar){
  cigar_char  = strsplit(cigar, split = c())
  cigar_char = cigar_char[[1]]
  is_letter = cigar_char %in% c("+","-","D","#")
  idx_letter = c(1:length(cigar_char))
  idx_letter = idx_letter[is_letter]
  
  cursor = 1
  LETTERS = cigar_char[idx_letter]
  NUM = c()
  for(i in idx_letter){
    num = as.integer(paste(cigar_char[c(cursor:(i-1))],collapse = ''))
    cursor = i + 1
    NUM = c(NUM,num) 
  }
  #print(LETTERS)
  #print(NUM)
  
  return(LETTERS)
}
read_cigars_number <- function(cigar){
  cigar_char  = strsplit(cigar, split = c())
  cigar_char = cigar_char[[1]]
  is_letter = cigar_char %in% c("+","-","D","#")
  idx_letter = c(1:length(cigar_char))
  idx_letter = idx_letter[is_letter]
  
  cursor = 1
  LETTERS = cigar_char[idx_letter]
  NUM = c()
  for(i in idx_letter){
    num = as.integer(paste(cigar_char[c(cursor:(i-1))],collapse = ''))
    cursor = i + 1
    NUM = c(NUM,num) 
  }
  #print(LETTERS)
  #print(NUM)
  
  return(NUM)
}
crop <- function(DATA,type,seq){
  if(nrow(DATA)==0){
    return(DATA)
  }else{
    L = strsplit(as.character(DATA[,seq]),split=c())
    CROP = rep(TRUE, nrow(DATA))
    for(i in c(1:length(L))){
      CROP[i] = type %in% L[[i]]
    }
    DATA_crop = DATA[!CROP,]
    return(DATA_crop)
  }
}
keep <- function(DATA,type,seq){
  if(nrow(DATA)==0){
    return(DATA)
  }else{
    L = strsplit(as.character(DATA[,seq]),split=c())
    KEEP = rep(FALSE, nrow(DATA))
    for(i in c(1:length(L))){
      KEEP[i] = type %in% L[[i]]
    }
    DATA_keep = DATA[KEEP,]
    return(DATA_keep)
  }
}
compute_min <- function(DATA, columns){
  count_min = rep(0,nrow(DATA))
  for(i in c(1:nrow(DATA))){
    count_min[i] = min(DATA[i,columns])
    #print(count_min[i])
  }
  return(count_min)
}
compute_max <- function(DATA, columns){
  count_max = rep(0,nrow(DATA))
  for(i in c(1:nrow(DATA))){
    count_max[i] = max(DATA[i,columns])
    #print(count_min[i])
  }
  return(count_max)
}
partition_data <- function(DATA, case){
  
  # cases:
  # 1: deletion in seed vs. deletion out
  print('-----------------------------------------------------------')
  print(paste(c(as.character(nrow(DATA)),' SIZ in DATA'),collapse=''))
  
  DATA_out = c()
  
  #CASE_1
  if (case == 1){
    print('>> partition: deletion in seed vs. deletion out')
    # count number of nt in seed:
    seed_length = sum(read_cigars_number(as.character(DATA[1,'cigar_seed'])))
    print(paste(c(as.character(seed_length),' nt in seed'),collapse=''))
    
    # extract CIZ with seed free
    DATA_part_1 = DATA[DATA$cigar_seed == paste(c(as.character(seed_length),'+'),collapse=''),]
    print(paste(c(as.character(nrow(DATA_part_1)),' SIZ seed intact'),collapse=''))
    # keep only the deletion ones
    DATA_part_1 = keep(DATA_part_1,'*','seq_zone')
    print(paste(c(as.character(nrow(DATA_part_1)),' SIZ seed intact and deletion only'),collapse=''))
    
    # extract CIZ with seed deleted
    DATA_part_2 = keep(DATA,'*','seq_seed')
    print(paste(c(as.character(nrow(DATA_part_2)),' SIZ seed deleted'),collapse=''))
    
    if(nrow(DATA_part_1)>1 && nrow(DATA_part_2)>1){
      DATA_part_1$partition = 'SIZ_deletion_seed_intact'
      DATA_part_2$partition = 'SIZ_deletion_seed_touched'
      
      DATA_out = rbind(DATA_part_1,DATA_part_2)
    }
    
  }
  return(DATA_out)
  
  
}
compute_mutation_percent <- function(DATA){
  
  mutation_percent = rep(0,nrow(DATA)) # contains percent of mutation on CIZ
  zone_length = sum(read_cigars_number(as.character(DATA[1,'cigar_zone']))) # required to compute percent
  print('**ZONE LENGTH**')
  print(zone_length)
  
  for (i in c(1:nrow(DATA))){
    
    current_cigar = as.character(DATA[i,'cigar_zone'])
    current_cigar_letter = read_cigars_letter(current_cigar)
    current_cigar_number = read_cigars_number(current_cigar)
    
    for (j in c(1:length(current_cigar_letter))){
      if (current_cigar_letter[j] == '-'){ # contains a mutation
        mutation_percent[i] = mutation_percent[i] + current_cigar_number[j] # add number of nt mutated
      }
    }
    mutation_percent[i] = mutation_percent[i]/zone_length
  }
  return(mutation_percent)
}
convert_to_XDX_seq <- function(sequences){
  
  if(identical(sequences, character(0))){
    return = c()
  }else{
    sequences_XDX = rep('',length(sequences))
    for(i in c(1:length(sequences))){
      seq_split = strsplit(sequences[i],c())
      seq_split = seq_split[[1]]
      seq_split_XDX = rep('',length(seq_split))
      for(j in c(1:length(seq_split))){
        if(seq_split[j] == '*'){
          seq_split_XDX[j] = '*'
        }else{
          seq_split_XDX[j] = 'X'
        }
      }
      sequences_XDX[i] = paste(seq_split_XDX,collapse='')
    }
    return(sequences_XDX)
  }
}
length_string <- function(my_string){
  my_string_split = strsplit(my_string,c())
  my_string_split = my_string_split[[1]]
  return(length(my_string_split))
}
string_to_char_array <- function(my_string){
  my_string_split = strsplit(my_string,c())
  my_string_split = my_string_split[[1]]
  return(my_string_split)
}
is_mutated <- function(vector_of_seq){
  
  is_m = rep(FALSE,length(vector_of_seq))
  
  for(i in c(1:length(vector_of_seq))){
    current_seq = string_to_char_array(vector_of_seq[i])
    is_m[i] = c('*') %in% current_seq
  }
  
  return(is_m)
}
compute_mutation_percent_from_cigar <- function(cigar){
  
  current_cigar_letter = read_cigars_letter(cigar)
  current_cigar_number = read_cigars_number(cigar)
  
  mutation_percent = 0
  
  for (j in c(1:length(current_cigar_letter))){
    if (current_cigar_letter[j] == '-'){ # contains a mutation
      mutation_percent = mutation_percent + current_cigar_number[j] # add number of nt mutated      
    }
  }
  mutation_percent = mutation_percent / sum(current_cigar_number)
  
  return(mutation_percent)
  
}

is_mre_only <- function(DATA, MRE){
  
  N = nrow(DATA)
  mre_only = rep(FALSE,N)
  
  for(i in c(1:N)){
    
    CIZ = string_to_char_array(as.character(DATA$seq_zone_XDX[i]))
    X = c(1:length(CIZ))
    
    X_MRE = X[X %in% MRE]
    X_MRE_out = X[!(X %in% MRE)]
  
    mre_only[i] = (sum(CIZ[X_MRE] == '*') >= 1) && (sum(CIZ[X_MRE_out] == '*') == 0)
  }
  return(mre_only)
}


display_deletion <- function(DATA,XDX,sort_score,alpha, type){
  
  DATA_plot = DATA[order(DATA[,c(sort_score)]),]
  
  print(DATA_plot[,c(sort_score)])
  
  plot_x = c()
  plot_y = c()
  plot_y_bis = c()
  plot_alpha = c()
  plot_type = c()
  
  for(i in c(1:nrow(DATA_plot))){
    seq = string_to_char_array(as.character(DATA_plot[i,c(XDX)]))
    seq_idx = c(1:length(seq))
    del_idx = seq_idx[seq == '*']
    
    plot_x = c(plot_x, del_idx)
    plot_y = c(plot_y, rep(i, length(del_idx))) #DATA[i,c(sort_score)]
    plot_y_bis = c(plot_y_bis, rep(DATA_plot[i,c(sort_score)], length(del_idx))) #DATA[i,c(sort_score)]
    plot_alpha = c(plot_alpha, rep(DATA_plot[i,c(alpha)], length(del_idx)))
    plot_type = c(plot_type, rep(DATA_plot[i,c(type)], length(del_idx)))
  }
  
  DATA_plot_2 = data.frame(x = plot_x, y = plot_y, norm_mirScore = plot_y_bis, a = plot_alpha, type = plot_type)
  g = ggplot(DATA_plot_2, aes(x = x, y = y, alpha = log10(norm_mirScore), color = factor(type))) + geom_point(aes(size = log10(a)))
  #g = ggplot(DATA_plot_2, aes(x = x, y = y, color = log10(norm_mirScore))) + geom_point(aes(size = log10(norm_mirScore), shape = factor(type)))
  g = g + annotate("rect", xmin = DATA_plot[1,'seed_start']-0.5, 
                   xmax = DATA_plot[1,'seed_end']+0.5, 
                   ymin = min(plot_y)-1, ymax = max(plot_y)+1,alpha = .2)
  return(g)
}

match_names <- function(DATA, name_of_col){
  
  MATCH = read.csv('../../META/OTHERS/name_match.csv')
  return(as.character(MATCH$Gene.Name[match(DATA[,name_of_col],MATCH$files)]))

}
