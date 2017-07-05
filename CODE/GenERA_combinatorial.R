# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||| GenERA PIPELINE

# NAME: GenERA_Combinatorial.R

# FUNCTION:
# - remove low depth UDPs
# - amend data with information regarding the length of the deletion and wether or not it is continuous
# - Establish list of CEP and filter out CEP which have less than N representative UDPs available (here we chose N = 10)
# - Provide code to compare the UNS distribution between UDP mapping to zone 1 (user defined), zone 2 (user defined) or both. With possibility for quantile filtering
# - Test the additive property of combinatorial editing by finding all CEP trios where CEP3 (parent) is a combination of CEP1 and CEP2 (children). Extracts corresponding UNS for all 3 CEPs.

# INPUT:
# - '../DATA/CSV/B_combinatorial/data_combinatorial_GenERA.csv', csv data file as produced by GenERA_6.R

# OUTPUT:
# - 'zone_boxplot.csv' contains the result of the two zone analysis
# - './export/cepComb_---x---_1110000_0000111_.csv' series of .csv files containing the result of the parent children combinatorial analysis

# DEPENDENCIES: 
# - mreCRISPR_functions.R

# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

setwd("~/Desktop/GenERA_GitHub/CODE") # sets path to working directory; /!\ changes based on system architecture

rm(list = ls()) # clear all variables
library("stringr") # load string handling library
library("ggplot2") # load graphics library
library("Hmisc") # load Hmisc library

source('./GenERA_functions.R') 
DATA = read.csv('../DATA/CSV/B_combinatorial/data_combinatorial_GenERA.csv', header = T)
DATA = DATA[DATA$count_min > 10,] # remove low depth UDPs
#/////////////////////////////////////////////////////////////// MRE coordinates
# provide coordinates of the regulatory elements
MREs = list(92:97,
             100:105,
             128:133,
             241:246,
             259:264,
             263:268,
             279:284)

#/////////////////////////////////////////////////////////////// add length of deletion and continuity to DATA

deletion.length = rep(0,nrow(DATA)) # will contain the length of the deletion
deletion.continuous = rep(NA, nrow(DATA)) # will contain a boolean assessing if the deletion is continuous
CEP = rep(NA, nrow(DATA)) # will contain the combinatorial editing pattern of each UDP

ROI.length = length_string(as.character(DATA$seq_zone_XDX[1]))
ROI = 1:ROI.length

for(i in 1:nrow(DATA)){
  
  current_udp = str_split(as.character(DATA$seq_zone_XDX[i]), pattern = "")[[1]]
  deletion.length[i] = sum(current_udp == '*')
  deletion.region = ROI[current_udp == '*']
  if(deletion.length[i] >= 1){
    if('X' %in% current_udp[min(deletion.region):max(deletion.region)]){
      deletion.continuous[i] = F
    }else{
      deletion.continuous[i] = T
    }
  }
  #/////////////////////////////
  cep = rep('x', length(MREs))
  for(j in 1:length(MREs)){
    if('*' %in% current_udp[MREs[[j]]]){
      cep[j] = '-'
    }
  }
  CEP[i] = paste(cep, collapse = '')
}

DATA$deletion.length = deletion.length
DATA$deletion.continuous = deletion.continuous
DATA$cep = CEP

#/////////////////////////////////////////////////////////////// Filter to CEP with minimum 
# CEP which have less than 10 UDPs are removed from downstream analysis

CEP.table = table(as.character(DATA$cep)) 
CEP.list = names(CEP.table[CEP.table >= 10])
DATA.short = DATA[DATA$cep %in% CEP.list,] 

# g = ggplot(DATA.short, aes(x = cep, y = log10(score_norm_divide))) + geom_jitter(aes(colour = cep, shape = deletion.continuous))
# g = g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# print(g)

#/////////////////////////////////////////////////////////////// zone analysis 1, 2 and 1+2
# Compare the UNS of UDPs in zone 1, 2 or 1+2 to check for combinatorial effects 

ZONES = list(48:180, 181:436) # set the zones
deletion.zone = rep("none",nrow(DATA)) # will contain z1 (contained in zone 1), z2 (contained in zone 2), or z1&2 (contained in zone 1&2)

for(i in 1:nrow(DATA)){
  
  current_udp = str_split(as.character(DATA$seq_zone_XDX[i]), pattern = "")[[1]]
  
  if(('*' %in% current_udp[ZONES[[1]]]) & !('*' %in% current_udp[ZONES[[2]]])){
    deletion.zone[i] = "z1"
  }
  if(('*' %in% current_udp[ZONES[[2]]]) & !('*' %in% current_udp[ZONES[[1]]])){
    deletion.zone[i] = "z2"
  }
  if(('*' %in% current_udp[ZONES[[1]]]) & ('*' %in% current_udp[ZONES[[2]]])){
    deletion.zone[i] = "z1&2"
  }

}

DATA$deletion.zone = deletion.zone

# g = ggplot(DATA, aes(x = deletion.zone, y = log10(score_norm_divide)))  + geom_jitter(aes(colour = deletion.zone)) + geom_boxplot(alpha = 0.5)
# print(g)

# add quantile filtering: remove extreme UNSs in each category
zone.names = unique(as.character(DATA$deletion.zone))

DATA.quantile = data.frame()

for(current.name in zone.names){
  data = DATA[DATA$deletion.zone == current.name,]
  P = quantile(data$score_norm_divide, c(0.25, 0.75))
  data = data[(data$score_norm_divide >= P[1]) & (data$score_norm_divide <= P[2]),]
  DATA.quantile = rbind(DATA.quantile, data)
}
write.csv(DATA.quantile[,c('deletion.zone', 'score_norm_divide')], file = '../DATA/CSV/B_combinatorial/zone_boxplot.csv')

#/////////////////////////////////////////////////////////////// plot CEP and corresponding UNS values


# PLOT.x = c()
# PLOT.y = c()
# 
# PLOT.x.mean = c()
# PLOT.y.mean = c()
# 
# g = ggplot()
# g = g + geom_point(aes(x = c(0,0,10,10), y = c(0,10,0,10)))
# 
# for(i in 1:length(CEP.list)){
# 
#   cep = str_split(CEP.list[[i]], pattern = '')[[1]]
#   print(cep)
# 
#   for(j in 1:length(cep)){
# 
#     if(cep[j] == '-'){
#       g = g + annotate("rect", xmin = j, xmax = j+0.5, ymin = i, ymax = (i-0.5),
#                        alpha = .2, fill = 'red')
# 
#     }else{
#       g = g + annotate("rect", xmin = j, xmax = j+0.5, ymin = i, ymax = (i-0.5),
#                        alpha = .2, fill = 'gray')
# 
#     }
# 
#   }
# 
# 
#   DATA.current = DATA[DATA$cep == CEP.list[[i]],]
# 
#   PLOT.x = c(PLOT.x, length(cep) + 1 - min(log10(DATA$score_norm_divide)) + log10(DATA.current$score_norm_divide))
#   PLOT.y = c(PLOT.y, rep(i, nrow(DATA.current)) + runif(nrow(DATA.current), min = -0.5, max = 0))
# 
#   PLOT.x.mean = c(PLOT.x.mean, length(cep) + 1 - min(log10(DATA$score_norm_divide)) + mean(log10(DATA.current$score_norm_divide)))
#   PLOT.y.mean = c(PLOT.y.mean, i-0.25)
# }
# 
# x.mean = length(cep) + 1 - min(log10(DATA$score_norm_divide)) + log10(1)
# g = g + annotate("segment", x = x.mean, xend = x.mean, y = 0, yend = length(CEP.list),
#                 colour = "blue")
# g = g + geom_point(aes(x = PLOT.x, y = PLOT.y), alpha = 0.5)
# g = g + geom_point(aes(x = PLOT.x.mean, y = PLOT.y.mean), alpha = 0.5, colour = 'red', size = 2)
# print(g)

#/////////////////////////////////////////////////////////////// 
#/////////////////////////////////////////////////////////////// test additive properties

fn_convert_01 <- function(input){
  output = rep("", length(input))
  for(i in 1:length(input)){
    o = str_split(input[i], pattern = "")[[1]]
    o01 = rep(0,length(o))
    o01[o == '-'] = 1
    output[i] = paste(as.character(o01), collapse = '')
  }
  return(output)
}


Qx.export = data.frame()
DATA$cep01 = fn_convert_01(as.character(DATA$cep)) # convert to binary number

for(cep in CEP.list){

  # check if the number of deleted elements is greater than 2
  nb.deleted.elements = sum(string_to_char_array(cep) == '-')

  if(nb.deleted.elements >= 2){

    # go over each pairs and see if the mix pattern is equal to the target pattern
    for(i in 1:(length(CEP.list)-1)){
      for(j in (i+1):length(CEP.list)){
        
        cep1 = as.numeric(string_to_char_array(fn_convert_01(CEP.list[i])))
        cep2 = as.numeric(string_to_char_array(fn_convert_01(CEP.list[j])))
        cep12 = paste(as.character(cep1+cep2), collapse = '')
        
        if((cep12 == fn_convert_01(cep)) & (sum(cep1)*sum(cep2) != 0)){
          
          print('match--------------------------------')
          print(cep)
          print(cep1)
          print(cep2)


          #//////////////////////////////////////////////////////////////////

          data.cep = DATA[DATA$cep == cep,]
          data.cep1 = DATA[DATA$cep == CEP.list[i],]
          data.cep2 = DATA[DATA$cep == CEP.list[j],]

          P = quantile(data.cep$score_norm_divide, c(0.25,0.75))
          data.cep = data.cep[(data.cep$score_norm_divide >= P[1]) & (data.cep$score_norm_divide <= P[2]),]

          P = quantile(data.cep1$score_norm_divide, c(0.25,0.75))
          data.cep1 = data.cep1[(data.cep1$score_norm_divide >= P[1]) & (data.cep1$score_norm_divide <= P[2]),]

          P = quantile(data.cep2$score_norm_divide, c(0.25,0.75))
          data.cep2 = data.cep2[(data.cep2$score_norm_divide >= P[1]) & (data.cep2$score_norm_divide <= P[2]),]

          score.sum = c()
          for(k in 1:nrow(data.cep1)){
            for(l in 1:nrow(data.cep2)){
              

              x = data.cep1$score_norm_divide[k] * data.cep2$score_norm_divide[l]
              score.sum = c(score.sum, x)
              
            }
          }
          
          data_plot = data.frame(class = c(rep('obs-parent', nrow(data.cep)), 
                                           rep('pred', length(score.sum)),
                                           rep('obs1-daughter', nrow(data.cep1)),
                                           rep('obs2-daughter', nrow(data.cep2))), 
                                 v = c(data.cep$score_norm_divide, 
                                       score.sum,
                                       data.cep1$score_norm_divide,
                                       data.cep2$score_norm_divide))
          
          # find the case A, B or C
          if((nrow(data.cep)>=10) & (nrow(data.cep1)>=10) & (nrow(data.cep2)>=10)){
            
            # note the that 10 UDPs requirement on both target and children CEPs are post quantile filtering
            data_plot$id = paste(c(cep, paste(cep1, collapse = ''), paste(cep2, collapse = '')), collapse = '|')
            
            write.csv(data_plot, file = paste(c('../DATA/CSV/B_combinatorial/cepComb', cep, paste(cep1, collapse = ''), paste(cep2, collapse = ''),'.csv'), collapse = '_'))
            Qx.export = rbind(Qx.export, data_plot)
            
            # g = ggplot(data_plot, aes(x = v)) + geom_density(aes(colour = class))
            # g = g + ggtitle(paste(c(cep, paste(cep1, collapse = ''), paste(cep2, collapse = ''), nrow(data.cep), nrow(data.cep1), nrow(data.cep2)), collapse = '|'))
            # g = g + scale_x_log10()
            # print(g)
          }
          
          
          
          
        }
      }
    }

  }


}

# list.cases = read.csv('list_cases.csv')
# Qx.export$case =  list.cases$class[match(Qx.export$id, list.cases$case)]
# write.csv(Qx.export, file = 'UNS_pred_obs_combin_multiply.csv')

