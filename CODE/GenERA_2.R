# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||| GenERA PIPELINE

# NAME: GenERA_2.R

# FUNCTION:
# - reads nucleotide deletion profile files and determines a Region Of Interest (ROI) to work with for each amplicon (deletion peak calling)

# INPUT:
# - "./META/nucleotide_deletion_pattern.csv": contains the deletion counts for gDNA and cDNA at single nucleotide resolution across the referecne genome
# (see output of mreCRISPR_1.py)

# OUTPUT:
# - "./META/ROI.csv": contains for each amplicon the amplicon ROI coordinates containing both gDNA and cDNA peaks.
# One row per amplicon, columns are as follow:
# - name: name of the amplicon
# - start: left coordinate of the zone (in the genome reference)
# - end: left coordinate of the zone (in the genome reference)
# - divergence: measure of overlap between gDNA and cDNA deletion peak

# DEPENDENCIES: none

# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

setwd("~/Desktop/GenERA_GitHub/CODE") # sets path to working directory; /!\ changes based on system architecture 

# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

rm(list = ls()); # cleans all variables
library('ggplot2') # loads graphics library; /!\ might require an install.packages('ggplot2')

# ||||||||||||||||||||||||||| functions

# smoothing functions
vector_smooth <- function(X){
  X_smooth = X
  for(i in c(2:length(X)-1)){
    X_smooth[i] = mean(X[c(i-1,i,i+1)])
  }
  return(X_smooth)
}
vector_smooth_n <- function(X,n){
  X_smooth = X
  for(i in c(1:n)){
    X_smooth= vector_smooth(X_smooth)
  }
  return(X_smooth)
}

# computes the derivative of a vector as discrete function
vector_derivative <- function(X){
  X_der = X
  for(i in c(2:length(X)-1)){
    der_left = X[i]-X[i-1]
    der_right = X[i+1]-X[i]
    X_der[i] = mean(c(der_left,der_right))
  }
  return(X_der)
}

# uses previous functions to find the deletion peak
find_peak <- function(data,smooth_factor, cutoff_factor, create_plot){
  
  output = c(0,0,0)
  
  x = as.integer(data[,2])
  y.smooth = vector_smooth_n(as.numeric(data[,3]),smooth_factor)
  y.max = max(y.smooth)
  y.base = 0.01*y.max 
  y.der = vector_derivative(y.smooth)
    
  BLOBs_count = c()
  BLOBs_min = c()
  BLOBs_max = c()
  blob.x = c()
  blob.y = c()
  cutoff = cutoff_factor*max(abs(y.der))
    
  for(i in c(2:(length(x)-1))){
    
    value = max(c(y.der[i-1]),abs(y.der[i]),abs(y.der[i+1]))
    if ((value >= cutoff) | (y.smooth[i] > y.base)){      
      blob.x = c(blob.x,x[i])
    }else{
      if(length(blob.x) > 0){
        BLOBs_count = c(BLOBs_count,length(blob.x))
        BLOBs_min = c(BLOBs_min,min(blob.x))
        BLOBs_max = c(BLOBs_max,max(blob.x))
        
        blob.x = c()
        blob.y = c()
      }
    }
  }

  if(length(BLOBs_count) > 0){
    
    idx = order(BLOBs_count, decreasing = TRUE)
    BLOBs_count = BLOBs_count[idx]
    BLOBs_min = BLOBs_min[idx]
    BLOBs_max = BLOBs_max[idx]
    
    output = c(BLOBs_min[1],BLOBs_max[1],BLOBs_count[1]) 
    
    if(create_plot){
      g = ggplot() 
      g = g + geom_point(aes(x = data$Position,y= y.der),color='blue')
      g = g + geom_line(aes(x = data$Position,y= y.der),color='blue')
      g = g + annotate("rect",xmin = BLOBs_min[1], xmax = BLOBs_max[1], ymin = -cutoff, ymax = +cutoff, alpha = 0.1, fill="blue")
      g = g + ggtitle(paste (col_names[k], length(BLOBs_count), sep = " ", collapse = NULL))
      print(g)
      
      g = ggplot() 
      g = g + geom_point(aes(x = data$Position,y= data[,3]),color='blue')
      g = g + geom_line(aes(x = data$Position,y= y.smooth),color='red')
      g = g + annotate("rect",xmin = BLOBs_min[1], xmax = BLOBs_max[1], ymin = min(data[,3]), ymax = max(data[,3]), alpha = 0.1, fill="blue")
      g = g + ggtitle(paste (col_names[k], length(BLOBs_count), sep = " ", collapse = NULL))
      print(g)
    }
  }
  
  return(output)
}


# /////////////////////////////////////////////////////////////// MAIN

DATA = read.csv("../META/nucleotide_deletion_pattern.csv",header = TRUE)

col_names = colnames(DATA)
print(col_names)

OUTPUT_names = c() # name of the amplicon
OUTPUT_x_start = c() # zone start
OUTPUT_x_end = c() # zone end
OUTPUT_seed_start = c() # seed start
OUTPUT_seed_end = c() # seed end
OUTPUT_overlap = c() # gDNA, cDNA peak overlap

for(k in c(3:ncol(DATA))){
  
  print(k) 
  
  # grabs data from a single amplicon
  data = DATA[,c(1,2,k)]
  data$Position = as.integer(data$Position)
  data[,3] = as.numeric(data[,3])
  
  data.gDNA = data[data$SampleType == 'gDNA',]
  data.cDNA = data[data$SampleType == 'cDNA',]
  
  factor_smooth = 10 # number of smoothing loop (arbitrary)
  factor_threshold = 0.005 # number of smoothing loop (arbitrary)
  
  output.gDNA = find_peak(data.gDNA,factor_smooth, factor_threshold, FALSE)
      
      # plots the gDNA nucleotide deletion profile with zone coordinates
      g = ggplot() 
      g = g + geom_point(aes(x = data.gDNA$Position,y= data.gDNA[,3]),color='blue',alpha = 0.5)
      #g = g + geom_line(aes(x = data.gDNA$Position,y= data.gDNA[,3]))
      g = g + geom_line(aes(x = data.gDNA$Position,y= vector_smooth_n(data.gDNA[,3],10)),color='blue')
      g = g + annotate("rect",xmin = output.gDNA[1], xmax = output.gDNA[2], ymin = min(data.gDNA[,3]), ymax = max(data.gDNA[,3]), alpha = 0.1, fill="blue")
  
  output.cDNA = find_peak(data.cDNA,factor_smooth, factor_threshold, FALSE)
      
      # plots the cDNA nucleotide deletion profile with zone coordinates
      g = g + geom_point(aes(x = data.cDNA$Position,y= data.cDNA[,3]),color='red',alpha = 0.5)
      # g = g + geom_line(aes(x = data.cDNA$Position,y= data.cDNA[,3]))
      g = g + geom_line(aes(x = data.cDNA$Position,y= vector_smooth_n(data.cDNA[,3],10)),color='red')
      g = g + annotate("rect",xmin = output.cDNA[1], xmax = output.cDNA[2], ymin = min(data.cDNA[,3]), ymax = max(data.cDNA[,3]), alpha = 0.1, fill="red")
      g = g + ggtitle(col_names[k])
      
      #print(g) # for display in IDE
      ggsave(file=paste("../FIGURES/ROI/",col_names[k],".pdf",sep = ""))
  
  OUTPUT_names = c(OUTPUT_names, col_names[k])
  print(col_names[k]) # shows current amplicon
  
  start = min(output.gDNA[1],output.cDNA[1])
  end = max(output.gDNA[2],output.cDNA[2])
  
  # displays zone informations
  print('gDNA zone')
  print(c(output.gDNA[1],output.gDNA[2]))
  print('cDNA zone')
  print(c(output.cDNA[1],output.cDNA[2]))
  
  OUTPUT_x_start = c(OUTPUT_x_start, start)
  OUTPUT_x_end = c(OUTPUT_x_end, end)
  
  diverge = length(setdiff(c(output.gDNA[1]:output.gDNA[2]),c(output.cDNA[1]:output.cDNA[2])))/max(output.gDNA[3],output.cDNA[3])
  OUTPUT_overlap = c(OUTPUT_overlap,diverge)
}

OUTPUT = data.frame(name = OUTPUT_names, start = OUTPUT_x_start, end = OUTPUT_x_end, divergence = OUTPUT_overlap)

write.csv(OUTPUT, file = "../META/ROI.csv", quote = FALSE)

