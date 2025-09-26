## 2021.9.26 Input the feature genes from each method to obtain the pathway gene-network pairs
rm(list=ls())

library(dplyr)        
library(tidyr)
library(tidyverse)    
library(igraph)


setwd('D:\\E\\博士\\R_程序\\SVM\\Data\\RTCGA')
net <- as.matrix(read.csv("allgene_comp_net.csv",header = T))
net[1,]



setwd('D:\\E\\博士\\R_程序\\SVM\\Data\\RTCGA\\result\\featurenew')

myfile = list.files("CsvdataFeatureSeleOnce")                 
dir = paste("./CsvdataFeatureSeleOnce/", myfile, sep="")      
n = length(dir)                                   


mynet <- function(i){
  # i <- 8
  node <- as.matrix(read.csv(file = dir[i],header = T))
  node[1]
  
  node_used <- node
  net_used <- net
  k1 <- which(net_used[,1] %in% node_used)   # 562
  k2 <- which(net_used[,2] %in% node_used)   # 929
  length(intersect(k1,k2))    # 4
  used <- net_used[intersect(k1,k2),]
  if (length(used) == 0){break;}
  
  PP <- graph_from_data_frame(used,directed = F)
  p1 <- simplify(PP)   
  ed <- as_edgelist(p1, names = TRUE)
  # union(ed[,1],ed[,2])
  # setdiff(as.character(node), union(ed[,1],ed[,2]))
  
  name <- dir[i]
  name1 <- str_split_fixed(name, "./", 2)
  name2 <- str_split_fixed(name1[2], "/", 2)
  name3 <- str_split_fixed(name2[2], ".", 2)
  name4 <- str_split_fixed(name3[2], " .", 2)
  name5 <- name4[1]
  name6 <- str_c("Feature_net_", name5)
  path <- paste("./Feature_net_sele_One/",paste(name6,".csv"))
  write.csv(ed, path, row.names = F, quote = F)
} 

mynet(8)
dir[i]


## -----------------------------------------------------------

for (i in 1:n) {
  # i <- 2
  node <- as.matrix(read.csv(file = dir[i],header = T))
  node[1]
  
  node_used <- node
  # net_used <- net[,c(2,4)]
  net_used <- net
  k1 <- which(net_used[,1] %in% node_used)   # 562
  k2 <- which(net_used[,2] %in% node_used)   # 929
  length(intersect(k1,k2))    # 4
  used <- net_used[intersect(k1,k2),]
  if (length(used) == 0){
    break;
  }
  
  PP <- graph_from_data_frame(used,directed = F)
  p1 <- simplify(PP)   
  ed <- as_edgelist(p1, names = TRUE)
  
  name <- dir[i]
  name1 <- str_split_fixed(name, "._", 2);
  name2 <- str_sub(name1[2], end=3)
  name3 <- str_c("Feature_net_", name2)
  path <- paste("./featurenet/",paste(name3,".csv"))
  write.csv(ed, path, row.names = F, quote = F)
  
}