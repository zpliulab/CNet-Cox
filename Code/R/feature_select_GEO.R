## 2023.5.10 LLY made, 2026.1.14 LLY used, adjust ‘Data/Independent_data’
## This code is used to select expression of six genes in PRS 
## (with time and srate information) from external GEO datasets 


##############################################################################
##############################################################################
## clear
rm(list = ls())

## set pathway
# path <- '/home/lly/R/'
path <- '/Users/lilingyu/E/PhD/R/'
setwd(paste(path,'CNetCox', sep=''))

library(dplyr)       
library(tidyr)
library(tidyverse)   

myfile = list.files("Data/Independent_data")                
dir = paste("./Data/Independent_data/", myfile, sep="")    
n = length(dir)                                  

Data1 = read.table(file = dir[1], header=T, check.names = FALSE) 
dim(Data1)     # 13703   159

coef_gene <- read.csv("Data/TCGA_NEW/UniMutVariate_markergene.csv", header = T, sep=',')
gene <- as.matrix(coef_gene[,1])
dim(gene)      # 6   1

colnames(gene) <- c('gene')
Data2 <- cbind(rownames(Data1), Data1)
colnames(Data2) <- c('gene', colnames(Data1))

genedata <- merge(gene, Data2, by = "gene")
genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
genedata2 <- rbind(Data1[c(1,2),],genedata1)   

name <- dir[1]
name
## Press ./Split, Extract
name1 <- str_split_fixed(name, "./", 2);
name2 <- str_split_fixed(name1[2], "/", 2)
name3 <- str_split_fixed(name2[2], "[.]", 2);
name4 <- substr(name3[1], 1, 9)
name5 <- substr(name3[1], 23, 27)
name6 <- str_c(name4, name5)

path <- paste("./Feature_data/Data_GEO/",paste(name6,".txt", sep=""), sep="")
# write.table(genedata2, path, quote = F)


myfile = list.files("Data/Independent_data")              
dir = paste("./Data/Independent_data/", myfile, sep="")   
n = length(dir)                                 


for (i in 2:n){

  # i <- 14
  Data1 = read.table(file = dir[i],header=T, check.names = FALSE) 
  dim(Data1)
  
  coef_gene <- read.csv("Data/TCGA_NEW/UniMutVariate_markergene.csv", header = T, sep=',')
  gene <- as.matrix(coef_gene[,1])

  colnames(gene) <- c('gene')
  Data2 <- cbind(rownames(Data1), Data1)
  colnames(Data2) <- c('gene', colnames(Data1))
  
  genedata <- merge(gene, Data2, by = "gene")
  genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
  genedata2 <- rbind(Data1[c(1,2),],genedata1)  
  
  name <- dir[i]
 
  name1 <- str_split_fixed(name, "./", 2);
  name2 <- str_split_fixed(name1[2], "/", 2)
  name3 <- str_split_fixed(name2[2], "[.]", 2);
  name4 <- substr(name3[1], 1, 9)
  name5 <- substr(name3[1], 23, 27)
  name6 <- str_c(name4, name5)
  
  
  path <- paste("./Feature_data/Data_GEO/",paste(name6,".txt", sep=""), sep="")
  # write.table(genedata2, path, quote = F)

  
  print("***i***")
  print(i)

}                                               






