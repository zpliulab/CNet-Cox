
## clear
rm(list = ls())

## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- '/home/lly/R/'

## load function
source(paste(path, 'CNetCox/R/myoverlap_separate2GroupsCox.R', sep=''))

setwd(paste(path,'CNetCox/Feature_data', sep=''))
library(tidyverse)

myfile = list.files("Data_GEO")               
dir = paste("./Data_GEO/", myfile, sep="")   
n = length(dir) 

data = read.table(file = dir[1],header=T, check.names = FALSE) 
datat <- as.data.frame(t(data))
Data <- as.matrix(data)

x <- model.matrix(status ~., datat)[,-c(1,2)]   
x_hat <- data.frame(x)
y <- as.matrix(datat[,c(1,2)])  
colnames(y) <- c("time", "status")
y_hat <- data.frame(y)

## coef
setwd(paste(path,'CNetCox/Data/TCGA_NEW/', sep=''))
coef_gene <- read.csv("UniMutVariate_markergene.csv", header = T, sep=',')
coef <- as.matrix(coef_gene[,2])
rownames(coef) <- coef_gene[,1]


# my_overlap_coef ------------------------------------------------------
my_overlap <- function(x, y){

  coefs.v <- x[,1] %>% { .[. != 0]}
  coefs.v %>% {
    data.frame(gene.name   = names(.),
               coefficient = .,
               stringsAsFactors = FALSE)
  } %>%
    arrange(gene.name) %>%
    knitr::kable()
  
  sele <- rownames(as.matrix(coefs.v))
  gene <- rownames(y)[-c(1,2)]
  overlap <- intersect(sele, gene)
  
  lab <- x[,1] %>% { .[. != 0]} %>% as.matrix
  coefs.v <- lab[overlap,]
  
  my <- list(coefs.v, overlap)
  return(my)
} 



# plot ----------------------------------------------------------------------

library(ggpubr)
library(magrittr)
library(survminer)
# library(glmSparseNet)
coef_test <- my_overlap(coef, Data)
plotp_Train <- separate2GroupsCox(as.vector(coef_test[[1]]), x_hat[, coef_test[[2]]], 
                                 plot.title = 'GSE1456_smfs', as.data.frame(y), 
                                 legend.outside = T)
plot_train <- plotp_Train$plot


## for Xtile
p_index <- cbind(y,plotp_Train$index)
colnames(p_index) <- c(colnames(y), "riskscore")


name <- dir[1]
name1 <- str_split_fixed(name, "./", 2);
name2 <- str_split_fixed(name1[2], "/", 2)
name3 <- str_split_fixed(name2[2], "[.]", 2);
name4 <- name3[1]
name5 <- str_c(name4,"_ex")


setwd(paste(path,'CNetCox/Feature_data', sep=''))
path <- paste("./Xtile/",paste(name4,".txt", sep=""), sep="")
write.table(p_index, path, quote = F, row.names = F, sep="\t")
path <- paste("./Figure/", paste(name5,".pdf", sep=""), sep="")
pdf(path, width = 4.5, height = 4.5, onefile = FALSE) 
plot_train
dev.off() 



###########################################################################

myfile_GEO = list.files("Data_GEO")               
dir_GEO = paste("./Data_GEO/", myfile_GEO, sep="")   
n = length(dir_GEO) 


for (i in 1:n){
  # i <- 18
  
  data = read.table(file = dir_GEO[i],header=T, check.names = FALSE)
  datat <- as.data.frame(t(data))
  Data <- as.matrix(data)
  dim(Data)

  x <- as.matrix(datat[,-c(1,2)])   
  colnames(x) <- rownames(Data)[-c(1,2)]
  x_hat <- data.frame(x)
  y <- as.matrix(datat[,c(1,2)])  
  colnames(y) <- c("time", "status")
  y_hat <- data.frame(y)
  
  ## set pathway
  path <- '/Users/lilingyu/E/PhD/R/'
  # path <- '/home/lly/R/'
  setwd(paste(path,'CNetCox/Data/TCGA_NEW/', sep=''))
  
  
  coef_gene <- read.csv("UniMutVariate_markergene.csv", header = T, sep=',')
  coef <- as.matrix(coef_gene[,2])
  rownames(coef) <- coef_gene[,1]
  
  coef_test <- my_overlap(coef, Data)
  plotp_Train <- separate2GroupsCox(as.vector(coef_test[[1]]), x_hat[, coef_test[[2]]], 
                                    plot.title = 'GSE1456_smfs', as.data.frame(y), 
                                    legend.outside = T)

  plot_train <- plotp_Train$plot
  plot_train
  
  ## for Xtile
  p_index <- cbind(y,plotp_Train$index)
  colnames(p_index) <- c(colnames(y), "riskscore")

  
  name <- dir[i]
  name1 <- str_split_fixed(name, "./", 2);
  name2 <- str_split_fixed(name1[2], "/", 2)
  name3 <- str_split_fixed(name2[2], "[.]", 2);
  name4 <- name3[1]
  name5 <- str_c(name4,"_ex")

  setwd(paste(path,'CNetCox/Feature_data', sep=''))
  path <- paste("./Xtile/",paste(name4,".txt", sep=""), sep="")
  write.table(p_index, path, quote = F, row.names = F, sep="\t")
  path <- paste("./Figure/", paste(name5,".pdf", sep=""), sep="")
  pdf(path, width = 4.5, height = 4.5, onefile = FALSE)
  print(plot_train)
  plot_train
  dev.off()

  print("***i***")
  print(i)

  
}




