## 2023.4.12 LLY create
## Split TCGA data into Training and Test data with clinical data


## clear
rm(list = ls())

## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- '/home/lly/R/'

## package
library(caret)
library(dplyr)

# Creat Files ----------------------------------------------------------------------
setwd(paste(path, 'CNetCox/Data/', sep=''))
dir.create("Data_train")
dir.create("Data_test")

# Load Data ----------------------------------------------------------------------
x <- read.table("TCGA_NEW/TCGA_BRCA_clin_546_1080_scale.txt", header = T, check.names = FALSE)
data <- data.frame(t(x))
dim(data)    # 1080  546
# View(data[,1:10])

# See mean and var
i <- 12
mean(data[,i])
var(data[,i])

# data split -----------------------------------------------------------

## set 20 seed
# condset <- c(3,13,14,17:20,23:25,28,32,33,36,38:41,44,46)
condset <- c(1:20)

## first
i <- 1
set.seed(123*i)
training_samples <- data$status %>% createDataPartition(p = 0.7, list = FALSE)
train_data  <- data[training_samples, ] 
test_data <- data[-training_samples, ] 
dim(train_data)    # 756 546
dim(test_data)     # 324 546

name <- as.character(i)
train <- as.matrix(train_data)
test <- as.matrix(test_data)

path <- paste("./Data_train/",paste(name,".txt",sep=""),sep="")
write.table(train,path,quote = F)

path <- paste("./Data_test/",paste(name,".txt",sep=""),sep="")
write.table(test,path,quote = F)

## second
for (i in condset) {
  print(paste("******* i= ********",i))
  set.seed(123*i)
  training_samples <- data$status %>% createDataPartition(p = 0.7, list = FALSE)
  ##  # sample*gene(col 1 -label)
  train_data  <- as.matrix(data[training_samples, ])  
  test_data <- as.matrix(data[-training_samples, ])  
  
  name <- as.character(i)
  train <- as.matrix(train_data)
  test <- as.matrix(test_data)
  
  path <- paste("./Data_train/",paste(name,".txt",sep=""),sep="")
  write.table(train,path,quote = F)
  
  path <- paste("./Data_test/",paste(name,".txt",sep=""),sep="")
  write.table(test,path,quote = F)
}
