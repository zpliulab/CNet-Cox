## 2023.4.12 LLY create
## Split TCGA data into Training and Test data with clinical data
## 2023.1.21 using five datasets and see the intersect

## clear
rm(list = ls())

## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- '/home/lly/R/'

## package
# library(caret)
# library(dplyr)

# Creat Files ------------------------------------------------------------------
setwd(paste(path, 'CNetCox/Data/', sep=''))



# load gene and net and cut-----------------------------------------------------
gene <- read.csv('TCGA_NEW/UNgene_component.csv',header = T)
net <- read.csv('TCGA_NEW/UNgene_comp_net.csv')
cut <- read.table("TCGA_NEW/cut_vector_UNgene.txt", header = T, check.names = FALSE, sep = "\t")

net[which(net[,1] == "CTCFL"),2]
net[which(net[,1] == "CTCFL"),2] == "EP300"

# load coefficient --------------------------------------------------------

coef <- read.csv("Result/Result1/theta1.csv", header = F, check.names = FALSE, sep = ",")
coefabs <- abs(coef)[-1,]

which(abs(cut[,1]) == 1)
value <- coefabs[which(abs(cut[,1]) == 1)]
coefsort <- sort(coefabs)

###############################################################################
## test once 
# value[1]
# num <- which(coefsort == value[1])
# lab <- c((num-2):(num+2)) 
###############################################################################

## The number of neigbourhold : l
l <- 2
lab <- c()
for (i in 1:5) {
  num <- which(coefsort == value[i])
  lab <- c(lab, c((num-l):(num+l)) )
}
unique(lab)
gene[which( abs(coef[-1,1]) %in% coefsort[lab]), ]
node <- gene[which( abs(coef[-1,1]) %in% coefsort[lab]), ]

###############################################################################
## test code, onle the cut node
# which(abs(coef[-1,1]) == value[1])    # 174
# cutoff <- max(coef[which(abs(cut[,1]) == 1)+1,1])
# which(abs(coef[-1,1]) <= cutoff)    # 174
# node <- gene[which(abs(coef[-1,1]) <= cutoff),]
###############################################################################

node_used <- node
net_used <- net
k1 <- which(net_used[,1] %in% node_used)   # 562
k2 <- which(net_used[,2] %in% node_used)   # 929
length(intersect(k1,k2))    # 4
used <- net_used[intersect(k1,k2),]
if (length(used) == 0){break;}

library(igraph)
PP <- graph_from_data_frame(used,directed = F)
# PP <- graph_from_data_frame(used,directed = T)
p1 <- simplify(PP)  
ed <- as_edgelist(p1, names = TRUE)
g <- p1
plot(g)

diameter(g)
get_diameter(g)

deg <- degree(g, mode="all")
V(g)$size <- deg * 6
plot(g, edge.color=edge.col, 
     edge.curved=.1, 
     edge.arrow.size = .3)


intersect(ed[,1], ed[,2])
marker <- union(ed[,1], ed[,2])
write.csv(marker, "Result/Result1/marker13.csv", row.names = F)




# my_CI function ----------------------------------------------------------
my_CI <- function(feature_plus, data){
  univ_formulas <- sapply(feature_plus, 
                          function(x) as.formula(paste('Surv(y1_hat$time, y1_hat$status)~', x)))
  univ_models <- lapply(univ_formulas, function(x){coxph(x, data = data)})
  CI <- univ_models[[1]]$concordance[6]
  return(CI)
}




# load test data and claculus the CI value --------------------------------
data <- read.table("Data_test/1.txt", header = T, sep = "")
Data_test <- t(data)
y1 <- t(Data_test[c(1,2),])
colnames(y1) <- c("time", "status")
y1_hat <- data.frame(y1)
x1 <- t(Data_test[-c(1,2),])
x1_hat <- data.frame(x1)


feature_plus <- paste(marker,collapse="+")

library(survival)
my_CI(feature_plus, x1_hat)
















# Load Data ----------------------------------------------------------------------
x <- read.table("Data_test/1.txt", header = T, check.names = FALSE)
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
