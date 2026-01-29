## 3023.4.21 Using other methods to select biomarkers
## copy the code from CoxReg Github TCGA_pro_clin_cox_1828_rep20_ridge.R


##############################################################################
# Compare CNet-Cox with other methods using same input dataset
##############################################################################

## clear
rm(list = ls())

## set pathway
# path <- '/home/lly/R/'
path <- '/Users/lilingyu/E/PhD/R/'
setwd(paste(path, 'CNetCox/Data/', sep=''))


## packages
library(caret) 
library(tidyverse)
library(survival)
library(ggpubr)
library(magrittr)
library(survminer)
library(glmSparseNet)
library(glmnet)
library(ncvreg)
# https://cran.r-project.org/src/contrib/Archive/APML0/
# install.packages("/Users/lilingyu/Downloads/APML0_0.10.tar.gz", repos=NULL, type="source")
library(APML0)
library(ncpen)


## revision v1 2026.01.16
# install.packages("timeROC")
library(timeROC)
# install.packages("pec")
library(pec)
# install.packages("rmda")
# library(rmda)
library(riskRegression)
library(rmda)

##############################################################################
# Some functions 
##############################################################################

## Calculate CI for a set of features of survival analysis using Cox models
my_CI <- function(feature_plus){
  # Create a list of univariate Cox regression formulas for each feature
  univ_formulas <- sapply(feature_plus, 
                          function(x) as.formula(paste('Surv(y1_hat$time, y1_hat$status)~', x)))
  
  # Fit Cox proportional hazards models for each formula using the data in x1_hat
  univ_models <- lapply(univ_formulas, function(x){coxph(x, data = x1_hat)})
  
  # Extract CI from the first fitted model(the 6th element of the concordance vector)
  CI <- univ_models[[1]]$concordance[6]

  return(CI)
}


my_overlap <- function(x, y){
  # Extract the first column from 'x' and keep only non-zero elements
  coefs.v <- x[,1] %>% { .[. != 0]}
  
  # Create a data frame with gene names and their corresponding coefficients,
  # then arrange by gene name and display as a table
  coefs.v %>% {
    data.frame(gene.name   = names(.),
               coefficient = .,
               stringsAsFactors = FALSE)
  } %>%
    arrange(gene.name) %>%
    knitr::kable()
  
  # Get the row names of the non-zero coefficients
  sele <- rownames(as.matrix(coefs.v))
  
  # Get the row names of 'y', excluding the first two rows
  gene <- rownames(y[-c(1,2),])
  
  # Find the intersection (overlap) between the selected gene and gene from 'y'
  overlap <- intersect(sele, gene)
  
  # Extract the coefficients corresponding to the overlapping genes
  lab <- x[,1] %>% { .[. != 0]} %>% as.matrix
  coefs.v <- lab[overlap,]
  
  # Return a list containing the coefficients of overlapping genes and itself
  my <- list(coefs.v, overlap)
  return(my)
}


##############################################################################
# Run other methods on-by-one
# ENet-Cox, L0-Cox, L1/2-Cox, Lasso-Cox, MCP-Cox, SCAD-Cox
##############################################################################
# Read data ------------------------------------------------------------------
data = read.table("TCGA_NEW/TCGA_BRCA_clin_546_1080.txt", header = T, check.names = FALSE)
mean(data[-c(1,2),1])

## scale
Data <- data.frame(cbind( t(data)[,c(1,2)], scale(t(data)[,-c(1,2)]) ))
# write.table(t(Data), file = "TCGA_BRCA_clin_1142_1080_scale.txt",quote=F,sep="\t")

l <- 50
mean(Data[,l])
var(Data[,l])

## From ‘day’ to ‘year’
# Data$time <- Data$time/365 
# View(Data[,1:10])

## Data-training set + test set
set.seed(123*3)
training.samples <- Data$status %>% createDataPartition(p = 0.7, list = FALSE)

## old, only splited once
train.data  <- Data[training.samples, ]
test.data <- Data[-training.samples, ]

## load train and test data, the same with before two lines
# train.data  <-read.table("Data_train/1.txt", header = T, check.names = FALSE)
# test.data  <-read.table("Data_test/1.txt", header = T, check.names = FALSE)


Data_train <- t(train.data)
x <- model.matrix(status ~., train.data)[,-c(1,2)]   
x_hat <- data.frame(x)
y <- as.matrix(train.data[,c(1,2)])  
colnames(y) <- c("time", "status")
y_hat <- data.frame(y)


Data_test <- t(test.data)
y1 <- t(Data_test[c(1,2),])
colnames(y1) <- c("time", "status")
y1_hat <- data.frame(y1)
x1 <- t(Data_test[-c(1,2),])
x1_hat <- data.frame(x1)


######################## glmnet -- ridge ########################
Coef_ridge <- matrix()    
Ci_ridge <- matrix()   
P_ridge <- matrix()  

# for(i in 1:20){
  i = 1
  set.seed(200*i)
  cv.ridge = cv.glmnet(x, y, family="cox", alpha = 0, nfolds = 10)
  lambda.min <- cv.ridge$lambda.min
  lambda.1se <- cv.ridge$lambda.1se
  print("***lambda.min??lambda.1se***")
  print(lambda.min)
  print(lambda.1se)
  plot(cv.ridge)

  ridge.model <- glmnet(x, y, family="cox", alpha = 0, lambda = cv.ridge$lambda.min)
  coef_ridge <- as.matrix(coef(ridge.model))
  which(coef_ridge[,1] != 0)
  feature_ridge <- colnames(x)[which(coef_ridge[,1] != 0)]
  # View(feature_ridge)
  ridge_plus <- paste(feature_ridge,collapse="+")


  CI_ridge <- my_CI(ridge_plus)
  CI_ridge

  coef_test <- my_overlap(coef_ridge, Data_test)
  plotp_ridge <- separate2GroupsCox(as.vector(coef_test[[1]]), x1_hat[, coef_test[[2]]],
                                    plot.title = 'Ridge', as.data.frame(y1), legend.outside = T)
  p_ridge <- plotp_ridge$pvalue

  
  predict_score <- predict(ridge.model, newx = as.matrix(x1_hat), type = "link")
  times <- c(1,3,5) * 365  
  time_roc <- timeROC(T = y1_hat$time, delta = y1_hat$status, marker = as.vector(predict_score), 
                      cause = 1, times = times, iid = TRUE)
  auc_res <- data.frame(Time=times/365, AUC=time_roc$AUC)
  print(auc_res)
  

  tcross <- rep(i, length(coef_ridge))                 
  step_ridge <- data.frame(cbind(coef_ridge, tcross))
  Coef_ridge <- cbind(Coef_ridge, step_ridge)   


  kcross <- rep(i, length(CI_ridge))
  temp_ridge <- data.frame(cbind(CI_ridge, kcross))
  Ci_ridge <- cbind(Ci_ridge, temp_ridge)    


  pcross <- rep(i, length(p_ridge))
  pemp_ridge <- data.frame(cbind(p_ridge, pcross))
  P_ridge <- cbind(P_ridge, pemp_ridge)   

  print(paste(i))

# }

## save files 
# write.csv(Coef_ridge, file = "Compare/coef_ridge.csv")
# write.csv(P_ridge, file = "Compare/P_ridge.csv")
# write.csv(Ci_ridge, file = "Compare/CI_ridge.csv")
# write.csv(auc_res, file="Compare/AUC_timeROC_ridge.csv", row.names = FALSE)
  
# pdf("Compare/timeROC_curve_ridge.pdf", width = 4, height = 4.5)
# plot(time_roc, time = times[1], col = "red", title = FALSE)
# plot(time_roc, time = times[2], add = TRUE, col = "blue")
# plot(time_roc, time = times[3], add = TRUE, col = "green")
# legend("bottomright", legend=c("1 year","3 years","5 years"), col=c("red","blue","green"), lwd=2)
# dev.off()

  
# pdf("Compare/km_curve_ridge.pdf", width = 4, height = 4.5)
# plotp_ridge
# dev.off()
  

# Set seed and reclassify data ------------------------------------------------------------
# set.seed(12345)
# training.samples <- Data$status %>% createDataPartition(p = 0.7, list = FALSE)
# train.data  <- Data[training.samples, ]
# test.data <- Data[-training.samples, ]
# 
# 
# Data_train <- t(train.data)
# x <- model.matrix(status ~., train.data)[,-c(1,2)]   
# x_hat <- data.frame(x)
# y <- as.matrix(train.data[,c(1,2)])  
# colnames(y) <- c("time", "status")
# y_hat <- data.frame(y)
# 
# 
# Data_test <- t(test.data)
# y1 <- t(Data_test[c(1,2),])
# colnames(y1) <- c("time", "status")
# y1_hat <- data.frame(y1)
# x1 <- t(Data_test[-c(1,2),])
# x1_hat <- data.frame(x1)


######################### glmnet -- Lasso  #############################
Coef_lasso <- matrix()      
Ci_lasso <- matrix()  
P_lasso <- matrix()

for(i in 1:1){
  # i = 1
  set.seed(200*i)
  
  cv.lasso = cv.glmnet(x, y, family="cox", alpha = 1, nfolds = 10)
  lambda.min <- cv.lasso$lambda.min
  lambda.1se <- cv.lasso$lambda.1se
  print("***lambda.min??lambda.1se***")
  print(lambda.min)
  print(lambda.1se)
  plot(cv.lasso)
  cv.lasso$cvm
  
  ## ????
  lasso.model <- glmnet(x, y, family="cox", alpha = 1, lambda = cv.lasso$lambda.min)
  coef_lasso <- as.matrix(coef(lasso.model))
  which(coef_lasso[,1] != 0)
  feature_lasso <- colnames(x)[which(coef_lasso[,1] != 0)]
  lasso_plus <- paste(feature_lasso,collapse="+")
  
  
  CI_lasso <- my_CI(lasso_plus)
  CI_lasso
  
  coef_test <- my_overlap(coef_lasso, Data_test)
  plotp_lasso <- separate2GroupsCox(as.vector(coef_test[[1]]), x1_hat[, coef_test[[2]]], 
                                    plot.title = 'Lasso', as.data.frame(y1), legend.outside = T)
  p_lasso <- plotp_lasso$pvalue
  
  
  predict_score <- predict(lasso.model, newx = as.matrix(x1_hat), type = "link")
  times <- c(1,3,5) * 365  
  time_roc <- timeROC(T = y1_hat$time, delta = y1_hat$status, marker = as.vector(predict_score), 
                      cause = 1, times = times, iid = TRUE)
  auc_res <- data.frame(Time=times/365, AUC=time_roc$AUC)
  print(auc_res)
  
  tcross <- rep(i, length(coef_lasso))               
  step_lasso <- data.frame(cbind(coef_lasso, tcross))
  Coef_lasso <- cbind(Coef_lasso, step_lasso)  
  
  
  kcross <- rep(i, length(CI_lasso)) 
  temp_lasso <- data.frame(cbind(CI_lasso, kcross))
  Ci_lasso <- cbind(Ci_lasso, temp_lasso)  
  
  
  pcross <- rep(i, length(p_lasso)) 
  pemp_lasso <- data.frame(cbind(p_lasso, pcross))
  P_lasso <- cbind(P_lasso, pemp_lasso)  
  
  print(paste(i)) 

}
# write.csv(Coef_lasso, file = "Compare/coef_lasso.csv")
# write.csv(P_lasso, file = "Compare/P_lasso.csv")
# write.csv(Ci_lasso, file = "Compare/CI_lasso.csv")
# write.csv(auc_res, file="Compare/AUC_timeROC_lasso.csv", row.names = FALSE)

# pdf("Compare/timeROC_curve_lasso.pdf", width = 4, height = 4.5)
# plot(time_roc, time = times[1], col = "red", title = FALSE)
# plot(time_roc, time = times[2], add = TRUE, col = "blue")
# plot(time_roc, time = times[3], add = TRUE, col = "green")
# legend("bottomright", legend=c("1 year","3 years","5 years"), col=c("red","blue","green"), lwd=2)
# dev.off()

# pdf("Compare/km_curve_lasso.pdf", width = 4, height = 4.5)
# plotp_lasso
# dev.off()

#################################### elastic net -------------------------------------------------------------

Coef_elastic <- matrix()       
Ci_elastic <- matrix()   
P_elastic <- matrix() 


for(i in 1:1){
  
  set.seed(200*i)
  
  cv.elastic = cv.glmnet(x, y, family="cox", alpha = 0.5, nfolds = 10)
  lambda.min <- cv.elastic$lambda.min
  lambda.1se <- cv.elastic$lambda.1se
  print("***lambda.min??lambda.1se***")
  print(lambda.min)
  print(lambda.1se)
  plot(cv.elastic)
  
  
  ## fit
  elastic.model <- glmnet(x, y, family="cox", alpha = 0.5, lambda = cv.elastic$lambda.min)
  coef_elastic <- as.matrix(coef(elastic.model))
  which(coef_elastic[,1] != 0)
  feature_elastic <- colnames(x)[which(coef_elastic[,1] != 0)]
  elastic_plus <- paste(feature_elastic,collapse="+")
  
  
  CI_elastic <- my_CI(elastic_plus)
  CI_elastic
  
  coef_test <- my_overlap(coef_elastic, Data_test)
  plotp_elastic <- separate2GroupsCox(as.vector(coef_test[[1]]), x1_hat[, coef_test[[2]]], 
                                    plot.title = 'Elastic', as.data.frame(y1), legend.outside = T)
  p_elastic <- plotp_elastic$pvalue
  
  
  predict_score <- predict(elastic.model, newx = as.matrix(x1_hat), type = "link")
  times <- c(1,3,5) * 365  
  time_roc <- timeROC(T = y1_hat$time, delta = y1_hat$status, marker = as.vector(predict_score), 
                      cause = 1, times = times, iid = TRUE)
  auc_res <- data.frame(Time=times/365, AUC=time_roc$AUC)
  print(auc_res)
  
  
  tcross <- rep(i, length(coef_elastic))               
  step_elastic <- data.frame(cbind(coef_elastic, tcross))
  Coef_elastic <- cbind(Coef_elastic, step_elastic)   
  
  
  kcross <- rep(i, length(CI_elastic)) 
  temp_elastic <- data.frame(cbind(CI_elastic, kcross))
  Ci_elastic <- cbind(Ci_elastic, temp_elastic)   

  pcross <- rep(i, length(p_elastic)) 
  pemp_elastic <- data.frame(cbind(p_elastic, pcross))
  P_elastic <- cbind(P_elastic, pemp_elastic)  
  
  print(paste(i)) 
  
}

# write.csv(Coef_elastic, file = "Compare/coef_elastic.csv")
# write.csv(P_elastic, file = "Compare/P_elastic.csv")
# write.csv(Ci_elastic, file = "Compare/CI_elastic.csv")
# write.csv(auc_res, file="Compare/AUC_timeROC_elastic.csv", row.names = FALSE)

# pdf("Compare/timeROC_curve_elastic.pdf", width = 4, height = 4.5)
# plot(time_roc, time = times[1], col = "red", title = FALSE)
# plot(time_roc, time = times[2], add = TRUE, col = "blue")
# plot(time_roc, time = times[3], add = TRUE, col = "green")
# legend("bottomright", legend=c("1 year","3 years","5 years"), col=c("red","blue","green"), lwd=2)
# dev.off()


# pdf("Compare/km_curve_elastic.pdf", width = 4, height = 4.5)
# plotp_elastic
# dev.off()

######################### ncvreg -- SCAD ?ͷ? #############################

X <- x
y <- y


Coef_SCAD <- matrix()      
Ci_SCAD <- matrix()   
P_SCAD <- matrix()


for(i in 1:1){
  
  set.seed(200*i)
  
  cv.SCAD <- cv.ncvsurv(X, y, family ="cox", penalty="SCAD")  
  lambda.min <- cv.SCAD$lambda.min 
  plot(cv.SCAD)
  print("*** lambda.min ***")
  print(lambda.min)
  
  
  SCAD.model <- ncvsurv(X, y, family ="cox", lambda = cv.SCAD$lambda.min, penalty="SCAD") 
  coef_SCAD <- as.matrix(coef(SCAD.model))
  which(coef_SCAD[,1] != 0)
  feature_SCAD <- colnames(x)[which(coef_SCAD[,1] != 0)]
  SCAD_plus <- paste(feature_SCAD,collapse="+")
  
  
  CI_SCAD <- my_CI(SCAD_plus)
  CI_SCAD
  
  coef_test <- my_overlap(coef_SCAD, Data_test)
  plotp_SCAD <- separate2GroupsCox(as.vector(coef_test[[1]]), x1_hat[, coef_test[[2]]], 
                                    plot.title = 'SCAD', as.data.frame(y1), legend.outside = T)
  p_SCAD <- plotp_SCAD$pvalue 

  
  predict_score <- predict(SCAD.model, as.matrix(x1_hat), type = "link")    # unique this line
  times <- c(1,3,5) * 365  
  time_roc <- timeROC(T = y1_hat$time, delta = y1_hat$status, marker = as.vector(predict_score), 
                      cause = 1, times = times, iid = TRUE)
  auc_res <- data.frame(Time=times/365, AUC=time_roc$AUC)
  print(auc_res)
  
  
  tcross <- rep(i, length(coef_SCAD))               
  step_SCAD <- data.frame(cbind(coef_SCAD, tcross))
  Coef_SCAD <- cbind(Coef_SCAD, step_SCAD)   
  
  
  kcross <- rep(i, length(CI_SCAD)) 
  temp_SCAD <- data.frame(cbind(CI_SCAD, kcross))
  Ci_SCAD <- cbind(Ci_SCAD, temp_SCAD)  
  
  
  pcross <- rep(i, length(p_SCAD)) 
  pemp_SCAD <- data.frame(cbind(p_SCAD, pcross))
  P_SCAD <- cbind(P_SCAD, pemp_SCAD)   
  
  print(paste(i)) 


}  

# write.csv(Coef_SCAD, file = "Compare/coef_SCAD.csv")
# write.csv(P_SCAD, file = "Compare/P_SCAD.csv")
# write.csv(Ci_SCAD, file = "Compare/CI_SCAD.csv")
# write.csv(auc_res, file="Compare/AUC_timeROC_SCAD.csv", row.names = FALSE)

# pdf("Compare/timeROC_curve_SCAD.pdf", width = 4, height = 4.5)
# plot(time_roc, time = times[1], col = "red", title = FALSE)
# plot(time_roc, time = times[2], add = TRUE, col = "blue")
# plot(time_roc, time = times[3], add = TRUE, col = "green")
# legend("bottomright", legend=c("1 year","3 years","5 years"), col=c("red","blue","green"), lwd=2)
# dev.off()


######################### ncvreg -- MCP #############################

X <- x
y <- y


Coef_MCP <- matrix()     
Ci_MCP <- matrix()   
P_MCP <- matrix()

# set.seed(123456)

for(i in 1:1){
  
  set.seed(31*i)
  
  cv.MCP <- cv.ncvsurv(X, y, family ="cox", penalty="MCP")  
  lambda.min <- cv.MCP$lambda.min 
  print("*** lambda.min ***")
  print(lambda.min)
  plot(cv.MCP)
  
  MCP.model <- ncvsurv(X, y, family ="cox", lambda = cv.MCP$lambda.min, penalty="MCP") 
  coef_MCP <- as.matrix(coef(MCP.model))
  which(coef_MCP[,1] != 0)
  feature_MCP <- colnames(x)[which(coef_MCP[,1] != 0)]
  MCP_plus <- paste(feature_MCP,collapse="+")
  if(MCP_plus == ""){
    set.seed(31)
    cv.MCP <- cv.ncvsurv(X, y, family ="cox", penalty="MCP")  
    lambda.min <- cv.MCP$lambda.min 
    print("*** lambda.min ***")
    print(lambda.min)
    

    MCP.model <- ncvsurv(X, y, family ="cox", lambda = cv.MCP$lambda.min, penalty="MCP") 
    coef_MCP <- as.matrix(coef(MCP.model))
    which(coef_MCP[,1] != 0)
    feature_MCP <- colnames(x)[which(coef_MCP[,1] != 0)]
    MCP_plus <- paste(feature_MCP,collapse="+")
  }

    
  CI_MCP <- my_CI(MCP_plus)
  CI_MCP
  
  coef_test <- my_overlap(coef_MCP, Data_test)
  plotp_MCP <- separate2GroupsCox(as.vector(coef_test[[1]]), x1_hat[, coef_test[[2]]], 
                                    plot.title = 'Ridge', as.data.frame(y1), legend.outside = T)
  p_MCP <- plotp_MCP$pvalue
  p_MCP
  
  
  predict_score <- predict(MCP.model, as.matrix(x1_hat), type = "link")    # unique this line
  times <- c(1,3,5) * 365  
  time_roc <- timeROC(T = y1_hat$time, delta = y1_hat$status, marker = as.vector(predict_score), 
                      cause = 1, times = times, iid = TRUE)
  auc_res <- data.frame(Time=times/365, AUC=time_roc$AUC)
  print(auc_res)
  
  
  tcross <- rep(i, length(coef_MCP))             
  step_MCP <- data.frame(cbind(coef_MCP, tcross))
  Coef_MCP <- cbind(Coef_MCP, step_MCP)   
  
  
  kcross <- rep(i, length(CI_MCP)) 
  temp_MCP <- data.frame(cbind(CI_MCP, kcross))
  Ci_MCP <- cbind(Ci_MCP, temp_MCP)  

  pcross <- rep(i, length(p_MCP)) 
  pemp_MCP <- data.frame(cbind(p_MCP, pcross))
  P_MCP <- cbind(P_MCP, pemp_MCP)   
  
  
  print(paste(i)) 
  
}

# write.csv(Coef_MCP, file = "Compare/coef_MCP.csv")
# write.csv(P_MCP, file = "Compare/P_MCP.csv")
# write.csv(Ci_MCP, file = "Compare/CI_MCP.csv")
# write.csv(auc_res, file="Compare/AUC_timeROC_MCP.csv", row.names = FALSE)

# pdf("Compare/timeROC_curve_MCP.pdf", width = 4, height = 4.5)
# plot(time_roc, time = times[1], col = "red", title = FALSE)
# plot(time_roc, time = times[2], add = TRUE, col = "blue")
# plot(time_roc, time = times[3], add = TRUE, col = "green")
# legend("bottomright", legend=c("1 year","3 years","5 years"), col=c("red","blue","green"), lwd=2)
# dev.off()


############################ APML0 L0 ###################################


y <- Data[,c(1,2)]
colnames(y) <- c("time", "status")
x <- as.matrix(Data[,-c(1,2)])
dim(x)    # 1080  544

library(APML0)
# Warning in install.packages :
#   package ‘APML0’ is not available for this version of R
# 
# A version of this package for your version of R might be available elsewhere,
# see the ideas at
# https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages


set.seed(123456)
l0.model = APML0(x, y, family="cox", penalty="Lasso", nfolds=10) # Lasso
print(l0.model)

# lambda.min <- l0.model$lambda.min
# print("*** lambda.min ***")
# print(lambda.min)

# coef_l0 <- as.matrix(l0.model$Beta)
# rownames(coef_l0) <- rownames(coef_lasso)

lambda.opt <- l0.model$lambda.opt
print("*** lambda.opt ***")
print(lambda.opt)

coef_l0 <- as.matrix(l0.model$Beta0)
rownames(coef_l0) <- rownames(data)[-c(1,2)]

which(coef_l0[,1] != 0)
feature_l0 <- colnames(x)[which(coef_l0[,1] != 0)]
l0_plus <- paste(feature_l0,collapse="+")


CI_l0 <- my_CI(l0_plus)
CI_l0


coef_test <- my_overlap(coef_l0, Data_test)
plotp_l0 <- separate2GroupsCox(as.vector(coef_test[[1]]), x1_hat[, coef_test[[2]]], 
                               plot.title = 'L0', as.data.frame(y1), legend.outside = T)

P_l0 <- plotp_l0$pvalue
P_l0


predict_score <- as.matrix(x1_hat) %*% as.matrix(coef_l0)     # unique this line
# predict_score <- predict(l0.model, new.x.mat = as.matrix(x1_hat), type = "y")# unique this line
times <- c(1,3,5) * 365  
time_roc <- timeROC(T = y1_hat$time, delta = y1_hat$status, marker = as.vector(predict_score), 
                    cause = 1, times = times, iid = TRUE)
auc_res <- data.frame(Time=times/365, AUC=time_roc$AUC)
print(auc_res)


# write.csv(coef_l0, file = "Compare/coef_l0.csv")
# write.csv(feature_l0, file = "Compare/feature_l0.csv", row.names = F)
# write.csv(CI_l0, file = "Compare/CI_l0.csv")
# write.csv(P_l0, file = "Compare/P_l0.csv")
# pdf(file = "Compare/l0.pdf",width = 4,height = 3)
# plotp_l0$plot
# dev.off()


# write.csv(auc_res, file="Compare/AUC_timeROC_l0.csv", row.names = FALSE)

# pdf("Compare/timeROC_curve_l0.pdf", width = 4, height = 4.5)
# plot(time_roc, time = times[1], col = "red", title = FALSE)
# plot(time_roc, time = times[2], add = TRUE, col = "blue")
# plot(time_roc, time = times[3], add = TRUE, col = "green")
# legend("bottomright", legend=c("1 year","3 years","5 years"), col=c("red","blue","green"), lwd=2)
# dev.off()



######################### ncpen  Bridge 1/2 #############################


x <- model.matrix(status ~., train.data)[,-c(1,2)]   # ɾ????time, statuss~.
x_hat <- data.frame(x)
y <- as.matrix(train.data[,c(1,2)])  
colnames(y) <- c("time", "status")
y_hat <- data.frame(y)


Data_test <- t(test.data)
y1 <- t(Data_test[c(1,2),])
colnames(y1) <- c("time", "status")
y1_hat <- data.frame(y1)
x1 <- t(Data_test[-c(1,2),])
x1_hat <- data.frame(x1)


x_1 <- as.matrix(cbind(x,y[,2]))
y_1 <- as.matrix(y[,1])


library(ncpen)

set.seed(12345)
Bridge.model <- ncpen(y.vec = y_1, x.mat = x_1, family="cox", penalty="mbridge")
opt.lambda <- gic.ncpen(Bridge.model, pch="*", type="b")$opt.lambda
print("*** optional lambda ***")
print(opt.lambda)

n <- which(Bridge.model$lambda == opt.lambda)


coef <- coef(Bridge.model)
coef_Bridge <- as.matrix(coef[, n])   

which(coef_Bridge[,1] != 0)
feature_Bridge <- colnames(x)[which(coef_Bridge[,1] != 0)]
Bridge_plus <- paste(feature_Bridge,collapse="+")


CI_Bridge <- my_CI(Bridge_plus)
CI_Bridge


coef_test <- my_overlap(coef_Bridge, Data_test)
plotp_Bridge <- separate2GroupsCox(as.vector(coef_test[[1]]), x1_hat[, coef_test[[2]]], 
                                   plot.title = 'Bridge', as.data.frame(y1), legend.outside = T)
plotp_Bridge$plot
plotp_Bridge$km
P_Bridge <- plotp_Bridge$pvalue
P_Bridge


# coef_Bridge是全部变量的系数（有的为0）
# x1_hat是全部变量
predict_score <- as.matrix(x1_hat) %*% as.matrix(coef_Bridge)     # unique this line
# predict_score <- predict(Bridge.model, new.x.mat = as.matrix(x1_hat), type = "y")# unique this line
times <- c(1,3,5) * 365  
time_roc <- timeROC(T = y1_hat$time, delta = y1_hat$status, marker = as.vector(predict_score), 
                    cause = 1, times = times, iid = TRUE)
auc_res <- data.frame(Time=times/365, AUC=time_roc$AUC)
print(auc_res)


# write.csv(coef, file = "Compare/coef_Bridge.csv")
# write.csv(coef_Bridge, file = "Compare/coef_Bridge.csv")
# write.csv(feature_Bridge, file = "Compare/feature_Bridge.csv", row.names = F)
# write.csv(CI_Bridge, file = "Compare/CI_Bridge.csv")
# write.csv(P_Bridge, file = "Compare/P_Bridge.csv")
# pdf(file = "Compare/Bridge.pdf",width = 4,height = 3)
plotp_Bridge$plot
# dev.off()


# write.csv(auc_res, file="Compare/AUC_timeROC_Bridge.csv", row.names = FALSE)
# pdf("Compare/timeROC_curve_Bridge.pdf", width = 4, height = 4.5)
# plot(time_roc, time = times[1], col = "red", title = FALSE)
# plot(time_roc, time = times[2], add = TRUE, col = "blue")
# plot(time_roc, time = times[3], add = TRUE, col = "green")
# legend("bottomright", legend=c("1 year","3 years","5 years"), col=c("red","blue","green"), lwd=2)
# dev.off()



# rm(list = ls())

###########################################################################
# All results analysis
###########################################################################
ridge <- read.table("Compare/coef_ridge.csv", header=TRUE, sep = ',')
lasso <- read.table("Compare/coef_lasso.csv", header=TRUE, sep = ',')
Elastic <- read.table("Compare/coef_elastic.csv", header=TRUE, sep = ',')
Bridge <- read.table("Compare/coef_Bridge.csv", header=TRUE, sep = ',')
l0 <- read.table("Compare/coef_l0.csv", header=TRUE, sep = ',')
SCAD <- read.table("Compare/coef_SCAD.csv", header=TRUE, sep = ',')
MCP <- read.table("Compare/coef_MCP.csv", header=TRUE, sep = ',')


# Selected gene -----------------------------------------------------------
coef_ridge2 <- as.matrix(ridge[which(ridge[,4] != 0),1])
coef_lasso2 <-  as.matrix(lasso[which(lasso[,3] != 0),1])
coef_Elastic2 <- as.matrix(Elastic[which(Elastic[,3] != 0),1])
coef_l02 <- as.matrix(l0[which(l0[,2] != 0),1])
coef_Bridge2 <- as.matrix(Bridge[which(Bridge[,2] != 0),1]) 
coef_SCAD2 <- as.matrix(SCAD[which(SCAD[,3] != 0),1])
coef_MCP2 <- as.matrix(MCP[which(MCP[,3] != 0),1])


summary(coef_ridge2)
summary(coef_lasso2)
summary(coef_Elastic2)
summary(coef_l02)
summary(coef_Bridge2)
summary(coef_SCAD2)
summary(coef_MCP2)


## see the independent marker genes selected by CNet-Cox methods
marker <- as.matrix(read.csv("Result/Result3/marker3.csv", header = T, sep=','))
uniongene <- Reduce(union, list(coef_lasso2, coef_Elastic2,
                                coef_l02, coef_Bridge2, coef_SCAD2,coef_MCP2))
indengene <- setdiff(marker, uniongene)
## see Resistance gene RAD21 and TNFSF11, whether indenpent gene by CNet-Cox
which(indengene == "RAD21")    # 20
which(indengene == "TNFSF11")  # 39


# Load CI file --------------------------------------------------------------
Ci_ridge <- read.table("Compare/CI_ridge.csv", header=TRUE, sep = ',')
Ci_lasso <- read.table("Compare/CI_lasso.csv", header=TRUE, sep = ',')
Ci_Elastic <- read.table("Compare/CI_elastic.csv", header=TRUE, sep = ',')
Ci_Bridge <- read.table("Compare/CI_Bridge.csv", header=TRUE, sep = ',')
Ci_l0 <- read.table("Compare/CI_l0.csv", header=TRUE, sep = ',')
Ci_SCAD <- read.table("Compare/CI_SCAD.csv", header=TRUE, sep = ',')
Ci_MCP <- read.table("Compare/CI_MCP.csv", header=TRUE, sep = ',')


# Extract my_ci from file ---------------------------------------------------
Ci_ridge1 <- as.matrix(Ci_ridge[,3])
Ci_lasso1 <- as.matrix(Ci_lasso[,3])
Ci_Elastic1 <- as.matrix(Ci_Elastic[,3])
Ci_l01 <- Ci_l0[,2]
Ci_Bridge1 <- Ci_Bridge[,2]
Ci_SCAD1 <- as.matrix(Ci_SCAD[,3])
Ci_MCP1 <- as.matrix(Ci_MCP[,3])


mean(Ci_ridge1)
mean(Ci_lasso1)
mean(Ci_Elastic1)
Ci_l01
Ci_Bridge1
mean(Ci_SCAD1)
mean(Ci_MCP1)


# Load P file---------------------------------------------------------------
P_ridge <- read.table("Compare/P_ridge.csv", header=TRUE, sep = ',')
P_lasso <- read.table("Compare/P_lasso.csv", header=TRUE, sep = ',')
P_Elastic <- read.table("Compare/P_elastic.csv", header=TRUE, sep = ',')
P_Bridge <- read.table("Compare/P_Bridge.csv", header=TRUE, sep = ',')
P_l0 <- read.table("Compare/P_l0.csv", header=TRUE, sep = ',')
P_SCAD <- read.table("Compare/P_SCAD.csv", header=TRUE, sep = ',')
P_MCP <- read.table("Compare/P_MCP.csv", header=TRUE, sep = ',')


P_ridge1 <- as.matrix(P_ridge[,3])
P_lasso1 <- as.matrix(P_lasso[,3])
P_Elastic1 <- as.matrix(P_Elastic[,3])
P_l01 <- P_l0[,2]
P_Bridge1 <- P_Bridge[,2]
P_SCAD1 <- as.matrix(P_SCAD[,3])
P_MCP1 <- as.matrix(P_MCP[,3])

mean(P_ridge1)
mean(P_lasso1)
mean(P_Elastic1)
P_l01
P_Bridge1
mean(P_SCAD1)
mean(P_MCP1)



# For AUROC ---------------------------------------------------------------
## Load at revision at 2026/01/16
ridge <- read.table("Compare/AUC_timeROC_ridge.csv", header=TRUE, sep = ',')
lasso <- read.table("Compare/AUC_timeROC_lasso.csv", header=TRUE, sep = ',')
Elastic <- read.table("Compare/AUC_timeROC_elastic.csv", header=TRUE, sep = ',')
Bridge <- read.table("Compare/AUC_timeROC_Bridge.csv", header=TRUE, sep = ',')
l0 <- read.table("Compare/AUC_timeROC_l0.csv", header=TRUE, sep = ',')
SCAD <- read.table("Compare/AUC_timeROC_SCAD.csv", header=TRUE, sep = ',')
MCP <- read.table("Compare/AUC_timeROC_MCP.csv", header=TRUE, sep = ',')
CNet <- read.table("Compare/AUC_timeROC_CNet.csv", header=TRUE, sep = ',')


library(ggplot2)
library(dplyr)
CNet$method <- "CNet-Cox"
Elastic$method <- "ENet-Cox"
l0$method <- "L0-Cox"
Bridge$method <- "L1/2-Cox"
lasso$method <- "Lasso-Cox"
MCP$method <- "MCP-Cox"
ridge$method <- "Ridge-Cox"
SCAD$method <- "SCAD-Cox"


auc_all <- bind_rows(CNet, Elastic, l0, Bridge, lasso, MCP, ridge, SCAD)
auc_all$method <- factor(auc_all$method, levels = c("CNet-Cox","ENet-Cox","L0-Cox","L1/2-Cox","Lasso-Cox","MCP-Cox","Ridge-Cox","SCAD-Cox"))
head(auc_all)

# pdf("AUC_comparison_nature.pdf", width=4, height=4, family="sans")  
ggplot(auc_all, aes(x = method, y = AUC, fill = method)) +
  geom_boxplot(width = 0.5, outlier.size = 1, outlier.stroke = 0.5, outlier.alpha = 0.7) +
  theme_minimal(base_size = 12, ) +
  labs(x = "Method", y = "AUC", title = "Comparison of AUC across Methods") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7, colour = "black"),
        axis.ticks = element_line(size = 0.7, colour = "black"),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5),
        plot.margin = margin(10, 10, 10, 10)
        )
# dev.off()



## median
# 计算均值和中位数
stats <- auc_all %>%
  group_by(method) %>%
  summarize(mean_auc = mean(AUC),
            median_auc = median(AUC))

library(RColorBrewer)

# pdf("Compare/AUC_comparison_8.pdf", width=8, height=6, family="sans")

ggplot(auc_all, aes(x = method, y = AUC, fill = method)) +
  geom_boxplot(width = 0.5, outlier.size = 1, outlier.stroke = 0.5, outlier.alpha = 0.7) +
  # Mean point
  stat_summary(fun = mean, geom = "point", shape = 21, size = 3, fill = "white", color = "red", stroke=1.2) +
  # Mean value label
  geom_text(
    data = stats,
    # aes(x = method, y = mean_auc + 0.02, label = sprintf("Mean: %.3f", mean_auc)),
    aes(x = method, y = mean_auc + 0.01, label = sprintf("%.3f", mean_auc)),
    inherit.aes = FALSE,
    color = "red",
    size = 3.6,
    fontface = "bold"
  ) +
  # Median value label (added on the horizontal line in the middle of the box).
  geom_text(
    data = stats,
    # aes(x = method, y = median_auc, label = sprintf("Median: %.3f", median_auc)),
    aes(x = method, y = median_auc - 0.02, label = sprintf("%.3f", median_auc)),
    inherit.aes = FALSE,
    size = 3.6,
    vjust = -1, 
    # color = "red",
    fontface = "bold"
  ) +
  scale_fill_brewer(palette = "Set2") +  
  theme_minimal(base_size = 12) +
  labs(x = "Method", y = "AUC", title = "Comparison of AUC across Methods") +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.7, colour = "black"),
    axis.ticks = element_line(size = 0.7, colour = "black"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, hjust = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  )

# dev.off()
