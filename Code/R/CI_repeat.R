## 2023.4.12 LLY create, 2026.1.16 LLY used
## Split TCGA data into Training and Test data with clinical data
## 2023.1.21 using five datasets and see the intersect
## 2023.4.23 using nine datasets and see each of them


##############################################################################
# Figure 2B. 68 prognostic bioomarkers
##############################################################################

## clear
rm(list = ls())

## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- '/home/lly/R/'
setwd(paste(path, 'CNetCox/Data/', sep=''))

## load function
source(paste(path, 'CNetCox/R/MyFunctions.R', sep=''))

## load gene and net and cut
gene <- read.csv('TCGA_NEW/UNgene_component.csv',header = T)
net <- read.csv('TCGA_NEW/UNgene_comp_net.csv')
cut <- read.table("TCGA_NEW/cut_vector_UNgene.txt", header = T, check.names = FALSE, sep = "\t")


# Marker genes ------------------------------------------------------------
ci <- c()
p <- c()

# for (k in c(4.9e-2,
#             # 5.0e-2, 5.1e-2, 5.2e-2, 5.3e-2, 
#             # 5.31e-2, 5.33e-2, 5.35e-2, 5.37e-2, 5.39e-2,
#             # 5.4e-2, 
#             # 5.41e-2, 5.43e-2, 5.45e-2, 5.47e-5, 5.49e-2,
#             5.5e-2, 5.6e-2, 5.7e-2, 5.8e-2, 5.9e-2, 6e-2,
#             6.1e-2, 6.2e-2, 6.3e-2, 6.4e-2, 6.5e-2,
#             # 1.5e-1, 1.6e-1, 1.7e-1, 1.8e-1, 1.9e-1,
#             2e-1)) {
#   marker3 <- markerselect(1,k)[[1]]
#   coef3 <- markerselect(1,k)[[2]]
  


# marker1 <- markerselect(1,1.5e-1)[[1]]
# marker2 <- markerselect(2,4e-2)[[1]]
marker3 <- markerselect(3,5.15e-2)[[1]]
# marker4 <- markerselect(4,9e-2)[[1]]
# marker5 <- markerselect(5,5e-2)[[1]]
# marker6 <- markerselect(6,5e-2)[[1]]
# # marker7 <- markerselect(7,2e-1)[[1]]
# marker8 <- markerselect(8,2e-1)[[1]]
# marker9 <- markerselect(9,2e-1)[[1]]
# marker10 <- markerselect(10,2e-1)[[1]]

# marker1 <- markerunion(1)[[1]]
# marker2 <- markerunion(2)[[1]]
# marker3 <- markerunion(3)[[1]]
# marker4 <- markerunion(4)[[1]]
# marker5 <- markerunion(5)[[1]]
# marker6 <- markerunion(6)[[1]]
# # marker7 <- markerunion(7)[[1]]
# # marker8 <- markerunion(8)[[1]]
# marker9 <- markerunion(9)[[1]]
# marker10 <- markerunion(10)[[1]]

# markeru <- Reduce(union,
                  # list(marker1, marker2, marker3, marker4, marker5))

# write.csv(marker3, "Result/Result3/marker3.csv", row.names = F)

# Marker coefs ------------------------------------------------------------
# coef1 <- markerselect(1,1.5e-1)[[2]]
# coef2 <- markerselect(2,4e-2)[[2]]
coef3 <- markerselect(3,5.15e-2)[[2]]
# coef4 <- markerselect(4,9e-2)[[2]]
# coef5 <- markerselect(5,5e-2)[[2]]
# coef6 <- markerselect(6,5e-2)[[2]]
# # coef7 <- markerunion(7,2e-1)[[2]]
# # coef8 <- markerunion(8,2e-1)[[2]]
# coef9 <- markerselect(9,2e-1)[[2]]
# coef10 <- markerselect(10,2e-1)[[2]]


# coef1 <- markerunion(1)[[2]]
# coef2 <- markerunion(2)[[2]]
# coef3 <- markerunion(3)[[2]]
# coef4 <- markerunion(4)[[2]]
# coef5 <- markerunion(5)[[2]]
# coef6 <- markerunion(6)[[2]]
# # coef7 <- markerunion(7)[[2]]
# # coef8 <- markerunion(8)[[2]]
# coef9 <- markerunion(9)[[2]]
# coef10 <- markerunion(10)[[2]]

# coefu <- cbind(coef1, coef2, coef3, coef4, coef5)
# coefmean <- as.matrix(apply(coefu,1,mean)) 
# rownames(coefmean) <- gene[,1]

# write.csv(markeru, "Result/markerunion5.csv", row.names = F)


# see the network ---------------------------------------------------------
filenum <- "3"
markerinput <- marker3
coefinput <- coef3

g <- markerPlot(markerinput, net)[[1]]

# deg <- degree(g, mode="all")
# V(g)$size <- deg * 2
# plot(g,
#      edge.curved=.1,
#      edge.arrow.size = .3)

markerNet <- markerPlot(markerinput, net)[[2]]
# write.csv(markerNet, "Result/Result3/markerNet3.csv", row.names = F, quote = F)


##############################################################################
# Figure 2G. Survival analysis results of 68 prognostic BRCA markers.
##############################################################################
# load test data and calculus the CI value --------------------------------
# filenum <- "1"
# readtestdata <- function(filenum){
data <- read.table(paste("Data_test/", filenum, ".txt", sep = ""), header = T, sep = "")
Data_test <- t(data)
y1 <- t(Data_test[c(1,2),])
colnames(y1) <- c("time", "status")
y1_hat <- data.frame(y1)
x1 <- t(Data_test[-c(1,2),])
x1_hat <- data.frame(x1)
# }

library(survival)
feature_plus <- paste(markerinput,collapse="+")
my_CI(feature_plus, x1_hat)

library(dplyr)   
library(glmSparseNet)
coef_test <- my_overlap(coefinput, Data_test)
plotp <- separate2GroupsCox(as.vector(coef_test[[1]]), x1_hat[, coef_test[[2]]], 
                                  # plot.title = 'Lasso', 
                            as.data.frame(y1), legend.outside = T)

# plotp$plot
plotp$pvalue
# plotp$km

# ci <- c(ci, my_CI(feature_plus, x1_hat))
# p <- c(p, plotp$pvalue)
# }
# ci 
# p


###########################################################################
###########################################################################
# 1. 计算风险得分
sel_feats <- coef_test[[2]]
sel_coefs <- as.numeric(coef_test[[1]])
X_sub <- as.matrix(x1_hat[, sel_feats])
predict_score <- as.vector(X_sub %*% sel_coefs)
# predict_score <- xtile3    # same

# 2. 计算time-dependent AUC
library(timeROC)
times <- c(1, 3, 5) * 365
time_roc <- timeROC(
  T = as.numeric(y1_hat$time),
  delta = as.numeric(y1_hat$status),
  marker = predict_score,
  cause = 1,
  times = times,
  iid = TRUE
)
auc_res <- data.frame(Time = times / 365, AUC = time_roc$AUC)
print(auc_res)


## Improve version
library(glmnet)
fit <- cv.glmnet(as.matrix(x1_hat), Surv(y1_hat$time, y1_hat$status), family="cox", nfolds=5)
predict_score <- predict(fit, newx=as.matrix(x1_hat), s="lambda.min", type="link")

library(timeROC)
times <- c(1, 3, 5) * 365
time_roc <- timeROC(
  T = as.numeric(y1_hat$time),
  delta = as.numeric(y1_hat$status),
  marker = predict_score,
  cause = 1,
  times = times,
  iid = TRUE
)
auc_res <- data.frame(Time = times / 365, AUC = time_roc$AUC)
print(auc_res)
## save 
# write.csv(auc_res, file="Compare/AUC_timeROC_CNet.csv", row.names = FALSE)

# pdf("Compare/timeROC_curve_CNet.pdf", width = 4, height = 4.5)
# plot(time_roc, time = times[1], col = "red", title = FALSE)
# plot(time_roc, time = times[2], add = TRUE, col = "blue")
# plot(time_roc, time = times[3], add = TRUE, col = "green")
# legend("bottomright", legend=c("1 year","3 years","5 years"), col=c("red","blue","green"), lwd=2)
# dev.off()
###########################################################################
###########################################################################


## reference from six methods
# predict_score <- predict(lasso.model, newx = as.matrix(x1_hat), type = "link")
# times <- c(1,3,5) * 365  
# time_roc <- timeROC(T = y1_hat$time, delta = y1_hat$status, marker = as.vector(predict_score), 
#                     cause = 1, times = times, iid = TRUE)
# auc_res <- data.frame(Time=times/365, AUC=time_roc$AUC)
# print(auc_res)


xtile3 <- as.matrix(cbind(y1, as.matrix(x1_hat[, coef_test[[2]]]) %*% as.matrix(coef_test[[1]])))
colnames(xtile3) <- c("time","status", "riskscore")
# write.table(xtile3, paste(path, "CNetCox/Data/Result/Result3/xtile3.txt", sep = ""),
#             row.names = F, sep = "\t", quote = F)


cut_off <- mean(as.matrix(x1_hat[, coef_test[[2]]]) %*% as.matrix(coef_test[[1]]))
cut_off <- 0.6
data$riskscore <- ifelse(as.matrix(x1_hat[, coef_test[[2]]]) %*% as.matrix(coef_test[[1]]) > cut_off,"High","Low")
# cox <- coxph(Surv(y1_hat$time, y1_hat$status) ~ risk, data)
# summary(cox)


# data$time <- data$time/12
# data$riskscore <- ifelse(risk_score > cut_off, 'high','low')
# table(data$risk_score)

fit <- survfit(Surv(time, status)~riskscore, data = data)
library(survminer)
library(survival)
library(litedown)
# library(ggsurvplot)
p <- ggsurvplot(fit, data = data, 
                conf.int = F,  
                # surv.median.line = "hv",   
                risk.table = TRUE,  
                tables.height= 0.25,  
                cumcensor = T,     
                legend = c(0.83,0.95),  
                ## P value
                pval = TRUE, 
                pval.size=6, 
                font.pval= c(14, "bold", "black"),
                pval.coord = c(0.00, 0.05), # 
                # legend
                # legend.title = '', # gene_name
                # legend.labs=c("High risk", "Low risk"), 
                # font.legend= c(14, "plain", "black"),  
                # # font.main = c(100, "bold", "black"),
                # # xlim = c(0,72), # present narrower X axis, but not affect
                # # survival estimates.
                # palette=c("red", "blue"),
                # font.x = c(14, "plain", "black"),
                # font.y = c(14, "plain", "black"),  
                # font.tickslab = c(14, "plain", "black"),  
                # xlab = "Time in years", # customize X axis label. year
                # break.time.by = 2
)  


# pdf("TCGA_NEW/TCGA_68_os.pdf", width = 5.0, height = 6, onefile = FALSE)
p
# dev.off()





##############################################################################
# Table S4. Univariate and multivariate Cox analysis
##############################################################################
# Extract marker data and label -------------------------------------------
Data1 = read.table("TCGA_NEW/TCGA_BRCA_clin_546_1080_scale.txt", header = T, 
                   check.names = FALSE)
# View(Data1[,1:10])

gene <- as.matrix(read.csv("Result/Result3/marker3.csv", header = T, sep=','))
colnames(gene) <- c('gene')

Data2 <- cbind(rownames(Data1), Data1)
colnames(Data2) <- c('gene', colnames(Data1))
# View(Data2[,1:10])

genedata <- merge(gene, Data2, by = "gene")    # 68
genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
# View(genedata1[,1:10])
genedata1[1,1]
genedata2 <- rbind(Data1[c(1,2),],genedata1)    # 70
# View(genedata2[,1:10])
# write.table(genedata2,"TCGA_NEW/TCGA_Marker_Expr.txt",quote=F,sep="\t")  


# Univariable and Multivalable Cox analysis -------------------------------
## Step 1：Univariate Cox ##############################
library(survival)
library(plyr)

markerdata = read.table("TCGA_NEW/TCGA_Marker_Expr.txt", header = T, check.names = FALSE)
Data <- data.frame(t(markerdata))
str(Data)

## 928 -- 0;  152  -- 1
Data$status<-factor(Data$status)
summary(Data$status)

## step 1:
y <- Surv(time=Data$time, event=Data$status==0)  

## step 2: Uni_cox_model
Uni_cox_model<- function(x){
  # x <- variable_names
  FML <- as.formula(paste0 ("y~",x))
  cox <- coxph(FML,data=Data)
  cox1 <-summary(cox)
  coef <- cox1$coefficients[,1]
  HR <- round(cox1$coefficients[,2],2)
  PValue <- round(cox1$coefficients[,5],3)
  CI5 <- round(cox1$conf.int[,3],2)
  CI95 <- round(cox1$conf.int[,4],2)
  Uni_cox_model<- data.frame('Characteristics' =x,
                             names <-rownames(cox1$conf.int),
                             'coef' = coef,
                             'HR' = HR,
                             'CI5' = CI5,
                             'CI95' = CI95,
                             'Uni_P' = PValue
  )
  return(Uni_cox_model)}  

## step 3: see marker genes's col number and select marker genes
names(Data)
variable_names <- colnames(Data)[c(3:70)]

## step 4: construct model
Uni_cox <- lapply(variable_names, Uni_cox_model)
Uni_cox <- ldply(Uni_cox, data.frame)

## step 5: HR+95% CI+P 
Uni_cox$HR.CI95 <- paste0(Uni_cox$HR," (",Uni_cox$CI5,'-',Uni_cox$CI95,")")
Uni_cox <- Uni_cox[,-4:-6] #HR (95% CI)+P
# View(Uni_cox)

which(Uni_cox[,4] <= 0.05)
Uni_cox[which(Uni_cox[,4] <= 0.05),c(1,3)]


## Univariate Cox p<0.05 and Multivariate Cox #################################
## step 1: p<0.05
Uni_cox$Characteristics[Uni_cox$Uni_P<0.05]

## positive or negtivate
Uni_cox$Characteristics[(Uni_cox$Uni_P<0.05) & (Uni_cox$coef <0)]    # 24
Uni_cox$Characteristics[(Uni_cox$Uni_P<0.05) & (Uni_cox$coef >0)]    # 15

## step 2: construct Multivariate model
mul_cox_model<- as.formula(paste0 ("y~",
                                   paste0(Uni_cox$Characteristics[Uni_cox$Uni_P<0.05],
                                          collapse = "+")))
mul_cox<-coxph(mul_cox_model,data=Data)
cox4 <- summary(mul_cox) 
coef <- cox4$coefficients[,c(1,5)]
mul_coef <- coef[which(coef[,2] <= 0.05),1]
# View(mul_coef)
# write.csv(mul_coef, file = "TCGA_NEW/UniMutVariate_markergene.csv")

## step 3:
mul_corf <- cox4$coefficients[,1] 
mul_HR <- round(cox4$coefficients[,2],2) 
mul_PValue<- round(cox4$coefficients[,5],4) 
mul_CI1<-round(cox4$conf.int[,3],2)
mul_CI2<-round(cox4$conf.int[,4],2)

## step 4: multivariate cox1
## 4-1 HR(95%CI)+P
mul_HR.CI95 <- paste(mul_HR,"(",mul_CI1,'-',mul_CI2,")")
mul_cox1 <- data.frame("coef"=mul_corf, "mul_HR.CI95"=mul_HR.CI95,"P"=mul_PValue)

## 4-2 HR+95%CI+P
#mul_CI<-paste(mul_CI1,'-',mul_CI2)
#mul_cox1<- data.frame("HR"=mul_HR,"mul_CI"=mul_CI, "P"=mul_PValue)


## Step 3：Univariate & multivariate data integration ######################
Uni_cox <- Uni_cox[,-1]
colnames(Uni_cox)[1] <- 'Characteristics'
mul_cox1 <- cbind(rownames(mul_cox1), mul_cox1, row.names=NULL); names(mul_cox1 )[1]<-"Characteristics"
table2 <- merge.data.frame(Uni_cox, mul_cox1, by="Characteristics", all = T, sort = T)
table3 <- table2[,c(1,2,4,3,5,6,7)]
colnames(table3) <- c("Feature", "Coef_Uni", "HR(95%CI)_Uni", "P_Uni",
                      "Coef_Mul", "HR(95%CI)_Mul", "P_Mul")
View(table3)
# write.csv(table3, file = "TCGA_NEW/UniMutVariate_Cox_Coef.csv",row.names = F)


library(stargazer)
stargazer(table3, summary=FALSE, rownames=FALSE) 
stargazer(table3[,], summary=FALSE, rownames=FALSE) 


##############################################################################
# Figure 3F,G: Anthracycline-sensitive validation. 
##############################################################################
# Extract Drug Resistance data  ------------------------------------------------
ResisDE <- read.csv("DrugRest/biomolecules-12-01834-s001/supplementary Table S2.csv", header = T, sep = ",")
which(ResisDE$Genename %in% table3$Feature)
markerDE <- ResisDE[which(ResisDE$Genename %in% table3$Feature),]
which(markerDE$log2FC > 0)
length(which(markerDE$log2FC > 0))    # 25
length(which(markerDE$log2FC < 0))    # 43

## plot bar plot of p value
markerFC <- markerDE[,c(1, 6)]
markerFC$Genename
which(markerFC$Genename == "RAD21")    # 65
which(markerFC$Genename == "TNFSF11")  # 68
markerP <- markerFC[c(which(markerFC$Genename == "RAD21"),which(markerFC$Genename == "TNFSF11")),]
colnames(markerP) <- c("Gene", "P")

library(ggplot2)
p4 <- ggplot(data = markerP, aes(x = P, y = Gene)) +
  # geom_bar(stat = "identity", fill = "steelblue") +
  geom_col(aes(fill = P > 0), color = "gray") +
  scale_fill_manual(values = c("#d5d9e0", "#a4afd5")) +
  # scale_x_continuous(limits = c(-0.2, 0.3), expand = c(0, 0)) +
  ## change x and y axies
  # coord_flip()+
  # theme_classic() +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 14),
        axis.line = element_line(colour = "black", size = 0.5),
        legend.position = "none",
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 16),
        plot.caption = element_text(size = 12, hjust = 1, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        plot.margin = margin(40, 40, 20, 20, "pt")) +
  labs(x = "P value")

# pdf("DrugRest/DistribP.pdf",width = 4, height = 5)
p4
# dev.off()


## plot bar plot of LogFC
markerDE <- ResisDE[which(ResisDE$Genename %in% table3$Feature),]
markerFC <- markerDE[,c(1, 3)]

which(markerFC[,2] < 0)
markerFC[which(markerFC[,2] < 0),3] <- "Down"    # 43
markerFC[which(markerFC[,2] > 0),3] <- "Up"      # 25
markerFC <- markerFC[,c(1,3,2)]
colnames(markerFC) <- c("Gene", "variable", "value" )
markerFC$value <- round(markerFC$value, 3)

library(ggplot2)
library(reshape2)
df <- markerFC
p<- ggplot(df, aes(
  x = factor(Gene,levels = unique(Gene)),   
  y = ifelse(variable == "variable", -value, value),  
  fill = variable)) +
  geom_bar(stat = 'identity')+   
  coord_flip()+
  geom_text(                                                  
    aes(label=value, 
        hjust = ifelse(variable == "variable", -0.4, 1.1)),size=2) +
  scale_fill_manual(values = c('#fec79e','#8ec4cb'))+
  labs(x='Gene',y='Log2FC') +
  theme_test(base_size = 10) 

# pdf("DrugRest/DistribLogFC.pdf", width = 5, height = 8)
p
# dev.off()


##############################################################################
# Figure 4B: Internal validation result using six-gene PRS
##############################################################################
# Extract PRS data and label -------------------------------------------
data = read.table("TCGA_NEW/TCGA_BRCA_clin_546_1080_scale.txt", header = T, check.names = FALSE)
# data = read.table("TCGA_NEW/TCGA_BRCA_clin_546_1080.txt", header = T, check.names = FALSE)
datat <- as.data.frame(t(data))
Data <- as.matrix(data)

x <- model.matrix(status ~., datat)[,-c(1,2)]  
x_hat <- data.frame(x)
y <- as.matrix(datat[,c(1,2)])  
colnames(y) <- c("time", "status")
y_hat <- data.frame(y)

# coef_gene <- read.csv("TCGA_NEW/UniMutVariate_Cox_Coef.csv", header = T, sep=',')
coef_gene <- read.csv("TCGA_NEW/UniMutVariate_markergene.csv", header = T, sep=',')
coef <- as.matrix(coef_gene[,2])
rownames(coef) <- coef_gene[,1]


##############################################################################
# Calculating riskscore using separate2GroupsCox
##############################################################################
# Load function -------------------------------------------------------------
source(paste(path, 'CNetCox/R/myoverlap_separate2GroupsCox.R', sep=''))
library(ggpubr)
library(magrittr)
library(survminer)

coef_test <- my_overlap(coef, Data)
plotp_Train <- separate2GroupsCox(as.vector(coef_test[[1]]), x_hat[, coef_test[[2]]], 
                                  plot.title = 'TCGA_OS_PRS_gene', as.data.frame(y), 
                                  legend.outside = T)
plot_train <- plotp_Train$plot
plot_train

## for Xtile
p_index <- cbind(y,plotp_Train$index)
colnames(p_index) <- c(colnames(y), "riskscore")
# write.table(p_index, file = "TCGA_NEW/TCGA_OS_PRS_gene_scale.txt", quote = F, row.names = F, sep="\t")
# write.table(p_index, file = "TCGA_NEW/TCGA_OS_PRS_gene.txt", quote = F, row.names = F, sep="\t")
# write.table(p_index, file = "TCGA_NEW/TCGA_OS_PRS_gene.csv", quote = F, row.names = F)
# write.table(p_index, file = "TCGA_NEW/TCGA_OS_PRS_gene_scale.csv", quote = F, row.names = F)


##############################################################################
# Figure 4B
##############################################################################
# Extract features on TCGA --------------------------------------------------
library(survminer)
library(survival)
library(ggplot2)
# data = read.table("TCGA_NEW/TCGA_OS_PRS_gene.txt", header = T, check.names = FALSE)
data = read.table("TCGA_NEW/TCGA_OS_PRS_gene_scale.txt", header = T, check.names = FALSE)

gene_name <- colnames(data)[3]
exprSet <- as.data.frame(t(data))

## Set cutoff value  
alpha <- 0.75 # (928/1080)
risk_score  <- t(as.matrix(exprSet[gene_name,]))
## This row: Error in xtfrm.data.frame(x) : cannot xtfrm data frames (26/01/16)
# cut_off <- rep(as.numeric(quantile(exprSet[gene_name,],alpha)), dim(exprSet)[2])
exprSet <- as.matrix(exprSet)
cut_off <- rep(as.numeric(quantile(exprSet[gene_name, ], alpha)), ncol(exprSet))

data$time <- data$time/365
data$riskscore <- ifelse(risk_score > cut_off, 'high','low')
table(data$risk_score)

fit <- survfit(Surv(time, status)~riskscore, data = data)

p <- ggsurvplot(fit, data = data, 
                conf.int = F, 
                # surv.median.line = "hv", 
                risk.table = TRUE, 
                tables.height= 0.25, 
                cumcensor = T,   
                legend = c(0.83,0.95),
                
                ## P value
                pval = TRUE, 
                pval.size=6, 
                font.pval= c(14, "bold", "black"),
                pval.coord = c(0.00, 0.05),
                
                ## legend
                legend.title = '', # gene_name
                legend.labs=c("High risk", "Low risk"), 
                font.legend= c(14, "plain", "black"), 
                # font.main = c(100, "bold", "black"),
                # xlim = c(0,72), # present narrower X axis, but not affect
                # survival estimates.
                palette=c("red", "blue"),
                font.x = c(14, "plain", "black"),
                font.y = c(14, "plain", "black"), 
                font.tickslab = c(14, "plain", "black"),
                xlab = "Time in years", 
                break.time.by = 6
)  
p

## add HR and CI
res_cox <- coxph(Surv(time, status) ~riskscore, data=data)
HR <- round(summary(res_cox)$conf.int[1],2)
ci_l <- round(summary(res_cox)$conf.int[3],2)
ci_r <- round(summary(res_cox)$conf.int[4],2)

p1 <- p
p1$plot = p1$plot + 
  ggplot2::annotate("text",x = 3.13, y = 0.3, label = paste("HR : ", HR), size = 5) + 
  ggplot2::annotate("text",x = 6.28, y = 0.2,
                    label = paste("(","95%CI : ", ci_l,"-",ci_r,")", sep = ""), size = 5)

# pdf("TCGA_NEW/TCGA_6_os.pdf", width = 4.8, height = 6, onefile = FALSE)
p1
# dev.off()


##############################################################################
# Try to plot ROC curves (DONT used in CNet-Cox v1)
##############################################################################
# Extract PRS gene expression and PRS value -------------------------------
library(dplyr)        
library(tidyr)
library(tidyverse)   

# Data1 = read.table("TCGA_NEW/TCGA_BRCA_clin_546_1080.txt", header = T, check.names = FALSE)
Data1 = read.table("TCGA_NEW/TCGA_BRCA_clin_546_1080_scale.txt", header = T, check.names = FALSE)
# View(Data1[,1:10])
coef_gene <- read.csv("TCGA_NEW/UniMutVariate_markergene.csv", header = T, sep=',')
gene <- as.matrix(coef_gene[,1])
colnames(gene) <- c('gene')
Data2 <- cbind(rownames(Data1), Data1)
colnames(Data2) <- c('gene', colnames(Data1))
# View(Data2[,1:10])

genedata <- merge(gene, Data2, by = "gene")
genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
# View(genedata1[,1:10])
genedata1[1,1]
genedata2 <- rbind(Data1[c(1,2),],genedata1)    # 55
# View(genedata2[,1:10])

# write.table(genedata2,"TCGA_NEW/TCGA_PRS_Expr.txt",quote=F,sep="\t")  
# write.table(genedata2,"TCGA_NEW/TCGA_PRS_Expr_scale.txt",quote=F,sep="\t")



##############################################################################
# Figure Sx in revision
##############################################################################
# Plot KM curves of each gene ---------------------------------------------
library(survival)
library(survminer)

# Convert data: genedata2 should be a matrix; transpose to get samples as rows
genedata2 <- as.matrix(genedata2)
genedata3 <- as.data.frame(t(genedata2))
genedata3$time <- as.numeric(genedata3$time) / 365   # Convert time to years
genedata3$status <- as.numeric(genedata3$status)     # Ensure status is numeric

# Set quantile cutoff
alpha <- 0.75 # (928/1080)
gene_list <- c("EGR1", "IGFBP5", "JUN", "MAFK", "MYC", "TCF7")
plot_list <- list()

# Loop through each gene and generate survival plots
for (gene in gene_list) {
  gene_expr <- as.numeric(genedata3[[gene]])
  cut_off <- as.numeric(quantile(gene_expr, alpha, na.rm = TRUE))
  
  # Group samples by expression level
  genedata3[[paste0(gene, "_group")]] <- ifelse(gene_expr > cut_off, "high", "low")
  
  # Fit survival curve
  fit <- survfit(Surv(time, status) ~ genedata3[[paste0(gene, "_group")]], data = genedata3)
  
  # Create survival plot (without risk table for compact display)
  p <- ggsurvplot(
    fit,
    data = genedata3,
    conf.int = FALSE,
    risk.table = TRUE,      # Hide risk table for multi-panel plot
    cumcensor = FALSE,
    legend.title = gene,
    legend.labs = c("High", "Low"),
    palette = c("red", "blue"),
    font.legend = c(12, "plain", "black"),
    font.x = c(12, "plain", "black"),
    font.y = c(12, "plain", "black"),
    font.tickslab = c(10, "plain", "black"),
    xlab = "Years",
    break.time.by = 6,
    pval = TRUE,
    pval.size = 4,
    font.pval = c(10, "bold", "black")
  )

  plot_list[[gene]] <- p
}

# pdf("TCGA_NEW/TCGA_6_os_gene.pdf", width = 14, height = 10, onefile = FALSE)
arrange_ggsurvplots(plot_list, nrow = 2, ncol = 3)
# dev.off()



# Extract PRS gene expression of P/M GSE147995 data-------------------------------
library(dplyr)       
library(tidyr)
library(tidyverse)   

Data1 = read.csv("DrugData/GSE147995_vst.csv", header = T, row.names= 1, check.names = FALSE, sep = ",")
# View(Data1[,1:10])
# colnames(Data1)
# status <- rep(c("Primary", "Metastases"), 13)
# Data1 <- rbind(status, Data0)
# Data0[1,1]
# rownames(Data1) <- c("status", rownames(Data0))

coef_gene <- read.csv("TCGA_NEW/UniMutVariate_markergene.csv", header = T, sep=',')
gene <- as.matrix(coef_gene[,1])
colnames(gene) <- c('gene')
Data2 <- cbind(rownames(Data1), Data1)
colnames(Data2) <- c('gene', colnames(Data1))
# View(Data2[,1:10])
colnames(Data2)

genedata <- merge(gene, Data2, by = "gene")
genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
# View(genedata1[,1:10])

# write.table(genedata1,"DrugData/GSE147995_PRS_Expr.txt",quote=F,sep="\t")  

## 1/9/11/14/15/17  为ER+ 对应行号为 1、/6、/7、/10、/11、/12
genedata1ER_plus <- genedata1[,c(c(1,11,13,19,21,23), c(1,11,13,19,21,23)+1)]
genedata1ER_mini <- genedata1[,-c(c(1,11,13,19,21,23), c(1,11,13,19,21,23)+1)]
# write.table(genedata1ER_plus,"DrugData/GSE147995_PRS_Expr_plus.txt",quote=F,sep="\t")
# write.table(genedata1ER_mini,"DrugData/GSE147995_PRS_Expr_mini.txt",quote=F,sep="\t")


# TCGA of PRS ----------------------------------------------------------------
# data = read.table("TCGA_NEW/TCGA_PRS_Expr.txt", header = T, check.names = FALSE)
data = read.table("TCGA_NEW/TCGA_PRS_Expr_scale.txt", header = T, check.names = FALSE)
datat <- as.data.frame(t(data))
prs = read.table("TCGA_NEW/TCGA_OS_PRS_gene.txt", header = T, check.names = FALSE)
datat$PRS <- prs$riskscore

data1 <- datat[,c(1,2,9,3,4,5,6,7,8)]
# write.table(data1, file = "TCGA_NEW/TCGA_prs_for_hiplot.txt", quote=F,sep="\t",row.names = F)  


# ROC: single gene and single time point ----------------------------------
## https://mp.weixin.qq.com/s/6TLqj6i-T3Ps9Rs7YKeL3A

library(survivalROC)
## add PRS score
td <- as.data.frame(cbind(colnames(data), data1))
td <- cbind(colnames(data), as.data.frame(apply(td[,2:dim(td)[2]], 2, as.numeric) ))
td$time <- td$time/365
colnames(td) <- c("id", colnames(data1)[-3], colnames(data1)[3])
# td <- datat[,c(2,1,3,4,5,6,7,8)]

## plot one gene as one time point
# par(mar= c(5,5,1,1),cex.lab=1.2,cex.axis= 1.2) #先设置一下图形的边界
# sROC <- survivalROC(Stime=td$time, # 生存时间
#                  status=td$status, # 生存状态
#                  marker = td$EGR1, #选择gene87
#                  predict.time =5,  # 看5年的时间段
#                  method="KM")
# sROC
# plot(sROC$FP, sROC$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="red",
#      xlab="False positive rate", ylab="True positive rate",
#      lwd = 2, cex.main=1.3, cex.lab=1.5, cex.axis=1.2, font=1.2)
# 
# abline(0,1)
# aucText=paste0("5 years"," (AUC=",sprintf("%.3f",sROC$AUC),")")
# legend("bottomright", aucText,
#        lwd=2,bty="n",col=c("red","green","blue"),cex=1.2)


## define function
library(timeROC)
library(survival)
tRocFuction=function(td=null,gene=null){
  tROC <-timeROC(T=td$time, delta = td$status, marker = gene,
                 cause = 1,times = c(1,3,5),ROC=T)
  par(mar= c(4.5,4.5,1,1), cex.lab=1.5, cex.axis= 1.2) #cex.lab=2横纵坐标的label变大，cex.axis=1.5坐标刻度数字变大
  plot(tROC,time=1,col="red",title=F, lwd=3) # lwd 线的粗线
  plot(tROC,time=3,col="green",add=T,title=F,lwd=3)
  plot(tROC,time=5,col="blue",add=T,title=F,lwd=3)
  legend(0.3,0.6, # 改变legend坐标位置
         # c(paste0("AUC at 1 years  ", round(tROC$AUC[1], 2)),
         #   paste0("AUC at 3 years  ", round(tROC$AUC[2], 2)),
         #   paste0("AUC at 5 years  ", round(tROC$AUC[3], 2))),
         c(paste0("AUC", round(tROC$AUC[1], 2)),
           paste0("AUC", round(tROC$AUC[2], 2)),
           paste0("AUC", round(tROC$AUC[3], 2))),
         col=c("red","green","blue"), lwd=2, cex=1.2, bty="n") # 改变legend粗细及大小 lwd=2,cex=1.5
}

## save plot
# pdf(file="gene87.sROC.pdf",width=6,height=5)  
tRocFuction(td=td,gene=td$IGFBP5)  
# dev.off()  


## save on a list
TimeROC <- list()
par(mfrow=c(3,3))
TimeROC[[1]] <- tRocFuction(td=td,gene=td$EGR1)  
TimeROC[[2]] <- tRocFuction(td=td,gene= - td$IGFBP5)  
TimeROC[[3]] <- tRocFuction(td=td,gene= td$JUN)  
TimeROC[[4]] <- tRocFuction(td=td,gene= - td$MAFK) 
TimeROC[[5]] <- tRocFuction(td=td,gene=td$MYC)  
TimeROC[[6]] <- tRocFuction(td=td,gene= - td$TCF7) 
TimeROC[[7]] <- tRocFuction(td=td,gene= -td$PRS)  
par(mfrow=c(1,1))

# # pdf("CompareFeature.pdf",width = 10, height = 5)
# cowplot::plot_grid(plotlist = list(TimeROC[[1]], TimeROC[[2]] , TimeROC[[3]] ,
#                                    TimeROC[[4]], TimeROC[[5]] , TimeROC[[6]],
#                                    TimeROC[[7]]), 
#                    # nrow = 3, rel_widths = c(5.5, 5.5, 5.5, 
#                    #                          5.5,5.5, 5.5, 5.5),
#                    labels = c('(a)','(b)','(c)','(d)','(e)','(f)', '(g)'),
#                    scale = c(0.95)) 
# # dev.off()


# ROC: Six gene at  single time point ----------------------------------
## https://mp.weixin.qq.com/s/eqw9ooOLcaFQUdh06jHQww
rbCol=rainbow(6)
mRocFuction=function(td=null,pt = null,of=null){ #td要使用的数据名，pt要看的时间段
  rbCol=rainbow(6)
  par(mar= c(4.5, 4.5, 1,1),cex.lab=1.2,cex.axis= 1.2)
  sROC=survivalROC(Stime=td$time, status = td$status, marker = td$EGR1,
                   predict.time = pt, #时间改为pt，方便改要看的时间段
                   method="KM")
  plot(sROC$FP, sROC$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rbCol[1], 
       xlab="False positive rate", ylab="True positive rate",
       lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
  abline(0,1)
  aucText=paste0("EGR1"," (AUC=",sprintf("%.3f",sROC$AUC),")")
  j=1
  for(i in colnames(td[,of])){
    sROC=survivalROC(Stime=td$time, status=td$status, marker = td[,i], predict.time = pt, method="KM")
    j=j+1
    lines(sROC$FP, sROC$TP, type="l", xlim=c(0,1), ylim=c(0,1),col = rbCol[j],lwd = 2)
    aucText=c(aucText,paste0(i," (AUC=",sprintf("%.3f",sROC$AUC),")"))
  }
  legend("bottomright", aucText,lwd=2,bty="n",cex = 1.2,col=rbCol) # cex改为1.2倍字体
  legend("topleft", paste(pt,"Years"),bty="n",cex = 1.8) # 加了一个时间点的legend
  
}

## The second 'td' should be replaced with the name of the data you want to view, 
## and 'pt' should be the time period you want to view, which is 3 years.
mRocFuction(td=td, pt=3, of=c("IGFBP5","JUN","MAFK","MYC","TCF7")) 

## save on a list
par(mfrow=c(2,2))
mRocFuction(td=td, pt=1, of=c("IGFBP5","JUN","MAFK","MYC","TCF7")) 
mRocFuction(td=td, pt=3, of=c("IGFBP5","JUN","MAFK","MYC","TCF7")) 
mRocFuction(td=td, pt=5, of=c("IGFBP5","JUN","MAFK","MYC","TCF7")) 
par(mfrow=c(1,1))



##############################################################################
# For revision v1  2026.01.16
##############################################################################
# Validation of oncotype21 and mamamprint70 -------------------------------

## Load data
DX21gene  <- read.csv("Prior_Infor/OncotypeDX21gene.csv", header = T, sep=',')
DX21gene <- DX21gene[,1]
mamaprint <- read.csv("Prior_Infor/70_mama.csv", header = T, sep=',')
mamaprint <- mamaprint[,1]
cnnetcox <- as.matrix(read.csv("Result/Result3/marker3.csv", header = T, sep=','))
cnnetcox <- cnnetcox[,1]
prs_gene <- read.csv("TCGA_NEW/UniMutVariate_markergene.csv", header = T, sep=',')
prs_gene <- prs_gene[,1]

## Venn plot
# install.packages("ggVennDiagram")
library(ggVennDiagram)

## All
gene_list <- list(
  `Oncotype DX 21` = DX21gene,
  `MammaPrint 70` = mamaprint,
  `CNet-Cox` = cnnetcox
)

p <- ggVennDiagram(gene_list, 
                   label_alpha = 0, 
                   label_size = 5,
                   edge_size = 0.8,
                   set_color = c("#377eb8", "#4daf4a", "#e41a1c"),
                   category.names = names(gene_list)) +
  scale_fill_gradient(low = "#FFFFFF", high = "#FFB6B6") +   
  theme_void(base_size = 15) +
  theme(legend.position = "none", 
        text = element_text(family = "Helvetica")) 
p
# ggsave("TCGA_NEW/venn_compare.pdf", p, width = 4, height = 4)


## PRS
gene_list <- list(
  `Oncotype DX 21` = DX21gene,
  `MammaPrint 70` = mamaprint,
  `PRS score` = prs_gene
)

p <- ggVennDiagram(gene_list, 
                   label_alpha = 0, 
                   label_size = 5,
                   edge_size = 0.8,
                   set_color = c("#377eb8", "#4daf4a", "#e41a1c"),
                   category.names = names(gene_list)) +
  scale_fill_gradient(low = "#FFFFFF", high = "#FFB6B6") +   
  theme_void(base_size = 15) +
  theme(legend.position = "none", 
        text = element_text(family = "Helvetica"))  
p
# ggsave("TCGA_NEW/venn_compare_prs.pdf", p, width = 4, height = 4)


intersect(DX21gene, mamaprint)
intersect(mamaprint, cnnetcox)     # "BBC3"   "EGLN1"  "IGFBP5" "RASSF7" "RECQL5"
intersect(mamaprint, prs_gene)     # IGFBP5
Reduce(intersect, list(DX21gene, mamaprint, cnnetcox)) 

