## 2026.01.19 LLY create, copy from R/BRCA/GSE42568_expr.R
## options(encoding = "UTF-8")

##############################################################################
# ER validation
##############################################################################


## clear
rm(list = ls())

## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- '/home/lly/R/'
setwd(paste(path, 'CNetCox/Data/RevisedData/', sep=''))
# setwd("/home/lingyu/ssd/Python/CNetCox/CNetCox_local/Data/RevisedData")

library(BiocGenerics)
library(parallel)
library(Biobase)
library(dplyr)        
library(tidyr)
library(tidyverse)    
library(fdrtool)      
library(data.table)   
library(fdrtool)      
library(data.table)   

## unzip the download data (1.85GB)
# (base) lilingyu@LiLingyudeMacBook-Pro RevisedData % 
# gunzip GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz

# Load gene expr data -----------------------------------------------------
exprSet <- read.table(
  "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz",
  header = TRUE,
  sep = ",",
  fill = TRUE,
  strip.white = TRUE
)
# View(exprSet[,1:2])
dim(exprSet)    #  30865  3410

## replace X using ID_REF, keep the same with other GEO dataset
exprSet <- exprSet %>%rename(GeneID = X)
exprSet$GeneID <- as.character(exprSet$GeneID)

## pro-processing
expset4 <- exprSet %>%
  mutate(rowIQR = apply(.[, -1], 1, IQR)) %>%    # -1, delete 'GeneID'
  arrange(desc(rowIQR)) %>%
  distinct(GeneID, .keep_all = TRUE) %>%
  select(-rowIQR) %>%
  column_to_rownames("GeneID")

View(expset4[1:10,1:10])
dim(expset4)   # 30865  3409
# write.table(expset4,"GSE96058_expr.txt",quote=F,sep="\t")


#############################################################################
#################### Load ER positive clinical information ##################
library(openxlsx)
# label <- read.xlsx("GSE96058-GPL18573_series_matrix.xlsx", sheet='Sheet1')
label <- read.csv("GSE96058-GPL18573_series_matrix.csv", header = T, sep=',')
dim(label)    # 340  24
# View(label)  
colnames(label)
## Selected ER=1
label_er <- label[label$er.status=='er status: 1',]
dim(label_er) # 289  24

## sample id
sample_er <- as.matrix(label_er$Sample_title)
colnames(sample_er) <- c("sample_er") 
# View(sample_er)

# e_os -------------------------------------------------------------------
e_os1 <- label_er$overall.survival.event 
e_os2 <- as.numeric(sub(".*: *", "", e_os1))  # my add
e_os2 <- as.matrix(e_os2) 
colnames(e_os2) <- "e_os" 
dim(e_os2)    # 289   1
# View(e_os2)

# t_os -------------------------------------------------------------------

t_os1 <- label_her2$overall.survival.days
t_os2 <- as.numeric(sub(".*: *", "", t_os1))  # my add
t_os2 <- as.matrix(t_os2)
colnames(t_os2) <- "t_os"
dim(t_os2)    # 289   1
# View(t_os2)

# os ---------------------------------------------------------------------
os <- cbind(sample_er, t_os2, e_os2)
# View(os)


# Add clinical information into expr data ---------------------------------
data <- read.table("GSE96058_expr.txt",header=T,sep='\t',
                   fill=TRUE,strip.white = T)
dim(data)    # 30865  3409
# View(data[,1:10])
# View(expset4[,1:10])
mean(data[,1])    # -0.2087899
mean(as.numeric(data[1,]))    # 9.432776


## selected samples in ER  [289] from 'data' [3409]
## F7 and F8 are removed
data_er <- data[, sample_er[sample_er %in% colnames(data)], drop = FALSE]
dim(data_er)    # 30865   289
dim(sample_er)  # 289   1

## use scale
expset6 <- t(scale(t(data_er)))
dim(expset6)    # 30865  3392
# View(t(expset6)[,1:10])

# ## another scale  
# expset6 <- scale(data)
# dim(expset6)    # 21835   104 
# # View(expset6[,1:10])

# data_os ---------------------------------------------------------------------
# View(as.matrix(t(os[,c(2,3)])))
data_os = rbind(as.matrix(t(os[,c(2,3)])), as.matrix(expset6))
colnames(data_os) <- c(colnames(expset6))
View(data_os[,1:10])
dim(data_os)    # 30867   289
# write.table(data_os,"GSE96058_scale_outcome_os.txt",quote=F,sep="\t")


##############################################################################
# ER validation, calculate PRS
# feature_select_GEO.R
##############################################################################
Data1 <- read.table("GSE96058_scale_outcome_os.txt", header=TRUE, 
                    sep="\t", check.names=FALSE)
dim(Data1)     # 30867   289

coef_gene <- read.csv(paste(path, 'CNetCox/Data/TCGA_NEW/UniMutVariate_markergene.csv', sep=''), 
                      header = T, sep=',')
gene <- as.matrix(coef_gene[,1])
dim(gene)      # 6   1

colnames(gene) <- c('gene')
Data2 <- cbind(rownames(Data1), Data1)
colnames(Data2) <- c('gene', colnames(Data1))

genedata <- merge(gene, Data2, by = "gene")
genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
genedata2 <- rbind(Data1[c(1,2),],genedata1)   

# write.table(genedata2, 'GSE96058_os.txt', quote = F)


##############################################################################
# Calculating riskscore using separate2GroupsCox
##############################################################################
# Extract PRS data and label -------------------------------------------
data = read.table("GSE96058_os.txt", header = T, check.names = FALSE)
datat <- as.data.frame(t(data))
Data <- as.matrix(data)
rownames(Data)

x <- model.matrix(e_os ~., datat)[,-c(1,2)]  
x_hat <- data.frame(x)
y <- as.matrix(datat[,c(1,2)])  
colnames(y) <- c("time", "status")
y_hat <- data.frame(y)

coef_gene <- read.csv(paste(path, 'CNetCox/Data/TCGA_NEW/UniMutVariate_markergene.csv', sep=''), 
         header = T, sep=',')

coef <- as.matrix(coef_gene[,2])
rownames(coef) <- coef_gene[,1]


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
# write.table(p_index, file = "GSE96058_OS_PRS_gene_scale.txt", quote = F, row.names = F, sep="\t")


##############################################################################
# ER validation, calculate PRS, plot KM curve
# feature_select_GEO.R, CI_repeat.R [Fig. 4B]
##############################################################################
library(survminer)
library(survival)
library(ggplot2)

data <- read.table("GSE96058_OS_PRS_gene_scale.txt", header = T, check.names = FALSE)
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
data$riskscore <- ifelse(risk_score > cut_off, 'high', 'low')
# table(data$risk_score)

fit <- survfit(Surv(time, status)~riskscore, data = data)

p <- ggsurvplot(fit, data = data, 
                conf.int = F, 
                # surv.median.line = "hv", 
                risk.table = TRUE, 
                tables.height= 0.25, 
                cumcensor = T,   
                legend = c(0.83,0.25),    # legend location
                # legend = "bottom",
                
                ## P value
                pval = TRUE, 
                pval.size=6, 
                font.pval= c(14, "bold", "black"),
                pval.coord = c(0.00, 0.05),
                
                ## legend
                legend.title = 'ER positive', # gene_name
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
                break.time.by = 1
)  
p

## add HR and CI
res_cox <- coxph(Surv(time, status) ~riskscore, data=data)
HR <- round(summary(res_cox)$conf.int[1],2)
ci_l <- round(summary(res_cox)$conf.int[3],2)
ci_r <- round(summary(res_cox)$conf.int[4],2)

p1 <- p
p1$plot = p1$plot + 
  ggplot2::annotate("text",x = 1.13, y = 0.3, label = paste("HR : ", HR), size = 5) + 
  ggplot2::annotate("text",x = 1.58, y = 0.2,
                    label = paste("(","95%CI : ", ci_l,"-",ci_r,")", sep = ""), size = 5)

# pdf("GSE96058_6_os.pdf", width = 4.8, height = 6, onefile = FALSE)
p1
# dev.off()

#############################################################################
#############################################################################


#############################################################################
#################### Load her2 positive clinical information ################
library(openxlsx)
label <- read.csv("GSE96058-GPL18573_series_matrix.csv", header = T, sep=',')
dim(label)    # 340  24
# View(label)  
colnames(label)
## Selected ER=1
label_her2 <- label[label$her2.status=='her2 status: 1',]
dim(label_her2) # 46 24

## sample id
sample_her2 <- as.matrix(label_her2$Sample_title)
colnames(sample_her2) <- c("sample_her2") 
# View(sample_her2)

# e_os -------------------------------------------------------------------
e_os1 <- label_her2$overall.survival.event 
e_os2 <- as.numeric(sub(".*: *", "", e_os1))  # my add
e_os2 <- as.matrix(e_os2) 
colnames(e_os2) <- "e_os" 
dim(e_os2)    # 46   1
# View(e_os2)

# t_os -------------------------------------------------------------------

t_os1 <- label_her2$overall.survival.days
t_os2 <- as.numeric(sub(".*: *", "", t_os1))  # my add
t_os2 <- as.matrix(t_os2)
colnames(t_os2) <- "t_os"
dim(t_os2)    # 46   1
# View(t_os2)

# os ---------------------------------------------------------------------
os <- cbind(sample_her2, t_os2, e_os2)
# View(os)


# Add clinical information into expr data ---------------------------------
data <- read.table("GSE96058_expr.txt",header=T,sep='\t',
                   fill=TRUE,strip.white = T)
dim(data)    # 30865  3400
# View(data[,1:10])
# View(expset4[,1:10])
mean(data[,1])    # -0.2087899
mean(as.numeric(data[1,]))    # 9.432776


## selected samples in HER2  [46] from 'data' [3400]
## F7 and F8 are removed
data_her2 <- data[, sample_her2[sample_her2 %in% colnames(data)], drop = FALSE]
dim(data_her2)    # 30865   46
dim(sample_her2)  # 46   1

## use scale
expset6 <- t(scale(t(data_her2)))
dim(expset6)    # 30865  46
# View(t(expset6)[,1:10])

# ## another scale  
# expset6 <- scale(data)
# dim(expset6)    # 21835   104 
# # View(expset6[,1:10])

# data_os ---------------------------------------------------------------------
# View(as.matrix(t(os[,c(2,3)])))
data_os = rbind(as.matrix(t(os[,c(2,3)])), as.matrix(expset6))
colnames(data_os) <- c(colnames(expset6))
# View(data_os[,1:10])
dim(data_os)    # 30867   46
# write.table(data_os,"GSE96058_scale_outcome_os_her2.txt",quote=F,sep="\t")


##############################################################################
# ER validation, calculate PRS
# feature_select_GEO.R
##############################################################################

Data1 <- read.table("GSE96058_scale_outcome_os_her2.txt", header=TRUE, 
                    sep="\t", check.names=FALSE)
dim(Data1)     # 30867   287

coef_gene <- read.csv(paste(path, 'CNetCox/Data/TCGA_NEW/UniMutVariate_markergene.csv', sep=''), 
                      header = T, sep=',')
gene <- as.matrix(coef_gene[,1])
dim(gene)      # 6   1

colnames(gene) <- c('gene')
Data2 <- cbind(rownames(Data1), Data1)
colnames(Data2) <- c('gene', colnames(Data1))

genedata <- merge(gene, Data2, by = "gene")
genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
genedata2 <- rbind(Data1[c(1,2),],genedata1)   

# write.table(genedata2, 'GSE96058_os_her2.txt', quote = F)



##############################################################################
# Calculating riskscore using separate2GroupsCox
##############################################################################
# Extract PRS data and label -------------------------------------------
data = read.table("GSE96058_os_her2.txt", header = T, check.names = FALSE)
datat <- as.data.frame(t(data))
Data <- as.matrix(data)
rownames(Data)

x <- model.matrix(e_os ~., datat)[,-c(1,2)]  
x_hat <- data.frame(x)
y <- as.matrix(datat[,c(1,2)])  
colnames(y) <- c("time", "status")
y_hat <- data.frame(y)

coef_gene <- read.csv(paste(path, 'CNetCox/Data/TCGA_NEW/UniMutVariate_markergene.csv', sep=''), 
                      header = T, sep=',')

coef <- as.matrix(coef_gene[,2])
rownames(coef) <- coef_gene[,1]


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
# write.table(p_index, file = "GSE96058_OS_PRS_gene_scale_her2.txt", quote = F, row.names = F, sep="\t")


##############################################################################
# ER validation, calculate PRS, plot KM curve
# feature_select_GEO.R, CI_repeat.R [Fig. 4B]
##############################################################################
library(survminer)
library(survival)
library(ggplot2)

data <- read.table("GSE96058_OS_PRS_gene_scale_her2.txt", header = T, check.names = FALSE)
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
data$riskscore <- ifelse(risk_score > cut_off, 'high', 'low')
table(data$risk_score)

fit <- survfit(Surv(time, status)~riskscore, data = data)

p <- ggsurvplot(fit, data = data, 
                conf.int = F, 
                # surv.median.line = "hv", 
                risk.table = TRUE, 
                tables.height= 0.25, 
                cumcensor = T,   
                legend = c(0.83,0.25),    # legend location
                # legend = "bottom",
                
                ## P value
                pval = TRUE, 
                pval.size=6, 
                font.pval= c(14, "bold", "black"),
                pval.coord = c(0.00, 0.05),
                
                ## legend
                legend.title = 'ER positive', # gene_name
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
                break.time.by = 1
)  
p

## add HR and CI
res_cox <- coxph(Surv(time, status) ~riskscore, data=data)
HR <- round(summary(res_cox)$conf.int[1],2)
ci_l <- round(summary(res_cox)$conf.int[3],2)
ci_r <- round(summary(res_cox)$conf.int[4],2)

p1 <- p
p1$plot = p1$plot + 
  ggplot2::annotate("text",x = 1.13, y = 0.3, label = paste("HR : ", HR), size = 5) + 
  ggplot2::annotate("text",x = 1.58, y = 0.2,
                    label = paste("(","95%CI : ", ci_l,"-",ci_r,")", sep = ""), size = 5)

# pdf("GSE96058_6_os.pdf", width = 4.8, height = 6, onefile = FALSE)
p1
# dev.off()


#############################################################################
#############################################################################


#############################################################################
#################### Load pgr positive clinical information ################
library(openxlsx)
label <- read.csv("GSE96058-GPL18573_series_matrix.csv", header = T, sep=',')
dim(label)    # 340  24
# View(label)  
colnames(label)
## Selected ER=1
label_pgr <- label[label$pgr.status=='pgr status: 1',]
dim(label_pgr) # 267 24

## sample id
sample_pgr <- as.matrix(label_pgr$Sample_title)
colnames(sample_pgr) <- c("sample_pgr") 
# View(sample_pgr)

# e_os -------------------------------------------------------------------
e_os1 <- label_pgr$overall.survival.event 
e_os2 <- as.numeric(sub(".*: *", "", e_os1))  # my add
e_os2 <- as.matrix(e_os2) 
colnames(e_os2) <- "e_os" 
dim(e_os2)    # 267   1
# View(e_os2)

# t_os -------------------------------------------------------------------

t_os1 <- label_pgr$overall.survival.days
t_os2 <- as.numeric(sub(".*: *", "", t_os1))  # my add
t_os2 <- as.matrix(t_os2)
colnames(t_os2) <- "t_os"
dim(t_os2)    # 267   1
# View(t_os2)

# os ---------------------------------------------------------------------
os <- cbind(sample_pgr, t_os2, e_os2)
# View(os)


# Add clinical information into expr data ---------------------------------
data <- read.table("GSE96058_expr.txt",header=T,sep='\t',
                   fill=TRUE,strip.white = T)
dim(data)    # 30865  3400
# View(data[,1:10])
# View(expset4[,1:10])
mean(data[,1])    # -0.2087899
mean(as.numeric(data[1,]))    # 9.432776


## selected samples in pgr  [46] from 'data' [3400]
## F7 and F8 are removed
data_pgr <- data[, sample_pgr[sample_pgr %in% colnames(data)], drop = FALSE]
dim(data_pgr)    # 30865   46
dim(sample_pgr)  # 267   1

## use scale
expset6 <- t(scale(t(data_pgr)))
dim(expset6)    # 30865  267
# View(t(expset6)[,1:10])

# ## another scale  
# expset6 <- scale(data)
# dim(expset6)    # 21835   104 
# # View(expset6[,1:10])

# data_os ---------------------------------------------------------------------
# View(as.matrix(t(os[,c(2,3)])))
data_os = rbind(as.matrix(t(os[,c(2,3)])), as.matrix(expset6))
colnames(data_os) <- c(colnames(expset6))
# View(data_os[,1:10])
dim(data_os)    # 30867   267
# write.table(data_os,"GSE96058_scale_outcome_os_pgr.txt",quote=F,sep="\t")


##############################################################################
# ER validation, calculate PRS
# feature_select_GEO.R
##############################################################################

Data1 <- read.table("GSE96058_scale_outcome_os_pgr.txt", header=TRUE, 
                    sep="\t", check.names=FALSE)
dim(Data1)     # 30867   267

coef_gene <- read.csv(paste(path, 'CNetCox/Data/TCGA_NEW/UniMutVariate_markergene.csv', sep=''), 
                      header = T, sep=',')
gene <- as.matrix(coef_gene[,1])
dim(gene)      # 6   1

colnames(gene) <- c('gene')
Data2 <- cbind(rownames(Data1), Data1)
colnames(Data2) <- c('gene', colnames(Data1))

genedata <- merge(gene, Data2, by = "gene")
genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
genedata2 <- rbind(Data1[c(1,2),],genedata1)   

# write.table(genedata2, 'GSE96058_os_pgr.txt', quote = F)


##############################################################################
# Calculating riskscore using separate2GroupsCox
##############################################################################
# Extract PRS data and label -------------------------------------------
data = read.table("GSE96058_os_pgr.txt", header = T, check.names = FALSE)
datat <- as.data.frame(t(data))
Data <- as.matrix(data)
rownames(Data)

x <- model.matrix(e_os ~., datat)[,-c(1,2)]  
x_hat <- data.frame(x)
y <- as.matrix(datat[,c(1,2)])  
colnames(y) <- c("time", "status")
y_hat <- data.frame(y)

coef_gene <- read.csv(paste(path, 'CNetCox/Data/TCGA_NEW/UniMutVariate_markergene.csv', sep=''), 
                      header = T, sep=',')

coef <- as.matrix(coef_gene[,2])
rownames(coef) <- coef_gene[,1]


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
# write.table(p_index, file = "GSE96058_OS_PRS_gene_scale_pgr.txt", quote = F, row.names = F, sep="\t")


##############################################################################
# ER validation, calculate PRS, plot KM curve
# feature_select_GEO.R, CI_repeat.R [Fig. 4B]
##############################################################################
library(survminer)
library(survival)
library(ggplot2)

data <- read.table("GSE96058_OS_PRS_gene_scale_pgr.txt", header = T, check.names = FALSE)
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
data$riskscore <- ifelse(risk_score > cut_off, 'high', 'low')
table(data$risk_score)

fit <- survfit(Surv(time, status)~riskscore, data = data)

p <- ggsurvplot(fit, data = data, 
                conf.int = F, 
                # surv.median.line = "hv", 
                risk.table = TRUE, 
                tables.height= 0.25, 
                cumcensor = T,   
                legend = c(0.83,0.25),    # legend location
                # legend = "bottom",
                
                ## P value
                pval = TRUE, 
                pval.size=6, 
                font.pval= c(14, "bold", "black"),
                pval.coord = c(0.00, 0.05),
                
                ## legend
                legend.title = 'ER positive', # gene_name
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
                break.time.by = 1
)  
p

## add HR and CI
res_cox <- coxph(Surv(time, status) ~riskscore, data=data)
HR <- round(summary(res_cox)$conf.int[1],2)
ci_l <- round(summary(res_cox)$conf.int[3],2)
ci_r <- round(summary(res_cox)$conf.int[4],2)

p1 <- p
p1$plot = p1$plot + 
  ggplot2::annotate("text",x = 1.13, y = 0.3, label = paste("HR : ", HR), size = 5) + 
  ggplot2::annotate("text",x = 1.58, y = 0.2,
                    label = paste("(","95%CI : ", ci_l,"-",ci_r,")", sep = ""), size = 5)

# pdf("GSE96058_6_os.pdf", width = 4.8, height = 6, onefile = FALSE)
p1
# dev.off()



#############################################################################
#################### Load TNBC clinical information ################
library(openxlsx)
label <- read.csv("GSE96058-GPL18573_series_matrix.csv", header = T, sep=',')
dim(label)    # 340  24
# View(label)  
colnames(label)
## Selected ER=1
label_tnbc <- label[
  label$pgr.status == 'pgr status: 0' &
    label$er.status == 'er status: 0' &
    label$her2.status == 'her2 status: 0',
]

dim(label_tnbc) # 18 24

## sample id
sample_tnbc <- as.matrix(label_tnbc$Sample_title)
colnames(sample_tnbc) <- c("sample_tnbc") 
# View(sample_tnbc)

# e_os -------------------------------------------------------------------
e_os1 <- label_tnbc$overall.survival.event 
e_os2 <- as.numeric(sub(".*: *", "", e_os1))  # my add
e_os2 <- as.matrix(e_os2) 
colnames(e_os2) <- "e_os" 
dim(e_os2)    # 18   1
# View(e_os2)

# t_os -------------------------------------------------------------------

t_os1 <- label_tnbc$overall.survival.days
t_os2 <- as.numeric(sub(".*: *", "", t_os1))  # my add
t_os2 <- as.matrix(t_os2)
colnames(t_os2) <- "t_os"
dim(t_os2)    # 18   1
# View(t_os2)

# os ---------------------------------------------------------------------
os <- cbind(sample_tnbc, t_os2, e_os2)
# View(os)


# Add clinical information into expr data ---------------------------------
data <- read.table("GSE96058_expr.txt",header=T,sep='\t',
                   fill=TRUE,strip.white = T)
dim(data)    # 30865  3400
# View(data[,1:10])
# View(expset4[,1:10])
mean(data[,1])    # -0.2087899
mean(as.numeric(data[1,]))    # 9.432776


## selected samples in tnbc  [18] from 'data' [3400]
## F7 and F8 are removed
data_tnbc <- data[, sample_tnbc[sample_tnbc %in% colnames(data)], drop = FALSE]
dim(data_tnbc)    # 30865   18
dim(sample_tnbc)  # 18   1

## use scale
expset6 <- t(scale(t(data_tnbc)))
dim(expset6)    # 30865  18
# View(t(expset6)[,1:10])

# ## another scale  
# expset6 <- scale(data)
# dim(expset6)    # 21835   104 
# # View(expset6[,1:10])

# data_os ---------------------------------------------------------------------
# View(as.matrix(t(os[,c(2,3)])))
data_os = rbind(as.matrix(t(os[,c(2,3)])), as.matrix(expset6))
colnames(data_os) <- c(colnames(expset6))
# View(data_os[,1:10])
dim(data_os)    # 30867   18
# write.table(data_os,"GSE96058_scale_outcome_os_tnbc.txt",quote=F,sep="\t")


##############################################################################
# ER validation, calculate PRS
# feature_select_GEO.R
##############################################################################

Data1 <- read.table("GSE96058_scale_outcome_os_tnbc.txt", header=TRUE,
                    sep="\t", check.names=FALSE)
dim(Data1)     # 30867   18

coef_gene <- read.csv(paste(path, 'CNetCox/Data/TCGA_NEW/UniMutVariate_markergene.csv', sep=''), 
                      header = T, sep=',')
gene <- as.matrix(coef_gene[,1])
dim(gene)      # 6   1

colnames(gene) <- c('gene')
Data2 <- cbind(rownames(Data1), Data1)
colnames(Data2) <- c('gene', colnames(Data1))

genedata <- merge(gene, Data2, by = "gene")
genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
genedata2 <- rbind(Data1[c(1,2),],genedata1)   

# write.table(genedata2, 'GSE96058_os_tnbc.txt', quote = F)


##############################################################################
# Calculating riskscore using separate2GroupsCox
##############################################################################
# Extract PRS data and label -------------------------------------------
data = read.table("GSE96058_os_tnbc.txt", header = T, check.names = FALSE)
datat <- as.data.frame(t(data))
Data <- as.matrix(data)
rownames(Data)

x <- model.matrix(e_os ~., datat)[,-c(1,2)]  
x_hat <- data.frame(x)
y <- as.matrix(datat[,c(1,2)])  
colnames(y) <- c("time", "status")
y_hat <- data.frame(y)

coef_gene <- read.csv(paste(path, 'CNetCox/Data/TCGA_NEW/UniMutVariate_markergene.csv', sep=''), 
                      header = T, sep=',')

coef <- as.matrix(coef_gene[,2])
rownames(coef) <- coef_gene[,1]


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
# write.table(p_index, file = "GSE96058_OS_PRS_gene_scale_tnbc.txt", quote = F, row.names = F, sep="\t")


##############################################################################
# ER validation, calculate PRS, plot KM curve
# feature_select_GEO.R, CI_repeat.R [Fig. 4B]
##############################################################################
library(survminer)
library(survival)
library(ggplot2)

data <- read.table("GSE96058_OS_PRS_gene_scale_tnbc.txt", header = T, check.names = FALSE)
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
data$riskscore <- ifelse(risk_score > cut_off, 'high', 'low')
table(data$risk_score)

fit <- survfit(Surv(time, status)~riskscore, data = data)

p <- ggsurvplot(fit, data = data, 
                conf.int = F, 
                # surv.median.line = "hv", 
                risk.table = TRUE, 
                tables.height= 0.25, 
                cumcensor = T,   
                legend = c(0.83,0.25),    # legend location
                # legend = "bottom",
                
                ## P value
                pval = TRUE, 
                pval.size=6, 
                font.pval= c(14, "bold", "black"),
                pval.coord = c(0.00, 0.05),
                
                ## legend
                legend.title = 'ER positive', # gene_name
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
                break.time.by = 1
)  
p

## add HR and CI
res_cox <- coxph(Surv(time, status) ~riskscore, data=data)
HR <- round(summary(res_cox)$conf.int[1],2)
ci_l <- round(summary(res_cox)$conf.int[3],2)
ci_r <- round(summary(res_cox)$conf.int[4],2)

p1 <- p
p1$plot = p1$plot + 
  ggplot2::annotate("text",x = 1.13, y = 0.3, label = paste("HR : ", HR), size = 5) + 
  ggplot2::annotate("text",x = 1.58, y = 0.2,
                    label = paste("(","95%CI : ", ci_l,"-",ci_r,")", sep = ""), size = 5)

# pdf("GSE96058_6_os_tnbc.pdf", width = 4.8, height = 6, onefile = FALSE)
p1
# dev.off()




#############################################################################
#################### Visualize clinical information #########################
library(openxlsx)
label <- read.csv("GSE96058-GPL18573_series_matrix.csv", header = T, sep=',')
dim(label)    # 340  24
# View(label)  
colnames(label)
## Selected ER=1
# label_tnbc <- label[
#   label$pgr.status == 'pgr status: 0' &
#     label$er.status == 'er status: 0' &
#     label$her2.status == 'her2 status: 0',
# ]



# Bar plot ----------------------------------------------------------------
## Create a statistics box
status_stats <- data.frame(
  status = c("er.status", "pgr.status", "her2.status"),
  zero = c(
    sum(label$er.status == 'er status: 0', na.rm=TRUE),
    sum(label$pgr.status == 'pgr status: 0', na.rm=TRUE),
    sum(label$her2.status == 'her2 status: 0', na.rm=TRUE)
  ),
  one = c(
    sum(label$er.status == 'er status: 1', na.rm=TRUE),
    sum(label$pgr.status == 'pgr status: 1', na.rm=TRUE),
    sum(label$her2.status == 'her2 status: 1', na.rm=TRUE)
  ),
  NA_count = c(
    sum(label$er.status == 'er status: NA', na.rm=TRUE),
    sum(label$pgr.status == 'pgr status: NA', na.rm=TRUE),
    sum(label$her2.status == 'her2 status: NA', na.rm=TRUE)
  )
)
print(status_stats)


## Replace group labels
library(reshape2)
df_long <- melt(status_stats, id.vars = "status", 
                variable.name = "value", value.name = "count")
df_long$value <- factor(df_long$value,
                        levels = c("zero", "one", "NA_count"),
                        labels = c("Status: -", "Status: +", "Status: NA"))

library(ggplot2)
ggplot(df_long, aes(x=status, y=count, fill=value)) +
  geom_bar(stat="identity", position="dodge") +
  labs(x="Status", y="Count", fill="Group", 
       title="ER/PGR/HER2 Status Distribution") +
  theme_minimal() +
  theme(text = element_text(size=14))



# stacked bar -------------------------------------------------------------
library(reshape2)
library(scales)
# 转长格式
df_long <- melt(status_stats, id.vars = "status", 
                variable.name = "value", value.name = "count")
df_long$value <- factor(df_long$value,
                        levels = c("zero", "one", "NA_count"),
                        labels = c("Negative", "Positive", "Missing"))

# 计算比例
df_long <- merge(
  df_long,
  aggregate(count ~ status, df_long, sum),
  by = "status",
  suffixes = c("", "_total")
)
df_long$prop <- df_long$count / df_long$count_total

# 画堆叠条形图
p2 <- ggplot(df_long, aes(x=status, y=prop, fill=value)) +
  geom_bar(stat="identity", position="stack", color="black", width=0.7) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("#4E79A7", "#F28E2B", "#E15759")) + # 蓝橙红，典型期刊配色
  labs(
    x = NULL,
    y = "Proportion (%)",
    fill = NULL,
    # title = "Distribution of ER, PR, and HER2 Status"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size=13, color="black"),
    axis.title = element_text(size=15),
    plot.title = element_text(size=17,hjust=0.5),
    legend.position = "right",
    legend.text = element_text(size=13)
  ) +
  geom_text(
    aes(label = ifelse(prop > 0.03, paste0(round(prop*100), "%"), "")),
    position = position_stack(vjust = 0.5), color="white", size=4.5
  )


# pdf("stack_bar_revision.pdf", width = 5.5, height = 6, onefile = FALSE)
p2
# dev.off()

