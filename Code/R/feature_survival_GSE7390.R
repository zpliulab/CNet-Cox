rm(list = ls())

## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- '/home/lly/R/'

setwd(paste(path, 'BRCA/Bulk_RNA-seq_data/GEO/GSE9195', sep=''))
Data = read.table("GSE9195_scale_outcome_dmfs.txt", header = T, check.names = FALSE)
# Data = read.table("GSE9195_scale_outcome_rfs.txt", header = T, check.names = FALSE)

dim(Data)    # 21837    88
# View(Data[,1:10])
class(Data)


###################################################################################
y <- t(Data[c(1,2),])
colnames(y) <- c("time", "status")
dim(y)    # 565   2
# View(y)

y1 <- cbind(as.numeric(y[,1]), as.data.frame(rownames(y)), as.numeric(y[,2]))
colnames(y1) <- c("time", "sample", "status")
dim(y1)    # 252   3
# View(y1)

x <- t(Data[-c(1,2),])
dim(x)    # 252 21835
# View(x[,1:10])

x1 <- cbind(rownames(x),x)
dim(x1)    # 252 21836
colnames(x1) <- c("sample", colnames(x))
# View(x1[,1:10])
class(x1)
x2 <- as.data.frame(x1)
class(x2)


library(dplyr)
exprSet <- x2 %>% 
  # then make it a tibble (nice printing while debugging)
  as_tibble() %>% 
  # then get just a few genes
  dplyr::select(sample, CEL, FLRT1, FUZ, NFE2, PDP1, PHACTR1, SARM1, SGCE) %>% 
  # then join back to clinical data
  inner_join(y1,by="sample") 

dim(exprSet)    # 88  11
# View(exprSet)   
class(exprSet)

# Survival -----------------------------------------------------------------

library(survival)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(survminer)


splots <- list()

# View(exprSet$CEL)
# n1 <- as.numeric(exprSet$CEL)
# median(n1)

group <- ifelse(exprSet$CEL>median(as.numeric(exprSet$CEL)),'high','low')
sfit <- survfit(Surv(time, status)~group, data=exprSet)
splots[[1]] <- ggsurvplot(sfit, title="CEL survival",  
                          conf.int=F, pval=TRUE)

group <- ifelse(exprSet$FLRT1>median(as.numeric(exprSet$FLRT1)),'high','low')
sfit <- survfit(Surv(time, status)~group, data=exprSet)
splots[[2]] <- ggsurvplot(sfit, title="FLRT1 survival",  
                          conf.int=F, pval=TRUE)

group <- ifelse(exprSet$FUZ>median(as.numeric(exprSet$FUZ)),'high','low')
sfit <- survfit(Surv(time, status)~group, data=exprSet)
splots[[3]] <- ggsurvplot(sfit, title="FUZ survival",  
                          conf.int=F, pval=TRUE)

group <- ifelse(exprSet$NFE2>median(as.numeric(exprSet$NFE2)),'high','low')
sfit <- survfit(Surv(time, status)~group, data=exprSet)
splots[[4]] <- ggsurvplot(sfit, title="NFE2 survival",  
                          conf.int=F, pval=TRUE)

group <- ifelse(exprSet$PDP1>median(as.numeric(exprSet$PDP1)),'high','low')
sfit <- survfit(Surv(time, status)~group, data=exprSet)
splots[[5]] <- ggsurvplot(sfit, title="PDP1 survival",  
                          conf.int=F, pval=TRUE)

group <- ifelse(exprSet$PHACTR1>median(as.numeric(exprSet$PHACTR1)),'high','low')
sfit <- survfit(Surv(time, status)~group, data=exprSet)
splots[[6]] <- ggsurvplot(sfit, title="PHACTR1 survival", 
                          conf.int=F, pval=TRUE)

group <- ifelse(exprSet$SARM1>median(as.numeric(exprSet$SARM1)),'high','low')
sfit <- survfit(Surv(time, status)~group, data=exprSet)
splots[[7]] <- ggsurvplot(sfit, title="SARM1 survival",  
                          conf.int=F, pval=TRUE)

group <- ifelse(exprSet$SGCE>median(as.numeric(exprSet$SGCE)),'high','low')
sfit <- survfit(Surv(time, status)~group, data=exprSet)
splots[[8]] <- ggsurvplot(sfit, title="SGCE survival",  
                          conf.int=F, pval=TRUE)


arrange_ggsurvplots(splots, print = TRUE,
                    nrow = 4, ncol = 3, risk.table.height = 0.4)


res <- arrange_ggsurvplots(splots, print = TRUE,  
                           nrow = 4, ncol = 3, risk.table.height = 0.4)

setwd("D:\\E\\??ʿ\\R_????\\BRCA\\Bulk_RNA-seq_data\\TCGA\\rep50\\result")
ggsave('TCGA_BRCA_clin_adjp_GSE9195_dmfs.pdf', res, width = 15, height = 21)
# ggsave('TCGA_BRCA_clin_adjp_GSE9195_rfs.pdf', res, width = 15, height = 21)

