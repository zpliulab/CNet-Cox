
rm(list = ls())

library(survminer)
library(survival)
library(ggplot2)


## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- '/home/lly/R/'
setwd(paste(path,'CNetCox/Feature_data', sep=''))

myfile = list.files("Xtile")                #list.files???input?ļ??????????ļ???????a
dir = paste("./Xtile/", myfile, sep="")     #??paste?????·????��dir
n = length(dir) 


# ????ȡcut-off,??????ͼ ---------------------------------------------------------
i <- 3

d9

data = read.table(file = dir[i],header=T, check.names = FALSE) #??????һ???ļ????ݣ????Բ????ȶ?һ????????Ϊ?˼򵥣?ʡȥ????data.frame??ʱ?䣬??ѡ???ȶ???һ???ļ?
dir[i]
gene_name <- colnames(data)[3]
exprSet <- data.frame(t(data))

## ???ýض?ֵ  https://www.jianshu.com/p/e72605df6348
alpha <- 0.328

risk_score  <- t(as.matrix(exprSet[gene_name,]))
cut_off <- rep(as.numeric(quantile(exprSet[gene_name,],alpha)), dim(exprSet)[2])


data$time <- data$time/12
data$riskscore <- ifelse(risk_score > cut_off, 'high','low')
table(data$risk_score)

fit <- survfit(Surv(time, status)~riskscore, data = data)

p <- ggsurvplot(fit, data = data, 
                conf.int = F, # ????????????
                # surv.median.line = "hv",  # ??????λ????ʱ??
                risk.table = TRUE, # ?????ۼƷ???????
                tables.height= 0.25, # ????Ϊ????��?????????ĸ߶?
                cumcensor = T,    # ?????ۻ?ɾʧ??
                legend = c(0.83,0.95), # ָ??ͼ??λ??
                
                # P value
                pval = TRUE, 
                pval.size=6, 
                font.pval= c(14, "bold", "black"),
                pval.coord = c(0.00, 0.05), #????Pval??λ??
                
                # legend
                legend.title = '', # gene_name
                legend.labs=c("High risk", "Low risk"), #??ǩ
                font.legend= c(12, "plain", "black"),  # ͼ??????
                # font.main = c(100, "bold", "black"),
                # xlim = c(0,72), # present narrower X axis, but not affect
                # survival estimates.
                palette=c("red", "blue"),
                font.x = c(12, "plain", "black"),
                font.y = c(12, "plain", "black"), # x?ᡢy??????
                font.tickslab = c(12, "plain", "black"), # ?̶ȱ?ǩ????
                xlab = "Time in years", # customize X axis label. year
                break.time.by = 2
) # break X axis in time intervals by 500.
p


# ????HR??CI https://www.zhihu.com/question/347864496/answer/836230749 -------

res_cox <- coxph(Surv(time, status) ~riskscore, data=data)
HR <- round(summary(res_cox)$conf.int[1],2)
ci_l <- round(summary(res_cox)$conf.int[3],2)
ci_r <- round(summary(res_cox)$conf.int[4],2)

p1 <- p
p1$plot = p1$plot + 
  ggplot2::annotate("text",x = 5.15, y = 0.3, label = paste("HR : ", HR), size = 5) + 
  ggplot2::annotate("text",x = 6.15, y =,
                    label = paste("(","95%CI : ", ci_l,"-",ci_r,")", sep = ""), size = 5)
p1


library(tidyverse)
name <- dir[i]
name1 <- str_split_fixed(name, "./", 2);
name2 <- str_split_fixed(name1[2], "/", 2)
name3 <- str_split_fixed(name2[2], "[.]", 2);
name4 <- name3[1]
name5 <- str_c(name4,"_ex")

path <- paste("./Xtilefigure/", paste(name5,".pdf", sep=""), sep="")
pdf(path, width = 4.8, height = 6, onefile = FALSE)  # ????onefile ????ΪFALSE ????ɢ??ͼ?Ḳ??ǰ???Ŀհ?
p1
dev.off() 


