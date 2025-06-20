
rm(list = ls())

# load R packages --------------------------------------------------------------------

# install.packages("glmnet")

library(glmSparseNet)
library(curatedTCGAData)
library(TCGAutils)
library(dplyr)


## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- '/home/lly/R/'
setwd(paste(path, 'CNetCox/Data/TCGA_NEW', sep=''))


## ????????ģ?͹???????Ҫ???????Ǵ???????????Ʒ???????????϶?Ӧ?????????˵???????Ϣ

# params <- list(seed = 29221)  

brca <- curatedTCGAData(diseaseCode = "BRCA", assays = "RNASeq2GeneNorm",
                        version = "1.1.38", dry.run = FALSE)
dim(brca)

brca.primary.solid.tumor <- TCGAutils::splitAssays(brca, '01')
xdata.raw <- t(assay(brca.primary.solid.tumor[[1]]))
dim(xdata.raw)    # 1093 20501
# View(xdata.raw[,1:10])

# Get survival information
ydata.raw <- colData(brca.primary.solid.tumor) %>% as.data.frame %>%
  # Keep only data relative to survival or samples
  select(patientID, vital_status, gender,
         Days.to.date.of.Death, Days.to.Date.of.Last.Contact,
         days_to_death,         days_to_last_followup,
         Vital.Status,    # ???????Ժ????Լ????ϵ?clinָ??
         years_to_birth, pathologic_stage, pathology_T_stage, pathology_N_stage, pathology_M_stage) %>%
  # Convert days to integer (??��???ּ???˫????)
  mutate(Days.to.date.of.Death = as.integer(Days.to.date.of.Death)) %>%
  mutate(Days.to.Last.Contact  = as.integer(Days.to.Date.of.Last.Contact)) %>%
  # Find max time between all days (ignoring missings)
  rowwise %>%
  mutate(time = max(days_to_last_followup,        Days.to.date.of.Death,
                    Days.to.Last.Contact, days_to_death, na.rm = TRUE)) %>%
  # Keep only survival variables and codes
  select(patientID, status = vital_status, time,   gender,  # ???????Ժ????Լ????ϵ?clinָ??
         years_to_birth, pathologic_stage, pathology_T_stage, pathology_N_stage, pathology_M_stage) %>%
  # Discard individuals with survival time less or equal to 0
  filter(!is.na(time) & time > 0) %>% as.data.frame


# Set index as the patientID
rownames(ydata.raw) <- ydata.raw$patientID

# Get matches between survival and assay data
xdata.raw_1 <- xdata.raw[TCGAbarcode(rownames(xdata.raw)) %in%
                           rownames(ydata.raw),]

dim(xdata.raw_1)    # 1080 20501
# View(xdata.raw_1[,1:10])
# Order ydata the same as assay
ydata.raw    <- ydata.raw[TCGAbarcode(rownames(xdata.raw_1)), ]
# View(ydata.raw)


xdata <- xdata.raw_1
ydata <- ydata.raw 
# ydata <- ydata.raw %>% select(time, status)
# View(xdata[,1:10])
# View(ydata)

data <- as.matrix(cbind(ydata, xdata))
rownames(data) <- rownames(xdata) 
dim(data)    # 1080 20221
# View(data[,1:10])

ydata.raw_2 <- data[,1:9]
ydata.raw_2[,1] <- rownames(ydata.raw_2)
View(ydata.raw_2)
# write.csv(ydata.raw_2, file = "TCGA_clin_information.csv", row.names = F)
# write.csv(ydata.raw_2, file = "TCGA_clin_information_gender.csv", row.names = F)

###########################################################################################
rm(list = ls())


## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- '/home/lly/R/'
setwd(paste(path, 'CNetCox/Data/TCGA_NEW', sep=''))
# clin <- read.csv("TCGA_clin_information.csv", header = T)
clin <- read.csv("TCGA_clin_information_gender.csv", header = T)
# View(clin)


######################### ???漸???ļ???????Sanger box ####################################

# clin_age <- clin[,c(3,2,4)]
# # write.csv(clin_age, file = "result\\Sanger_age_TCGA.csv", row.names = F)
# 
# # pathologic_stage      8 ??NA
# clin[which(clin[,5] == "stage i"), 5] <- c(1)   
# clin[which(clin[,5] == "stage ia"), 5] <- c(1)  
# clin[which(clin[,5] == "stage ib"), 5] <- c(1)
# clin[which(clin[,5] == "stage ii"), 5] <- c(2)  
# clin[which(clin[,5] == "stage iia"), 5] <- c(2)
# clin[which(clin[,5] == "stage iib"), 5] <- c(2)  
# clin[which(clin[,5] == "stage iii"), 5] <- c(3)
# clin[which(clin[,5] == "stage iiia"), 5] <- c(3)
# clin[which(clin[,5] == "stage iiib"), 5] <- c(3)
# clin[which(clin[,5] == "stage iiic"), 5] <- c(3)
# clin[which(clin[,5] == "stage iv"), 5] <- c(4)
# clin[which(clin[,5] == "stage x"), 5] <- c(5)
# 
# clin_stage <- clin[,c(3,2,5)]
# # write.csv(clin_stage, file = "result\\Sanger_stage_TCGA.csv", row.names = F)


#########################################################################################

# # years_to_birth
# 
# mean(na.omit(clin[,4]))    # ƽ??????
# 
# clin[which(clin[,4] <= "58"), 4] <- c(0) 
# clin[which(clin[,4] > "58"), 4] <- c(1) 

clin$pathologic_stage
# pathologic_stage      8 NA
# clin[which(clin[,5] == "stage i"), 5] <- c(0)
# clin[which(clin[,5] == "stage ia"), 5] <- c(0)
# clin[which(clin[,5] == "stage ib"), 5] <- c(0)
# clin[which(clin[,5] == "stage ii"), 5] <- c(0)
# clin[which(clin[,5] == "stage iia"), 5] <- c(0)
# clin[which(clin[,5] == "stage iib"), 5] <- c(0)
# clin[which(clin[,5] == "stage iii"), 5] <- c(1)
# clin[which(clin[,5] == "stage iiia"), 5] <- c(1)
# clin[which(clin[,5] == "stage iiib"), 5] <- c(1)
# clin[which(clin[,5] == "stage iiic"), 5] <- c(1)
# clin[which(clin[,5] == "stage iv"), 5] <- c(1)
# clin[which(clin[,5] == "stage x"), 5] <- c(1)


## 2023.5.9 re-split (I 期 0； others 1)
k <- 6
clin[which(clin[,k] == "stage i"), k] <- c(0)
clin[which(clin[,k] == "stage ia"), k] <- c(0)
clin[which(clin[,k] == "stage ib"), k] <- c(0)
clin[which(clin[,k] == "stage ii"), k] <- c(1)
clin[which(clin[,k] == "stage iia"), k] <- c(1)
clin[which(clin[,k] == "stage iib"), k] <- c(1)
clin[which(clin[,k] == "stage iii"), k] <- c(1)
clin[which(clin[,k] == "stage iiia"), k] <- c(1)
clin[which(clin[,k] == "stage iiib"), k] <- c(1)
clin[which(clin[,k] == "stage iiic"), k] <- c(1)
clin[which(clin[,k] == "stage iv"), k] <- c(1)
clin[which(clin[,k] == "stage x"), k] <- c(0)



# pathology_T_stage
## https://zhuanlan.zhihu.com/p/79308806
kt <- 7
clin[which(clin[,kt] == "t1"), kt] <- c(0)   
clin[which(clin[,kt] == "t1a"), kt] <- c(0)  
clin[which(clin[,kt] == "t1b"), kt] <- c(0)
clin[which(clin[,kt] == "t1c"), kt] <- c(0)  
clin[which(clin[,kt] == "t2"), kt] <- c(0)
clin[which(clin[,kt] == "t2a"), kt] <- c(0)
clin[which(clin[,kt] == "t2b"), kt] <- c(0)
clin[which(clin[,kt] == "t3"), kt] <- c(1)
clin[which(clin[,kt] == "t3a"), kt] <- c(1)
clin[which(clin[,kt] == "t4"), kt] <- c(1)
clin[which(clin[,kt] == "t4b"), kt] <- c(1) 
clin[which(clin[,kt] == "t4d"), kt] <- c(1)
clin[which(clin[,kt] == "tx"), kt] <- c(0) 


# pathology_N_stage
## N0表示淋巴结未受影响，N1-N3依次表示淋巴结受影响程度和范围的增加，
kn <- 8
## Nx表示淋巴结受影响状况无法评估。
clin[!duplicated(clin[,kn]),kn]
clin[which(clin[,kn] == "n0"), kn] <- c(0)   
clin[which(clin[,kn] == "n0 (i-)"), kn] <- c(0) 
clin[which(clin[,kn] == "n0 (i+)"), kn] <- c(0)
clin[which(clin[,kn] == "n0 (mol+)"), kn] <- c(0)
clin[which(clin[,kn] == "n1"), kn] <- c(1) 
clin[which(clin[,kn] == "n1a"), kn] <- c(1)
clin[which(clin[,kn] == "n1b"), kn] <- c(1) 
clin[which(clin[,kn] == "n1c"), kn] <- c(1) 
clin[which(clin[,kn] == "n1mi"), kn] <- c(1)
clin[which(clin[,kn] == "n2"), kn] <- c(1)
clin[which(clin[,kn] == "n2a"), kn] <- c(1)
clin[which(clin[,kn] == "n3"), kn] <- c(1) 
clin[which(clin[,kn] == "n3a"), kn] <- c(1)
clin[which(clin[,kn] == "n3b"), kn] <- c(1) 
clin[which(clin[,kn] == "n3c"), kn] <- c(1)
clin[which(clin[,kn] == "nx"), kn] <- c(0)


# pathology_M_stage
## M0表示没有转移，M1表示有远处转移。
kk <- 9
clin[!duplicated(clin[,kk]),kk]
clin[which(clin[,kk] == "m0"), kk] <- c(0)
clin[which(clin[,kk] == "cm0 (i+)"), kk] <- c(0)
clin[which(clin[,kk] == "m1"), kk] <- c(1) 
clin[which(clin[,kk] == "mx"), kk] <- c(0)


# Gender
kg <- 4
clin[!duplicated(clin[,4]),4]
clin[which(clin[,kg] == "female"), kg] <- c(1)
clin[which(clin[,kg] == "male"), kg] <- c(0)


# clin$years_to_birth <- factor(clin$years_to_birth,levels = c(0,1),labels = c("<=58",">58"))
# clin$pathologic_stage <- factor(clin$pathologic_stage,levels = c(0,1),labels = c("I/II","III/IV/V"))
# clin$pathology_T_stage <- factor(clin$pathology_T_stage,levels = c(0,1),labels = c("T1/T2","T3/T4/TX"))
# clin$pathology_N_stage <- factor(clin$pathology_N_stage,levels = c(0,1),labels = c("N0/N1","N2/N3/NX"))
# clin$pathology_M_stage <- factor(clin$pathology_M_stage,levels = c(0,1),labels = c("M0/CM0","M1/Mx"))


class(clin)
clinn <- as.matrix(clin[,-1])


# prs <- read.table("TCGA_OS_PRS_gene.txt",header = T, check.names = FALSE)
prs <- read.table("TCGA_OS_PRS_gene_scale.txt", header = T, check.names = FALSE)
clin1 <- cbind(clinn[,1], clinn[,2], prs[,3], clinn[,3], clinn[,4], clinn[,5], 
               clinn[,6], clinn[,7], clinn[,8])
# colnames(clin)
colnames(clin1) <- c("status", "time", "PRS", "Gender", "Age", "Pathologic stage", 
                     "Pathology T stage", "Pathology N stage", "Pathology M stage")

clin <- as.data.frame(clin1)
clin[3,3]
class(clin)

# as.data.frame(lapply(clin,as.numeric))
clin_n <- as.data.frame(lapply(clin[,c(1,2,3,4,5,6,7,8,9)],as.numeric))
clin_n[3,3]
# write.csv(clin_n, file = "TCGA_clin_information_hitplot_NEW.csv", row.names = F)


# install.packages("rms")
library(survival)
library(rms) 

## Error in summary.rms(f) : adjustment values not defined here or with datadist for years_to_birth

dd <- datadist(clin_n)
options(datadist="dd")

# ??��COX?ع鷽??
f <- cph(Surv(time, status) ~ PRS+Gender+Age+Pathologic.stage+Pathology.T.stage+Pathology.N.stage+Pathology.M.stage,
         x=T, y=T, surv=T, data=clin_n, time.inc=365)
clin_n$Pathology.T.stage

# ??COX?ع鷽?̽??н???
summary(f)


# ????????ͼ
surv <- Survival(f)
nom <- nomogram(f, fun=list(function(x) surv(365, x), 
                            function(x) surv(1095, x),   # 36-3year, 60-5year
                            function(x) surv(1825, x)), 
                lp=T, 
                funlabel=c("1 year survival", "3 year survival", "5 year survival"),
                maxscale=100,
                fun.at = c(0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))
plot(nom)

plot(nom, lplabel="Linear Predictor",
     xfrac=.35, varname.label=TRUE, varname.label.sep="=", ia.space=.2, 
     tck=NA, tcl=-0.20, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE, 
     cap.labels=FALSE, cex.var = 1.2, cex.axis = 1.05, lwd=5,
     label.every = 1, col.grid = gray(c(0.8, 0.95)))


pdf(file = "nomogram_prs.pdf",width = 12,height = 8)
# plot(nom)
plot(nom, lplabel="Linear Predictor",
     xfrac=.35, varname.label=TRUE, varname.label.sep="=", ia.space=.2, 
     tck=NA, tcl=-0.20, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE, 
     cap.labels=FALSE, cex.var = 1.2, cex.axis = 1.05, lwd=5,
     label.every = 1, col.grid = gray(c(0.8, 0.95)))
dev.off()



# ??У?????? -------------------------------------------------------------------
# https://www.iikx.com/news/statistics/1746.html
library(rms) 
s <- Surv(clin_n$time, clin_n$status, type="right")

## 365??
f <- cph(s~PRS+Gender+Age+Pathologic.stage+Pathology.T.stage+Pathology.N.stage+Pathology.M.stage, 
         x=TRUE, y=TRUE, surv = TRUE, time.inc=365, data=clin_n)

cal <- calibrate(f,u=365,cmethod='KM',m=180)  #??????uӦ??time.inc????һ??  216

# pdf(file = "calibration_1_NEW.pdf",width = 4.5,height = 5.0)
plot(cal, xlim = c(0.88,1), ylim= c(0.88,1),  #xlim??ylim?޶?x??y???????䣬???Ե???
     errbar.col=c(rgb(0,0,0,maxColorValue=255)),
     col=c(rgb(255,0,0,maxColorValue=255)))
abline(0,1,lty=3,lwd=2,col=c(rgb(0,0,255,maxColorValue= 255)))
# dev.off()


## 365*3   1095
f <- cph(s~PRS+Gender+Age+Pathologic.stage+Pathology.T.stage+Pathology.N.stage+Pathology.M.stage, 
         x=TRUE, y=TRUE, surv = TRUE, time.inc=1095, data=clin_n)

cal <- calibrate(f,u=1095,cmethod='KM',m=180)  #??????uӦ??time.inc????һ??  216

# pdf(file = "calibration_3_NEW.pdf",width = 4.5,height = 5.0)
plot(cal, xlim = c(0.68,1), ylim= c(0.68,1),  #xlim??ylim?޶?x??y???????䣬???Ե???
     errbar.col=c(rgb(0,0,0,maxColorValue=255)),
     col=c(rgb(255,0,0,maxColorValue=255)))
abline(0,1,lty=3,lwd=2,col=c(rgb(0,0,255,maxColorValue= 255)))
# dev.off()


## 365*3   1825
f <- cph(s~PRS+Gender+Age+Pathologic.stage+Pathology.T.stage+Pathology.N.stage+Pathology.M.stage, 
         x=TRUE, y=TRUE, surv = TRUE, time.inc=1825, data=clin_n)

cal <- calibrate(f,u=1825,cmethod='KM',m=180)  #??????uӦ??time.inc????һ??  216

# pdf(file = "calibration_5_NEW.pdf",width = 4.5,height = 5.0)
plot(cal, xlim = c(0.46,1), ylim= c(0.46,1),  #xlim??ylim?޶?x??y???????䣬???Ե???
     errbar.col=c(rgb(0,0,0,maxColorValue=255)),
     col=c(rgb(255,0,0,maxColorValue=255)))
abline(0,1,lty=3,lwd=2,col=c(rgb(0,0,255,maxColorValue= 255))) 
# dev.off()
