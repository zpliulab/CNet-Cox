
## clear
rm(list = ls())

## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- '/home/lly/R/'

setwd(paste(path, 'CNetCox/', sep=''))


# Real Data1 --------------------------------------------------------------


# df <- readxl::read_excel("Pathologicstage.xlsx")   # , sheet="AUROC"
df <- read.csv("Pathologicstage.csv", header = T, sep = ",")   # , sheet="AUROC"

library(stargazer)
# options(digits = 3)	

df0 <- as.matrix(df) 
stargazer(df0)
