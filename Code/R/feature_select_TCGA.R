
## 2021.4.1 ??4?址?????union??62??gene ????????


library(dplyr)       # ??>?? ?艿篮????牡??茫?????
library(tidyr)
library(tidyverse)   # tibble ?牡???


rm(list = ls())

setwd("D:\\E\\??士\\R_????\\BRCA\\Data\\TCGA_NEW")
Data1 = read.table("TCGA_BRCA_clin_1142_1080_scale.txt", header = T, check.names = FALSE)
# View(Data1[,1:10])

# setwd("D:\\E\\??士\\R_????\\BRCA\\Data\\TCGA_NEW\\result")
# gene = read.csv("union_gene.csv", header=TRUE, sep = ',')
# dim(Data1)    # 21836   121

setwd("D:\\E\\??士\\R_????\\BRCA\\Data\\TCGA_NEW\\result")
coef_gene <- read.csv("univariate_cox_coef.csv", header = T, sep=',')
gene <- as.matrix(coef_gene[,1])

col
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

setwd("D:\\E\\??士\\R_????\\BRCA\\Data\\TCGA_NEW\\result")
write.table(genedata2,"GSE_TCGA_3_ori.txt",quote=F,sep="\t")  # 2020.6.29# ??????TCGA?玫???PRS指?? ----------------------------------------------------------------

rm(list = ls())

setwd("D:\\E\\??士\\R_????\\BRCA\\Data\\TCGA_NEW\\result")

data = read.table("GSE_TCGA_3_ori.txt", header = T, check.names = FALSE)
datat <- as.data.frame(t(data))
prs = read.table("TCGA_OS_3gene.txt", header = T, check.names = FALSE)

datat$PRS <- prs$riskscore

data1 <- datat[,c(1,2,6,3,4,5)]

write.table(data1, file = "prs_TCGA_for_hiplot.txt", quote=F,sep="\t",row.names = F)  












