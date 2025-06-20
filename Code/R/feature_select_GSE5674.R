rm(list = ls())

## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- '/home/lly/R/'

# setwd(paste(path, 'CNetCox/Result/ComMethd', sep=''))


## data set
# GSE42568  GSE5764  

setwd(paste(path, 'BRCA/Bulk_RNA-seq_data/GEO_noclin/GSE5764', sep=''))


# load package ------------------------------------------------------------
library(dplyr)       
library(tidyr)
library(tidyverse)    

# load data ---------------------------------------------------------------
# Data1 = read.table("GSE5764_outcome_scale_IDC.txt", header = T, check.names = FALSE)

Data1 = read.table("GSE5764_outcome_scale_all.txt", header = T, check.names = FALSE)


# View(Data1[,1:10])
dim(Data1)    # 21836   121


setwd(paste(path, 'CNetCox/Data/', sep=''))
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
genedata1 <- apply(genedata1, 2, as.numeric)

# ## ?????? -- ??Ч
# genedata1 <- scale(genedata1)
# rownames(genedata1) <- genedata[,1]
# # View(genedata1[1:10,1:10])
# genedata1[3,3]

## ?Ի???-- ??Ч
genedata1 <-  t(scale(t(genedata1)))
rownames(genedata1) <- genedata[,1]
# View(genedata1[1:10,1:10])
genedata1[3,3]


genedata2 <- rbind(Data1[1,],genedata1)
# View(genedata2[,1:10])
write.table(genedata2,"Bulk_RNA-seq_data/GEO_noclin/GSE5764/GSE5764_6_ori.txt",quote=F,sep="\t")


# ????PRS??ϵ??,??ȡ??Ҫ??ϵ?? ----------------------------------------------------------------
phi_coef <- read.csv("TCGA_NEW/UniMutVariate_markergene.csv", header = T, sep=',')
genedata3 <- read.table("Bulk_RNA-seq_data/GEO_noclin/GSE5764/GSE5764_6_ori.txt", header = T, check.names = FALSE)
# View(genedata3[,1:10])

colnames(phi_coef) <- c('gene',"coef")
Data3 <- cbind(rownames(genedata3[-1,]), genedata3[-1,])
colnames(Data3) <- c('gene', colnames(genedata3))
# View(Data3[,1:10])

genedata <- merge(phi_coef, Data3, by = "gene")
# View(genedata[,1:10])
genedata[2,3]

## ????
A <- t(as.matrix(genedata[,3:dim(genedata)[2]]))
A <- apply(A, 2, as.numeric)
dim(A)    # 121  67
# View(A[,1:10])
A[2,1]
b <- as.matrix(genedata[,2])
# View(b)
dim(b)
b[2]
phi <- A %*% b
row.names(phi) <- colnames(genedata3)
phi[1,1]

phi1 <- as.data.frame(cbind(phi, t(genedata3[1,])))
colnames(phi1) <- c("PRS", "Label")
phi1$PRS <- as.numeric(phi1$PRS)
phi1[1,1]
write.table(phi1, file = "Bulk_RNA-seq_data/GEO_noclin/GSE5764/prs_GSE5764_NEW.txt", quote=F,sep="\t")


# ??????ͼ --------------------------------------------------------------------

pred_log = read.table("Bulk_RNA-seq_data/GEO_noclin/GSE5764/prs_GSE5764_NEW.txt", header = T, check.names = FALSE)
pred_log[which(pred_log[,2] == "Cancer"),2]  <- c("Tumor")
dim(pred_log)


# ??pֵ ---------------------------------------------------------------------


library(ggplot2)
library(ggpubr)
library(Rcpp)

my_comparisons <- list( c("Tumor", "Normal") )
plot <- ggboxplot(pred_log, x = "Label", y = "PRS",
                  color = "Label", palette = "npg"#,  "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty"
) +  # add = "jitter" , fill = c("#90CAF9", "#F48FB1")
  stat_boxplot(geom = "errorbar", width=0.30, size=0.6) +   # ʹ??????????ĩ?˶̺???
  geom_boxplot(size=0.6, fill=c("#099052", "#EF3122")) + 
  #size????????ͼ?ı߿??ߺͺ??????߿??ȣ?fill??????????ɫ
  stat_compare_means(comparisons=my_comparisons, method = "t.test") + # Add p-value
  theme(legend.position="none", #????Ҫͼ??
        axis.text.x=element_text(colour="black",size=14), #????x???̶ȱ?ǩ??????????
        axis.text.y=element_text(size=14,face="plain"), #????x???̶ȱ?ǩ??????????
        axis.title.y=element_text(size = 14,face="plain"), #????y???ı?????????????
        axis.title.x=element_text(size = 14,face="plain"), #????x???ı?????????????
        plot.title = element_text(size=15,face="bold",hjust = 0.5), #?????ܱ?????????????
        panel.grid.major = element_blank(), #????ʾ??????
        panel.grid.minor = element_blank())+
  ylab("Risk score")+xlab("") #????x????y???ı???
# stat_summary(fun = mean, geom = "point", 
# shape = 18, size = 4, color = "blue")  # "#C5E1A5" Add global p-value


plot

# pdf(file = "Box_GSE5764_NEW_color.pdf",width = 3,height = 4)
# plot
# dev.off()



# С????ͼ --------------------------------------------------------------------

ggplot(pred_log, aes(x=Label, y=PRS), color = "Label") + # , fill = c("#90CAF9", "#F48FB1")
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + # Add pairwise comparisons p-value
  geom_violin(aes(fill=Label))  # С????ͼ
