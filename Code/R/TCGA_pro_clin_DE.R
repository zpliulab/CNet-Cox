## 2023.4.20 LLY 
## This code is used to select DEGs, from two ways
## output DEGs and volcano plot


##############################################################################
##############################################################################
# DEGs with normal and tumor information ----------------------------------

rm(list = ls())

## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- 'lly/home/R'
setwd(paste(path, 'CNetCox/Data/', sep=''))


# install R package --------------------------------------------------------------------

# if (!require("BiocManager"))
#   install.packages("BiocManager")
# BiocManager::install("curatedTCGAData")

# if (!require("BiocManager"))
#   install.packages("BiocManager")
# BiocManager::install("glmSparseNet")


# load R package --------------------------------------------------------------------
library(glmSparseNet)
library(curatedTCGAData)
library(TCGAutils)
library(dplyr)


# get data from TCGA----------------------------------------------------
# brca  <- curatedTCGAData(diseaseCode = "BRCA", assays = "RNASeq2GeneNorm", FALSE) 
brca <- curatedTCGAData(diseaseCode = "BRCA", assays = "RNASeq2GeneNorm",
                        version = "1.1.38", dry.run = FALSE)
# keep only solid tumour (code: 01)
brca <- TCGAutils::TCGAsplitAssays(brca, c('01','11'))
xdata.raw <- t(cbind(assay(brca[[1]]), assay(brca[[2]])))
dim(xdata.raw)    # 1205 20501
# View(xdata.raw[,1:10])


# Get matches between survival and assay data
class.v        <- TCGAbiospec(rownames(xdata.raw))$sample_definition
# class.v        <- TCGAbiospec(rownames(xdata.raw))$sample_definition %>% factor
names(class.v) <- rownames(xdata.raw)
# View(class.v)

# View(rownames(xdata.raw))
# View( TCGAbiospec(rownames(xdata.raw))$sample_definition %>% factor )


# keep features with standard deviation > 0
xdata_raw <- xdata.raw %>%
{ (apply(., 2, sd) != 0) } %>%
{ xdata.raw[, .] }

dim(xdata_raw)    # 1205 20222
View(xdata_raw[,1:10])

## select some specific genes
# small.subset <- c('CD5', 'CSF2RB', 'HSF1', 'IRGC', 'LRRC37A6P', 'NEUROG2',
#                   'NLRC4', 'PDE11A', 'PIK3CB', 'QARS', 'RPGRIP1L', 'SDC1',
#                   'TMEM31', 'YME1L1', 'ZBTB11',
#                   sample(colnames(xdata.raw), 100))
# xdata <- xdata.raw[, small.subset[small.subset %in% colnames(xdata.raw)]]

small_subset <- colnames(xdata.raw)
xdata <- xdata_raw[, small_subset[small_subset %in% colnames(xdata_raw)]]
dim(xdata)    # 1205 20222
# View(xdata[,1:10])
xdatat <- t(xdata)

ydata <- class.v
# View(ydata)
# class(ydata)


# 112 Normal and 112 Tumor -----------------------------------------
## Normal
sample_N <- rownames(xdata)[which(ydata == "Solid Tissue Normal")]
sample_N1 <- sample_N %>%
  as_tibble() %>%
  mutate(sample_N = substr(sample_N, 1, 12))
## Tumor
sample <- rownames(xdata)
sample1 <- sample %>%
  as_tibble() %>%
  mutate(sample = substr(sample, 1, 12))

## match
lab <- which(as.matrix(sample1[,2]) %in% as.matrix(sample_N1[,2]))
xdata_TN <- xdata[lab,]
dim(xdata_TN)    # 224 20222
# View(xdata_TN[,1:10])

ydata_TN <- rbind(as.matrix(rep(c("Tumor"), 112)), as.matrix(rep(c("Normal"), 112)))
colnames(ydata_TN) <- c("outcome")
# View(ydata_TN)

data_TN <- cbind(ydata_TN, xdata_TN)
dim(data_TN)    # 224 20223
# write.table(t(data_TN), "TCGA_pro_outcome_TN.txt",quote=F,sep="\t")


##############################################################################
##############################################################################
# Performe DEG analysis ---------------------------------------------------
rm(list = ls())

## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- 'lly/home/R'
setwd(paste(path, 'CNetCox/Data/', sep=''))


## load R package -------------------------------------------------------------------
# BiocManager::install("DESeq2",force = TRUE)
# BiocManager::install("pasilla")
# install.packages("rlang")   ## need install in R not Rstudio
library(DESeq2)
library(pasilla)


## load data
data <- read.table("TCGA_pro_outcome_TN.txt",header=T,sep='\t', check.names = F)
dim(data)   # 20223   224
# View(data[,1:10])
xdatat <- t(as.matrix(apply(as.matrix(data[-1,]),1,function(x) as.numeric(x))))
colnames(xdatat) <- colnames(data)
# View(xdatat[,1:10])
xdatat[1,2]
ydata <- data[1,]
class(xdatat)


# DESeq RNA-seq -------------------------------------------------
exprSet <- round(xdatat)
dim(exprSet)    # 20222   224
# View(exprSet[,1:10])

ydata_TN <- rbind(as.matrix(rep(c("Tumor"), 112)), as.matrix(rep(c("Normal"), 112)))
colnames(ydata_TN) <- c("outcome")
# View(ydata_TN)

group_list <- as.factor(ydata_TN)
# group_list <- xydata[,1]
colData <- data.frame(row.names=colnames(exprSet), group_list=group_list)
colnames(colData) <- c("outcome")
dim(colData)

dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ outcome)

##Step2 Using DESeq function directly
dds2 <- DESeq(dds)  
resultsNames(dds2)
suppressMessages(dds)


# Extract the DGE results, compare Normal with Tumor 
res <-  results(dds2, contrast=c("outcome","Tumor","Normal"))
summary(res) 
## 查看fdr校正后的P<0.05的个数
table(res$padj<0.05)     # 14034
# plotMA(res)
# # BiocManager::install("apeglm")
# library(apeglm)
# resLFC <- lfcShrink(dds2, coef="outcome_Tumor_vs_Normal", type="apeglm") #????lfcShrink ????log2 fold change
# plotMA(resLFC) 

res_order <- res[order(res$padj),]
res_order <- as.data.frame(res_order)
# write.csv(res_order,file= "DEG_res_order_TN.csv")


# set thrshold -----------------------------------------------------------
## FDR adjusted --Defult FDR is 0.1
res1 <-  results(dds2, alpha = 0.01)  
# write.csv(res1,file= "DEG_res.csv")

## log2(10) = 3.321928
diff_gene_deseq2 <- subset(res1, padj < 0.01 & abs(log2FoldChange) > 3.321928)
diff_gene_deseq2 <- as.data.frame(diff_gene_deseq2)
dim(diff_gene_deseq2)
## save the DEGs information
# write.csv(diff_gene_deseq2,file= "DEG_Tumor_vs_Normal_10.csv")


# Plot volcano ------------------------------------------------------------
logFC <- res1$log2FoldChange
## pvalue is adjusted p value
pvalue <- res1$padj
genes <- rownames(res1)
## LogFC is Log2FC 
data <- data.frame(logFC, pvalue, genes)

# Example dataset
# logFC <- c(2.1, 1.9, 0.3, -0.5, -1.2, -1.8)
# pvalue <- c(0.0001, 0.0005, 0.01, 0.05, 0.1, 0.2)
# genes <- c("A", "B", "C", "D", "E", "F")
# data <- data.frame(logFC, pvalue, genes)

# Volcano plot
library(ggplot2)

# Define colors
color_significant <- "darkred"
color_nonsignificant <- "black"

# Define significance threshold
sig_threshold <- -log10(0.01)
# sig_threshold <- -log2(10)

# Create plot
ggplot(data, aes(x=logFC, y=-log10(pvalue))) +
  geom_point(aes(color=ifelse(abs(logFC)>log2(10) & -log10(pvalue)>2, "significant", "nonsignificant")), size=2) +
  scale_color_manual(values=c("significant"=color_significant, "nonsignificant"=color_nonsignificant)) +
  geom_vline(xintercept=c(-1, 1), linetype="dashed", color="gray50", size=0.5) +
  geom_hline(yintercept=sig_threshold, linetype="dashed", color="gray50", size=0.5) +
  labs(x="Log2 Fold Change", y="-Log10 P.adj-value", color="Significance") +
  theme_classic() +
  theme(axis.line=element_line(color="black", size=0.5),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10),
        legend.position="bottom",
        legend.key.size=unit(0.5,"cm"),
        panel.background = element_rect(fill="white"),
        plot.title=element_text(hjust=0.5, size=14, face="bold"),
        plot.subtitle=element_text(hjust=0.5, size=12, face="italic"),
        plot.caption=element_text(hjust=1, size=10),
        plot.margin=unit(c(1,1,1,1),"cm"))
        # plot.margin=unit(c(1,1,1,1),"cm")) +
  # ggtitle("Volcano Plot of Gene Expression Data",
  #         subtitle="Genes with absolute log2 fold change >3.321928 and -log10 p-value >2 are considered significant") +
  # labs(caption="Source: BRCA dataset")


# Data normalized -------------------------------------------------------------
# normalized_counts <- counts(dds2, normalized=TRUE)
# View(normalized_counts[,1:10])
vst_dat <- vst(dds2, blind = TRUE)
dat111 <- assay(vst_dat)
dim(dat111)    # 20222   224
# View(dat111[,1:10])
dat111[4,3]
# write.csv(dat111,file= "deseq_nor.csv")
data0 <- rbind(t(ydata_TN), dat111)
dim(data0)    # 20223   224
# View(data0[,1:10])
data0[3,3]
# View(data0[,1:10])
## the data with ""
# write.table(data0,"TCGA_pro_outcome_TN_log.txt",quote=F,sep="\t") 


##############################################################################
##############################################################################
# DEGs with clinical information ------------------------------------------
rm(list = ls())

## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- 'lly/home/R'
setwd(paste(path, 'CNetCox/Data/', sep=''))


# load R package --------------------------------------------------------------------
library(glmSparseNet)
library(curatedTCGAData)
library(TCGAutils)
library(dplyr)


# xtract the expression matrix of TNBC subtypes from the BRCA dataset of the TCGA database --------------------------------- --------
# The data required for constructing a survival analysis model is a pure tumor sample expression matrix plus the corresponding tumor patient survival information.
# params <- list(seed = 29221)  
# brca <- curatedTCGAData(diseaseCode = "BRCA", assays = "RNASeq2GeneNorm", FALSE)
brca <- curatedTCGAData(diseaseCode = "BRCA", assays = "RNASeq2GeneNorm",
                        version = "1.1.38", dry.run = FALSE)
head(getSubtypeMap(brca))
head(getClinicalNames("BRCA"))


brca.primary.solid.tumor <- TCGAutils::splitAssays(brca, '01')
xdata.raw <- t(assay(brca.primary.solid.tumor[[1]]))
# dim(xdata.raw)    # 1093 20501
# View(xdata.raw[,1:10])


# keep features with standard deviation > 0
xdata.raw <- xdata.raw %>%
  { (apply(., 2, sd) != 0) } %>%
  { xdata.raw[, .] }
dim(xdata.raw)    # 1093 20220
# View(xdata.raw[,1:10])


# Get survival information
ydata.raw <- colData(brca.primary.solid.tumor) %>% as.data.frame %>%
  # Keep only data relative to survival or samples
  select(patientID, vital_status,
         Days.to.date.of.Death, Days.to.Date.of.Last.Contact,
         days_to_death,         days_to_last_followup,
         Vital.Status) %>%
  # Convert days to integer (本来数字加了双引号)
  mutate(Days.to.date.of.Death = as.integer(Days.to.date.of.Death)) %>%
  mutate(Days.to.Last.Contact  = as.integer(Days.to.Date.of.Last.Contact)) %>%
  # Find max time between all days (ignoring missings)
  rowwise %>%
  mutate(time = max(days_to_last_followup,        Days.to.date.of.Death,
                    Days.to.Last.Contact, days_to_death, na.rm = TRUE)) %>%
  # Keep only survival variables and codes
  select(patientID, status = vital_status, time) %>%
  # Discard individuals with survival time less or equal to 0
  filter(!is.na(time) & time > 0) %>% as.data.frame

# View(ydata.raw)
# dim(ydata.raw)    # 1080    3


# Set index as the patientID
rownames(ydata.raw) <- ydata.raw$patientID

# Get matches between survival and assay data
xdata.raw_1 <- xdata.raw[TCGAbarcode(rownames(xdata.raw)) %in%
                           rownames(ydata.raw),]

dim(xdata.raw_1)    # 1080 20220
# View(xdata.raw_1[,1:10])


# Order ydata the same as assay
ydata.raw    <- ydata.raw[TCGAbarcode(rownames(xdata.raw_1)), ]
# View(ydata.raw)

xdata <- xdata.raw_1
ydata <- ydata.raw %>% select(status)
# ydata <- ydata.raw %>% select(time, status)
# View(xdata[,1:10])
# View(ydata)

data <- as.matrix(cbind(ydata, xdata))
rownames(data) <- rownames(xdata) 
dim(data)    # 1080 20221
# View(data[,1:10])
data1 <- t(data)
# write.table(data1,"TCGA_pro_outcome.txt",quote=F,sep="\t") 


##############################################################################
##############################################################################
# Performe DEG analysis ---------------------------------------------------

## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- 'lly/home/R'
setwd(paste(path, 'CNetCox/Data/', sep=''))


# Perform difference analysis according to survival status and normalize -------------------------------------------
## load data
data_outcome <- read.table("TCGA_NEW/TCGA_pro_outcome.txt",header=T,sep='\t', check.names = F)
dim(data_outcome)   # 20221  1080
# View(data_outcome[,1:10])
data_outcome[2,1]

## expression data
xdatat <- data_outcome[-1,]
dim(xdatat)   # 20220  1080
# View(xdatat[,1:10])

## Survival status
ydata <- t(data_outcome[1,])
sum(ydata[,1])    # 152 dead

## Label grouping
ydata[which(ydata[,1] == "1"), 1] <- c("Dead")
ydata[which(ydata[,1] == "0"), 1] <- c("Alive")


# DESeq2 package for differential analysis of RNA-seq data -------------------------------------------------

exprSet <- round(xdatat)
dim(exprSet)    # 20220  1080
# View(exprSet[,1:10])
group_list <- as.factor(ydata)
colData <- data.frame(row.names=colnames(exprSet), group_list=group_list)
colnames(colData) <- c("outcome")
dim(colData)


# load R package --------------------------------------------------------------------

# BiocManager::install("DESeq2")
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ outcome)

dds2 <- DESeq(dds)  ##The second step is to use the DESeq function directly
resultsNames(dds2)


# # Extract the DGE results, compare Dead with Alive 
res <-  results(dds2, contrast=c("outcome","Dead","Alive"))
summary(res) 
# plotMA(res)

res_order <- res[order(res$padj),]
res_order <- as.data.frame(res_order)
# write.csv(res_order,file= "DEG_res_order_DA.csv")


## ZIshaung‘s code, both count the same
res1 <-  results(dds2, alpha = 0.01)
# write.csv(res1,file= "DEG_res.csv")

# diff_gene_deseq2 <- subset(res1, padj < 0.01)    # 501
diff_gene_deseq2 <- subset(res1, padj < 0.01 & abs(log2FoldChange) > log2(2))
diff_gene_deseq2 <- as.data.frame(diff_gene_deseq2)
dim(diff_gene_deseq2)    # 501   6    196    6
# write.csv(diff_gene_deseq2,file= "DEG_Dead_vs_Alive_1.csv")


# Volcano plot ------------------------------------------------------------
logFC <- res$log2FoldChange
## pvalue is adjusted p value
pvalue <- res$padj
genes <- rownames(res)
## LogFC is Log2FC 
data <- data.frame(logFC, pvalue, genes)
# write.csv(data, "TCGA_NEW/Volcano_DE.csv", row.names = F)

# Volcano plot
library(ggplot2)

# Define colors
color_significant <- "#e7897d"
color_nonsignificant <- "#9fa0b5"

# Define significance threshold
sig_threshold <- -log10(0.01)
# sig_threshold <- -log2(10)

# Create plot
ggplot(data, aes(x=logFC, y=-log10(pvalue))) +
  geom_point(aes(color=ifelse(abs(logFC)>log2(2) & -log10(pvalue)>2, "significant", "nonsignificant")), size=2) +
  scale_color_manual(values=c("significant"=color_significant, "nonsignificant"=color_nonsignificant)) +
  geom_vline(xintercept=c(-1, 1), linetype="dashed", color="gray50", size=0.5) +
  geom_hline(yintercept=sig_threshold, linetype="dashed", color="gray50", size=0.5) +
  labs(x="Log2 Fold Change", y="-Log10 P.adj-value", color="Significance") +
  theme_classic() +
  theme(axis.line=element_line(color="black", size=0.5),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10),
        legend.position="bottom",
        legend.key.size=unit(0.5,"cm"),
        panel.background = element_rect(fill="white"),
        plot.title=element_text(hjust=0.5, size=14, face="bold"),
        plot.subtitle=element_text(hjust=0.5, size=12, face="italic"),
        plot.caption=element_text(hjust=1, size=10),
        plot.margin=unit(c(1,1,1,1),"cm"))
# plot.margin=unit(c(1,1,1,1),"cm")) +
# ggtitle("Volcano Plot of Gene Expression Data",
#         subtitle="Genes with absolute log2 fold change >1 and -log10 p-value >2 are considered significant") +
# labs(caption="Source: BRCA dataset")



# Volcano plot example  ---------------------------------------------------
library(ggplot2)
library(ggrepel)   

# Read volcano map data file
data1 = read.delim("https://www.bioladder.cn/shiny/zyp/bioladder2/demoData/Volcano/Volcano.txt",
                  header = T     
)

FC = 1.5 
PValue = 0.05  

data$sig[(-1*log10(data$PValue) < -1*log10(PValue)|data$PValue=="NA")|(log2(data$FC) < log2(FC))& log2(data$FC) > -log2(FC)] <- "NotSig"
data$sig[-1*log10(data$PValue) >= -1*log10(PValue) & log2(data$FC) >= log2(FC)] <- "Up"
data$sig[-1*log10(data$PValue) >= -1*log10(PValue) & log2(data$FC) <= -log2(FC)] <- "Down"

data$label=ifelse(data$Marker == 1, as.character(data$Name), '')

## plot
ggplot(data,aes(log2(data$FC),-1*log10(data$PValue))) +     
  geom_point(aes(color = sig)) +                            
  labs(title="volcanoplot",                                 
       x="log[2](FC)", 
       y="-log[10](PValue)") + 
  # scale_color_manual(values = c("red","green","blue")) +  
  geom_hline(yintercept=-log10(PValue),linetype=2)+         
  geom_vline(xintercept=c(-log2(FC),log2(FC)),linetype=2)+  
  geom_text_repel(aes(x = log2(data$FC),                    
                      y = -1*log10(data$PValue),          
                      label=label),                       
                  max.overlaps = 10000,                     
                  size=3,                                   
                  box.padding=unit(0.5,'lines'),            
                  point.padding=unit(0.1, 'lines'), 
                  segment.color='black',                   
                  show.legend=FALSE)


# Scale-------------------------------------------------------------------

vst_dat <- vst(dds2, blind = TRUE)
dat111 <- assay(vst_dat)
dim(dat111)    # 20220  1080
# View(dat111[,1:10])


data_norm_clin <- t(dat111)
dim(data_norm_clin)    # 1080 20220
# View(data_norm_clin[,1:10])

# 带有生存时间的 -----------------------------------------------------------------

ydata_clin <- ydata.raw %>% select(time, status)
# View(xdata[,1:10])
# View(ydata_clin)


data_clin <- as.matrix(cbind(ydata_clin, data_norm_clin))
rownames(data_clin) <- rownames(data_norm_clin) 
dim(data_clin)    # 1080 20222
View(data_clin[,1:10])

data_clin1 <- t(data_clin)
# write.table(data_clin1,"TCGA_pro_norm_clin.txt",quote=F,sep="\t") 


# Intergather with prior knowledges ---------------------------------------
rm(list = ls())

## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- 'lly/home/R'
setwd(paste(path, 'CNetCox/Data/', sep=''))


## load data
Data = read.table("TCGA_NEW/TCGA_pro_norm_clin.txt", header = T, check.names = FALSE)
dim(Data)    # 20222  1080
gene <- as.vector(rownames(Data)[-c(1,2)])


# DE_outcome <- read.csv("DEG_Tumor_vs_Normal_10.csv", header = T, sep=',')    # 489
DE_clin <- read.csv("TCGA_NEW/DEG_Dead_vs_Alive_1.csv", header = T, sep=',')    # 196
# Degene <- as.vector(union(DE_outcome[,1], DE_clin[,1]))    # 956
Degene <- as.vector(DE_clin[,1])    # 196


## see non DEGs
sub <- setdiff(gene, Degene)   


# biomarker ---------------------------------------------------------------
mark <- as.matrix(read.csv("Prior_infor/128 biomarkers.csv", header = T, sep=','))
scmark <- as.matrix(c("KRT15","UBE2C","TOP2A","KRT6B","MKI67","HMGB2","ASPM","CDC20","KIF20A","CDK"))
mama_70 <- as.matrix(read.csv("Prior_infor/70_mama.csv", header = T, sep=',') )
KEGG_147 <- as.matrix(read.csv("Prior_infor/KEGG147.csv", header = T, sep=',') )
tf <- as.matrix(read.csv("Prior_infor/tfgene119.csv", header = T, sep = ','))


## input marker genes
marker <- as.matrix(read.csv("Result/Result3/marker3.csv", header = T, sep=','))


## 68 marker genes and prior knowledge's overlap
dim(as.matrix(which(marker %in% mark)))    # 10
marker_128 <- as.matrix(marker[which(marker %in% mark),])
colnames(marker_128) <- "OSbrca"
# write.csv(marker_128,"Result/Result3/marker_128.csv", row.names = F)

dim(as.matrix(which(marker %in% mama_70)))    # 5
marker_70 <- as.matrix(marker[which(marker %in% mama_70),])
colnames(marker_70) <- "Mama"
# write.csv(marker_70,"Result/Result3/marker_70.csv", row.names = F)

dim(as.matrix(which(marker %in% KEGG_147)))    # 27
marker_147 <- as.matrix(marker[which(marker %in% KEGG_147),])
colnames(marker_147) <- "KEGG"
# write.csv(marker_147,"Result/Result3/marker_147.csv", row.names = F)

dim(as.matrix(which(marker %in% tf)))    # 32
marker_tf <- as.matrix(marker[which(marker %in% tf),])
colnames(marker_tf) <- "tf"
# write.csv(marker_tf,"Result/Result3/marker_tf.csv", row.names = F)

dim(as.matrix(which(marker %in% as.matrix(Degene))))    # 32
marker_DE <- as.matrix(marker[which(marker %in% Degene),])
colnames(marker_DE) <- "DE"
# write.csv(marker_DE,"Result/Result3/marker_DE.csv", row.names = F)

DE <- as.matrix(read.csv("TCGA_NEW/DEG_Dead_vs_Alive_1.csv", header = T,  sep = ','))
dim(as.matrix(which(as.matrix(DE[,1]) %in% marker)))    # 32


# 110+9+56+145+117  gene Integration ------------------------------------------------------------
add_gene <- as.matrix( union(union(union(union(intersect(sub, mark), 
                                                     intersect(sub, scmark)), 
                                               intersect(sub, mama_70)), 
                                         intersect(sub, KEGG_147)), 
                                  intersect(sub, tf)) )    # 807

allgene <- as.matrix(rbind(add_gene, as.matrix(Degene)))    #616
dim(allgene)    # 1003   ----    616
rownames(allgene) <- allgene[,1] 



# Somponent selection and node cut set ------------------------------------

net <- as.matrix(read.csv("Prior_infor/Regnetwork_hum.csv",header = T))
net[1,]
node_used <- allgene   
dim(node_used)    # 616
# net_used <- net[,c(2,4)]
net_used <- net
k1 <- which(net_used[,1] %in% node_used)   #  
k2 <- which(net_used[,2] %in% node_used)   #  
length(intersect(k1,k2))    # 10414  ----   5980
# length(union(k1,k2))
used <- net_used[intersect(k1,k2),]
dim(used)     # 5980    2

# install.packages('igraph')
library(igraph)
PP <- graph_from_data_frame(used,directed = F)
p1 <- simplify(PP, remove.loops = T, remove.multiple = T)   
# ed <- as_edgelist(p1, names = TRUE)


# Compute the largest (weakly or strongly) connected component of a graph ---------------------------------------------------------
# g <- sample_gnp(20, 1/20)
clu <- components(p1)
groups(clu)    # 544


component <- groups(clu)$'1'
# write.csv(component, file = "TCGA_NEW/UNgene_component.csv", row.names = F)
node_used <- component
dim(as.matrix(node_used))    # 544
# net_used <- net[,c(2,4)]
k1 <- which(net_used[,1] %in% node_used)   #  
k2 <- which(net_used[,2] %in% node_used)   #  
length(intersect(k1,k2))    # 58980
used <- net_used[intersect(k1,k2),]

# install.packages('igraph')
library(igraph)
PP <- graph_from_data_frame(used,directed = F)
p1 <- simplify(PP, remove.loops = T, remove.multiple = T)  # 最终的数对
ed <- as_edgelist(p1, names = TRUE)
# write.csv(ed,"TCGA_NEW/UNgene_comp_net.csv",row.names = F,quote = F)


# Extract DEgene_component data ------------------------------------------------------------
data <- Data[component,]
dim(data)    # 544 1080
# View(data[,1:10])

all_data <- rbind(Data[c(1,2),], data)
dim(all_data)    # 546 1080
# View(all_data[,1:10])
# write.table(all_data, file = "TCGA_NEW/TCGA_BRCA_clin_546_1080.txt",quote = F, sep = "\t")


# scale -------------------------------------------------------------------
all_data_scale <- rbind(all_data[c(1,2),], t(scale(t(all_data[-c(1,2),]))))
dim(all_data_scale)    # 546 1080
View(all_data_scale[,1:10])
# write.table(all_data_scale, file = "TCGA_NEW/TCGA_BRCA_clin_546_1080_scale.txt",quote = F, sep = "\t")


##############################################################################
##############################################################################
# Node cut set ------------------------------------------------------------

## clear
rm(list = ls())

## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- 'lly/home/R'
setwd(paste(path, 'CNetCox/Data/', sep=''))


library(igraph)
gene <- read.csv('TCGA_NEW/UNgene_component.csv')
genes <- data.frame(gene[,1])
net <- read.csv('TCGA_NEW/UNgene_comp_net.csv')


net[which(net[,1] == "EP300"),2] == "CTCFL"

g <- graph_from_data_frame(net, directed=F)


plot(g, layout=layout.fruchterman.reingold,  
     vertex.size=4,   
     vertex.label = V(g)$name,  
     vertex.label.cex=0.7,  
     vertex.label.dist=0.4,  
     vertex.label.color = "black"   
)

# Calculate the maximum distance between two points ------------------------------------------------------------------

## diameter, breadth-first search
diameter(g) 

## TRUE, the diameters of the connected components
diameter(g, unconnected=TRUE)

## FALSE, the number of vertices
diameter(g, unconnected=FALSE)

## Weighted diameter
# E(g)$weight <- sample(seq_len(ecount(g)))   

## returns a path with the actual diameter
get_diameter(g) 

## returns two vertex ids, connected by the diameter path.
farthest_vertices(g) 
# diameter(g, weights=NA)

# Find the smallest cut  ----------------------------------------------------------------

## FALSE, the edges in the cut and a the two (or more) partitions are also returned.
min_cut(g, source = "BRF1", target = "EGLN1", value.only = FALSE)
# min_cut(g, source = 2, target = 5, value.only = FALSE)

## TRUE, only the minumum cut value is returned
min_cut(g, source = "BRF1", target = "EGLN1", value.only = TRUE)
# min_cut(g, source = 2, target = 5, value.only = TRUE)

# Convert to a directed graph, use stmincut ------------------------------------------------------------------
## https://stackoverflow.com/questions/29375138/calculating-minimum-s-t-cuts-is-not-implemented-yet-in-igraph

dg <- as.directed(g)
# diameter(dg)
# get_diameter(dg)

min_cut(dg, value.only = FALSE)

# st_min_cuts(dg, source=2, target=5)
# cut <- st_min_cuts(dg, source=2, target=5)

st_min_cuts(dg, source = "BRF1", target = "EGLN1")
cut <- st_min_cuts(dg, source = "BRF1", target = "EGLN1")

# Number of edges removed (unweighted)
cut$value
# The minus side
E(dg)[cut$cuts[[1]]]
## Vertices in the first partition
V(dg)[cut$partition1s[[1]]]
## If there are multiple cuts
cut$cuts[[2]] 
cut$partition1s[[2]]


##  Minimum size vertex separators ------------------------------------------
min_separators(g)


## Assembling matrices and vectors -----------------------------------------------------------------
library(igraph)

gene <- read.csv('TCGA_NEW/UNgene_component.csv',header = T)
net <- read.csv('TCGA_NEW/UNgene_comp_net.csv')
g <- graph_from_data_frame(net, directed=F)
node <- as.matrix(get_diameter(g))
lab1 <- as.matrix(rownames(node))
colnames(lab1) <- c("node")

# 2020.7.20 Output 1 vector, the vector is the coefficient of the fixed point and the cut point ----------------------------------------------------
my_vector <- function(gene, lab1){
  k <- length(as.matrix(gene))
  vec1 <- matrix(0,1,k)
  l <- length(lab1)
  for (i in 2:l-1) {
    # i <- 4
    # which(gene[,1]==lab1[i])
    vec1[which(gene[,1]==lab1[i])] <- c("-1")
    # View(t(vec1))
  }
  ## This line appears before the loop: The first endpoint is -1
  x_first <- which(gene[,1]==lab1[1])
  x_end <- which(gene[,1]==lab1[l])
  vec1[,c(x_first,x_end)] <- c("1")
  return(vec1)
}


vec2 <- my_vector(gene, lab1)
vec3 <- t(vec2)
View(vec3)
# write.table(vec3,file = "TCGA_NEW/cut_vector_UNgene.txt", quote=F, sep="\t", row.names = F)
node
get_diameter(g) 
gene[c(40,455),]


##############################################################################
##############################################################################
# Adj matrix --------------------------------------------------------------
rm(list = ls())

## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- 'lly/home/R'
setwd(paste(path, 'CNetCox/Data/', sep=''))

## load data
gene <- read.csv('TCGA_NEW/UNgene_component.csv')
genes <- data.frame(gene[,1])
net <- read.csv('TCGA_NEW/UNgene_comp_net.csv')

G <- graph_from_data_frame(net, directed=F, vertices=genes)
print(G, e=TRUE, v=TRUE)
# plot(G)


# Convert graph to adjacency matrix ---------------------------------------------------------------
# adj <- as_adjcaency_matrix(G,sparse=FALSE)  # Same function as get.adjacency
adj <- get.adjacency(G,sparse=FALSE) 
View(adj[1:10,1:10])
# write.csv(adj, 'adjmatrix_comp_UNG.csv')


# Laplacian Matrix ------------------------------------------------------------------
# Non-Normalized Laplacian Matrix from adjacency matrix
Non.NormalizedLaplacianMatrix = function(adj){
  diag(adj) <- 0
  deg <- apply(adj,1,sum)
  D = diag(deg)
  L = D - adj             
  return(L)
}

L <- Non.NormalizedLaplacianMatrix(adj)

# Eigenvalue ---------------------------------------------------------------------
a.e <- eigen(L,symmetric=T)

Vector <- a.e$vectors
eigvalue <- a.e$values
# View(a.e$vectors)
a.e$vectors
# View(a.e$values)
# write.csv(Vector, 'TCGA_NEW/Vector_R.csv')
# write.csv(eigvalue, 'TCGA_NEW/eigvalue_R.csv')


# Normalized Laplacian matrix --------------------------------------------------------------

# Normalized Laplacian Matrix from adjacency matrix
laplacianMatrix = function(adj){
  diag(adj) <- 0                    
  # Degree matrix element (diagonal)--the sum of the absolute values ​​of each row of the adjacency matrix
  deg <- apply(abs(adj),1,sum)     
  p <- ncol(adj)
  L <- matrix(0,p,p)                
  nonzero <- which(deg!=0)          
  for (i in nonzero){
    for (j in nonzero){
      L[i,j] <- -adj[i,j]/sqrt(deg[i]*deg[j])   
    }
  }
  diag(L) <- 1                                  
  return(L)
}

L_norm <- laplacianMatrix(adj)
# View(L_norm[1:20,1:20])
# Eigenvalue ---------------------------------------------------------------------
a.e <- eigen(L_norm,symmetric=T)
Vector <- a.e$vectors
eigvalue <- a.e$values
# View(a.e$vectors)
a.e$vectors
View(a.e$values)
# write.csv(Vector, 'TCGA_NEW/adj_vector_norm.csv')
# write.csv(eigvalue, 'TCGA_NEW/adj_eigvalue_norm.csv')