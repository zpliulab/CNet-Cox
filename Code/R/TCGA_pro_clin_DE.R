## 2023.4.20 LLY made, 2026.1.14 LLY used
## This code is used to select DEGs (surviving vs deceased)


##############################################################################
##############################################################################
# Data download: gene expression data + clinical information -------------------

## clear
rm(list = ls())

## set pathway
# path <- 'lly/home/R'
path <- '/Users/lilingyu/E/PhD/R/'
setwd(paste(path, 'CNetCox/Data/TCGA_NEW/', sep=''))

## load R package 
library(glmSparseNet)
library(curatedTCGAData)
library(TCGAutils)
library(dplyr)

# Extract expr matrix of TNBC from the BRCA in TCGA ----------------------------
brca <- curatedTCGAData(diseaseCode = "BRCA", assays = "RNASeq2GeneNorm",
                        version = "1.1.38", dry.run = FALSE)
head(getSubtypeMap(brca))
head(getClinicalNames("BRCA"))
# [1] "years_to_birth"  "vital_status"  "days_to_death"  "days_to_last_followup"
# [5] "tumor_tissue_site"  "pathologic_stage" 

# brca.primary.solid.tumor <- TCGAutils::splitAssays(brca, '01')
brca.primary.solid.tumor <- TCGAutils::TCGAsplitAssays(brca, '01')    #2026.1.14
xdata.raw <- t(assay(brca.primary.solid.tumor[[1]]))
dim(xdata.raw)    # 1093 20501
# View(xdata.raw[,1:10])

## keep features with standard deviation > 0
xdata.raw <- xdata.raw %>%
  { (apply(., 2, sd) != 0) } %>%
  { xdata.raw[, .] }
dim(xdata.raw)    # 1093 20220
# View(xdata.raw[,1:10])

## Get survival information
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
dim(ydata.raw)    # 1080    3

## Set index as the patientID
rownames(ydata.raw) <- ydata.raw$patientID    
## Get matches between survival and assay data
xdata.raw_1 <- xdata.raw[TCGAbarcode(rownames(xdata.raw)) %in%
                           rownames(ydata.raw),]
dim(xdata.raw_1)    # 1080 20220
# View(xdata.raw_1[,1:10])
## Order ydata the same as assay
ydata.raw    <- ydata.raw[TCGAbarcode(rownames(xdata.raw_1)), ]
# View(ydata.raw)

xdata <- xdata.raw_1
ydata <- ydata.raw %>% select(status)
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
# DEG analysis based on survival status and normalize. --------------------------

## set pathway
# path <- 'lly/home/R'
path <- '/Users/lilingyu/E/PhD/R/'
setwd(paste(path, 'CNetCox/Data/', sep=''))

## load data
data_outcome <- read.table("TCGA_NEW/TCGA_pro_outcome.txt",header=T,sep='\t', check.names = F)
dim(data_outcome)   # 20221  1080
# View(data_outcome[,1:10])
data_outcome[2,1]

## expression data
xdatat <- data_outcome[-1,]
dim(xdatat)   # 20220  1080
# View(xdatat[,1:10])

## survival status
ydata <- t(data_outcome[1,])
sum(ydata[,1])    # 152 dead

## label and group
ydata[which(ydata[,1] == "1"), 1] <- c("Dead")
ydata[which(ydata[,1] == "0"), 1] <- c("Alive")


# DESeq2 makes DEG -------------------------------------------------
exprSet <- round(xdatat)
dim(exprSet)    # 20220  1080
# View(exprSet[,1:10])
group_list <- as.factor(ydata)
colData <- data.frame(row.names=colnames(exprSet), group_list=group_list)
colnames(colData) <- c("outcome")
dim(colData)    # 1080    1

# BiocManager::install("DESeq2")
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ outcome)
dds2 <- DESeq(dds)  
resultsNames(dds2)

## Extract the DGE results (Dead vs Alive )
res <-  results(dds2, contrast=c("outcome", "Dead", "Alive"))
summary(res) 
# plotMA(res)

res_order <- res[order(res$padj),]
res_order <- as.data.frame(res_order)
# write.csv(res_order,file= "DEG_res_order_DA.csv")

## From Zishuang Zhang, both calculate the same number
res1 <-  results(dds2, alpha = 0.01)
# write.csv(res1,file= "DEG_res.csv")

# diff_gene_deseq2 <- subset(res1, padj < 0.01)    # 501
diff_gene_deseq2 <- subset(res1, padj < 0.01 & abs(log2FoldChange) > log2(2))
diff_gene_deseq2 <- as.data.frame(diff_gene_deseq2)
dim(diff_gene_deseq2)    #196    6
# write.csv(diff_gene_deseq2,file= "DEG_Dead_vs_Alive_1.csv")


##############################################################################
# Fig. 2A
##############################################################################
# Volcano plot ------------------------------------------------------------
logFC <- res$log2FoldChange    # LogFC is Log2FC 
pvalue <- res$padj    # pvalue is adjusted p value
genes <- rownames(res)
data <- data.frame(logFC, pvalue, genes)
# write.csv(data, "TCGA_NEW/Volcano_DE.csv", row.names = F)

library(ggplot2)
## Define colors and significance threshold
color_significant <- "#e7897d"
color_nonsignificant <- "#9fa0b5"
sig_threshold <- -log10(0.01)
# sig_threshold <- -log2(10)

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


# Normalized -------------------------------------------------------------------
vst_dat <- vst(dds2, blind = TRUE)
dat111 <- assay(vst_dat)
dim(dat111)    # 20220  1080
# View(dat111[,1:10])

data_norm_clin <- t(dat111)
dim(data_norm_clin)    # 1080 20220
# View(data_norm_clin[,1:10])

# with survival time -----------------------------------------------------------
ydata_clin <- ydata.raw %>% select(time, status)
# View(xdata[,1:10])
# View(ydata_clin)

data_clin <- as.matrix(cbind(ydata_clin, data_norm_clin))
rownames(data_clin) <- rownames(data_norm_clin) 
dim(data_clin)    # 1080 20222
View(data_clin[,1:10])

data_clin1 <- t(data_clin)
# write.table(data_clin1,"TCGA_pro_norm_clin.txt",quote=F,sep="\t") 


##############################################################################
##############################################################################
# Integrate DEGs with prior knowledges -----------------------------------------

## clear
rm(list = ls())

## set pathway
# path <- 'lly/home/R'
path <- '/Users/lilingyu/E/PhD/R/'
setwd(paste(path, 'CNetCox/Data/', sep=''))

## load data
Data = read.table("TCGA_NEW/TCGA_pro_norm_clin.txt", header = T, check.names = FALSE)
dim(Data)    # 20222  1080
gene <- as.vector(rownames(Data)[-c(1,2)])

DE_clin <- read.csv("TCGA_NEW/DEG_Dead_vs_Alive_1.csv", header = T, sep=',')    # 196
Degene <- as.vector(DE_clin[,1])    # 196

## see non DEGs
sub <- setdiff(gene, Degene)   


# Prior biomarker --------------------------------------------------------------
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

dim(as.matrix(which(marker %in% as.matrix(Degene))))    # 0
marker_DE <- as.matrix(marker[which(marker %in% Degene),])
colnames(marker_DE) <- "DE"
# write.csv(marker_DE,"Result/Result3/marker_DE.csv", row.names = F)

DE <- as.matrix(read.csv("TCGA_NEW/DEG_Dead_vs_Alive_1.csv", header = T,  sep = ','))
dim(as.matrix(which(as.matrix(DE[,1]) %in% marker)))    # 0


# Gene integration 110+9+56+145+117 --------------------------------------------
add_gene <- as.matrix( union(union(union(union(intersect(sub, mark), 
                                                     intersect(sub, scmark)), 
                                               intersect(sub, mama_70)), 
                                         intersect(sub, KEGG_147)), 
                                  intersect(sub, tf)) )     # 420

allgene <- as.matrix(rbind(add_gene, as.matrix(Degene)))    # 616
dim(allgene)    # 616
rownames(allgene) <- allgene[,1] 


##############################################################################
##############################################################################
# Select component and node-cut set --------------------------------------------
net <- as.matrix(read.csv("Prior_infor/Regnetwork_hum.csv",header = T))
net[1,]
node_used <- allgene   
dim(node_used)    # 616
# net_used <- net[,c(2,4)]
net_used <- net
k1 <- which(net_used[,1] %in% node_used)   
k2 <- which(net_used[,2] %in% node_used)  
length(intersect(k1,k2))    # 5980
# length(union(k1,k2))
used <- net_used[intersect(k1,k2),]
dim(used)     # 5980    2

# install.packages('igraph')
library(igraph)
PP <- graph_from_data_frame(used,directed = F)
p1 <- simplify(PP, remove.loops = T, remove.multiple = T)  # final pairs
# ed <- as_edgelist(p1, names = TRUE)


# Compute maximum connected components -----------------------------------------
clu <- components(p1)
groups(clu)    # 544
  
## plot graph 
# g <- p1
# plot(g, layout=layout.fruchterman.reingold,  
#      vertex.size=4,   
#      vertex.label = V(g)$name,  
#      vertex.label.cex=0.7,   
#      vertex.label.dist=0.4,  
#      vertex.label.color = "black"   
# )

component <- groups(clu)$'1'
# write.csv(component, file = "TCGA_NEW/UNgene_component.csv", row.names = F)
node_used <- component
dim(as.matrix(node_used))    # 544
# net_used <- net[,c(2,4)]
k1 <- which(net_used[,1] %in% node_used)   
k2 <- which(net_used[,2] %in% node_used)     
length(intersect(k1,k2))    # 5980
used <- net_used[intersect(k1,k2),]

# install.packages('igraph')
library(igraph)
PP <- graph_from_data_frame(used,directed = F)
p1 <- simplify(PP, remove.loops = T, remove.multiple = T)  # final pairs
ed <- as_edgelist(p1, names = TRUE)
# write.csv(ed,"TCGA_NEW/UNgene_comp_net.csv",row.names = F,quote = F)


# Extract DEgene_component data ------------------------------------------------
data <- Data[component,]
dim(data)    # 544 1080
# View(data[,1:10])

all_data <- rbind(Data[c(1,2),], data)
dim(all_data)    # 546 1080
# View(all_data[,1:10])
# write.table(all_data, file = "TCGA_NEW/TCGA_BRCA_clin_546_1080.txt",quote = F, sep = "\t")


# scale ------------------------------------------------------------------------
all_data_scale <- rbind(all_data[c(1,2),], t(scale(t(all_data[-c(1,2),]))))
dim(all_data_scale)    # 546 1080
View(all_data_scale[,1:10])
# write.table(all_data_scale, file = "TCGA_NEW/TCGA_BRCA_clin_546_1080_scale.txt",quote = F, sep = "\t")


##############################################################################
##############################################################################
# Node cut set -----------------------------------------------------------------

## clear
rm(list = ls())

## set pathway
# path <- 'lly/home/R'
path <- '/Users/lilingyu/E/PhD/R/'
setwd(paste(path, 'CNetCox/Data/', sep=''))


library(igraph)
gene <- read.csv('TCGA_NEW/UNgene_component.csv')
genes <- data.frame(gene[,1])
net <- read.csv('TCGA_NEW/UNgene_comp_net.csv')


net[which(net[,1] == "EP300"),2] == "CTCFL"

g <- graph_from_data_frame(net, directed=F)


# Calculate two modes with the largest distance.--------------------------------
## diameter，breadth-first search
diameter(g) 

## TRUE, the diameters of the connected components
diameter(g, unconnected=TRUE)

## FALSE, the number of vertices
diameter(g, unconnected=FALSE)

## returns a path with the actual diameter
get_diameter(g)     # BRF1  EGLN1

## returns two vertex ids, connected by the diameter path.
farthest_vertices(g)     # 4
# diameter(g, weights=NA)

# Find the smallest cut  -------------------------------------------------------
## FALSE, the edges in the cut and a the two (or more) partitions are also returned.
min_cut(g, source = "BRF1", target = "EGLN1", value.only = FALSE)
# min_cut(g, source = 2, target = 5, value.only = FALSE)

## TRUE, only the minumum cut value is returned
min_cut(g, source = "BRF1", target = "EGLN1", value.only = TRUE)
# min_cut(g, source = 2, target = 5, value.only = TRUE)

# Convert to a directed graph using stmincut -----------------------------------
## https://stackoverflow.com/questions/29375138/calculating-minimum-s-t-cuts-is-not-implemented-yet-in-igraph

dg <- as.directed(g)
# diameter(dg)
# get_diameter(dg)

min_cut(dg, value.only = FALSE)

# st_min_cuts(dg, source=2, target=5)
# cut <- st_min_cuts(dg, source=2, target=5)

st_min_cuts(dg, source = "BRF1", target = "EGLN1")
cut <- st_min_cuts(dg, source = "BRF1", target = "EGLN1")

## Number of edges removed (unweighted)
cut$value
E(dg)[cut$cuts[[1]]]   # The side that was cut off
V(dg)[cut$partition1s[[1]]]    # Vertex in the first partition
cut$cuts[[2]] 
cut$partition1s[[2]]

# Minimum size vertex separators -----------------------------------------------
min_separators(g)

# Assemble matrices and vectors ------------------------------------------------
library(igraph)

gene <- read.csv('TCGA_NEW/UNgene_component.csv',header = T)
net <- read.csv('TCGA_NEW/UNgene_comp_net.csv')
g <- graph_from_data_frame(net, directed=F)
node <- as.matrix(get_diameter(g))
lab1 <- as.matrix(rownames(node))
colnames(lab1) <- c("node")

# 2020.7.20 Output vector, coefficients of fixed points and cut points ---------
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
  ## This line, placed before the loop, will result in the first endpoint being -1.
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
# Adj matrix -------------------------------------------------------------------

rm(list = ls())

## set pathway
# path <- 'lly/home/R'
path <- '/Users/lilingyu/E/PhD/R/'
setwd(paste(path, 'CNetCox/Data/', sep=''))

## load data
gene <- read.csv('TCGA_NEW/UNgene_component.csv')
genes <- data.frame(gene[,1])
net <- read.csv('TCGA_NEW/UNgene_comp_net.csv')

G <- graph_from_data_frame(net, directed=F, vertices=genes)
print(G, e=TRUE, v=TRUE)
# plot(G)


# Convert graph into an adjacency graph ----------------------------------------
adj <- as_adjcaency_matrix(G,sparse=FALSE)  
# adj <- get.adjacency(G,sparse=FALSE) 
View(adj[1:10,1:10])
# write.csv(adj, 'adjmatrix_comp_UNG.csv')


# laplacian -------------------------------------------------------------------
# Non-Normalized Laplacian Matrix from adjacency matrix
Non.NormalizedLaplacianMatrix = function(adj){
  diag(adj) <- 0
  deg <- apply(adj,1,sum)
  D = diag(deg)
  L = D - adj            
  return(L)
}
L <- Non.NormalizedLaplacianMatrix(adj)

# Eigenvalues ------------------------------------------------------------------
a.e <- eigen(L,symmetric=T)
Vector <- a.e$vectors
eigvalue <- a.e$values
# View(a.e$vectors)
a.e$vectors
# View(a.e$values)
# write.csv(Vector, 'TCGA_NEW/Vector_R.csv')
# write.csv(eigvalue, 'TCGA_NEW/eigvalue_R.csv')


# Normalized Laplace matrix ----------------------------------------------------
laplacianMatrix = function(adj){
  diag(adj) <- 0                    
  
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
a.e <- eigen(L_norm,symmetric=T)
Vector <- a.e$vectors
eigvalue <- a.e$values
# View(a.e$vectors)
a.e$vectors
View(a.e$values)
# write.csv(Vector, 'TCGA_NEW/adj_vector_norm.csv')
# write.csv(eigvalue, 'TCGA_NEW/adj_eigvalue_norm.csv')
