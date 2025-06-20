## 2023.5.12 make functional analysis for CNet-RCPH 68 markers


## clear
rm(list = ls())

## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- '/home/lly/R/'

# Creat Files ------------------------------------------------------------------
setwd(paste(path, 'CNetCox/Data/', sep=''))

# load marker genes
gene <- as.matrix(read.csv("Result/Result3/marker3.csv", header = T, sep=','))
colnames(gene) <- c('gene')


# R packages --------------------------------------------------------------------


# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::available()
# 
# library(BiocManager)

# BiocManager::install(version = "3.9")   #228???????°?װ
# BiocManager::install("clusterProfiler")


# install.packages("BiocInstaller",
#                  repos="https://bioconductor.org/packages/3.8/bioc")
# install.packages("BiocInstaller",
#                  repos="http://bioconductor.org/packages/3.8/bioc")

# 
# if(!require(devtools)) install.packages("devtools")
# devtools::install_github('GuangchuangYu/clusterProfiler')


# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("clusterProfiler")


# load R package --------------------------------------------------------------------
library(org.Hs.eg.db)
library(clusterProfiler)


# ????genes symbol ----------------------------------------------------------

genelist <- as.character(gene)
eg <- bitr(genelist, fromType="SYMBOL", toType=c("ENTREZID","GENENAME"), OrgDb="org.Hs.eg.db"); 
head(eg)


# go ----------------------------------------------------------------------

geneList <- eg$ENTREZID
## 17??gene,????У??
go_BP <- enrichGO(gene = geneList,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
head(go_BP)[,1:6]
# View(go_BP[1:50,])
class(go_BP)

View(go_BP@result)
setwd(paste(path, 'CNetCox/Result/cluster', sep=''))
# write.csv(go_BP@result, file = "cluster_GO_cluster.csv", row.names = F)

# ???м򵥵Ŀ??ӻ? ----------------------------------------------------------------

barplot(go_BP, y =go_BP@result$Description, showCategory=23,drop=T)    # showCategory=10,title="Enrichment GO"
go_BP@result$Description
library(ggplot2)
library(stringr)
dotplot(go_BP, font.size=12, showCategory=23,title="Enrichment GO Top 23") + 
  scale_size(rang=c(5.20)) + 
  scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))



# 选择要显示的术语 ----------------------------------------------------------------
setwd('/Users/lilingyu/E/PhD/R/CNetCox/Result/cluster/')
## load go list
set <- c(1,2,9:16, 19:23, 26,27, 29:31, 37, 43, 47)
aa0 <- read.csv("cluster_GO_cluster.csv",header = T, sep = ",")[set,c(1,2)]
selected_terms <- aa0$ID

# 根据选定的术语对数据进行子集化
go_BPsub <- go_BP[go_BP@result$ID %in% selected_terms,]
# 然后自己计算Fold Enrichment，并按照Fold Enrichment升序排序
library(stringr)
gr1 <- as.numeric(str_split(go_BPsub$GeneRatio,"/",simplify = T)[,1])
gr2 <- as.numeric(str_split(go_BPsub$GeneRatio,"/",simplify = T)[,2])
bg1 <- as.numeric(str_split(go_BPsub$BgRatio,"/",simplify = T)[,1])
bg2 <- as.numeric(str_split(go_BPsub$BgRatio,"/",simplify = T)[,2])
go_BPsub$fold <- (gr1/gr2)/(bg1/bg2)
go_BPsub <- arrange(go_BPsub,fold)


# 使用 dotplot() 绘制图形
dotplot(go_BPsub)
go_BPsub$GeneRatio


go_BPsub$Description = factor(go_BPsub$Description,
                              levels = go_BPsub$Description,
                              ordered = F)

ggplot(go_BPsub, aes(x = fold,y = Description)) +
  geom_point(aes(color = p.adjust,
                 size = Count)) +
  scale_color_gradient(low = "red", high = "blue") +
  xlab("Fold Enrichment") +
  theme_bw() +
  #edit legends
  guides(
    #reverse color order (higher value on top)
    color = guide_colorbar(reverse = TRUE)) +
#reverse size order (higher diameter on top) 
#size = guide_legend(reverse = TRUE))
  # theme_test(base_size = 12) +
  theme(panel.border = element_rect(size=1,fill = 'transparent'),
        axis.text = element_text(color='black'),
        axis.text.x = element_text(color='black',  
                                   hjust = 0.5, vjust = 0.5, size = 11)) 


# plotGOgraph(go_BP) 	#GOͼ?????????????Գ??????Ͻ?????Ϊpdf
library(ggnewscale)
cnetplot(go_BP)

bb <- go_BP@result[set,]



cnetplot(go_BP,
        showCategory = 5,
        foldChange = NULL,
        layout = "kk",
        colorEdge = T,
        circular = T,
        node_label = "all"
)
  

enrichGO = DOSE::setReadable(go_BP, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
enrichGO
library(ggnewscale)
cnetplot(enrichGO)



# KEGG --------------------------------------------------------------------

enrichKK <- enrichKEGG(gene         =  geneList,
                       organism     = 'hsa',
                       #universe     = gene_all,
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
head(enrichKK)[,1:6] 
dotplot(enrichKK)
# enrichKK = DOSE::setReadable(enrichKK, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
# enrichKK 

barplot(enrichKK,showCategory=20)
dotplot(enrichKK)

#(3)չʾtop5ͨ·?Ĺ?ͬ??????Ҫ?Ŵ󿴡?
#Gene-Concept Network 
# install.packages("ggnewscale")
library(ggnewscale)
cnetplot(enrichKK)



# ???? ----------------------------------------------------------------------

# go_BP_result <- as.data.frame(go_BP@result)    # ?????ɿɶ?
# p_value <- which(go_BP_result$pvalue < 0.05)   # ?ҳ?p<0.05??go
# go_BP_pvalue <- go_BP_result[p_value,]
# write.csv(as.data.frame(go_BP@result), file="go_BP_result_pvalue.csv")  # ????p<0.05??go

## ????????????????
# go_BP <- enrichGO(OrgDb="org.Hs.eg.db",
#                    gene = geneList,
#                    pvalueCutoff = 0.05,
#                    ont = "BP",
#                    readable=TRUE)
# head(go_BP)



# ?????Ի???GO????????ϵͼ??????ֵ??ע????????????????ֻ???Ǹ???һ??GOͨ·??BP??CC??MF???????? ---------------------


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rgraphviz")

# install.packages('Rgraphviz')

library(topGO)
library(Rgraphviz)


go.BP <- enrichGO(geneList, 
                  OrgDb = org.Hs.eg.db, 
                  ont='BP',
                  pAdjustMethod = 'BH', 
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.05, 
                  keyType = 'ENTREZID')
plotGOgraph(go.BP)



# 123 ---------------------------------------------------------------------

library(enrichplot)
library(DOSE)
# data(geneList)
# de <- names(geneList)[1:100]
# x <- enrichDO(de)
# x2 <- pairwise_termsim(x)

x2 <- go_BP

cnetplot(x2)
# use `layout` to change the layout of map
cnetplot(x2, layout = "star")
# use `showCategory` to select the displayed terms. It can be a number of a vector of terms.
cnetplot(x2, showCategory = 10)
# categorys <- c("pre-malignant neoplasm", "intestinal disease",
#                "breast ductal carcinoma", "non-small cell lung carcinoma")

categorys <- c("cell fate commitment", "canonical Wnt signaling pathway",
               "cardiac chamber development", "mitral valve morphogenesis",
               "Notch signaling pathway", "regulation of Wnt signaling pathway")

# cnetplot(x2, showCategory = categorys,
#          colorEdge = T,
#          circular = F,
#          node_label = "all"
# )

plot <-  cnetplot(go_BP,
         showCategory = categorys,
         foldChange = NULL,
         layout = "kk",    # kk, gem
         colorEdge = T,
         circular = F,
         node_label = "gene",
         cex_category = 1.0,
         cex_gene = 1.0,
         cex_label_category = 1.0,
         cex_label_gene = 0.6,
         shadowtext = "none")

# setwd("D:\\E\\??ʿ\\DR_paper\\Paper4\\coxͼ?ͱ?\\???ܷ???")
# pdf(file = "cnetplot_cluster_NEW.pdf", width = 6.5, height = 4.5)
# plot
# dev.off()


## ??star??, ??circle??, ??gem??, ??dh??, ??graphopt??, ??grid??, ??mds??, ??randomly??, ??fr??, ??kk??, ??drl??, ??lgl??

# It can also graph compareClusterResult
data(gcSample)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
data(gcSample)
# gcSample <- go_BP
xx <- compareCluster(gcSample, fun="enrichGO", OrgDb="org.Hs.eg.db")
xx2 <- pairwise_termsim(xx)
cnetplot(xx2)


