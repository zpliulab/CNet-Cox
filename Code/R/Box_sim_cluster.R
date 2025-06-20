

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("GOSemSim")

## 2023.5.12 make ss-measure for CNet-RCPH 68 markers


## clear
rm(list = ls())

## set pathway
path <- '/Users/lilingyu/E/PhD/'
# path <- '/home/lly/R/'

# Creat Files ------------------------------------------------------------------
# setwd(paste(path, 'CNetCox/Data/', sep=''))


library(GOSim)
library(GOSemSim)


setOntology(ont = "BP", loadIC=TRUE, DIR=NULL)
H <-getAncestors()
HH <-getOffsprings()
HHH <-unlist(HH)
d <- godata('org.Hs.eg.db', ont="BP", computeIC=FALSE)

# Creat Files ------------------------------------------------------------------
setwd(paste(path, 'R/CNetCox/Result/cluster/', sep=''))

## load go list
set <- c(1,2,9:16, 19:23, 26,27, 29:31, 37, 43, 47)
b <- as.matrix(read.csv("cluster_GO_cluster.csv",header = T, sep = ",")[set,1])

a <- b
## ??ѡ23??go
# a <- b[c(1,2,3,12,13,14,15,17,19,20,22,23,24,25,26,28,30,34,35,37,38,46,48),]
# View(a)
# write.csv(a, "cluster_GO_NEW_select23.csv", row.names = F)

# setwd("D:\\E\\??ʿ\\DR_paper\\Paper4\\coxͼ?ͱ?\\???ܷ???")
# a <-as.matrix(read.csv("NOA.csv",header = T))


setwd(paste(path, 'R/BRCA/Data', sep=''))
aa <-as.matrix(read.csv("bc shishi.csv"))[,2]
amlsim <- as.matrix(mgoSim(a,aa,semData=d,measure = "Wang",combine=NULL))

# setwd(paste(path, 'CNetCox/Result/cluster', sep=''))
# write.csv(amlsim,"TCGA\\687_37\\result\\consimi.csv", row.names = F)



# ?ҵ? go terms --------------------------------------------------------------------

library(corrplot)
corrplot(amlsim)  #[c(1:10),]

dasimr <- as.matrix(apply(amlsim,1,max))
dasimc <- as.matrix(apply(amlsim,2,max))
congoterm <-as.matrix(rbind(dasimr, dasimc))

l <- length(dasimr)
n <- length(dasimr) + length(dasimc)

type1 <- as.matrix(rep("Selected",n))
congoterm1 <- data.frame(cbind(congoterm,type1))
colnames(congoterm1) <- c("sim","type")
congoterm1$sim <-as.numeric(as.vector(congoterm1$sim))
congoterm1hou <-as.matrix(congoterm1)

dasimrmean <-mean(dasimr)
dasimcmean <-mean(dasimc)



# selected go terms -------------------------------------------------------------

Boxdata  <- matrix() 
set <- c(1,2,4,5,6,7,9,11,14,19)
length(set)


for(i in set){
  
  # i = 19
  set.seed(i)
  # set.seed(9531*i)

  # set.seed(1)
  rand <- as.matrix(sample(HHH, size=l, replace = FALSE))
  row.names(rand) <- NULL
  is.na(rand)
  randsim <-as.matrix(mgoSim(rand,aa,semData=d,measure = "Wang",combine=NULL))
  dim(randsim)
  
  randsimr <- as.matrix(apply(randsim,1,max))
  randsimc <- as.matrix(apply(randsim,2,max))
  randomterm <- as.matrix(rbind(randsimr,randsimc))
  type2 <- as.matrix(rep("Random",n))
  randomterm1 <- data.frame(cbind(randomterm,type2))
  randomterm1$X1 <- as.numeric(as.vector(randomterm1$X1))
  randomterm1hou <- as.matrix(randomterm1)
  
  
  # һ?µ?go term???????????? --------------------------------------------------------
  
  boxdata <- data.frame(rbind(congoterm1hou,randomterm1hou))
  boxdata$sim <- as.numeric(as.vector(boxdata$sim))
  
  randsimrmean <- mean(randsimr)
  randsimcmean <- mean(randsimc)
  
  
  tcross <- rep(i, length(boxdata))                # i?ǵڼ???ѭ?????棬??K??
  step <- data.frame(cbind(boxdata, tcross))
  Boxdata <- cbind(Boxdata, step)                  #temp???к?pred?ϲ?
  
  
  # kcross <- rep(i, length(randsimrmean)) 
  # temp <- data.frame(cbind(randsimrmean, kcross))
  # Ci_ridge <- cbind(Ci_ridge, temp)   #temp???к?pred?ϲ?
  # 
  # 
  # pcross <- rep(i, length(p_ridge)) 
  # pemp <- data.frame(cbind(p_ridge, pcross))
  # P_ridge <- cbind(P_ridge, pemp)   #temp???к?pred?ϲ?
  
  
  print(paste(i)) 
  
}


# my_cbind ??ȡ20?ε?coef -----------------------------------------------------

my_cbind <- function(x){
  x <- Boxdata
  x1 <- matrix()
  x1 <- as.character(x[,3])
  for (i in 0:9){
    x1 <- as.matrix(cbind(x1, x[, 2+3*i]))
  }
  return(x1)
}

coef <- my_cbind(Boxdata)
coef2 <- apply(coef[,-1], 2, as.numeric)

class(boxdata)

boxdata1 <- as.data.frame(cbind(rowMeans(coef2), boxdata[,2]))
colnames(boxdata1) <- c("SSmeasure", "Type")
rownames(boxdata1) <- rownames(boxdata)
boxdata1$SSmeasure <- as.numeric(as.vector(boxdata1$SSmeasure))
class(boxdata[2,1])
class(boxdata1[2,1])

# ??????ͼ --------------------------------------------------------------------
library(ggpubr)
library(ggplot2)

# Set theme and colors
theme_set(theme_cowplot())
colors <- c("#E69F00", "#56B4E9")

P1 <- ggplot(boxdata1, aes(x = Type, y = SSmeasure, fill = Type)) +
  geom_boxplot(alpha = 0.8, color = "#404040") +
  scale_fill_manual(values = colors) +
  labs(title = "", x = "", y = "SS-measure") +
  theme_bw() +
  coord_flip() +
  # theme_classic() +
  # theme_minimal() +
  theme(plot.title = element_text(size = rel(1.2)), # , face = "bold"
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.position = "none") +
  stat_compare_means(label = "p.format", method = "t.test", 
                     size = 6, label.x = 1.5, angle = -90, 
                     label.y = max(boxdata1$SSmeasure)) 


P1

# ???? ----------------------------------------------------------------------

# pdf(file = "TCGA_NEW\\result\\box_sim_cluster_NEW.pdf",width = 4.5,height = 4.5)
# p
# dev.off()

# ??ѡcluster GO ---------------------------------------------------------------

pr <- 0.6
which(dasimr[,1] >= pr)


setwd(paste(path, 'R/CNetCox/Result/cluster/', sep=''))

## load go list
set <- c(1,2,9:16, 19:23, 26,27, 29:31, 37, 43, 47)
aa0 <- as.matrix(read.csv("cluster_GO_cluster.csv",header = T, sep = ",")[set,c(1,2)])
goselect <- aa0[which(dasimr[,1] >= pr),]
goselect[,2]
# write.csv(goselect, "go_6_select_NEW.csv", row.names = F)




# functional enrichment analysis ------------------------------------------


## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- '/home/lly/R/'

# Creat Files ------------------------------------------------------------------
setwd(paste(path, 'CNetCox/Data/', sep=''))

# load marker genes
gene <- as.matrix(read.csv("Result/Result3/marker3.csv", header = T, sep=','))
colnames(gene) <- c('gene')

# load R package --------------------------------------------------------------------
library(org.Hs.eg.db)
library(clusterProfiler)


# load genes symbol ----------------------------------------------------------

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
# setwd(paste(path, 'CNetCox/Result/cluster', sep=''))
# write.csv(go_BP@result, file = "cluster_GO_cluster.csv", row.names = F)

# plot bar and plot pot  ----------------------------------------------------------------

barplot(go_BP, y =go_BP@result$Description, showCategory=23,drop=T)    
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


GO23 <- as.matrix(read.csv("cluster_GO_cluster.csv",header = T, sep = ",")[set,c(1,2,5,6)])
dim(GO23)
library(stargazer)
# options(scipen = 3)
GO23p <- format(as.numeric(GO23[,3]), scientific = T)
GO23[,3] <-  format(as.numeric(GO23[,3]), scientific = T) 
GO23[,4] <-  format(as.numeric(GO23[,4]), scientific = T) 
options(digits = 3)
df0 <- GO23 
stargazer(df0)
stargazer(GO23p)

# 然后自己计算Fold Enrichment，并按照Fold Enrichment升序排序
library(stringr)
gr1 <- as.numeric(str_split(go_BPsub$GeneRatio,"/",simplify = T)[,1])
gr2 <- as.numeric(str_split(go_BPsub$GeneRatio,"/",simplify = T)[,2])
bg1 <- as.numeric(str_split(go_BPsub$BgRatio,"/",simplify = T)[,1])
bg2 <- as.numeric(str_split(go_BPsub$BgRatio,"/",simplify = T)[,2])
go_BPsub$fold <- (gr1/gr2)/(bg1/bg2)
go_BPsub <- arrange(go_BPsub,fold)

## rank description
go_BPsub$Description = factor(go_BPsub$Description,
                              levels = go_BPsub$Description,
                              ordered = F)

P2 <- ggplot(go_BPsub, aes(x = fold,y = Description)) +
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
  theme(panel.border = element_rect(size=1, fill = 'transparent'),
        axis.text = element_text(color='black'),
        axis.text.x = element_text(color='black',  
                                   hjust = 0.5, vjust = 0.5, size = 11)) 

P2



# 基因—通路图 	 ----------------------------------------------------

library(ggnewscale)
cnetplot(go_BP)

library(enrichplot)
library(DOSE)

x2 <- go_BP
cnetplot(x2)
# use `layout` to change the layout of map
cnetplot(x2, layout = "star")
# use `showCategory` to select the displayed terms. It can be a number of a vector of terms.
cnetplot(x2, showCategory = 10)

## selected term
categorys <- c("canonical Wnt signaling pathway", "gland development","protein kinase B signaling" )
P3 <-  cnetplot(go_BP,
                  showCategory = categorys,
                  foldChange = NULL,
                  layout = "kk",    # kk, gem
                  colorEdge = T,
                  circular = F,
                  node_label = "all",  ## gene
                  cex_category = 1.0,
                  cex_gene = 1.0,
                  cex_label_category = 1.0,
                  cex_label_gene = 0.6,
                  shadowtext = "none")

P3 


# KEGG pathway ------------------------------------------------------------

enrichKK <- enrichKEGG(gene         =  geneList,
                       organism     = 'hsa',
                       keyType = 'kegg', 
                       #universe     = gene_all,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = 'BH', 
                       minGSSize = 3,
                       maxGSSize = 500,
                       qvalueCutoff = 0.05,
                       use_internal_data = FALSE)


head(enrichKK)[,1:6] 
dotplot(enrichKK)

library(ggplot2)
library(stringr)
P4 <- dotplot(enrichKK,showCategory=5) + 
  scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))

P4

browseKEGG(enrichKK,'hsa05224') 
barplot(enrichKK,showCategory=10)
enrichKK@result$Description[1:10]

#Gene-Concept Network 
# install.packages("ggnewscale")
library(ggnewscale)
cnetplot(enrichKK)


library(patchwork)
P1 + P3






P1 + P3 + plot_layout(nrow = 1, widths = c(1.5, 3)) 
#不想图片都是一样大的，使用plot_layout来设置高度/宽度的分配
#这边为a分配了四分之三的高，而b只有四分之一。

(P1 + P3 + plot_layout(nrow = 1, widths = c(1.5, 3)) )

PP <- (P2 | ((P1/P3) + plot_layout(height=c(2,3)))) + 
  plot_layout(widths = c(2, 4))

PP
ggsave(filename="beautiful.pdf", plot=PP, width = 14, height=8, units="in")







 # ??ѡGO term ---------------------------------------------------------------

which(dasimr[,1] >= pr)
go <- as.matrix(read.csv("TCGA\\687_37\\result\\chart_96AD.csv",header = T)[-c(10,13,25),])

library(tidyverse)
goterm <- str_split_fixed(go[,2], "[~]", 22)

goterm1 <- as.matrix(cbind(goterm[,1], goterm[,2], go[,c(3,5,6)]))
colnames(goterm1) <- c("Term", "Describe", colnames(go)[c(3,5,6)])

goselect <- goterm1[which(dasimr[,1] >= pr),c(1,2)]
goselect[,2]


library(stargazer)
stargazer(goselect, summary=FALSE, rownames=FALSE) #????stargazer????
# write.csv(goselect, "TCGA\\687_37\\result\\chart_5goselect.csv", row.names = F)
