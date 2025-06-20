

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
setwd(paste(path, 'Paper/Paper13Cox/Metascape/CNet/all.tk57fofgx/Enrichment_GO/', sep=''))

## load go list
b <- read.csv("GO_AllLists.csv",header = T, sep = ",")
b1 <- as.matrix(b[which(b$Category == "GO Biological Processes"),3])

a <- b1[1:23,]
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
## ??ÿһ??????ֵ
boxdata1 <- as.data.frame(cbind(rowMeans(coef2), boxdata[,2]))
colnames(boxdata1) <- c("SSmeasure", "Type")
rownames(boxdata1) <- rownames(boxdata)
boxdata1$SSmeasure <- as.numeric(as.vector(boxdata1$SSmeasure))

# write.csv(boxdata1, "TCGA\\687_37\\result\\boxdata_cluster.csv", row.names = F)
class(boxdata[2,1])
class(boxdata1[2,1])

# ??????ͼ --------------------------------------------------------------------
library(ggpubr)
library(ggplot2)
# my_comparisons <- list( c("Selected", "Random") )
# p <- ggboxplot(boxdata1, x = "Type", y = "SSmeasure",
#           color = "Type", palette = "lancet", # "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty"
#           add = "jitter", fill = c("#90CAF9", "#F48FB1")) + 
#   stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
#   stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "#C5E1A5")# Add global p-value
# p


# Set theme and colors
theme_set(theme_cowplot())
colors <- c("#E69F00", "#56B4E9")

ggplot(boxdata1, aes(x = Type, y = SSmeasure, fill = Type)) +
  geom_boxplot(alpha = 0.8, color = "#404040") +
  scale_fill_manual(values = colors) +
  labs(title = "", x = "", y = "SS-measure") +
  theme(plot.title = element_text(size = rel(1.2)), # , face = "bold"
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.position = "none") +
  stat_compare_means(label = "p.format", method = "t.test", size = 6, label.x = 1.5, label.y = max(boxdata1$SSmeasure)) 


# ???? ----------------------------------------------------------------------

# pdf(file = "TCGA_NEW\\result\\box_sim_cluster_NEW.pdf",width = 4.5,height = 4.5)
# p
# dev.off()

# ??ѡcluster GO ---------------------------------------------------------------

pr <- 0.7
which(dasimr[,1] >= pr)

setwd("D:\\E\\??ʿ\\DR_paper\\Paper4\\coxͼ?ͱ?\\???ܷ???")
bb0 <- as.matrix(read.csv("cluster_GO_NEW.csv",header = T)[,c(1,2)]) 
aa0 <- as.matrix(bb0[c(1,2,3,12,13,14,15,17,19,20,22,23,24,25,26,28,30,34,35,37,38,46,48),])
aa

goselect <- aa0[which(dasimr[,1] >= pr),]
goselect[,2]
# write.csv(goselect, "go_6_select_NEW.csv", row.names = F)



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
