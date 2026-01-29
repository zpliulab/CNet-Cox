library(BiocGenerics)
library(parallel)
library(Biobase)
library(dplyr)       # ％>％ 管道函数的调用，传参
library(tidyr)
library(tidyverse)   # tibble 的调用
library(fdrtool)     # fdr校正
library(data.table)  # 用 fread 读入.soft 文件
library(fdrtool)     # 处理 p 值
library(data.table)  # 使用 fread 读入数据

setwd("D:\\E\\博士\\R_程序\\BRCA\\Bulk_RNA-seq_data\\GEO\\GSE20711")


# 读入数据 --------------------------------------------------------------------

## 读入基因表达数据
exprSet <- read.table("GSE20711_series_matrix.txt",header=T,sep='\t',fill=TRUE,strip.white = T)
# View(exprSet[,1:2])
dim(exprSet)    # 54675    91
exprSet$ID_REF <- as.character(exprSet$ID_REF)
# exprSet[1:10,1]
# exprSet$ID_REF[1:10]


## 读入芯片探针注释文件
anno <- read.table("GPL570.txt",header=T,sep='\t',fill=TRUE,strip.white = T,quote = "")  # 探针和基因ID的对应文件
# anno1 <- read.table("GPL570.annot",header=T,sep='\t',fill=TRUE,strip.white = T,quote = "")  # 探针和基因ID的对应文件
# （如上，其实行数一样）
# View(anno)

library(dplyr)
library(tidyr)
anno2 <- anno %>%     # GSE73685_family.soft 文件的内容 
  select(ID,ENTREZ_GENE_ID) %>%             # 提取文件中 ID,ENTREZ_GENE_ID 两列 
  filter(ENTREZ_GENE_ID != '') 
# anno2 <- anno[,c('ID','GeneID')] 
colnames(anno2) <- c('ID_REF','EntrezID')    # 将这两列的标签替换, 'GeneSymbol''Gene.Symbol'
anno2$ID_REF <- as.character(anno2$ID_REF)   # 将ID_REF列全部转换成符号,为了同expSet合并
# View(anno2)

## 将基因表达数据与芯片注释文件的探针名进行对应
exprset2 <- exprSet %>%                      # ％>％来自dplyr包的管道函数，
  inner_join(anno2,by='ID_REF') %>%        # 作用是将前一步的结果直接传参给下一步的函数，省略中间的赋值步骤
  select(ID_REF,EntrezID, everything())    # %>%# 重新排列
# View(exprset2[,1:2])                            # expset2 与 expset 相比，第2列多了基因号


## 整理芯片注释文件，把其中一个探针对应多个基因的拆分开
exprset3 <- exprset2
a <- tibble(exprset3[,1:2])
# a <- exprset3[,1:2]
# View(a)# 把第 1 和第 2 列提取出来，放在 a 中  

test1 <- apply(a,1, function(x){
  str_split(x[2],'///',simplify=T)       # 把 a 的第2列提取出，给test1
} )
# View(test1)

test2 <- apply(a, 1, function(x){          # 将探针号和基因号，进行---的链接
  paste(x[1],str_split(x[2],'///', simplify=T), sep = "---")
})
# View(test2)

unlist(test2)                              # 将 list 数据变成字符串向量或者数字向量的形式
# View(unlist(test2) )
x <- tibble(unlist(test2))                 # tibble，取代传统data.frame，读取并自动添加列名：unlist(test2)
colnames(x) <- "lala"                      # 改变 x 的列名：将 unlist(test2) 定义为 lala
# View(x)


x2 <- separate(x,lala,c("id","entrezID"),sep = '---')     # 识别 lala 中的 ---，将数据分离，单独成列并附新标签
# View(x2)
x2[1:10,1]
exprset3[1:10,1]
x3 <- merge(x2,exprset3,by.x = "id",by.y="ID_REF",all=FALSE)  #  将两个文件按顺序合并为一个，
# View(x3[,1:3]) 

x4<-x3[,-c(1,3)]                           # 将 第1 和第3 两列删除, 剩下的数据还是字符型的，带着“ "
View(x4[,1:3])
# dim(x4)
# x4[1,1]
x4[1:6,1]

zz <- as.matrix(apply(as.matrix(x4[,1]),1,function(x) as.numeric(x)))
View(zz)
# zz[1:6,1]
# dim(zz)

XX <- x4[,-1]
colnames(XX)[1:3]
XX1 <- cbind(zz,XX)
colnames(XX1) <- c("entrezID",colnames(XX))
# XX1[1:6,1]

## 用基因id对整理好的芯片注释文件进行基因名的更新
homo<-read.table("homo.txt",header=T,sep='\t')
# dim(homo)
#homo[5,]
homo[1:6,1]
x5 <- merge(homo,XX1,by.x="GeneID",by.y = "entrezID",all=FALSE) 
# 合并， x5 从54613，降为46886，基因号与基因名进行匹配
# View(x5[,1:10])
# dim(x5)    # 46896    92


## 探针名匹配基因名，取出多个探针对应一个基因的数据计算IQR，保留IQR最大的探针数据
expset4 <- x5 %>%
  dplyr::select(-GeneID) %>%              # 去掉多余信息
  mutate(rowIQR =apply(.[,-1],1,IQR)) %>% # 计算每行的IQR
  arrange(desc(rowIQR)) %>%               # 把表达量的平均值按从大到小排序
  distinct(Symbol,.keep_all = T) %>%      # symbol留下第一个
  dplyr::select(-rowIQR) %>%                 # 反向选择去除rowIQR这一列
  tibble::column_to_rownames(colnames(.)[1]) # 把第一列变成行名并删除
# View(expset4[1:10,1:10])
dim(expset4)   # 21835    90
# write.table(expset4,"GSE20711_expr.txt",quote=F,sep="\t")  


#################### 贴标签 从266个样本--252个样本 #######################################

lable = read.csv("GSE20711_all.csv", header = T, sep=',')
View(lable[,1:10])

sample <- as.matrix(lable$Accession)
sample1 <- matrix(NA)
for (i in 1:90) {
  sample1[i] <- sample[1+6*(i-1),]
}
sample2 <- as.matrix(sample1)
colnames(sample2) <- c("sample") 
View(sample2)

character <- as.matrix(lable$Characteristics)
View(character)



###################################################################################

# e_rfs -------------------------------------------------------------------

e_rfs1 <- matrix(NA)
for (i in 1:90) {
  e_rfs1[i] <- character[3+6*(i-1),]
}

e_rfs2 <- as.matrix(e_rfs1)
colnames(e_rfs2) <- c("e_rfs")
# View(e_rfs2)

e_rfs3 <- e_rfs2 %>% 
  # then make it a tibble (nice printing while debugging)
  as_tibble() %>% 
  # then trim the barcode (see head(clin), and substr)
  mutate(e_rfs = substr(e_rfs, 7, 8)) # 取19-25个字符(可以>25)
# View(e_rfs3)


# t_rfs -------------------------------------------------------------------

t_rfs1 <- matrix(NA)
for (i in 1:90) {
  t_rfs1[i] <- character[4+6*(i-1),]
}

t_rfs2 <- as.matrix(t_rfs1)
colnames(t_rfs2) <- c("t_rfs")
# View(t_rfs2)

t_rfs3 <- t_rfs2 %>% 
  # then make it a tibble (nice printing while debugging)
  as_tibble() %>% 
  # then trim the barcode (see head(clin), and substr)
  mutate(t_rfs = substr(t_rfs, 7, 11)) # 取19-25个字符(可以>25)
# View(t_rfs3)


# rfs ---------------------------------------------------------------------

rfs <- cbind(sample1, t_rfs3, e_rfs3)[-c(89,90),]
View(rfs)

###################################################################################



# e_os -------------------------------------------------------------------

e_os1 <- matrix(NA)
for (i in 1:90) {
  e_os1[i] <- character[5+6*(i-1),]
}

e_os2 <- as.matrix(e_os1)
colnames(e_os2) <- c("e_os")
# View(e_os2)

e_os3 <- e_os2 %>% 
  # then make it a tibble (nice printing while debugging)
  as_tibble() %>% 
  # then trim the barcode (see head(clin), and substr)
  mutate(e_os = substr(e_os, 6, 7)) # 取19-25个字符(可以>25)
# View(e_os3)


# t_os -------------------------------------------------------------------

t_os1 <- matrix(NA)
for (i in 1:90) {
  t_os1[i] <- character[6+6*(i-1),]
}

t_os2 <- as.matrix(t_os1)
colnames(t_os2) <- c("t_os")
# View(t_os2)

t_os3 <- t_os2 %>% 
  # then make it a tibble (nice printing while debugging)
  as_tibble() %>% 
  # then trim the barcode (see head(clin), and substr)
  mutate(t_os = substr(t_os, 6, 10)) # 取19-25个字符(可以>25)
# View(t_os3)


# os ---------------------------------------------------------------------

os <- cbind(sample1, t_os3, e_os3)[-c(89,90),]
View(os)

###################################################################################


data <- read.table("GSE20711_expr.txt",header=T,sep='\t',fill=TRUE,strip.white = T)
dim(data)    # 21835    90
View(data[,1:10])
View(expset4[,1:10])

expset6 <- t(scale(t(data)))
dim(expset6)    # 21835    90  
# View(t(expset6)[,1:10])


# data_rfs ---------------------------------------------------------------------

View(as.matrix(expset6[,-c(89,90)])[,1:10])

data_rfs = rbind(as.matrix(t(rfs[,c(2,3)])), as.matrix(expset6[,-c(89,90)]))
colnames(data_rfs) <- c(colnames(expset6[,-c(89,90)]))
View(data_rfs[,1:10])
dim(data_rfs)    # 21837    88
# write.table(data_rfs,"GSE20711_scale_outcome_rfs.txt",quote=F,sep="\t")  

# data_os ---------------------------------------------------------------------

data_os = rbind(as.matrix(t(os[,c(2,3)])), as.matrix(expset6[,-c(89,90)]))
colnames(data_os) <- c(colnames(expset6[,-c(89,90)]))
View(data_os[,1:10])
dim(data_os)    # 21837    88
# write.table(data_os,"GSE20711_scale_outcome_os.txt",quote=F,sep="\t")  





