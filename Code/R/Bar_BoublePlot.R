## 2023.4.24 Plot the double Y Figure
## try to cut-off Y axis https://blog.csdn.net/woodcorpse/article/details/111148666


rm(list = ls())



library(tidyverse) #ggplot2包等
library(reshape2) #用到melt()函数将宽数据转换成长数据
# install.packages("gg.gap")
library(gg.gap)

# color setting -----------------------------------------------------------

colors8 <- c("#d79e9e", "#8ccddc", "#5dbdaf", "#c6cdd7", "#efbaa7", "#cab2d6", "#b3dad1", "#c5c48b")



# figure save -------------------------------------------------------------

Figure <- list()


# Compare with seven methods --------------------------------------------------------------

## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- '/home/lly/R/'
setwd(paste(path, 'CNetCox/Result/ComMethd', sep=''))

 
## feature number and Ci value
data <- as.data.frame(read.csv("CompareFeature.csv", header = T, sep = ","))
data[,c(1:3)]
data$Method <- factor(data$Method)
data$Number.of.features



p1 <- ggplot(data)+
  geom_col(aes(Method,Number.of.features,fill=Method),
           position = 'dodge',width = 0.6)+
  scale_fill_manual(values = colors8)+
  labs(x='Method',y='Number of features') +
  theme_test(base_size = 15) +
  theme(legend.position = 'none',
        panel.border = element_rect(size=2,fill = 'transparent'),
        axis.text = element_text(color='black'),
        axis.text.x = element_text(color='black', angle = 45, 
                                   hjust = 0.5, vjust = 0.5, size = 11)) 

p2 = gg.gap(plot = p1,
            segments = c(80, 500),
            tick_width = 20,
            rel_heights = c(0.25, 0, 0.1),# 设置分隔为的三个部分的宽度
            ylim = c(0, 550)
)

p2
Figure[[2]] <- p2


# col1 <- 80
# col2 <- 1.2
# k1 <- 20
# k2 <- 0.2
# 
# p2 + scale_y_continuous(limits = c(0.0, col1),
#                         breaks = seq(0.0, col1, k1),
#                         expand = c(0,0),
#                         sec.axis = sec_axis(~./(col1/col2),
#                                             name = 'C-index',
#                                             breaks = seq(0, col2, k2))) +
#   geom_point(data = data,
#              aes(Method, CI*(col1/col2)+0.01,
#                  color=variable),
#              size=4) +
#   geom_line(data = data,
#             aes(Method, CI*(col1/col2)+0.01,
#                 color=variable, group=variable),
#             cex=1.3) +
#   scale_color_manual(values = '#1e8b9b') 


col1 <- 80
col2 <- 1.05
k1 <- 20
k2 <- 0.2
# col1 <- 550
# col2 <- 1.2
# k1 <- 50
# k2 <- 0.2

variable <- rep("method", 8)
p3 <- ggplot(data)+
  geom_col(aes(Method,Number.of.features,fill=Method),
           position = 'dodge',width = 0.6)+
  scale_fill_manual(values = colors8)+
  labs(x='Method',y='Number of features') +
  theme_test(base_size = 15)+
  theme(legend.position = 'none',
        panel.border = element_rect(size=2,fill = 'transparent'),
        axis.text = element_text(color='black'),
        axis.text.x = element_text(color='black', angle = 45, 
                                   hjust = 0.5, vjust = 0.5, size = 11)) +
  scale_y_continuous(limits = c(0.0, col1),
                     breaks = seq(0.0, col1, k1),
                     expand = c(0,0),
                     sec.axis = sec_axis(~./(col1/col2),
                                         name = 'C-index',
                                         breaks = seq(0, col2, k2))) +
  geom_point(data = data,
             aes(Method, CI*(col1/col2)+0.01,
                 color=variable),
             size=3) +
  geom_line(data = data,
            aes(Method, CI*(col1/col2)+0.01,
                color=variable, group=variable),
            cex=1.0) +
  scale_color_manual(values = '#1e8b9b') 




# pdf('FigurescBioRevised/RealdataAUROC_Data2.pdf', width = 5.5, height = 5.5)
p3
# dev.off()

Figure[[3]] <- p3


# pdf("CompareFeature.pdf",width = 10, height = 5)
cowplot::plot_grid(plotlist = list(Figure[[2]], Figure[[3]] ), 
                   nrow = 1, rel_widths = c(5.5, 5.5),
                   labels = c('(a)','(b)'),
                   scale = c(0.95)) 
# dev.off()




# Plot the bar of coef ----------------------------------------------------
## Example
# create example dataset
df <- data.frame(
  variable = c("A", "B", "C", "D", "E"),
  weight = c(-0.5, 0.8, -1.2, 0.4, 0.7)
)

# create bar graph
library(ggplot2)
ggplot(df, aes(x = variable, y = weight)) +
  geom_col(aes(fill = weight > 0), color = "black") +
  scale_fill_manual(values = c("#FF6B6B", "#77DD77")) +
  coord_flip() +
  labs(x = "Variable Name", y = "Weight Coefficient") +
  theme_minimal() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



## set path
setwd(paste(path, 'CNetCox/Data/TCGA_NEW/', sep=''))

marker6 <- read.csv("UniMutVariate_markergene.csv", header = T, sep = ",")
colnames(marker6) <- c("Gene", "Coefficient")


library(ggplot2)
p4 <- ggplot(data = marker6, aes(x = Coefficient, y = Gene)) +
  # geom_bar(stat = "identity", fill = "steelblue") +
  geom_col(aes(fill = Coefficient > 0), color = "gray") +
  scale_fill_manual(values = c("#e6c7c2", "#a4afd5")) +
  scale_x_continuous(limits = c(-0.2, 0.3), expand = c(0, 0)) +
  theme_classic() +
  # theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 14, hjust = ifelse(df$Coefficient < 0, 0, 1)),
        axis.line = element_line(colour = "black", size = 0.5),
        legend.position = "none",
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 16),
        plot.caption = element_text(size = 12, hjust = 1, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        plot.margin = margin(40, 40, 20, 20, "pt")) +
  labs(# title = "Prognostic biomarkers",
       # subtitle = "Prognostic biomarkers",
       # caption = "Source: Generated by R"
       x = "Coefficient")



Figure[[4]] <- p4 
setwd(paste(path, 'CNetCox/Result/ComMethd', sep=''))
# pdf("CompareFeature_marker.pdf",width = 10, height = 5)
cowplot::plot_grid(plotlist = list(Figure[[2]], Figure[[4]] ), 
                   nrow = 1, rel_widths = c(5.5, 5.5),
                   labels = c('(a)','(b)'),
                   scale = c(0.95)) 
# dev.off()


# DEGs plot of six genes --------------------------------------------------

setwd(paste(path, 'CNetCox/Data', sep=''))
data = read.table("TCGA_NEW/TCGA_PRS_Expr.txt", header = T, check.names = FALSE)
datat <- as.data.frame(t(data)[,-1])

# datat$status <- "Dead"
which(datat$status == 1)
datat[which(datat$status == 1),1] <- "Dead"
datat[which(datat$status == 0),1] <- "Alive"

# Set theme and colors
theme_set(theme_cowplot())
colors <- c("#E69F00", "#56B4E9")

# Function to calculate p values for t-test between two groups
calc_p <- function(x) {
  t.test(x ~ datat$status)$p.value
}

library(ggplot2)
library(ggpubr)
# Plot differential expression box plots for each variable
p1 <- ggplot(datat, aes(x = status, y = EGR1, fill = status)) +
  geom_boxplot(alpha = 0.8, color = "#404040") +
  scale_fill_manual(values = colors) +
  labs(title = "EGR1", x = "", y = "Expression value") +
  theme(plot.title = element_text(size = rel(1.2)), # , face = "bold"
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.position = "none") +
  stat_compare_means(label = "p.format", method = "t.test", size = 6, label.x = 1.5, label.y = max(datat$EGR1)) 

p2 <- ggplot(datat, aes(x = status, y = IGFBP5, fill = status)) +
  geom_boxplot(alpha = 0.8, color = "#404040") +
  scale_fill_manual(values = colors) +
  labs(title = "IGFBP5", x = "", y = "Expression value") +
  theme(plot.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.3)),
        legend.position = "none") +
  stat_compare_means(label = "p.format", method = "t.test", size = 6, label.x = 1.5, label.y = max(datat$IGFBP5)) 

p3 <- ggplot(datat, aes(x = status, y = JUN, fill = status)) +
  geom_boxplot(alpha = 0.8, color = "#404040") +
  scale_fill_manual(values = colors) +
  labs(title = "JUN", x = "", y = "Expression value") +
  theme(plot.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.3)),
        legend.position ="none") +
  stat_compare_means(label = "p.format", method = "t.test", size = 6, label.x = 1.5, label.y = max(datat$JUN))

p4 <- ggplot(datat, aes(x = status, y = MAFK, fill = status)) +
  geom_boxplot(alpha = 0.8, color = "#404040") +
  scale_fill_manual(values = colors) +
  labs(title = "MAFK", x = "", y = "Expression value") +
  theme(plot.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.3)),
        legend.position = "none") +
  stat_compare_means(label = "p.format", method = "t.test", size = 6, label.x = 1.5, label.y = max(datat$MAFK))

p5 <- ggplot(datat, aes(x = status, y = MYC, fill = status)) +
  geom_boxplot(alpha = 0.8, color = "#404040") +
  scale_fill_manual(values = colors) +
  labs(title = "MYC", x = "", y = "Expression value") +
  theme(plot.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.3)),
        legend.position = "none") +
  stat_compare_means(label = "p.format", method = "t.test", size = 6, label.x = 1.5, label.y = max(datat$MYC))

p6 <- ggplot(datat, aes(x = status, y = TCF7, fill = status)) +
  geom_boxplot(alpha = 0.8, color = "#404040") +
  scale_fill_manual(values = colors) +
  labs(title = "TCF7", x = "", y = "Expression value") +
  theme(plot.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.3)),
        legend.position = "none") +
  stat_compare_means(label = "p.format", method = "t.test", size = 6, label.x = 1.5, label.y = max(datat$TCF7)) 

# Arrange the plots on a grid
library(gridExtra)
grid.arrange(p1,p2, p3, p4, p5, p6, ncol = 3, widths = c(1.0, 1.0, 1.0))


# setwd(paste(path, 'CNetCox/Result/ComMethd', sep=''))
# pdf("DEG_TCGA.pdf",width = 10, height = 10)
# grid.arrange(p1,p2, p3, p4, p5, p6, ncol = 3, widths = c(1.0, 1.0, 1.0))
# dev.off()



# P/M GSE147995 DEGs plot of six genes --------------------------------------------------

# DEGs plot of six genes --------------------------------------------------

setwd(paste(path, 'CNetCox/Data', sep=''))
## orignal
data = read.table("DrugData/GSE147995_PRS_Expr.txt", header = T, check.names = FALSE)
datat <- as.data.frame(t(data))
datat$status <-  rep(c("Primary", "Metastases"), 13)
## ER plus 
# data = read.table("DrugData/GSE147995_PRS_Expr_plus.txt", header = T, check.names = FALSE)
# datat <- as.data.frame(t(data))
# datat$status <-  rep(c("Primary", "Metastases"), 6)
## ER plus or mini
# data = read.table("DrugData/GSE147995_PRS_Expr_mini.txt", header = T, check.names = FALSE)
# datat <- as.data.frame(t(data))
# datat$status <-  rep(c("Primary", "Metastases"), 7)

# Set theme and colors
theme_set(theme_cowplot())
colors <- c("#E69F00", "#56B4E9")

# Function to calculate p values for t-test between two groups
calc_p <- function(x) {
  t.test(x ~ datat$status)$p.value
}

# Plot differential expression box plots for each variable
p1 <- ggplot(datat, aes(x = status, y = EGR1, fill = status)) +
  geom_boxplot(alpha = 0.8, color = "#404040") +
  scale_fill_manual(values = colors) +
  labs(title = "EGR1", x = "", y = "Expression value") +
  theme(plot.title = element_text(size = rel(1.2)), # , face = "bold"
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.position = "none") +
  stat_compare_means(label = "p.format", method = "t.test", size = 6, label.x = 1.5, label.y = max(datat$EGR1)) 

p2 <- ggplot(datat, aes(x = status, y = IGFBP5, fill = status)) +
  geom_boxplot(alpha = 0.8, color = "#404040") +
  scale_fill_manual(values = colors) +
  labs(title = "IGFBP5", x = "", y = "Expression value") +
  theme(plot.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.3)),
        legend.position = "none") +
  stat_compare_means(label = "p.format", method = "t.test", size = 6, label.x = 1.5, label.y = max(datat$IGFBP5)) 

p3 <- ggplot(datat, aes(x = status, y = JUN, fill = status)) +
  geom_boxplot(alpha = 0.8, color = "#404040") +
  scale_fill_manual(values = colors) +
  labs(title = "JUN", x = "", y = "Expression value") +
  theme(plot.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.3)),
        legend.position ="none") +
  stat_compare_means(label = "p.format", method = "t.test", size = 6, label.x = 1.5, label.y = max(datat$JUN))

p4 <- ggplot(datat, aes(x = status, y = MAFK, fill = status)) +
  geom_boxplot(alpha = 0.8, color = "#404040") +
  scale_fill_manual(values = colors) +
  labs(title = "MAFK", x = "", y = "Expression value") +
  theme(plot.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.3)),
        legend.position = "none") +
  stat_compare_means(label = "p.format", method = "t.test", size = 6, label.x = 1.5, label.y = max(datat$MAFK))

p5 <- ggplot(datat, aes(x = status, y = MYC, fill = status)) +
  geom_boxplot(alpha = 0.8, color = "#404040") +
  scale_fill_manual(values = colors) +
  labs(title = "MYC", x = "", y = "Expression value") +
  theme(plot.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.3)),
        legend.position = "none") +
  stat_compare_means(label = "p.format", method = "t.test", size = 6, label.x = 1.5, label.y = max(datat$MYC))

p6 <- ggplot(datat, aes(x = status, y = TCF7, fill = status)) +
  geom_boxplot(alpha = 0.8, color = "#404040") +
  scale_fill_manual(values = colors) +
  labs(title = "TCF7", x = "", y = "Expression value") +
  theme(plot.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.3)),
        legend.position = "none") +
  stat_compare_means(label = "p.format", method = "t.test", size = 6, label.x = 1.5, label.y = max(datat$TCF7)) 

# Arrange the plots on a grid
grid.arrange(p1,p2, p3, p4, p5, p6, ncol = 3, widths = c(1.0, 1.0, 1.0))


# setwd(paste(path, 'CNetCox/Result/ComMethd', sep=''))
# pdf("DEG_GSE147995.pdf",width = 10, height = 10)
# grid.arrange(p1,p2, p3, p4, p5, p6, ncol = 3, widths = c(1.0, 1.0, 1.0))
# dev.off()

# Only show the last gene TCF7 --------------------------------------------

setwd(paste(path, 'CNetCox/Data', sep=''))
## orignal
data = read.table("DrugData/GSE147995_PRS_Expr.txt", header = T, check.names = FALSE)
datat <- as.data.frame(t(data))
datat$Status <- rep(c("Primary", "Metastases"), 13)
library(stringr)
datat$Patient <- rep(str_c("P", c(1,2,5:7,9,11:15,17,8)), each = 2)


# create box plots with scatter points and connecting lines
library(ggplot2)
library(ggpubr)
# ggplot(datat, aes(x = Status y = TCF7, color = Status)) +
## 修改横轴次序
p7 <- ggplot(datat, aes(x = factor(Status, levels = unique(Status)), y = TCF7, color = Status)) +
  geom_boxplot(alpha =0.5,size=1,outlier.shape = NA) +
  scale_color_manual(limits=c("Metastases","Primary"), 
                     # values=c("#E29827","#922927")) +   
                     values=c("#E69F00", "#56B4E9")) + 
  stat_compare_means(method = "t.test",paired = TRUE,
                     comparisons=list(c("Metastases","Primary"))) +
  geom_jitter(alpha = 0.3,size=4, aes(fill=Patient), shape=21, position = position_dodge(0.5)) +
  geom_line(aes(group = Patient), 
            color = 'grey40', lwd = 0.5, position = position_dodge(0.5))+ #添加连线
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  theme_bw() +
  # theme_classic() +
  # theme_minimal() +
  theme(panel.grid =element_blank(),
        panel.border = element_rect(size=1.5,fill = 'transparent'),
        axis.text = element_text(size = 12, colour = "black"),
        legend.position = 'right') +
  theme(plot.title = element_text(hjust =0.5)) +    # 标题居中
  labs(title = "TCF7", x = "", y = "Gene Expression") 

p7
Figure[[7]] <- p7

# create heatmap plots among primary and metastates 
library(pheatmap)

selected <- cbind(t(datat[seq(1,25,2),-c(7,8)]), t(datat[-seq(1,25,2),-c(7,8)]))
# View(selected[1:5,1:5])
Label = read.table("DrugData/phe_PM.txt", header = TRUE, sep = "\t")
Label <- factor(Label[,"class"])
Label <- data.frame(Label)
rownames(Label) = colnames(selected)
selected[1,1]

selectednum <- apply(selected, 2, as.numeric)
# t = pheatmap(selectednum, scale="row", 
#              show_colnames=F,
#              cluster_row=F)

# ## 关键一步，进行平滑，每一行都给他平滑一下子
matrix.smooth <- apply(selectednum,1,function(x){smooth.spline(x)$y})
## 然后进行标准化一下子
normalized_matrix <- apply(matrix.smooth,1,function(x){(x-mean(x))/sd(x)})
## 图c,再来看就有那种平滑的感觉了
pheatmap (normalized_matrix, show_colnames=F,cluster_row=F, cluster_col=F)
#图d，加上monocle3独家冠名的颜色，看着是不是就有那么一点儿意思了。如果是比较混乱的可能更好看
bks <- seq(-3.1, 3.1, by = 0.1)

colnames(normalized_matrix) <- colnames(selectednum) <- colnames(selected)
rownames(normalized_matrix) <- rownames(selectednum) <- rownames(selected)
p8 <- pheatmap(normalized_matrix,    # normalized_matrix    selectednum
               annotation_col = Label, 
               # color = colorRampPalette(c("blue", "white","red"))(100),
               fontsize_row = 8,
               # scale = "row", 
               cutree_cols = 2,
               cluster_col = F,
               border_color = NA,
               cluster_row = F,
               show_colnames = T,
               show_rownames = T)

p8

Figure[[8]] <- p8

# p <- pheatmap(selectednum,annotation_col = Label, 
#               color = colorRampPalette(c("blue", "white","red"))(100),
#               fontsize_row = 8,scale = "row", cutree_cols = 2, border_color = NA)


setwd(paste(path, 'CNetCox/Result/ComMethd', sep=''))
# pdf("DrugAnalysis.pdf",width = 10, height = 5)
cowplot::plot_grid(plotlist = list(Figure[[7]], Figure[[7]] ), 
                   nrow = 1, rel_widths = c(5.5, 5.5),
                   labels = c('(a)','(b)'),
                   scale = c(0.95)) 
# dev.off()



# pdf("DrugAnalysis_heatmap.pdf",width = 5.5, height = 5.5)
p8
# dev.off()



# Line Dot Plot -- Example -----------------------------------------------------------

## https://zhuanlan.zhihu.com/p/609173487

# create data.frame
df <- data.frame(
  Status = rep(c("Dead", "Alive"), each = 10),
  Patient = rep(c("p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "p10"), 2),
  Var1 = rnorm(20, 10, 2),
  Var2 = rnorm(20, 5, 1),
  Var3 = rnorm(20, 15, 3),
  Var4 = rnorm(20, 8, 1.5),
  Var5 = rnorm(20, 20, 4),
  Var6 = rnorm(20, 12, 2.5)
)

df <- rbind(df, df)
df$sample <- rep(c("Gene1", "Gene2"), each = 20)

# create box plots with scatter points and connecting lines
library(ggplot2)

ggplot(df, aes(x = Status, y = Var1, color = Status)) +
  geom_boxplot(alpha =0.5,size=1,outlier.shape = NA) +
  scale_color_manual(limits=c("Alive","Dead"), 
                     values=c("#E29827","#922927")) +
  stat_compare_means(method = "t.test",paired = TRUE,
                     comparisons=list(c("Alive","Dead"))) +
  geom_jitter(alpha = 0.3,size=3, aes(fill=Patient), shape=21, position = position_dodge(0.5)) +
  geom_line(aes(group = Patient), 
            color = 'grey40', lwd = 0.5, position = position_dodge(0.5))+ #添加连线
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  facet_wrap(~sample, scales = "free_y") +
  theme_bw() +
  theme(panel.grid =element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'top')+
  labs(y='Expression')




# Volcano plot ------------------------------------------------------------


## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- '/home/lly/R/'
setwd(paste(path, 'CNetCox/Data/TCGA_NEW', sep=''))

volcanoDE <- read.csv("Volcano_DE.csv", header = T, sep = ",")

# Volcano plot
library(ggplot2)

# Define colors
color_significant <- "#e7897d"
color_nonsignificant <- "#9fa0b5"

# Define significance threshold
sig_threshold <- -log10(0.01)
# sig_threshold <- -log2(10)

# Create plot
p9 <- ggplot(volcanoDE, aes(x=logFC, y=-log10(pvalue))) +
  geom_point(aes(color=ifelse(abs(logFC)>log2(2) & -log10(pvalue)>2, "significant", "nonsignificant")), size=2) +
  scale_color_manual(values=c("significant"=color_significant, "nonsignificant"=color_nonsignificant)) +
  geom_vline(xintercept=c(-1, 1), linetype="dashed", color="gray50", size=0.5) +
  geom_hline(yintercept=sig_threshold, linetype="dashed", color="gray50", size=0.5) +
  labs(x="Log2 Fold Change", y="-Log10(P.adj-value)", color="Significance") +
  theme_classic() +
  theme(axis.line=element_line(color="black", size=0.5),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),# ,face="bold"
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        legend.position="bottom",
        legend.key.size=unit(0.5,"cm"),
        panel.background = element_rect(fill="white"),
        plot.title=element_text(hjust=0.5, size=14, face="bold"),
        plot.subtitle=element_text(hjust=0.5, size=12, face="italic"),
        plot.caption=element_text(hjust=1, size=12),
        plot.margin=unit(c(1,1,1,1),"cm"))
# plot.margin=unit(c(1,1,1,1),"cm")) +
# ggtitle("Volcano Plot of Gene Expression Data",
#         subtitle="Genes with absolute log2 fold change >1 and -log10 p-value >2 are considered significant") +
# labs(caption="Source: BRCA dataset")


Figure[[9]] <- p9

# p <- pheatmap(selectednum,annotation_col = Label, 
#               color = colorRampPalette(c("blue", "white","red"))(100),
#               fontsize_row = 8,scale = "row", cutree_cols = 2, border_color = NA)


setwd(paste(path, 'CNetCox/Result/ComMethd', sep=''))
# pdf("CompareFeature_volcano_coef.pdf",width = 15, height = 5)
cowplot::plot_grid(plotlist = list(Figure[[9]], Figure[[3]], Figure[[4]]), 
                   nrow = 1, rel_widths = c(5.5, 5.5, 5.5),
                   labels = c('(a)','(b)','(c)'),
                   scale = c(0.95)) 
# dev.off()

# Functional CNet and ENet compare -----------------------------------------------------------

t.test(1:44.77, 1:27.90)

datafuction <- cbind(seq(1, 44.77, 0.5), seq(1, 27.90, 0.5))

P <- rbind(as.matrix(seq(1, 44.77, 0.5)), as.matrix(seq(1, 27.90, 0.5)))
Method <- c(rep("CNet", length(seq(1, 44.77, 0.5))), 
                        rep("ENet",length(seq(1, 27.90, 0.5)) ))
Sample <- str_c("P", seq(1, length(seq(1, 44.77, 0.5)) + length(seq(1, 27.90, 0.5))))
datafuction <- as.data.frame(cbind(Method, P))
colnames(datafuction) <- c("Method", "P")
datafuction$P <- as.numeric(datafuction$P)


## Bar plot

# Calculate p-value
t_test <- t.test(P ~ Method, data = datafuction)
p_value <- signif(t_test$p.value, 4)
# p_value <- round(t_test$p.value, 3)

# Plot bar chart in Cell academic journal style
plot <- ggplot(datafuction, aes(x = Method, y = P, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("CNet" = "darkblue", "ENet" = "lightblue")) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.position = "none") +
  labs(x = "Method", y = "P value")

# Add p-value to the bar chart
pfunction <- plot + 
  geom_text(aes(x = 1.5, y = max(P) + 1, label = paste0("p = ", p_value)), size = 5, colour = "black") +
  geom_segment(aes(x = 1, y = max(P), xend = 2, yend = max(P)), color = "black") +
  geom_segment(aes(x = 1, y = max(P), xend = 1, yend = max(P) - 1), color = "black") +
  geom_segment(aes(x = 2, y = max(P), xend = 2, yend = max(P), color = "black"))

# setwd(paste(path, 'CNetCox/Result/ComMethd', sep=''))
# pdf("FunctionP2method_bar.pdf",width = 5.5, height = 5.5)
plot
# dev.off()


## Bar plot example 
# # Create data
# method <- c("Method1", "Method1", "Method1", "Method1", "Method1",
#             "Method2", "Method2", "Method2", "Method2", "Method2")
# value <- c(10, 15, 12, 17, 20, 22, 18, 25, 24, 19)
# 
# data <- data.frame(method, value)
# 
# # Change the order of method names on the x-axis
# data$method <- factor(data$method, levels = c("Method2", "Method1"))
# 
# # Calculate p-value
# t_test <- t.test(value ~ method, data = data)
# p_value <- round(t_test$p.value, 3)
# 
# # Plot bar chart
# plot <- ggplot(data, aes(x = method, y = value, fill = method)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   theme_minimal() +
#   labs(x = "Method Name", y = "Value", title = "Comparison of Methods Bar Chart")
# 
# # Add p-value to the bar chart
# plot + 
#   geom_text(aes(x = 1.5, y = max(data$value) + 1, label = paste0("p = ", p_value)), size = 5, colour = "black") 



## Box plot
library(ggplot2)
library(ggpubr)
# Plot differential expression box plots for each variable
pfunction <- ggplot(datafuction, aes(x = Method, y = P, fill = Method)) +
  geom_boxplot(alpha = 0.8, color = "#404040") +
  scale_fill_manual(values = colors) +
  labs(title = " ", x = "", y = "-Log10(P)") +
  theme(plot.title = element_text(size = rel(1.2)), # , face = "bold"
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.position = "none") +
  theme_bw() +
  # theme_classic() +
  # theme_minimal() +
  stat_compare_means(label = "p.format", method = "t.test", size = 6, label.x = 1.5, label.y = max(datafuction$P)) 
pfunction

setwd(paste(path, 'CNetCox/Result/ComMethd', sep=''))
# pdf("FunctionP2method.pdf",width = 5.5, height = 5.5)
pfunction
# dev.off()




# nomogram ----------------------------------------------------------------

# Load required libraries
library(dplyr)

# Set seed for reproducibility
set.seed(42)

# Generate data
sample_names <- paste("Sample", 1:50, sep = "_")
survival_status <- sample(c(0, 1), 50, replace = TRUE)
survival_time <- round(runif(50, 1, 100))
age <- round(runif(50, 18, 90))
pathological_stage <- sample(c(0, 1), 50, replace = TRUE)

# Create data frame
sample_data <- data.frame(SampleName = sample_names,
                          SurvivalStatus = survival_status,
                          SurvivalTime = survival_time,
                          Age = age,
                          PathologicalStage = pathological_stage)

# Add three additional columns with random data
sample_data <- sample_data %>%
  mutate(Column6 = round(runif(50, 0, 1)),
         Column7 = round(runif(50, 0, 1)),
         Column8 = round(runif(50, 0, 1)))

# View the data frame
head(sample_data)


library(survival)
library(rms)

# Prepare data
sample_data$SurvivalStatus <- as.factor(sample_data$SurvivalStatus)
sample_data$PathologicalStage <- as.factor(sample_data$PathologicalStage)

# Create a survival object
surv_obj <- Surv(sample_data$SurvivalTime, sample_data$SurvivalStatus)

# Fit a Cox regression model
cox_model <- cph(surv_obj ~ Age + PathologicalStage, data = sample_data, x = TRUE, y = TRUE)

# Create a nomogram
nom_cox <- nomogram(cox_model, fun = NULL, lp = FALSE)

# Plot the nomogram
plot(nom_cox)