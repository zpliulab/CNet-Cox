## 2023.4.12 LLY create, 2026.1.16 LLY used
## Split TCGA data into Training and Test data with clinical data
## 2023.1.21 using five datasets and see the intersect
## 2023.4.23 using nine datasets and see each of them
## 2026.1.18 See C-index ofCNet on different cv, from  CNetCoxComResults.xlsx


##############################################################################
# Figure 2B. 68 prognostic bioomarkers
##############################################################################

## clear
rm(list = ls())

## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- '/home/lly/R/'
setwd(paste(path, 'CNetCox/Data/', sep=''))

## load function
source(paste(path, 'CNetCox/R/MyFunctions.R', sep=''))

## load gene and net and cut
gene <- read.csv('TCGA_NEW/UNgene_component.csv',header = T)
net <- read.csv('TCGA_NEW/UNgene_comp_net.csv')
cut <- read.table("TCGA_NEW/cut_vector_UNgene.txt", header = T, check.names = FALSE, sep = "\t")


# Marker genes ------------------------------------------------------------
ci <- c()
p <- c()

# for (k in c(4.9e-2,
#             # 5.0e-2, 5.1e-2, 5.2e-2, 5.3e-2, 
#             # 5.31e-2, 5.33e-2, 5.35e-2, 5.37e-2, 5.39e-2,
#             # 5.4e-2, 
#             # 5.41e-2, 5.43e-2, 5.45e-2, 5.47e-5, 5.49e-2,
#             5.5e-2, 5.6e-2, 5.7e-2, 5.8e-2, 5.9e-2, 6e-2,
#             6.1e-2, 6.2e-2, 6.3e-2, 6.4e-2, 6.5e-2,
#             # 1.5e-1, 1.6e-1, 1.7e-1, 1.8e-1, 1.9e-1,
#             2e-1)) {
#   marker3 <- markerselect(1,k)[[1]]
#   coef3 <- markerselect(1,k)[[2]]
  


# marker1 <- markerselect(1,1.4e-1)[[1]]
# marker2 <- markerselect(2,4e-2)[[1]]
# marker3 <- markerselect(3,5.15e-2)[[1]]
# marker4 <- markerselect(4,9.5e-2)[[1]]
# marker5 <- markerselect(5,5e-2)[[1]]
# marker6 <- markerselect(6,4.7e-2)[[1]]
# # marker7 <- markerselect(7,2e-1)[[1]]
# marker8 <- markerselect(8,2e-1)[[1]]
# marker9 <- markerselect(9,2e-1)[[1]]
# marker10 <- markerselect(10,2e-1)[[1]]

# marker1 <- markerunion(1)[[1]]
# marker2 <- markerunion(2)[[1]]
# marker3 <- markerunion(3)[[1]]
# marker4 <- markerunion(4)[[1]]
# marker5 <- markerunion(5)[[1]]
# marker6 <- markerunion(6)[[1]]
# # marker7 <- markerunion(7)[[1]]
# # marker8 <- markerunion(8)[[1]]
# marker9 <- markerunion(9)[[1]]
# marker10 <- markerunion(10)[[1]]

# markeru <- Reduce(union,
                  # list(marker1, marker2, marker3, marker4, marker5))

# write.csv(marker4, "Result/Result4/marker4_revision.csv", row.names = F)

# Marker coefs ------------------------------------------------------------
# coef1 <- markerselect(1,1.4e-1)[[2]]
# coef2 <- markerselect(2,4e-2)[[2]]
# coef3 <- markerselect(3,5.15e-2)[[2]]
# coef4 <- markerselect(4,9.5e-2)[[2]]
# coef5 <- markerselect(5,5e-2)[[2]]
coef6 <- markerselect(6,4.7e-2)[[2]]
# # coef7 <- markerunion(7,2e-1)[[2]]
# # coef8 <- markerunion(8,2e-1)[[2]]
# coef9 <- markerselect(9,2e-1)[[2]]
# coef10 <- markerselect(10,2e-1)[[2]]


# coef1 <- markerunion(1)[[2]]
# coef2 <- markerunion(2)[[2]]
# coef3 <- markerunion(3)[[2]]
# coef4 <- markerunion(4)[[2]]
# coef5 <- markerunion(5)[[2]]
# coef6 <- markerunion(6)[[2]]
# # coef7 <- markerunion(7)[[2]]
# # coef8 <- markerunion(8)[[2]]
# coef9 <- markerunion(9)[[2]]
# coef10 <- markerunion(10)[[2]]

# coefu <- cbind(coef1, coef2, coef3, coef4, coef5)
# coefmean <- as.matrix(apply(coefu,1,mean)) 
# rownames(coefmean) <- gene[,1]

# write.csv(markeru, "Result/markerunion5.csv", row.names = F)


# see the network ---------------------------------------------------------
filenum <- "6"
markerinput <- marker6
coefinput <- coef6

g <- markerPlot(markerinput, net)[[1]]

# deg <- degree(g, mode="all")
# V(g)$size <- deg * 2
# plot(g,
#      edge.curved=.1,
#      edge.arrow.size = .3)

markerNet <- markerPlot(markerinput, net)[[2]]
# write.csv(markerNet, "Result/Result3/markerNet3.csv", row.names = F, quote = F)


##############################################################################
# Figure 2G. Survival analysis results of 68 prognostic BRCA markers.
##############################################################################
# load test data and calculus the CI value --------------------------------
# filenum <- "3"
# readtestdata <- function(filenum){
data <- read.table(paste("Data_test/", filenum, ".txt", sep = ""), header = T, sep = "")
Data_test <- t(data)
y1 <- t(Data_test[c(1,2),])
colnames(y1) <- c("time", "status")
y1_hat <- data.frame(y1)
x1 <- t(Data_test[-c(1,2),])
x1_hat <- data.frame(x1)
# }

library(survival)
feature_plus <- paste(markerinput,collapse="+")
my_CI(feature_plus, x1_hat)



library(dplyr)   
library(glmSparseNet)
coef_test <- my_overlap(coefinput, Data_test)
plotp <- separate2GroupsCox(as.vector(coef_test[[1]]), x1_hat[, coef_test[[2]]], 
                                  # plot.title = 'Lasso', 
                            as.data.frame(y1), legend.outside = T)

# plotp$plot
plotp$pvalue
# plotp$km

# ci <- c(ci, my_CI(feature_plus, x1_hat))
# p <- c(p, plotp$pvalue)
# }
# ci 
# p


##############################################################################
# Figure 2G. 
##############################################################################
rm(list=ls())
library(openxlsx)
com_data <- read.xlsx("Compare/CNetCoxComResults.xlsx", 
                      sheet = "Plot_revision")

library(dplyr)
library(ggplot2)

## delete Ridge
# plot_dat <- com_data %>%
#   filter(Methods != "Ridge")


library(dplyr)
library(tidyr)
library(ggplot2)

## 1. First, fill down the NA values ​​in the Datasets.
com_data2 <- com_data %>%
  fill(Datasets, .direction = "down")

## 2. Remove Ridge.
plot_dat <- com_data2 %>%
  filter(Methods != "Ridge")

## 3. Set the factor order (optional)
plot_dat$Methods  <- factor(plot_dat$Methods,
                            levels = c("CNetCox","ENet","L0","L1/2","Lasso","MCP","SCAD"))
plot_dat$Datasets <- factor(plot_dat$Datasets,
                            levels = c("Data 1","Data 3","Data 4","Data 6"))



# pdf("Compare/C-index_comparison_7.pdf", width=8, height=5, family="sans")

## 4. Plotting: Each dataset uses a single color, 
## and all methods use the same color for points on that dataset.
ggplot(plot_dat, aes(x = Methods, y = CI)) +
  geom_boxplot(
    outlier.shape = NA,
    width   = 0.55,
    fill    = "grey92",
    color   = "black",
    linewidth = 0.6
  ) +
  geom_jitter(
    aes(color = Datasets),
    width = 0.15,
    size  = 2.5,
    alpha = 0.9
  ) +
  ## set the name of the x-axis display.
  scale_x_discrete(
    labels = c(
      "CNetCox" = "CNet-Cox",
      "ENet"    = "ENet-Cox",
      "L0"      = "L0-Cox",
      "L1/2"    = "L1/2-Cox",
      "Lasso"   = "Lasso-Cox",
      "MCP"     = "MCP-Cox",
      "SCAD"    = "SCAD-Cox"
    )
  ) +
  scale_color_brewer(palette = "Set1", name = "Dataset") +
  ylab("C-index") +
  xlab(NULL) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 0, hjust = 1, vjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

# dev.off()



library(ggplot2)
## Factor order (maintains the same as before)
plot_dat$Methods  <- factor(
  plot_dat$Methods,
  levels = c("CNetCox","ENet","L0","L1/2","Lasso","MCP","SCAD")
)
plot_dat$Datasets <- factor(
  plot_dat$Datasets,
  levels = c("Data 1","Data 3","Data 4","Data 6")
)

### 1）A Num line chart showing all methods on the same graph.
ggplot(plot_dat, aes(x = Datasets, y = Num,
                     group = Methods, color = Methods)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_color_brewer(palette = "Dark2", name = "Methods") +
  xlab("Dataset") +
  ylab("Num") +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.text.x  = element_text(angle = 0, hjust = 0.5),
    legend.position = "right"
  )

### 2）Num line chart with a separate small plot for each method
ggplot(plot_dat, aes(x = Datasets, y = Num, group = 1)) +
  geom_line(linewidth = 0.8, color = "black") +
  geom_point(size = 2.2, color = "black") +
  xlab("Dataset") +
  ylab("Num") +
  facet_wrap(~ Methods, nrow = 2) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.text.x  = element_text(angle = 0, hjust = 0.5)
  )



##############################################################################
# Figure S2. Overlap
##############################################################################

# Overlap between four times for CNet -------------------------------------
marker1 <- read.csv("Result/Result1/marker1_revision.csv",header = T)
marker3 <- read.csv("Result/Result3/marker3.csv",header = T)
marker4 <- read.csv("Result/Result4/marker4_revision.csv",header = T)
marker6 <- read.csv("Result/Result6/marker6_revision.csv",header = T)


genes1 <- unique(marker1$x)
genes3 <- unique(marker3$x)
genes4 <- unique(marker4$x)
genes6 <- unique(marker6$x)

venn_list <- list(
  "Data 1" = genes1,
  "Data 3" = genes3,
  "Data 4" = genes4,
  "Data 5" = genes6
)
head(marker1)
colnames(marker1)

# install.packages("ggvenn") 
library(ggvenn)

pdf("Compare/Feature_Overlap_CNet_venn.pdf", width=5, height=5, family="sans")
ggvenn(
  venn_list,
  fill_color   = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
  stroke_color = "black",
  stroke_size  = 0.5,
  set_name_size = 4,
  text_size    = 3
)
dev.off()


