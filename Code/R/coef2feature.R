## 2025.06.01 Copy from SVM and used  for CNet-RCPH

###########  theta_selection (select features from the coefficients theta: matlab)  ############

rm(list=ls())


## Load features of each method
setwd('D:\\OneDrive - The University Of Hong Kong\\R\\CNetCox\\Result\\ComMethd')  
feature_list <- read.csv('CompareMethod.csv', header = T, sep=',')[,1:7]
colnames(feature_list) <- c('Lasso-RCPH', 'ENet-RCPH', 'L0-RCPH', 'L1/2-RCPH', 'SCAD-RCPH', 'MCP-RCPH', 'CNet-RCPH')
feature <-feature_list['CNet-RCPH']  


# feature network ------------------------------------------------------------
setwd('D:\\OneDrive - The University Of Hong Kong\\R\\CNetCox\\Data\\TCGA_NEW')  
net_used <- as.matrix(read.csv("UNgene_comp_net.csv",header = T))
node_used <- as.matrix(feature)

k1 <- which(net_used[,1] %in% node_used) 
k2 <- which(net_used[,2] %in% node_used) 
length(intersect(k1,k2))   
used <- net_used[intersect(k1,k2),]


library(igraph)
PP <- graph_from_data_frame(used,directed = F)
p1 <- simplify(PP)
# p1 <- simplify(PP,remove.multiple = TRUE, remove.loops = TRUE) 
ed <- as_edgelist(p1, names = TRUE)


g <- p1
plot(g, layout=layout.fruchterman.reingold,  
     vertex.size = 8,   
     vertex.label = V(g)$name, 
     vertex.label.cex = 0.8, 
     vertex.label.dist = 0.1,  
     vertex.label.color = "black"   
)


# The gene in the figure is used as a biomarker ----------------------------------------------------
nodename <- get.vertex.attribute(g)
feature_new <- as.matrix(nodename$name)
colnames(feature_new) <- c('biomarker')
# View(feature_new)

setwd('D:\\OneDrive - The University Of Hong Kong\\R\\CNetCox\\Result\\ComMethd') 
write.csv(ed,"feature_net_CNet.csv",row.names = F,quote = F)
write.csv(feature_new,"feature_CNet.csv",row.names = F,quote = F)



# Extracting the largest network ------------------------------------------------------------------
clu <- components(p1)
groups(clu)

component <- groups(clu)$'1'
# write.csv(component, file = "result\\feature\\feature_component.csv", row.names = F)

node_used <- as.matrix(component)
k1 <- which(net_used[,1] %in% node_used) 
k2 <- which(net_used[,2] %in% node_used) 
length(intersect(k1,k2))   
used <- net_used[intersect(k1,k2),]


library(igraph)
PP <- graph_from_data_frame(used,directed = F)
p1 <- simplify(PP)
# p1 <- simplify(PP,remove.multiple = TRUE, remove.loops = TRUE) 
ed <- as_edgelist(p1, names = TRUE)

g <- p1
plot(g, layout=layout.fruchterman.reingold,  
     vertex.size = 8,  
     vertex.label = V(g)$name, 
     vertex.label.cex = 0.8,  
     vertex.label.dist = 0.1,  
     vertex.label.color = "black"   
)