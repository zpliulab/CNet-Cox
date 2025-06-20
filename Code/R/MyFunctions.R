## 2023.4.23 LLY create
## Some functions used in CI and P values of methods



# load coefficient --------------------------------------------------------

markerselect <- function(filenum,cutoff){
  
  # filenum <- 1
  # cutoff <- 2e-1
  
  coef <- read.csv(paste("Result/Result",filenum,"/theta",filenum,".csv", sep = ""),
                   header = F, check.names = FALSE, sep = ",")
  theta0 = coef[1,]  
  theta = as.data.frame(coef[-1,])
  ## Give thro and select feature
  label <- which(abs(theta)/max(abs(theta)) >= cutoff)
  node <- gene[label,1]
  
  ## plot sub figure
  node_used <- node
  net_used <- net
  k1 <- which(net_used[,1] %in% node_used)   # 562
  k2 <- which(net_used[,2] %in% node_used)   # 929
  length(intersect(k1,k2))    # 4
  used <- net_used[intersect(k1,k2),]
  if (length(used) == 0){break;}
  
  library(igraph)
  PP <- graph_from_data_frame(used,directed = F)
  # PP <- graph_from_data_frame(used,directed = T)
  p1 <- simplify(PP)
  ed <- as_edgelist(p1, names = TRUE)
  
  ## save marker
  marker <- union(ed[,1], ed[,2])
  ## save coef
  coefno1 <- coef
  coefno1[-(which( gene[,1] %in% marker)+1), ] <- 0
  coefno1 <- as.matrix(coefno1[-1,])
  rownames(coefno1) <- gene[,1]
  
  return(list(marker,coefno1))
}


# load coefficient --------------------------------------------------------

markerunion <- function(filenum){

  # filenum <- 1
  l <- 4
  
  coef <- read.csv(paste("Result/Result",filenum,"/theta",filenum,".csv", sep = ""),
                   header = F, check.names = FALSE, sep = ",")
  coefabs <- abs(coef)[-1,]
  coefsort <- sort(coefabs)
  value <- coefabs[which(abs(cut[,1]) == 1)]
  
  ## The number of neigbourhold : l
  # l <- 2
  lab <- c()
  for (i in 1:5) {
    num <- which(coefsort == value[i])
    lab <- c(lab, c((num-l):(num+l)) )
  }
  node <- gene[which( abs(coef[-1,1]) %in% coefsort[lab]), ]

  ## plot sub figure
  node_used <- node
  net_used <- net
  k1 <- which(net_used[,1] %in% node_used)   # 562
  k2 <- which(net_used[,2] %in% node_used)   # 929
  length(intersect(k1,k2))    # 4
  used <- net_used[intersect(k1,k2),]
  if (length(used) == 0){break;}
  
  library(igraph)
  PP <- graph_from_data_frame(used,directed = F)
  # PP <- graph_from_data_frame(used,directed = T)
  p1 <- simplify(PP)
  ed <- as_edgelist(p1, names = TRUE)
  
  ## save marker
  marker <- union(ed[,1], ed[,2])
  ## save coef
  coefno1 <- coef
  coefno1[-(which( gene[,1] %in% marker)+1), ] <- 0
  coefno1 <- as.matrix(coefno1[-1,])
  rownames(coefno1) <- gene[,1]
  
  return(list(marker,coefno1))
}



# Plot the selected markers -----------------------------------------------
markerPlot <- function(node_used,net){
  
  # node_used <- marker6
  
  net_used <- net
  k1 <- which(net_used[,1] %in% node_used)   # 562
  k2 <- which(net_used[,2] %in% node_used)   # 929
  length(intersect(k1,k2))    # 4
  used <- net_used[intersect(k1,k2),]
  if (length(used) == 0){break;}
  
  library(igraph)
  PP <- graph_from_data_frame(used,directed = F)
  # PP <- graph_from_data_frame(used,directed = T)
  p1 <- simplify(PP)  
  ed <- as_edgelist(p1, names = TRUE)
  return(list(p1,ed))
  
}


# my_CI function ----------------------------------------------------------
my_CI <- function(feature_plus, data){
  univ_formulas <- sapply(feature_plus, 
                          function(x) as.formula(paste('Surv(y1_hat$time, y1_hat$status)~', x)))
  univ_models <- lapply(univ_formulas, function(x){coxph(x, data = data)})
  CI <- univ_models[[1]]$concordance[6]
  std <- univ_models[[1]]$concordance[7]
  return(CI)
}


# my_overlap_coef ------------------------------------------------------
my_overlap <- function(x, y){
  coefs.v <- x[,1] %>% { .[. != 0]}
  coefs.v %>% {
    data.frame(gene.name   = names(.),
               coefficient = .,
               stringsAsFactors = FALSE)
  } %>%
    arrange(gene.name) %>%
    knitr::kable()
  
  sele <- rownames(as.matrix(coefs.v))
  gene <- rownames(y[-c(1,2),])
  overlap <- intersect(sele, gene)
  
  lab <- x[,1] %>% { .[. != 0]} %>% as.matrix
  coefs.v <- lab[overlap,]
  
  my <- list(coefs.v, overlap)
  return(my)
} 



