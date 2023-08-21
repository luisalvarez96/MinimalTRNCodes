# Author: Ian Leifer <ianleifer93@gmail.com> and Luis Alvarez <luisalvarez.10.96@gmail.com>


library(tidyr)
library(dplyr)
library(stringr)
library(foreach)
library(igraph)
library(gridExtra)
library(doParallel)
library(RColorBrewer)

getCircuits <- function(adjMatrix, graph){
  # run code to find circuits
  adjMatrix <- graph_from_adjacency_matrix(adjMatrix, mode = "directed", weighted = T, diag = T)
  circuits <- subgraph_isomorphisms(adjMatrix, graph, method = "lad", induced = T)
  if(length(circuits) == 0){return(NULL)}
  
  circuits <- foreach(i = 1:length(circuits), .combine = rbind) %do% {paste(names(unlist(circuits[i])), collapse = " ")}
  
  # reshape circuits for output
  circuits <- as.data.frame(circuits[, 1], stringsAsFactors = F)
  colnames(circuits)[1] <- "Nodes"
  return(circuits)
}

arrangeLine <- function(lineToSort) {
  lineToSort <- unlist(strsplit(lineToSort, split = " "))
  return(paste(sort(lineToSort), collapse = " "))
}

classify_ARs <- function(circuits_AR, edges){
  classes <- NULL
  for (circuit in circuits_AR$Nodes) {
    nodes <- unlist(strsplit(circuit, split = " "))
    type1 <- (edges %>% filter(Source == nodes[1] & Target == nodes[2]))$Type
    type2 <- (edges %>% filter(Source == nodes[2] & Target == nodes[1]))$Type
    
#    if(length(type1) > 1 | length(type2) > 1){print(paste('circuit:', circuit, ', types:', type1, type2))}
    
    allclasses <- NULL
    for (i in type1) {
      for (j in type2) {
        #onyl works for same labels of edge type 'positive', 'negative' and 'dual'
        if(i != j){class <- 'oscillator'}
        if(i == 'positive' & j == 'positive'){class <- 'lock-on'}
        if(i == 'negative' & j == 'negative'){class <- 'toggle-switch'}
        
        if(i == 'dual' | j == 'dual'){
          if(i == 'negative' | j == 'negative'){class <- c('toggle-switch', 'oscillator')}
          if(i == 'positive' | j == 'positive'){class <- c('lock-on', 'oscillator')}
          if(i == j){class <- 'everything'}
        }
        
        if(length(allclasses)){
          for (c in class) {
            if(!grepl(c, allclasses)){allclasses <- paste(allclasses, c, sep = '/', collapse = '/')}
          }
        } else {
          allclasses <- paste(class, collapse = '/')
        }
      }
    }
    
    classes <- rbind(classes, allclasses)
  }
  circuits_AR$Type <- as.character(classes)

  return(circuits_AR)
}

findCircuits <- function(edges) {
  circuit_AR <- matrix(
    data = c(0, 1,
             1, 0), ncol = 2)
  circuit_FFF <- matrix(
    data = c(0, 0, 0, 0, 0,
             1, 0, 0, 0, 0,
             1, 0, 0, 0, 0,
             0, 1, 0, 0, 1,
             0, 0, 1, 1, 0), ncol = 5)
  circuit_JK <- matrix(
    data = c(0, 0, 0, 0, 0,
             1, 0, 0, 0, 1,
             1, 0, 0, 1, 0,
             0, 1, 0, 0, 1,
             0, 0, 1, 1, 0), ncol = 5)
  
  # returns number of circuits found and writes circuits to files
#  nodes <- unique(c(edges$Source, edges$Target))
  edges <- edges %>% filter(Source != Target)
#  edges <- edges[!duplicated(edges[, 1:2]), ]
#  graph <- graph.data.frame(edges)
  graph <- graph.data.frame(edges[!duplicated(edges[, 1:2]), ])
  
  circuits_AR <- getCircuits(circuit_AR, graph)
  if(length(circuits_AR) != 0){
#    circuits_AR$Type <- 'AR'
    #class as osc or toggle-switch
    circuits_AR <- classify_ARs(circuits_AR, edges)
    }
  
  circuits_FFF <- getCircuits(circuit_FFF, graph)
  if(length(circuits_FFF) != 0){circuits_FFF$Type <- 'FFF'}
  
  circuits_JK <- getCircuits(circuit_JK, graph)
  if(length(circuits_JK) != 0){circuits_JK$Type <- 'JK'}
  
  circuits <- rbind(circuits_AR, circuits_FFF, circuits_JK)

  if(length(circuits) != 0){
    # remove circuits, which are rearrangement of one another
    # duplicationTable <- sapply(circuits$Nodes, function(x) arrangeLine(x))
    # circuits <- circuits[!duplicated(duplicationTable), ]
    circuits <- circuits[!duplicated(sapply(circuits$Nodes, function(x) arrangeLine(x))), ]
    circuits <- as.data.frame(circuits, stringsAsFactors = F)
    colnames(circuits)[1] <- "Nodes"
  }
  
  
  # if(nrow(circuits) != 0 & !pvalues) {
  #   write(unlist(unite(circuits, Output, sep = "\t")),
  #         paste("circuits3/", dataset, "_", name, ".txt", sep = ""))
  # }
  # 
  # if(!pvalues) {
  #   makePngs(dataset, name, circuits$Nodes, edges)
  # }
  return(circuits)
}

findCircuits(edges_e_operon)










