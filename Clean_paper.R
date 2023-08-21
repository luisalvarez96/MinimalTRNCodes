#Author Luis Alvarez luisalvarez.10.96@gmail.com


library(tidyr)
library(dplyr)
library(igraph)

#### Directories ####
###changed directory names ###

## bac2
#Full w/sigma spacing="\t"
Nodes <- "~/Documents/Fibers/Bacillus/bacilus/Sigma[bac2]/Results_nodes.csv"
Edges <- "~/Documents/Fibers/Bacillus/bacilus/Sigma[bac2]/bacilus_activationRepressionFull.txt"

#spacing = "\t
Nodes <- "~/Documents/Fibers/Bacillus/bacilus/Sigma[bac2]/bacilus_Results_nodes.csv"
Edges <- "~/Documents/Fibers/Bacillus/bacilus/Sigma[bac2]/bacilus.txt"
##
   
#Full_noSigma ##  spacing="\t"  bac
Nodes <- "~/Documents/Fibers/Bacillus/bacilus/NoSigma[bac]/bacilus_Results_nodes.csv"
Edges <- "~/Documents/Fibers/Bacillus/bacilus/NoSigma[bac]/bacilus_activationRepressionFull_noSigma.txt"

#Sigma spacing=" "   bac1
Nodes <- "~/Documents/Fibers/Bacillus/bac_Mishael[bac1]/bacilus_Results_nodes.csv"
Edges <- "~/Documents/Fibers/databases/B_subtilis_quantitative_transcription_network.txt"


#
#ecoli spacing=" "
Nodes <- "~/Dropbox (Graduate Center)/Fibers/Bacillus/ReproducciónEcoli/Ecoli_Results_nodes.csv"
Edges <- "~/Dropbox (Graduate Center)/Fibers/Bacillus/ReproducciónEcoli/Ecoli.txt"

Nodes <- "~/Documents/Fibers/Bacillus/ReproducciónEcoli/nodual_Results_nodes.csv"
Edges <- "~/Documents/Fibers/Bacillus/ReproducciónEcoli/ecoli_nodual.txt"
#

#regulondb_ian spacing="\t"
Nodes <- "~/Documents/Fibers/Bacillus/regulondb/ian/ian_results_nodes.csv"
Edges <- "~/Documents/Fibers/Bacillus/regulondb/ian/fullnetwork_ian.txt"


#### Plan ####
# 1) Collapse to fibers -> collapseFibers(...)
# 2) kout-0 core (outdegree & Core > 0 | 1)?   ->   nodes_plt <- filter(Out & Core > ?)   ->   getEdges_plt(...)
# 3) plot in gephi

#### Functions ####
read_network <- function(nodes, edges, spacing){
  #reads from file if given file path, if not then give nodes data frame
  if(is.character(nodes)){nodes <- read.csv(nodes, stringsAsFactors = F, sep = "\t")}
  nodes$Id = 1:nrow(nodes)
  
  ### Carefull with format ###
  #same as reading nodes
  if(is.character(edges)){edges <- read.csv(edges, stringsAsFactors = F, sep = spacing, header = F)}
  colnames(edges)[1:2] <- c("Source", "Target")
  if(ncol(edges) == 3){colnames(edges)[3] <- "Type"}
  
  #doesn't work if there is a node 0  -> use .txt not _Results.csv
  graph <- graph_from_edgelist(as.matrix(edges[, 1:2]), directed = TRUE)

  nodes <- nodes %>%
    group_by(FiberId) %>%
    mutate(FiberSize = n()) %>%
    ungroup()
 
  nodes <- nodes_info(nodes, edges)
  
  return(list(nodes, edges))
}

#collapsing on fibers: first only looks at FIberId and SCCSize, after collapsing adds Out,In & Coreness

collapseFibers <- function(nodes, edges){

  nodes_col <- c() 
  for(i in unlist(nodes %>% group_by(FiberId) %>% summarise())){
    fiber <- nodes %>% filter(FiberId == i)
    SCC <- fiber %>% filter(SCCSize > 1)
    if(nrow(SCC) > 1){
      print(paste("SCCs in fiber ", i, ": ", SCC$Label))
    }
    
    ##For including all SCCs of same fiber: 
  #   if(nrow(SCC) == 1){
  #     nodes_col <- rbind(nodes_col, SCC[1,])
  #   } else if(nrow(SCC) > 1){
  #     # if various SCCs per fiber, take all SCCs
  #     nodes_col <- rbind(nodes_col, SCC)
  #   } else {
  #     nodes_col <- rbind(nodes_col, fiber[1,])
  #   }
  # }
    
  if(nrow(SCC)){
    nodes_col <- rbind(nodes_col, SCC[1,])
  } else {
    nodes_col <- rbind(nodes_col, fiber[1,])
  }
}
  nodes_col <- nodes_col[complete.cases(nodes_col),]
  # colnames(nodes_col)[2] <- "Id"
  
  edges_col <- edges
  for(i in 1:nrow(edges)){
    
    if(!edges_col$Source[i] %in% nodes_col$Label){
      edges_col$Source[i] <- nodes_col$Label[ nodes_col$FiberId == nodes$FiberId[ nodes$Label == edges_col$Source[i] ] ]
    }
    if(!edges_col$Target[i] %in% nodes_col$Label){
      edges_col$Target[i] <- nodes_col$Label[ nodes_col$FiberId == nodes$FiberId[ nodes$Label == edges_col$Target[i] ] ]
    }
    

  }
  edges_col <- unique(edges_col)

  nodes_col <- nodes_info(nodes_col, edges_col)
  
  return(list(nodes_col, edges_col))
}

nodes_info <- function(nodes, edges){
  #doesn't work if there is a node 0  -> use .txt not _Results.csv
  #edges need to be label -> label, not ids
  graph <- graph_from_edgelist(as.matrix(edges[, 1:2]), directed = TRUE)
  
  SCCIds <- data.frame(SCCId = components(graph, mode = "strong")$membership)
  OutDegree <- data.frame(Out = degree(graph, mode = "out", loops = F, normalized = F))
  InDegree <- data.frame(In = degree(graph, mode = "in", loops = F, normalized = F))
  Coreness <- data.frame(Core = coreness(graph, mode = "out"))
  WeakIds <- data.frame(WeakId = components(graph, mode = "weak")$membership)
  
  for(i in 1:nrow(nodes)){
    nodes$SCCId[i] <- SCCIds[nodes$Label[i],]
    nodes$OutDegree[i] <- OutDegree[nodes$Label[i],]
    nodes$InDegree[i] <- InDegree[nodes$Label[i],]
    nodes$Coreness[i] <- Coreness[nodes$Label[i],]
    nodes$WeakId[i] <- WeakIds[nodes$Label[i],]
  }
  
  nodes <- nodes %>%
    group_by(SCCId) %>%
    mutate(SCCSize = n()) %>%
    ungroup() 

  nodes <- nodes %>%
    group_by(WeakId) %>%
    mutate(WeakSize = n()) %>%
    ungroup() 
  
  return(nodes)  
}



getEdges_plt <- function(edges, nodes_plt){
  edges_plt <- edges[edges$Target %in% nodes_plt$Label & edges$Source %in% nodes_plt$Label,]
  return(edges_plt)
}

gephi_format <- function(nodes, edges){
  #format to gephi
  for(i in 1:nrow(edges)){
    edges$Source[i] <- nodes$Id[nodes$Label == edges$Source[i]]
    edges$Target[i] <- nodes$Id[nodes$Label == edges$Target[i]]
  }
  return(edges)
}


filter_kshell <- function(nodes, edges){
  #takes nodes after taking Core > 0 & Out > 0  || edges need to be label -> label, not ids
  
  nodes_plt <- nodes_info(nodes, edges)
  nodes_plt <- nodes_plt %>% filter(OutDegree > 0) %>% filter(Coreness > 0)
  # nodes_plt <- nodes_plt %>% filter(Coreness > 0)
  edges_plt <- getEdges_plt(edges, nodes_plt)
  # for (node in nodes_plt$Label) {
  #   if(!node %in% edges_plt$Target & !node %in% edges_plt$Source){
  #     nodes_plt <- nodes_plt[nodes_plt$Label != node,]
  #   }
  # }
  
  return(list(nodes_plt, edges_plt))
}



find_core <- function(nodes, edges){
  print(paste("initial = ",nrow(nodes)))
  kshell <- filter_kshell(nodes, edges)
  print(paste("initial core=", nrow(kshell[[1]])))
  while(nrow(nodes) != nrow(kshell[[1]])){
    nodes <- kshell[[1]]
    edges <- kshell[[2]]
    kshell <- filter_kshell(nodes, edges)
    print(paste("core = ",nrow(kshell[[1]])))
  }
  return(list(kshell[[1]], kshell[[2]]))
}

#### Exe ####

minimal <- function(Nodes, Edges, spacing, suffix, order = 1, prnt = F) {
  library(tidyr)
  library(dplyr)
  library(igraph)
  
  if(order){
    #first collapse fibers and then kout
    
    print('Reading network')
    network_list <- read_network(Nodes, Edges, spacing)
    
    print('Collapsing network')
    network_col <- collapseFibers(network_list[[1]], network_list[[2]])
    
    print('Finding the minimal network')
    minimal_network <- find_core(network_col[[1]], network_col[[2]])
    
    
  } else{
    #first  kout and then collapse fibers
    
    network_list <- read_network(Nodes, Edges, spacing)
    
    network_kout <- find_core(network_list[[1]], network_list[[2]])
    
    minimal_network <- collapseFibers(network_kout[[1]], network_kout[[2]])
    
    
  }
  
  if(prnt){
    write.table(minimal_network[[1]], file = paste("nodes_", suffix, ".csv", sep = ""), quote = F, row.names = F, sep = "\t")
    write.table(gephi_format(minimal_network[[1]],  minimal_network[[2]]), file = paste("edges_", suffix, ".csv", sep = ""), quote = F, row.names = F, sep = "\t")
  }
  
  list <- NULL
  list$Nodes <- minimal_network[[1]]
  list$Edges <- minimal_network[[2]]
  return(list)
  #return(list(minimal_network[[1]], minimal_network[[2]]))
}





