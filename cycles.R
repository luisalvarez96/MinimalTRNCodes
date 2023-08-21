#Author Luis Alvarez luisalvarez.10.96@gmail.com

library(tidyr)
library(dplyr)
library(igraph)


#Do same for all SCCs
Cycles <- NULL
for (sccid in SCCids) {
  #Focus one by one in each SCC
#  sccg <- induced.subgraph(g, nodes$SCCId == 23)
  sccg <- induced.subgraph(g, nodes$SCCId == sccid)
  ids <- 1:vcount(sccg)
  
  #Cycles starting in every node in ascending order of vertex id , eliminating each node as taken into account
#  Cycles <- NULL
  while(vcount(sccg) != 1){
    print(paste('nodes = ', vcount(sccg)))
    #All cycles of smallest (vertexid) node 
    nodeid <- ids[1]
    temp_cycles <- NULL
    for (j in neighbors(sccg, nodeid, mode = 'out')) {
      print(paste('node =', get.vertex.attribute(sccg, "name", nodeid), 'neighbor = ', get.vertex.attribute(sccg, "name", j)))
      temp_cycles <- c(temp_cycles, lapply(all_simple_paths(sccg, j, nodeid, mode="out"), function(p) c(nodeid,p)))
      print(paste('temp_cycles =', length(temp_cycles)))
    }
      
    #add the new cycles to list of all cycles
    cycles <- NULL
    if(length(temp_cycles) != 0){
      for (j in 1:length(temp_cycles)) {
        cycles$SCCId <- sccid
        cycles$Path[j] <- paste(data.frame(c = get.vertex.attribute(sccg, "name", temp_cycles[[j]]))$c, collapse = ' ')
        cycles$Length[j] <- length(temp_cycles[[j]]) - 1
      }
    }
    Cycles <- rbind(Cycles, data.frame(cycles))
    print(paste('Cycles =', nrow(Cycles)))
    
    #Remove node i just taken into account
    sccg <- induced.subgraph(sccg, ids[-1])
    ids <- 1:vcount(sccg)
  }
  
  #View(Cycles)
  
}

getcycles_fromSCC <- function(graph, nodes, sccid){
  sccg <- induced.subgraph(graph, nodes)
  ids <- 1:vcount(sccg)
  Cycles <- NULL
  
  print(paste(vcount(sccg), 'nodes in SCC'))
  
  #Cycles starting in every node in ascending order of vertex id , eliminating each node as taken into account
  #  Cycles <- NULL
  while(vcount(sccg) != 1){
    #print(paste('nodes = ', vcount(sccg)))

    if(vcount(sccg) > 20 & vcount(sccg) %% 5 == 0){print(paste(vcount(sccg), 'nodes left in the SCC'))}
    
    #All cycles of smallest (vertexid) node 
    nodeid <- ids[1]
    
    neighbors <- neighbors(sccg, nodeid, mode = 'out')
    if(length(neighbors) > 10){print(paste(length(neighbors), 'neighbors'))}
    
    temp_cycles <- NULL
    for (j in neighbors) {
      #print(paste('node =', get.vertex.attribute(sccg, "name", nodeid), 'neighbor = ', get.vertex.attribute(sccg, "name", j)))
      temp_cycles <- c(temp_cycles, lapply(all_simple_paths(sccg, j, nodeid, mode="out"), function(p) c(nodeid,p)))
      #print(paste('temp_cycles =', length(temp_cycles)))
    }
    
    #add the new cycles to list of all cycles
    cycles <- NULL
    if(length(temp_cycles) != 0){
      for (j in 1:length(temp_cycles)) {
        cycles$SCCId <- sccid
        cycles$Length[j] <- length(temp_cycles[[j]]) - 1
        cycles$Path[j] <- paste(data.frame(c = get.vertex.attribute(sccg, "name", temp_cycles[[j]]))$c, collapse = ' ')
      }
    }
    Cycles <- rbind(Cycles, data.frame(cycles, stringsAsFactors = F))
    #print(paste('Cycles =', nrow(Cycles)))
    
    #Remove node i just taken into account
    sccg <- induced.subgraph(sccg, ids[-1])
    ids <- 1:vcount(sccg)
  }
  return(Cycles)
}

getcycles_fromEdges <- function(edges){
  #edges <- edges %>% filter(Source != Target)
  nodes <- as.data.frame(sort(unique(c(edges$Source, edges$Target))))
  colnames(nodes)[1] <- "Label"
  g <- graph_from_data_frame(d = edges, vertices = nodes, directed = 1)

  sc <- clusters(g, mode = "strong")
  if(length(sc$csize)){
    nodes$SCCId <- sc$membership
    SCCids <- NULL
    for (i in 1:length(sc$csize)) {
      if(sc$csize[i] > 1){SCCids <- c(SCCids, i)}
    }
    
    Cycles <- NULL
    #Find cycles in each SCC 
    for (sccid in SCCids){
      Cycles <- rbind(Cycles, getcycles_fromSCC(g, nodes$SCCId == sccid, sccid))
    }
    Cycles <- Cycles[!duplicated(Cycles),]
    return(Cycles)
  } else {
    return(NULL)
  }
  
}

##### Ecoli
dir <- "/home/luis/Dropbox (Graduate Center)/Fibers/Bacillus/ReproducciónEcoli"

edges <- read.csv("Ecoli.txt", stringsAsFactors = F, sep = ' ', header = F)
colnames(edges) <- c('Source', 'Target', 'Color')

edges <- read.csv("edges_col_kout.csv", stringsAsFactors = F, header = T, sep = '\t')

nodes <- read.csv("nodes_col_kout.csv", stringsAsFactors = F, header = T, sep = '\t')
for(i in 1:nrow(edges)){
  edges$Source[i] <- nodes$Label[nodes$Id == edges$Source[i]]
  edges$Target[i] <- nodes$Label[nodes$Id == edges$Target[i]]
}

# edges$Source <- edges$Source + 1
# edges$Target <- edges$Target + 1

Cycles <- getcycles_fromEdges(edges)
Cycles <- Cycles[!duplicated(Cycles),]

#### Bacillus
dir <- "/home/luis/Dropbox (Graduate Center)/Fibers/Bacillus/bacilus/Sigma[bac2]"
edges <- read.csv("edges_col_kout.csv", stringsAsFactors = F, header = T, sep = '\t')
nodes <- read.csv("nodes_col_kout.csv", stringsAsFactors = F, header = T, sep = '\t')
for(i in 1:nrow(edges)){
  edges$Source[i] <- nodes$Label[nodes$Id == edges$Source[i]]
  edges$Target[i] <- nodes$Label[nodes$Id == edges$Target[i]]
}

Cycles_bac <- getcycles_fromEdges(edges)

write.table(Cycles_bac, file = '/home/luis/Dropbox (Graduate Center)/Fibers/Bacillus/Cycles/bac_cycles.csv', quote = F, row.names = F, sep = "\t", col.names = T)

####

for (i in 1:nrow(files)) {
  print(paste('File = ', files$Input[i]))
  main(files$Input[i], ' ', OutputFile = files$Output[i])
}

#### idea for faster:
#### start with node with less out-degree then from neighbor with less out-degree also. 

