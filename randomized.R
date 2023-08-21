#Author Luis Alvarez luisalvarez.10.96@gmail.com

library(tidyr)
library(dplyr)
library(igraph)


current <- "/home/luis/Dropbox (Graduate Center)/Fibers/Bacillus"
setwd("/home/luis/Dropbox (Graduate Center)/Fibers/Bacillus")

#Viewing original data
hist(((nodes_e_operon %>% filter(FiberSize > 1) %>% group_by(FiberId)) %>% summarise(n = n()))$n, main = 'Original', xlab = 'Fiber Size')
View(nodes_e_operon_col_kout %>% group_by(SCCSize) %>% summarise(n = n()))

generate_random <- function(edges, iterations = 10000){
  colnames(edges)[1:2] <- c('Source','Target')
  rand_edges <- edges %>% filter(Source != Target)
  
  N <- iterations
  while (N > 0) {
    if(N%%500 == 0){print(paste('N =', N))}
    entry1 <- sample(nrow(rand_edges),1)
    entry2 <- sample(nrow(rand_edges),1)
    edge1 <- rand_edges[entry1,]
    edge2 <- rand_edges[entry2,]
    
    #avoid converting into self-loops
    while (edge1[1] == edge2[2] | edge2[1] == edge1[2]) {
      entry2 <- sample(nrow(rand_edges),1)
      edge2 <- rand_edges[entry2,]
      print('Creats self-loop, drawing again')
    }
    
    
    rand_edges[entry1,1] <- edge2[1]
    rand_edges[entry2,1] <- edge1[1]
    
    
    N <- N-1
  }
  rand_edges <- rbind(rand_edges, edges %>% filter(Source == Target))
  
  return(rand_edges)
}

setwd("~/Dropbox (Graduate Center)/Fibers")
write.table(rand_edges, file = 'randomized.txt', quote = F, row.names = F, sep = "\t", col.names = F)

analyze_graph <- function(file_name){
  fibs <- main(file_name, '\t')
  rand_nodes <- read_network(fibs$Nodes, file_name, '\t')[[1]]
#  write.table(rand_nodes, file = 'randomized_nodes.csv', quote = F, row.names = F, sep = "\t", col.names = T)
  
  rand_core <- minimal(rand_nodes, file_name, '\t')
  rand_nodes_core <- rand_core[[1]]
  rand_edges_core <- rand_core[[2]]
  rand_nodes$Min <- rand_nodes$Label %in% rand_nodes_core$Label
  rand_nodes$FiberId[rand_nodes$FiberSize == 1] <- ''
#  write.table(rand_nodes, file = 'randomized_nodes.csv', quote = F, row.names = F, sep = "\t", col.names = T)

  list <- NULL
  list$Nodes <- rand_nodes
  list$Nodes_Core <- rand_nodes_core
  list$Edges_Core <- rand_edges_core
  list$Blocks <- fibs$Blocks
  list$Circuits <- findCircuits(rand_edges_core)
  list$Cycles <- getcycles_fromEdges(rand_edges_core)

  return(list)
}

analyze_set <- function(path){
  # created distributions aggregate
  start_time <- Sys.time()
  
  files <- list.files(path, full.names = T)
#  print(files)  
  FiberSizes <- NULL
  BBClasses <- NULL
  SCCs <- NULL
  Circuits <- NULL
  Cycles <- NULL
  for (file in files) {
    id <- as.integer(gsub('[^0-9]', '', file))
    if(is.na(id)){id <- 0}
    cat(paste('\n', file, '\n'))
    analysis <- analyze_graph(file)
    
    FiberSizes <- c(FiberSizes, (analysis$Nodes %>% filter(FiberSize > 1) %>% group_by(FiberId) %>% summarise(n = n()))$n)
    
    analysis$Blocks$GraphId <- id
    BBClasses <- rbind(BBClasses, analysis$Blocks)
    
    sccs <- analysis$Nodes_Core %>% group_by(SCCSize) %>% summarise(NumberSCCs = n()) %>% mutate(NumberSCCs = NumberSCCs/SCCSize)
    sccs$GraphId <- id
    SCCs <- rbind(SCCs, sccs)

    if(!is.null(analysis$Circuits)){
      analysis$Circuits$GraphId <- id
      Circuits <- rbind(Circuits, analysis$Circuits)
    } # else { rbind(Circuits, data.frame(zeros,id))}
    
    if(!is.null(analysis$Cycles)){
      analysis$Cycles$GraphId <- id
      Cycles <- rbind(Cycles, analysis$Cycles)
    } # else { rbind(Circuits, data.frame(zeros,id))}
    
    if(id %% 5 == 0){
      print('Updating files')
      print(Sys.time() - start_time)
      
      write.table(FiberSizes, file = 'Randomized_FiberSizes.csv', quote = F, row.names = F, sep = "\t")
      write.table(BBClasses, file = 'Randomized_BlocksClass.csv', quote = F, row.names = F, sep = "\t")
      write.table(SCCs, file = 'Randomized_SCCs.csv', quote = F, row.names = F, sep = "\t")
      write.table(Circuits, file = 'Randomized_Circuits.csv', quote = F, row.names = F, sep = "\t")
      write.table(Cycles, file = 'Randomized_Cycles.csv', quote = F, row.names = F, sep = "\t")
    }
  }
  
  set <- NULL
  set$FiberSizes <- FiberSizes
  set$BBClasses <- BBClasses
  set$SCCs <- SCCs
  set$Circuits <- Circuits
  set$Cycles <- Cycles
  
  print('Finished')
  print(Sys.time() - start_time)
  
  return(set)
}


#### Ecoli ####

#generate
for (i in 1:99) {
  file <- paste('randomized/randomized', i,'.txt', sep = '')
  write.table(generate_random(edges_e_operon), file = file, quote = F, row.names = F, sep = "\t", col.names = F)
}

set <- analyze_set('/home/luis/Dropbox (Graduate Center)/Fibers/randomized')

#fiber sizes
hist(set$FiberSizes, main = 'Randomized', xlab = 'Fiber Size')

#Fibers classes
View(set$BBClasses %>% group_by(GraphId, nl) %>% summarise(N = n()) %>% group_by(nl) %>% summarise(sum(N)))

Blocks <- set$BBClasses %>% group_by(GraphId, nl) %>% summarise(N = n())
sd(c(Blocks$N[Blocks$nl == 'n = 0, l = 3'], integer(92)))
sd(c(Blocks$N[Blocks$nl == 'n = 1, l = 1'], integer(9)))
sd(c(Blocks$N[Blocks$nl == 'n = 1, l = 2'], integer(96)))
sd(c(Blocks$N[Blocks$nl == 'Fibonacci'], integer(86)))
sd(c(Blocks$N[Blocks$nl == 'Multi-layered Fiber'], integer(26)))

ones <- set$BBClasses %>% filter(nl == 'n = 0, l = 1') %>% select(Fiber)
size <- NULL
for (fiber in ones) {
  size <- sum(size, length(unlist(strsplit(fiber, split = ', '))))
}
size/nrow(ones) #5.17 ecoli=5.3


#SCCs breakdown
set$SCCs <- set$SCCs %>% filter(SCCSize > 1)
set$SCCs <- set$SCCs %>% group_by(GraphId) %>% mutate(N= sum(NumberSCCs))

View(1:99 %in% set$SCCs$GraphId)

SCCs <- set$SCCs
#SCCs[158,] <- data.frame(SCCSize = 0, NumberSCCs = 0, GraphId = 88, N = 0)
SCCs %>% group_by(GraphId) %>% summarise(N=sum(NumberSCCs)) %>% summarise(sum(N))
SCCs %>% group_by(GraphId) %>% summarise(N=sum(NumberSCCs)) %>% summarise(sd(c(N,0))) #graph 88 no SCC/cycles

mean(SCCs$SCCSize*SCCs$NumberSCCs)
sd(SCCs$SCCSize*SCCs$NumberSCCs)

mean((set$SCCs %>% filter(N == 1))$SCCSize)
mean((set$SCCs %>% filter(N != 1))$SCCSize)

set$SCCs %>% group_by(GraphId) %>% summarise(N=sum(NumberSCCs)) %>% group_by(N) %>% summarise(n())
SCCs %>% filter(N != 1) %>% summarise(S =sum(SCCSize*NumberSCCs)) %>% summarise(sum(S)/107)

#Circuits 
View(1:99 %in% set$Circuits$GraphId)

Circuits <- set$Circuits
Circuits$Type[Circuits$Type == 'oscillator/lock-on'] <- 'lock-on/oscillator'
Circuits$Type[Circuits$Type == 'oscillator/toggle-switch'] <- 'toggle-switch/oscillator'

set$Circuits %>% group_by(GraphId) %>% summarise(N = n()) %>% summarise(sum(N)/100)
set$Circuits %>% group_by(GraphId) %>% summarise(N = n()) %>% summarise(sd(c(N, integer(14)))) #14 graphs w/o circuits

Types <- Circuits %>% group_by(GraphId, Type) %>% summarise(N = n())
Types %>% group_by(Type) %>% summarise(sum(N))

sd(c(Types$N[Types$Type == 'toggle-switch'], integer(67)))
sd(c(Types$N[Types$Type == 'toggle-switch/oscillator'], integer(87)))
sd(c(Types$N[Types$Type == 'oscillator'], integer(47)))
sd(c(Types$N[Types$Type == 'lock-on'], integer(62)))
sd(c(Types$N[Types$Type == 'lock-on/oscillator'], integer(83)))
sd(c(Types$N[Types$Type == 'FFF'], integer(87)))


#Cycles
set$Cycles %>% group_by(GraphId) %>% summarise(N = n()) %>% summarise(sum(N))
set$Cycles %>% group_by(GraphId) %>% summarise(N = n()) %>% summarise(sd(c(N,0))) #graph 88 no SCC/cycles
                                                                      
mean(set$Cycles$Length)
sd(set$Cycles$Length)

set$Cycles %>% group_by(GraphId, SCCId) %>% summarise(N = n()) %>% ungroup() %>% summarise(mean(N))



##### Bacillus ####

#generate
for (i in 0:99) {
  file <- paste('randomized_bac/randomized_bac', i,'.txt', sep = '')
  write.table(generate_random(edges_bac2), file = file, quote = F, row.names = F, sep = "\t", col.names = F)
}

tic()
set_bac <- analyze_set('/home/luis/Dropbox (Graduate Center)/Fibers/Bacillus/randomized_bac')
toc()
#2h 10min

#fiber sizes
hist(set_bac$FiberSizes, main = 'Randomized', xlab = 'Fiber Size')

#Fibers classes
View(set_bac$BBClasses %>% group_by(GraphId, nl) %>% summarise(N = n()) %>% group_by(nl) %>% summarise(sum(N)))

Blocks <- set_bac$BBClasses %>% group_by(GraphId, nl) %>% summarise(N = n())
View(Blocks %>% group_by(nl) %>% summarise(sum(N)))
View(Blocks %>% group_by(nl) %>% summarise(sd(N)))
sd(c(Blocks$N[Blocks$nl == 'n = 0, l = 3'], integer(44)))
sd(c(Blocks$N[Blocks$nl == 'n = 1, l = 2'], integer(89)))
sd(c(Blocks$N[Blocks$nl == 'n > 1'], integer(99)))
sd(c(Blocks$N[Blocks$nl == 'Fibonacci'], integer(45)))
sd(c(Blocks$N[Blocks$nl == 'Multi-layered Fiber'], integer(1)))




#SCCs breakdown
set_bac$SCCs <- set_bac$SCCs %>% filter(SCCSize > 1)
set_bac$SCCs <- set_bac$SCCs %>% group_by(GraphId) %>% mutate(N= sum(NumberSCCs))

View(0:99 %in% set_bac$SCCs$GraphId) # all with SCC/cycles

SCCs <- set_bac$SCCs
#SCCs[158,] <- data.frame(SCCSize = 0, NumberSCCs = 0, GraphId = 88, N = 0)
SCCs %>% group_by(GraphId) %>% summarise(N=sum(NumberSCCs)) %>% summarise(sum(N))
SCCs %>% group_by(GraphId) %>% summarise(N=sum(NumberSCCs)) %>% summarise(sd(N)) 

mean(SCCs$SCCSize*SCCs$NumberSCCs)
sd(SCCs$SCCSize*SCCs$NumberSCCs)

mean((set_bac$SCCs %>% filter(N == 1))$SCCSize)
mean((set_bac$SCCs %>% filter(N != 1))$SCCSize)

set_bac$SCCs %>% group_by(GraphId) %>% summarise(N=sum(NumberSCCs)) %>% group_by(N) %>% summarise(n())

#Circuits 
View(0:99 %in% set_bac$Circuits$GraphId) # graphs 12 and 33 no circuits

set_bac$Circuits %>% group_by(GraphId) %>% summarise(N = n()) %>% summarise(sum(N)/100)
set_bac$Circuits %>% group_by(GraphId) %>% summarise(N = n()) %>% summarise(sd(c(N, integer(2)))) #2 graphs w/o circuits

Circuits <- set_bac$Circuits
Circuits$Type[Circuits$Type == 'oscillator/lock-on'] <- 'lock-on/oscillator'
Circuits$Type[Circuits$Type == 'oscillator/toggle-switch'] <- 'toggle-switch/oscillator'

Types <- Circuits %>% group_by(GraphId, Type) %>% summarise(N = n())
Types %>% group_by(Type) %>% summarise(sum(N))

sd(c(Types$N[Types$Type == 'toggle-switch/oscillator'], integer(63)))
sd(c(Types$N[Types$Type == 'oscillator'], integer(11)))
sd(c(Types$N[Types$Type == 'FFF'], integer(73)))


#Cycles
View(0:99 %in% set_bac$Cycles$GraphId) # all cycles

set_bac$Cycles %>% group_by(GraphId) %>% summarise(N = n()) %>% summarise(sum(N))
set_bac$Cycles %>% group_by(GraphId) %>% summarise(N = n()) %>% summarise(sd(N)) 

mean(set_bac$Cycles$Length)
sd(set_bac$Cycles$Length)


#### stuff ####
three_sccs <-  unlist(set$SCCs %>% filter(N > 2) %>% group_by(GraphId) %>% summarise())

library(tictoc)
tic()
for (g in three_sccs) {
  file <- paste('/home/luis/Dropbox (Graduate Center)/Fibers/randomized/randomized', g,'.txt', sep = '')
#  print(file)
  analysis <- analyze_graph(file)
  file <- paste('randomized', g, sep = '')
#  print(file)
#  print(paste(file, '_nodes.csv', sep = ''))
  write.table(analysis$Nodes_Core, file = paste(file, '_nodes.csv', sep = ''), quote = F, row.names = F, sep = "\t", col.names = T)
  write.table(analysis$Edges_Core, file = paste(file, '_core.txt', sep = ''), quote = F, row.names = F, sep = "\t", col.names = F)
}
toc()

####


SizeScatter <- SCCs %>% group_by(GraphId) %>% summarise(Size = max(SCCSize))
SizeScatter <- cbind(SizeScatter, set$Cycles %>% group_by(GraphId, SCCId) %>% summarise(N = n(), length = mean(Length)) %>% summarise(Cycles = max(N)) %>% select(Cycles)) #pierdo length
View(SizeScatter %>% filter(9 < Size) %>% filter(Size < 16))
SizeScatter %>% filter(9 < Size) %>% filter(Size < 16) %>% summarise(mean(Cycles))

# ecoli
sizeVcycles <- set$Cycles %>% group_by(GraphId, SCCId) %>% summarise(Cycles = n(), Length = mean(Length))
sizeVcycles <- sizeVcycles[-38,]
sizeVcycles <- sizeVcycles[-64,]
sizeVcycles <- sizeVcycles[-101,]
sizeVcycles <- sizeVcycles[-153,]
sizeVcycles <- sizeVcycles %>% group_by(GraphId) %>% slice_max(Cycles) %>% ungroup() #tira repetidos
View(sizeVcycles %>% group_by(GraphId) %>% filter(n() > 1))
sizeVcycles <- sizeVcycles[-c(30,41,68,77,78),]
sizeVcycles <- cbind(sizeVcycles, size = SCCs %>% group_by(GraphId) %>% summarise(Size = max(SCCSize)) %>% select(Size))

plot(sizeVcycles$Size, sizeVcycles$Cycles, title('ecoli'))
plot(sizeVcycles$Size, sizeVcycles$Cycles, title('ecoli'), xlim = c(10,15), ylim = c(0,100))
points(11, 32, col='red')

plot(sizeVcycles$Size, sizeVcycles$Length, title('ecoli'))
plot(sizeVcycles$Size, sizeVcycles$Length, title('ecoli'), xlim = c(10,15), ylim = c(0,8))
points(11, 4.03, col='red')

cycles_e_operon %>% filter(SCCId == 23) %>% summarise(mean(Length)) #4.03


# bacilus
sizeVcycles_bac <- set_bac$Cycles %>% group_by(GraphId, SCCId) %>% summarise(Cycles = n(), Length = mean(Length))
sizeVcycles_bac <- sizeVcycles_bac %>% group_by(GraphId) %>% slice_max(Cycles) %>% ungroup() #tira repetidos
View(sizeVcycles_bac %>% group_by(GraphId) %>% filter(n() > 1))
sizeVcycles_bac <- cbind(sizeVcycles_bac, set_bac$SCCs %>% group_by(GraphId) %>% summarise(Size = max(SCCSize)) %>% select(Size))

plot(sizeVcycles_bac$Size, sizeVcycles_bac$Cycles, title('Bac'))
plot(sizeVcycles_bac$Size, sizeVcycles_bac$Cycles, title('Bac'), xlim = c(5,20), ylim = c(0,170))
points(13, 45, col='red')

plot(sizeVcycles_bac$Size, sizeVcycles_bac$Length, title('Bac'))
plot(sizeVcycles_bac$Size, sizeVcycles_bac$Length, title('Bac'), xlim = c(5,20), ylim = c(0,10))
points(13, 5.07, col='red')

botcycles_bac2 %>% filter(SCCId == 3) %>% summarise(mean(Length)) #5.07

