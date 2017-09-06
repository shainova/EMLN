# NOTE!!! This file accompanies the following publication and can
# only be understood by reading the details in the manuscript and its
# SI. Please cite the original publication if using this code.
# 
# Pilosof S, Porter MA, Pascual M, Kefi S.
# The multilayer nature of ecological networks.
# Nature Ecology & Evolution (2017).

# Load libraries and functions -------------------------------------------------------------------------
library(bipartite)
library(igraph)
library(ggplot2)
library(reshape2)
library(stringr)
library(plyr)
library(RColorBrewer)


# Remove all objects, including functions
rm(list=ls())
# Remove all objects, except functions
rm(list = setdiff(ls(), lsf.str()))

df2matrix <- function(df,binary=F){
  rownames(df) <- df[,1]
  df <- df[,-1]
  df <- data.matrix(df)
  if (binary){df[df>0] <- 1}
  return(df)
}

ggplotToBrowser <- function(p, w=35, h=16.875) {
  ggsave(filename = tf_img <- tempfile(fileext = ".svg"), plot = p, width=w, height=h, units = 'cm')
  html <- sprintf('<html><body><img src="%s"></body></html>', paste0("file:///", tf_img))
  cat(html, file = tf_html <- tempfile(fileext = ".html"))
  if(Sys.info()[1]=="Linux"){
    options(browser="google-chrome")
    browseURL(tf_html)  
  } else {
    system(sprintf("open %s", tf_html))  
  }
}

theme_Publication <- function(base_size=14, base_family="helvetica", legend.position='bottom') {
  require(grid)
  require(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            #panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = legend.position,
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
} 

list2Array <- function(networkList){
  networkArray <- array(0, dim = list(nrow(networkList[[1]]),ncol(networkList[[1]]),length(networkList)),
                        dimnames = list(rownames(networkList[[1]]),colnames(networkList[[1]]),names(networkList)))
  for (l in 1:length(networkList)){
    networkArray[,,l] <- data.matrix(networkList[[l]])
  }
  return(networkArray)
}



#################################################################################################
# Analyses of robustness of main text
# The algorithms and rational for the removal scenarios are explained in detail in the paper
#################################################################################################
edgelist2Matrix <- function(elist){
  elist$lower.taxon <- as.character(elist$lower.taxon)
  elist$upper.taxon <- as.character(elist$upper.taxon)
  rows <- unique(elist$lower.taxon)
  cols <- unique(elist$upper.taxon)
  network <- matrix(0, nrow = length(rows), ncol = length(cols), dimnames = list(rows,cols))
  for (r in 1:nrow(elist)){
    lowerTaxon <- elist$lower.taxon[r]
    upperTaxon <- elist$upper.taxon[r]
    network[lowerTaxon,upperTaxon] <- elist$estimated.interaction.strength[r]
  }
  return(network)
}

binarize <- function(mat){
  m <- mat 
  m[which(m > 0)] <- 1
  return(m)
}

createSupraAdjMatrix <- function(layer1,layer2){
  s1 <- nrow(layer1);s2 <- ncol(layer1);s3 <- ncol(layer2)
  N=2*s1+s2+s3
  N_L1=s1+s2
  supraAdjMatrix <- matrix(0,N,N)
  # Intralayer edges layer 1
  supraAdjMatrix[1:s1,(s1+1):N_L1] <- layer1
  supraAdjMatrix[(s1+1):N_L1,1:s1] <- t(layer1)
  #Intralayer edges layer 2
  supraAdjMatrix[(N_L1+1):(N_L1+s1),(N_L1+s1+1):N] <- layer2
  supraAdjMatrix[(N_L1+s1+1):N,(N_L1+1):(N_L1+s1)] <- t(layer2)
  # Interlayer edges
  supraAdjMatrix[1:s1,(N_L1+1):(N_L1+s1)] <- diag(rep(1,s1))
  supraAdjMatrix[(N_L1+1):(N_L1+s1),1:s1] <- diag(rep(1,s1))
  colnames(supraAdjMatrix) <- rownames(supraAdjMatrix) <- make.unique(c(rownames(layer1),colnames(layer1),rownames(layer1),colnames(layer2)))
  if(!isSymmetric(supraAdjMatrix)){stop('Supra adjacency matrix is not symmetric')}
  return(supraAdjMatrix)
}

extractLayer <- function(network,s1,s2,s3){
  lyr1=network[1:s1,(s1+1):(s1+s2)]
  lyr2=network[(s1+s2+1):(s1+s2+s1),(s1+s2+s1+1):(2*s1+s2+s3)]
  return(list(lyr1,lyr2))
}

removeSpeciesMultilayer <- function(network,species){
  # Removing a species is akin to setting all its interactions to 0.
  network[which(str_detect(rownames(network),species)==T),] <- 0
  network[,which(str_detect(colnames(network),species)==T)] <- 0
  return(network)
}

extantSpecies <- function(network,s1,s2,s3,summary=F){ # Check how many species are extant in each layer
  # An extant species is a one for which the sum of interactions is different than 0.
  lyr1 <- extractLayer(network,s1,s2,s3)[[1]]
  lyr2 <- extractLayer(network,s1,s2,s3)[[2]]
  extantSet1_1 <- rowSums(lyr1)
  extantSet1_2 <- rowSums(lyr2)
  extantSet2 <- colSums(lyr1)
  extantSet3 <- colSums(lyr2)
  if (summary){
    extantSet1_1  <- sum(extantSet1_1>0)
    extantSet1_2 <- sum(extantSet1_2>0)
    extantSet2 <- sum( extantSet2 >0)
    extantSet3 <- sum( extantSet3 >0)
  }
  return(list(extantSet1_1=extantSet1_1,extantSet2=extantSet2,extantSet1_2=extantSet1_2,extantSet3=extantSet3))  
} 

extantSpeciesMonolayer <- function(network,summary=F){ # Check how many species are extant in each layer
  # An extant species is a one for which the sum of interactions is different than 0.
  extantSet1 <- rowSums(network)
  extantSet2 <- colSums(network)
  if (summary){
    extantSet1  <- sum(extantSet1>0)
    extantSet2 <- sum( extantSet2 >0)
  }
  return(list(extantSet1=extantSet1,extantSet2=extantSet2))  
} 

dataMain <- read.csv('Pocock2012_norwood.csv', header = T) # This data set can be downloaded from Dryad: http://dx.doi.org/10.5061/dryad.3s36r118
# These 3 species appear both as flower visitors and parasitoids. So I omitted them for simplicity
dataMain <- dataMain[dataMain$upper.taxon != 'Alysiinae sp',]
dataMain <- dataMain[dataMain$upper.taxon != 'Miscogaster maculata',]
dataMain <- dataMain[dataMain$upper.taxon != 'Bracon sp',]


#Create two bipartite matrices from the data
dataVisitors <- dataMain[dataMain$upper.guild=='flower visitor',]
dataParasitoids <- dataMain[dataMain$upper.guild=='leaf-miner parasitoid',]
visitorsLayer <- edgelist2Matrix(dataVisitors)
parasitoidsLayer <- edgelist2Matrix(dataParasitoids)

commonplants <- intersect(rownames(visitorsLayer), rownames(parasitoidsLayer))
visitorsLayer <- visitorsLayer[commonplants,]
parasitoidsLayer <- parasitoidsLayer[commonplants,]

visitorsLayer <- bipartite::empty(visitorsLayer)
parasitoidsLayer <- bipartite::empty(parasitoidsLayer)

# Work with binary layers because of the disparate scales of weights and the different methods of quantification
visitorsLayer <- binarize(visitorsLayer)
parasitoidsLayer <- binarize(parasitoidsLayer)

# Create the empirical multilayer network
supraAdjMatrix <- createSupraAdjMatrix(visitorsLayer,parasitoidsLayer)
s1 <- nrow(visitorsLayer);s2 <- ncol(visitorsLayer);s3 <- ncol(parasitoidsLayer)
N <- nrow(supraAdjMatrix)

multilayerNetwork <- supraAdjMatrix  
multilayerNetworkOrig <- supraAdjMatrix

## Scenario 1 - Monolayer (remove plants and follow parasitoids’ secondary extinction)
monolayerNetwork <- parasitoidsLayer
extinctionProcess <- matrix(unlist(extantSpeciesMonolayer(monolayerNetwork)),1,s1+s3)
colnames(extinctionProcess) <- unlist(dimnames(monolayerNetwork))
#Removal is by degree, starting with the least connected nodes.
extinctionOrder <- extantSpeciesMonolayer(monolayerNetwork)$extantSet1
extinctionOrder <- names(sort(extinctionOrder, decreasing = T))
for (x in 1:length(extinctionOrder)){
  print(x)
  #1. Remove a node from rows
  monolayerNetwork <- removeSpeciesMultilayer(monolayerNetwork,extinctionOrder[x]) 
  #2. Remove nodes from columns which are disconnected
  extinctNodes <- which(extantSpeciesMonolayer(monolayerNetwork,F)$extantSet2==0)
  if (length(extinctNodes)>0){
    for (E in names(extinctNodes)){monolayerNetwork <- removeSpeciesMultilayer(monolayerNetwork,E)} 
  }
  #3. Record
  extinctionProcess <- rbind(extinctionProcess,unlist(extantSpeciesMonolayer(monolayerNetwork)))
}
rownames(extinctionProcess) <- c('Initial',extinctionOrder)
rowsExtinction <- extinctionProcess[,1:s1]
rowsExtinction <- apply(rowsExtinction, 1, function(x) sum(x!=0))/s1
colsExtinction <- extinctionProcess[,(s1+1):ncol(extinctionProcess)]
colsExtinction <- apply(colsExtinction, 1, function(x) sum(x!=0))/s3
dScenario1 <- data.frame(order = rep(1:nrow(extinctionProcess),2),
                         extant = c(rowsExtinction,colsExtinction),
                         Extinction=c(rep('Plants primary extinction',nrow(extinctionProcess)),rep('Parasitoids secondary extinction',nrow(extinctionProcess))),
                         Guild=c(rep('Plants',nrow(extinctionProcess)),rep('Parasitoids',nrow(extinctionProcess))))
dScenario1 <- dScenario1[with(dScenario1, order(order,Extinction)), ]
dScenario1$extant <- as.numeric(dScenario1$extant)
dScenario1$order <- as.numeric(dScenario1$order)
dScenario1$order <- (dScenario1$order-1)/max(dScenario1$order)
# Plot it
plt=ggplot(dScenario1, aes(x=order, y=extant, group=Extinction, color=Extinction))+geom_point(size=4)+geom_line()+facet_wrap(~Extinction,scales="free_y")+
  labs(title='Removal in the Scenario1 by degree',x='Primary extinction', y='Number of species remaining')+ scale_colour_manual(values=c('orange','blue'))
# ggplotToBrowser(plt, w = 50,20)

## Scenario 1 - Monolayer RANDOM (remove plants and follow parasitoids’ secondary extinction)
extinctionProcessList=list()
for (run in 1:100){
  monolayerNetwork <- parasitoidsLayer
  extinctionProcess <- matrix(unlist(extantSpeciesMonolayer(monolayerNetwork)),1,s1+s3)
  colnames(extinctionProcess) <- unlist(dimnames(monolayerNetwork))
  #Removal is random
  extinctionOrder <- extantSpeciesMonolayer(monolayerNetwork)$extantSet1
  extinctionOrder <- sample(names(extinctionOrder), size = length(names(extinctionOrder)), replace = F)
  for (x in 1:length(extinctionOrder)){
    print(x)
    #1. Remove a node from rows
    monolayerNetwork <- removeSpeciesMultilayer(monolayerNetwork,extinctionOrder[x]) 
    #2. Remove nodes from columns which are disconnected
    extinctNodes <- which(extantSpeciesMonolayer(monolayerNetwork,F)$extantSet2==0)
    if (length(extinctNodes)>0){
      for (E in names(extinctNodes)){monolayerNetwork <- removeSpeciesMultilayer(monolayerNetwork,E)} 
    }
    #3. Record
    extinctionProcess <- rbind(extinctionProcess,unlist(extantSpeciesMonolayer(monolayerNetwork)))
  }
  rownames(extinctionProcess) <- c('Initial',extinctionOrder)
  rowsExtinction <- extinctionProcess[,1:s1]
  rowsExtinction <- apply(rowsExtinction, 1, function(x) sum(x!=0))/s1
  colsExtinction <- extinctionProcess[,(s1+1):ncol(extinctionProcess)]
  colsExtinction <- apply(colsExtinction, 1, function(x) sum(x!=0))/s3
  dScenario1_random <- data.frame(order = rep(1:nrow(extinctionProcess),2),
                           extant = c(rowsExtinction,colsExtinction),
                           Extinction=c(rep('Plants primary extinction',nrow(extinctionProcess)),rep('Parasitoids secondary extinction',nrow(extinctionProcess))),
                           Guild=c(rep('Plants',nrow(extinctionProcess)),rep('Parasitoids',nrow(extinctionProcess))))
  dScenario1_random <- dScenario1_random[with(dScenario1_random, order(order,Extinction)), ]
  dScenario1_random$extant <- as.numeric(dScenario1_random$extant)
  dScenario1_random$order <- as.numeric(dScenario1_random$order)
  dScenario1_random$order <- (dScenario1_random$order-1)/max(dScenario1_random$order)
  extinctionProcessList[[run]] <- dScenario1_random
}
dScenario1_random <- do.call(rbind.data.frame, extinctionProcessList)
library(dplyr)
dScenario1_random <- dScenario1_random %>% group_by(order,Extinction,Guild) %>% summarise(extant=mean(extant))
  # Plot it
plt=ggplot(dScenario1_random, aes(x=order, y=extant, group=Extinction, color=Extinction))+geom_point(size=4)+geom_line()+facet_wrap(~Extinction,scales="free_y")+
  labs(title='Removal in the Scenario1 by degree',x='Primary extinction', y='Number of species remaining')+ scale_colour_manual(values=c('orange','blue'))
# ggplotToBrowser(plt, w = 50,20)


## Scenario 2 -- Multilayer (remove flower visitors and follow parasitoids’ tertiary extinction)
multilayerNetwork <- multilayerNetworkOrig
extinctionProcess=matrix(unlist(extantSpecies(multilayerNetwork,s1,s2,s3,F)),1,length(colnames(multilayerNetwork)))  
colnames(extinctionProcess) <- colnames(multilayerNetwork)
extinctionOrder <- extantSpecies(multilayerNetwork,s1,s2,s3,F)$extantSet2 # Extinction occurs by removing nodes from Set2 in layer 1
extinctionOrder <- names(sort(extinctionOrder, decreasing = T))
for (x in 1:length(extinctionOrder)){
  print(x)
  #1. Remove a node from Set2 in layer 1
  multilayerNetwork <- removeSpeciesMultilayer(multilayerNetwork,extinctionOrder[x]) 
  #2. see if any nodes from Set1 in layer 1 have a degree of 0
  extinctSet1_1 <- which(extantSpecies(multilayerNetwork,s1,s2,s3,F)$extantSet1_1==0)
  if (length(extinctSet1_1)>0){
    #3. If so, remove these nodes from layer 2 as well
    for (E in names(extinctSet1_1)){multilayerNetwork <- removeSpeciesMultilayer(multilayerNetwork,E)}
  }
  #4. Remove nodes from Set3 in layer 2 which are disconnected
  extinctSet3 <- which(extantSpecies(multilayerNetwork,s1,s2,s3,F)$extantSet3==0)
  if (length(extinctSet3)>0){
    for (E in names(extinctSet3)){multilayerNetwork <- removeSpeciesMultilayer(multilayerNetwork,E)} 
  }
  #5. Record
  extinctionProcess <- rbind(extinctionProcess,unlist(extantSpecies(multilayerNetwork,s1,s2,s3,F)))
}
rownames(extinctionProcess) <- c('Initial',extinctionOrder)

# Plot
Set1Extinction <- extinctionProcess[,1:s1]
Set1Extinction <- apply(Set1Extinction, 1, function(x) sum(x!=0))/s1
Set3Extinction <- extinctionProcess[,(2*s1+s2+1):N]
Set3Extinction <- apply(Set3Extinction, 1, function(x) sum(x!=0))/s3
dScenario2 <- data.frame(order = rep(1:nrow(extinctionProcess),2),
                         extant = c(Set1Extinction,Set3Extinction),
                         Extinction=c(rep('Plants secondary extinction',nrow(extinctionProcess)),rep('Parasitoids tertiary extinction',nrow(extinctionProcess))),
                         Guild=c(rep('Plants',nrow(extinctionProcess)),rep('Parasitoids',nrow(extinctionProcess))))
dScenario2 <- dScenario2[with(dScenario2, order(order,Extinction)), ]
dScenario2$extant <- as.numeric(dScenario2$extant)
dScenario2$order <- as.numeric(dScenario2$order)
dScenario2$order <- dScenario2$order/max(dScenario2$order)
plt=ggplot(dScenario2, aes(x=order, y=extant, group=Extinction, color=Extinction))+geom_point(size=4)+geom_line()+facet_wrap(~Extinction,scales="free_y")+
  labs(title='Scenario 2--Flower visitors removed by degree',x='Primary extinction', y='Number of species remaining')+ scale_colour_manual(values=c('orange','blue'))
# ggplotToBrowser(plt, w = 50,20)


## Scenario 2 RANDOM -- Multilayer (remove flower visitors and follow parasitoids’ tertiary extinction)
extinctionProcessList=list()
for (r in 1:20){
  multilayerNetwork <- multilayerNetworkOrig
  extinctionProcess=matrix(unlist(extantSpecies(multilayerNetwork,s1,s2,s3,F)),1,length(colnames(multilayerNetwork)))  
  colnames(extinctionProcess) <- colnames(multilayerNetwork)
  extinctionOrder <- extantSpecies(multilayerNetwork,s1,s2,s3,F)$extantSet2 # Extinction occurs by removing nodes from Set2 in layer 1
  extinctionOrder <- sample(names(extinctionOrder), size = length(names(extinctionOrder)), replace = F)
  for (x in 1:length(extinctionOrder)){
    print(x)
    #1. Remove a node from Set2 in layer 1
    multilayerNetwork <- removeSpeciesMultilayer(multilayerNetwork,extinctionOrder[x]) 
    #2. see if any nodes from Set1 in layer 1 have a degree of 0
    extinctSet1_1 <- which(extantSpecies(multilayerNetwork,s1,s2,s3,F)$extantSet1_1==0)
    if (length(extinctSet1_1)>0){
      #3. If so, remove these nodes from layer 2 as well
      for (E in names(extinctSet1_1)){multilayerNetwork <- removeSpeciesMultilayer(multilayerNetwork,E)}
    }
    #4. Remove nodes from Set3 in layer 2 which are disconnected
    extinctSet3 <- which(extantSpecies(multilayerNetwork,s1,s2,s3,F)$extantSet3==0)
    if (length(extinctSet3)>0){
      for (E in names(extinctSet3)){multilayerNetwork <- removeSpeciesMultilayer(multilayerNetwork,E)} 
    }
    #5. Record
    extinctionProcess <- rbind(extinctionProcess,unlist(extantSpecies(multilayerNetwork,s1,s2,s3,F)))
  }
  rownames(extinctionProcess) <- c('Initial',extinctionOrder)
  
  # Plot
  Set1Extinction <- extinctionProcess[,1:s1]
  Set1Extinction <- apply(Set1Extinction, 1, function(x) sum(x!=0))/s1
  Set3Extinction <- extinctionProcess[,(2*s1+s2+1):N]
  Set3Extinction <- apply(Set3Extinction, 1, function(x) sum(x!=0))/s3
  dScenario2_random <- data.frame(order = rep(1:nrow(extinctionProcess),2),
                           extant = c(Set1Extinction,Set3Extinction),
                           Extinction=c(rep('Plants secondary extinction',nrow(extinctionProcess)),rep('Parasitoids tertiary extinction',nrow(extinctionProcess))),
                           Guild=c(rep('Plants',nrow(extinctionProcess)),rep('Parasitoids',nrow(extinctionProcess))))
  dScenario2_random <- dScenario2_random[with(dScenario2_random, order(order,Extinction)), ]
  dScenario2_random$extant <- as.numeric(dScenario2_random$extant)
  dScenario2_random$order <- as.numeric(dScenario2_random$order)
  dScenario2_random$order <- dScenario2_random$order/max(dScenario2_random$order)
  
  extinctionProcessList[[run]] <- dScenario2_random
}

dScenario2_random <- do.call(rbind.data.frame, extinctionProcessList)
dScenario2_random <- dScenario2_random %>% group_by(order,Extinction,Guild) %>% summarise(extant=mean(extant))


plt=ggplot(dScenario2_random, aes(x=order, y=extant, group=Extinction, color=Extinction))+geom_point(size=4)+geom_line()+facet_wrap(~Extinction,scales="free_y")+
  labs(title='Scenario 2--Flower visitors removed by degree',x='Primary extinction', y='Number of species remaining')+ scale_colour_manual(values=c('orange','blue'))
# ggplotToBrowser(plt, w = 50,20)


dScenario1$Scenario <- 'Monolayer'
dScenario1_random$Scenario <- 'Monolayer'
dScenario2$Scenario <- 'Multilayer'
dScenario2_random$Scenario <- 'Multilayer'
dScenario1$Type <- 'Degree'
dScenario1_random$Type <- 'Random'
dScenario2$Type <- 'Degree'
dScenario2_random$Type <- 'Random'


## Plot final plots
d_Plot <- rbind(dScenario1,as.data.frame(dScenario1_random),dScenario2,as.data.frame(dScenario2_random))

plt=ggplot(d_Plot, aes(x=order, y=extant, color=Scenario))+geom_point(size=3)+geom_line(aes(linetype=Type))+facet_wrap(~Guild)+
    labs(title='Removal by degree',x='Proportion of species removed', y='Proportion of leafminer parasitoids remaining')+
  scale_colour_manual(values=c('dark orange','#914200','mediumpurple1','forestgreen'))+theme_Publication()
# ggplotToBrowser(plt, w = 40,20)


# Plot parasitoid extinction as a function of plant removal, instead of species removal
d_Plot <- d_Plot[with(d_Plot, order(Scenario, order)),]

x1 <- subset(d_Plot, Guild=='Plants')
x2 <- subset(d_Plot, Guild=='Parasitoids')
x2$order <- 1-x1$extant # order is the number of plants removed
png('~/Dropbox/RESEARCH/Meetings and Workshops/2017 ESA/robustness_plant_removal.png', width = 2000, height = 780, units='px')
ggplot(x2, aes(x=order, y=extant, color=Scenario))+
  geom_line(aes(linetype=Type), size=3)+
  geom_point(size=7)+
  facet_wrap(~Type)+
  scale_colour_manual(values=c('dark orange','forestgreen'))+
  # labs(x='Proportion of plants removed', y='Proportion of leafminer parasitoids remaining')+
  labs(x='', y='')+
  theme_Publication(base_size = 60, legend.position = 'none')+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())
dev.off()
