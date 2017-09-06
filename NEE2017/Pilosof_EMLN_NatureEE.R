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

theme_Publication <- function(base_size=14, base_family="helvetica") {
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
            legend.position = "bottom",
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
# 1. prepare the data
#################################################################################################

# Each row in the data set is a single host specimen (individual). Matrix
# entries are the number of parasites recovered from each host specimen
dat <- read.csv('SupplementaryData1.csv')

attach(dat)
hostAbundYear <- as.matrix(table(Host, YearCollected)) # abundance of hosts in different years
parasiteAbundanceYear <- aggregate(.~dat$YearCollected,data = dat[,3:ncol(dat)], FUN = sum) # abundance of parasites in each year (across hosts)
names(parasiteAbundanceYear)[1] <- 'YearCollected'
parasiteAbundanceYear <- df2matrix(parasiteAbundanceYear)


#### Create and write network layers
data.list <- list()
years <- 1982:1987
for (y in years){
  idx <- which(years==y)
  d <- dat[dat$YearCollected==y,]
  d <- aggregate(.~d$Host, data=d[,2:ncol(d)], sum) # The total number of parasites found on a given host
  d <- d[-2]
  d <- df2matrix(d)
  d <- sweep(d, 1, hostAbundYear[rownames(hostAbundYear)%in%rownames(d), idx], '/') # Average parasite abundance per host
  missingHosts <- setdiff(rownames(hostAbundYear),rownames(d)) # All hosts have to appear in all matrices even if they were not present
  d <- rbind(d,matrix(0,length(missingHosts),ncol(d),dimnames =  list(missingHosts,colnames(d)))) # Add missing host species so all matrices will have the same size
  d <- d[sort(rownames(d)),] # sort by host so all matrices will have the same order
  data.list[[idx]] <- d
  write.table(d, paste('host_parasite_abundance_weighted_layer_',idx,'.csv',sep=''), row.names = F, col.names = F,sep=',')
}
names(data.list) <- years
sapply(data.list,dim)

### Create aggregated networks
#1. Aggregated network is the sum of all interactions across time.
aggregated_network_sum <- Reduce('+', data.list)
dim(aggregated_network_sum)
any(colSums(aggregated_network_sum)==0) #should be False
any(rowSums(aggregated_network_sum)==0) #should be False
write.table(aggregated_network_sum, 'aggregated_network_sum.csv', row.names = F, col.names = F,sep=',')
#2. Aggregated network is the mean of all interactions across time.
aggregated_network_mean <- apply(simplify2array(data.list), 1:2, mean)
dim(aggregated_network_mean)
any(colSums(aggregated_network_mean)==0) #should be False
any(rowSums(aggregated_network_mean)==0) #should be False
write.table(aggregated_network_mean, 'aggregated_network_mean.csv', row.names = F, col.names = F,sep=',')



### Use relative changes in species abundance as interlayer edge weights
mat <- hostAbundYear
tot <- ncol(mat)
interlayerEdgesHost <- matrix(0,nrow(mat),tot-1,dimnames=list(rownames(mat),colnames(mat)[-1]))
for (x in 1:(tot-1)){
  interlayerEdgesHost[,x] <- mat[,x+1]/mat[,x]
}
mat <- t(parasiteAbundanceYear)
tot <- ncol(mat)
interlayerEdgesParas <- matrix(0,nrow(mat),tot-1,dimnames=list(rownames(mat),colnames(mat)[-1]))
for (x in 1:(tot-1)){
  interlayerEdgesParas[,x] <- mat[,x+1]/mat[,x]
}
interlayerEdges <- rbind(interlayerEdgesHost,interlayerEdgesParas)
interlayerEdges[is.nan(interlayerEdges)] <- 0
interlayerEdges[is.infinite(interlayerEdges)] <- 0

ile <- c()
for (i in 1:5){
  ile <- c(ile, interlayerEdges[,i])
}

N=nrow(data.list[[1]])+ncol(data.list[[1]])
L=length(data.list)
mat=matrix(0,N*L,N*L)
delta <- row(mat) - col(mat)
mat[delta == N] <- ile
mat[delta == -N] <- ile
isSymmetric(mat)
# This file will be used in the analysis of modular strucutre (done in Matlab).
# It is actually a supra-adjacency matrix where all cells have a value of zero
# besides the off-block diagonals, which contain the interlayer edges.
write.table(mat,'interlayer_relative_abundance_matrix.csv', row.names = F, col.names = F,sep=',')


### Reshuffle interactions within each layer -- this is for the first null model in the analysis of modularity

# The idea here is to reshuffle the networks when
# they contain only the species that appeared in a given year and then paste
# them back to a matrix with all the speices/interactions
dir.create('reshuffled_networks_abundance')
realizations=1000
baseMatrix <- data.list[[1]]
baseMatrix[baseMatrix>0]=0
for (l in 1:6){
  y <- years[l]
  d <- dat[dat$YearCollected==y,]
  d <- aggregate(.~d$Host, data=d[,2:ncol(d)], sum) # The total number of parasites found on a given host
  d <- d[-2]
  d <- df2matrix(d)
  # see commsim help for details on null model. r0_both conserves the number of parasite species that infect a given host and their mean abundance
  nullnets <- vegan::nullmodel(bipartite::empty(d), method = 'r0_ind')  
  nullnets <- simulate(nullnets, nsim = realizations)
  dim(nullnets)
  for (i in 1:realizations){
    print(i)
    x=baseMatrix
    x[rownames(nullnets[,,i]),colnames(nullnets[,,i])] <- sweep(nullnets[,,i], 1, hostAbundYear[rownames(hostAbundYear)%in%rownames(nullnets[,,i]), l], '/')
    write.table(x, paste('reshuffled_networks_abundance/network_',i,'_layer_',l,'.csv',sep=''), row.names = F,col.names =F,sep=',')
  }
}

#################################################################################################
# 2. Analyze modular structure. This is done in Matlab
#################################################################################################


#################################################################################################
# 3. Post analysis of the modular structure
#################################################################################################

## Give each species an index
hosts <- unique(as.vector(sapply(data.list, rownames)))
parasites <- unique(as.vector(sapply(data.list, colnames)))
species <- c(hosts,parasites)
any(duplicated(species))
nodes.df <- data.frame(nodeID=1:length(species),nodeLabel=species)


output_folder='output'
# Analyze the values of the modularity function
Q.multilayerlObs <- unname(data.matrix(read.table(paste(output_folder,'/Q_obs.csv',sep=''),he=F,sep=',')))
Q.multilayerlNull1 <- unname(data.matrix(read.table(paste(output_folder,'/Q_null1.csv',sep=''),he=F,sep=',')))
Q.multilayerNull2_hosts <- unname(data.matrix(read.table(paste(output_folder,'/Q_null2_hosts.csv',sep=''),he=F,sep=',')))
Q.multilayerNull2_parasites <- unname(data.matrix(read.table(paste(output_folder,'/Q_null2_parasites.csv',sep=''),he=F,sep=',')))
Q.multilayerNull3 <- unname(data.matrix(read.table(paste(output_folder,'/Q_null3.csv',sep=''),he=F,sep=',')))

mean(Q.multilayerlObs)
mean(Q.multilayerlNull1)
mean(Q.multilayerNull2_hosts)
mean(Q.multilayerNull2_parasites)
mean(Q.multilayerNull3)

sum(Q.multilayerlNull1>mean(Q.multilayerlObs))/nrow(Q.multilayerlNull1)
sum(Q.multilayerNull2_hosts>mean(Q.multilayerlObs))/nrow(Q.multilayerNull2_hosts)
sum(Q.multilayerNull2_parasites>mean(Q.multilayerlObs))/nrow(Q.multilayerNull2_parasites)
sum(Q.multilayerNull3>mean(Q.multilayerlObs))/nrow(Q.multilayerNull3)

# Function to calculate the average number of modules and also the average propotion of species that switch modules.
# Input is the file that contains the module affiliations (memberships) that was produced in Matlab.
analyzeModuleAffiliation <- function(membership){
  runs=ncol(membership)
  membershipRuns <- matrix(0,length(species),runs,dimnames = list(species,1:runs)) # to store the module affiliation of species across runs
  modulesRuns <- rep(0,runs) #to store the number of modules across runs
  speciesFlexibility <- matrix(0,2,runs,dimnames = list(c('hosts','parasites'),1:runs))
  for (r in 1:runs){
    membershipRun <- matrix(membership[,r], nrow=78,ncol=6) # This reshapes the vector to a matrix of dimensions nodes X layers
    colnames(membershipRun) <- c('1982','1983','1984','1985','1986','1987')
    membershipRun <- as.data.frame(membershipRun)
    membershipRun <- cbind(species,c(rep('host',length(hosts)),rep('paras',length(parasites))),membershipRun)
    names(membershipRun)[2] <- 'type'
    # remove instances where a species was not present in a given year. This is
    # necessary because the original matrices contain all the species, regardless
    # if they were present. In the analysis of modular structure, these are placed
    # in their own modules which obviously do not really exist. This inflates the
    # number of real modules in any given year.
    speciesPresence <- rbind(hostAbundYear,t(parasiteAbundanceYear))
    speciesPresence[speciesPresence>0] <- 1
    if(!all(rownames(speciesPresence)==membershipRun$species)){stop('error')}
    membershipRun[,3:8] <- membershipRun[,3:8]*speciesPresence # make zero whenever a species is not in a layer
    # To how many communitites were species classified across years?
    membershipRuns[,r] <- apply(membershipRun[,3:8], 1, function(x) length(unique(x[x!=0])))
    # Number of hosts that were assigned to more than 1 module
    speciesFlexibility['hosts',r] <- sum(membershipRuns[1:length(hosts),r]>1)
    # Number of parasites that were assigned to more than 1 module
    speciesFlexibility['parasites',r] <- sum(membershipRuns[(length(hosts)+1):length(species),r]>1)
    # Number of modules
    modulesRuns[r] <- length(unique(unlist(membershipRun[,3:8])))-1 #reduce one to ommit the zero values which indicate species absence and not a module
  }
  return(list(modulesRuns=modulesRuns,speciesFlexibility=speciesFlexibility))
}

moduleAnalysisObs <- analyzeModuleAffiliation(read.csv('output/S_obs.csv', he=F))
moduleAnalysisNull1 <- analyzeModuleAffiliation(read.csv('output/S_null1.csv', he=F))
moduleAnalysisNull2_hosts <- analyzeModuleAffiliation(read.csv('output/S_null2_hosts.csv', he=F))
moduleAnalysisNull2_parasites <- analyzeModuleAffiliation(read.csv('output/S_null2_parasites.csv', he=F))
moduleAnalysisNull3 <- analyzeModuleAffiliation(read.csv('output/S_null3.csv', he=F))

# Average number of modules
mean(moduleAnalysisObs$modulesRuns)
mean(moduleAnalysisNull1$modulesRuns)
mean(moduleAnalysisNull2_hosts$modulesRuns)
mean(moduleAnalysisNull2_parasites$modulesRuns)
mean(moduleAnalysisNull3$modulesRuns)

# t-test to compare observed number of modules to rehsuffled networks
t.test(moduleAnalysisNull1$modulesRuns-mean(moduleAnalysisObs$modulesRuns),mu = 0)
t.test(moduleAnalysisNull2_hosts$modulesRuns-mean(moduleAnalysisObs$modulesRuns),mu = 0)
t.test(moduleAnalysisNull2_parasites$modulesRuns-mean(moduleAnalysisObs$modulesRuns),mu = 0)
t.test(moduleAnalysisNull3$modulesRuns-mean(moduleAnalysisObs$modulesRuns),mu = 0)

# Average host flexibility divided by the total number of hosts
mean(moduleAnalysisObs$speciesFlexibility[1,]) / 22
mean(moduleAnalysisNull1$speciesFlexibility[1,]) / 22
mean(moduleAnalysisNull2_hosts$speciesFlexibility[1,]) / 22
mean(moduleAnalysisNull2_parasites$speciesFlexibility[1,]) / 22
mean(moduleAnalysisNull3$speciesFlexibility[1,]) / 22

# Average parasite flexibility divided by the total number of parasites
mean(moduleAnalysisObs$speciesFlexibility[2,]) / 56
mean(moduleAnalysisNull1$speciesFlexibility[2,]) / 56
mean(moduleAnalysisNull2_hosts$speciesFlexibility[2,]) / 56
mean(moduleAnalysisNull2_parasites$speciesFlexibility[2,]) / 56
mean(moduleAnalysisNull3$speciesFlexibility[2,]) / 56

# t-test to compare observed host flexibility to that of rehsuffled networks
t.test(moduleAnalysisNull1$speciesFlexibility[1,]-rowMeans(moduleAnalysisObs$speciesFlexibility)[1],mu = 0)
t.test(moduleAnalysisNull2_hosts$speciesFlexibility[1,]-rowMeans(moduleAnalysisObs$speciesFlexibility)[1],mu = 0)
t.test(moduleAnalysisNull2_parasites$speciesFlexibility[1,]-rowMeans(moduleAnalysisObs$speciesFlexibility)[1],mu = 0)
t.test(moduleAnalysisNull3$speciesFlexibility[1,]-rowMeans(moduleAnalysisObs$speciesFlexibility)[1],mu = 0)
# t-test to compare observed parasite flexibility to that of rehsuffled networks
t.test(moduleAnalysisNull1$speciesFlexibility[2,]-rowMeans(moduleAnalysisObs$speciesFlexibility)[2],mu = 0)
t.test(moduleAnalysisNull2_hosts$speciesFlexibility[2,]-rowMeans(moduleAnalysisObs$speciesFlexibility)[2],mu = 0)
t.test(moduleAnalysisNull2_parasites$speciesFlexibility[2,]-rowMeans(moduleAnalysisObs$speciesFlexibility)[2],mu = 0)
t.test(moduleAnalysisNull3$speciesFlexibility[2,]-rowMeans(moduleAnalysisObs$speciesFlexibility)[2],mu = 0)


#### Plot the modularity membership
# To plot this first one should find out which is the run with the maximum Q_B
maxrun <- which(Q.multilayerlObs==max(Q.multilayerlObs))
# This is now the same code as in the for loop of the analyzeModuleAffiliation function
membership <- read.csv('output/S_obs.csv', he=F)
membershipRun <- matrix(membership[,maxrun], nrow=78,ncol=6) # This reshapes the vector to a matrix of dimensions nodes X layers
colnames(membershipRun) <- c('1982','1983','1984','1985','1986','1987')
membershipRun <- as.data.frame(membershipRun)
membershipRun <- cbind(species,c(rep('host',length(hosts)),rep('paras',length(parasites))),membershipRun)
names(membershipRun)[2] <- 'type'
speciesPresence <- rbind(hostAbundYear,t(parasiteAbundanceYear))
speciesPresence[speciesPresence>0] <- 1
if(!all(rownames(speciesPresence)==membershipRun$species)){stop('error')}
membershipRun[,3:8] <- membershipRun[,3:8]*speciesPresence # make zero whenever a species is not in a layer

# Plot module affiliation across years
membership_plot <- melt(membershipRun, id.vars = c('species','type'))
names(membership_plot) <- c('species','type','year','module')
membership_plot <- membership_plot[membership_plot$module!=0,]
unique(membership_plot$module)
moduleNumbers <- sort(unique(membership_plot$module)) # Make the module numbers consecutive
membership_plot$module <- match(membership_plot$module,moduleNumbers)
membership_plot$module <- factor(membership_plot$module)
membership_plot$species <- factor(membership_plot$species, levels=rev(levels(membership_plot$species)))

membership_plot$species <- str_replace(membership_plot$species,'_',' ')
p=ggplot(data = membership_plot, aes(x=year,y=species, color=module))+geom_point(size=3, shape=15)+
  facet_grid(type~.,space='free', scales='free')+theme_bw()+  scale_color_brewer(palette="Set1")
ggplotToBrowser(p, 20,40)

# svg('Fig_temporal_modualrity_assignment.svg', width = 15, height = 30)
# ggplot(data = membership_plot, aes(x=year,y=species, color=module))+geom_point(size=8, shape=15)+
#   facet_grid(type~.,space='free', scales='free')+theme_bw()+  scale_color_brewer(palette="Set1")
# dev.off()

# Now plot the module size across years
new <- merge(membership_plot, as.data.frame(with(membership_plot, table(module, year))), by = c('module','year')) 
p=ggplot(data = new, aes(x=year,y=module))+geom_point(aes(size=Freq),color='blue')+theme_bw()+ scale_size(range=c(3,15))
ggplotToBrowser(p)

# svg('Fig_temporal_modualrity_module_size.svg', width = 7, height = 7)
# ggplot(data = new, aes(x=year,y=module))+geom_point(aes(size=Freq))+theme_bw()+
#   scale_size_continuous(breaks=c(0,5,10,15,20,25,30),range=c(3,15))
# dev.off()


# Plot 6 layers and color STATE NODES by module
# Make the node color the same as in the membership_plot figure
# Track the species Myodes glareolus to see how it changes modules with time
moduleColordDF <- data.frame(module=unique(membership_plot$module), color=brewer.pal(6,"Set1"))
membership_plot$modColor <- moduleColordDF$color[match(membership_plot$module,moduleColordDF$module)]
shape <- c("circle", "square") # TO distinguish hosts and parasites
guildColor <- c("blue", "red") # TO distinguish hosts and parasites
for (l in 1:6){
  m <- data.list[[l]]
  m <- bipartite::empty(m)
  g <- graph.incidence(m)
  nodes <- membership_plot[membership_plot$year==c(1982:1987)[l],]
  V(g)$name <- str_replace(V(g)$name,'_',' ')
  nodes <- nodes[nodes$species %in% V(g)$name,]
  print(all(V(g)$name==nodes$species))
  V(g)$modColor <- as.character(nodes$modColor)
  displayNames <- nodes$species
  displayNames[displayNames!='Myodes glareolus'] <- NA
  # This plots node colors by their guild
  plot(g, layout=layout.gem, vertex.label=NA, vertex.size=7,vertex.color = guildColor[as.numeric(V(g)$type)+1], vertex.shape = shape[as.numeric(V(g)$type)+1])
  svg(paste('layer_',l,'.svg',sep=''),width = 6,height = 6)
  plot(g, layout=layout.gem, vertex.label=NA, vertex.size=7,vertex.color = guildColor[as.numeric(V(g)$type)+1], vertex.shape = shape[as.numeric(V(g)$type)+1])
  dev.off()
  # This plots node colors by their module
  plot(g, layout=layout.gem, vertex.label=displayNames, vertex.size=7,vertex.color = V(g)$modColor, vertex.shape = shape[as.numeric(V(g)$type)+1])
  svg(paste('layer_',l,'_modules.svg',sep=''),width = 6,height = 6)
  plot(g, layout=layout.gem, vertex.label=displayNames, vertex.size=7,vertex.color = V(g)$modColor, vertex.shape = shape[as.numeric(V(g)$type)+1])
  dev.off()
}

#################################################################################################
# A comparison of module affiliations between multilayer and monolayer networks
#################################################################################################

MI <- function (N) {
  S <- sum(N)
  CA = dim(N)[1]; CB = dim(N)[2]
  
  # Danon et al.'s 2005 method from http://arxiv.org/pdf/cond-mat/0505245.pdf
  Iup=0
  for (i in 1:CA){
    for (j in 1:CB){
      if (N[i,j] != 0){
        Ni.=sum(N[i,])
        N.j=sum(N[,j])
        Iup = Iup+N[i,j]*log((N[i,j]*S)/(Ni.*N.j))
      }
    }
  }
  Idown1=0;Idown2=0
  for (i in 1:CA){
    Ni.=sum(N[i,])
    Idown1=Idown1+Ni.*log(Ni./S)
  }
  for (j in 1:CB){
    N.j=sum(N[,j])
    Idown2=Idown2+N.j*log(N.j/S)
  }
  I=-2*Iup/(Idown1+Idown2)
  
  return(I)
}


# MULTILAYER: Get the module affiliations of each species-layer tuple
maxrun <- which(Q.multilayerlObs==max(Q.multilayerlObs))
# This is now the same code as in the for loop of the analyzeModuleAffiliation function
membership <- read.csv('output/S_obs.csv', he=F)
membershipRun <- matrix(membership[,maxrun], nrow=78,ncol=6) # This reshapes the vector to a matrix of dimensions nodes X layers
colnames(membershipRun) <- c('1982','1983','1984','1985','1986','1987')
membershipRun <- as.data.frame(membershipRun)
membershipRun <- cbind(species,c(rep('host',length(hosts)),rep('paras',length(parasites))),membershipRun)
names(membershipRun)[2] <- 'type'
speciesPresence <- rbind(hostAbundYear,t(parasiteAbundanceYear))
speciesPresence[speciesPresence>0] <- 1
if(!all(rownames(speciesPresence)==membershipRun$species)){stop('error')}
membershipRun[,3:8] <- membershipRun[,3:8]*speciesPresence # make zero whenever a species is not in a layer

# MONOLAYER: get the module affiliation of each species and compare to the multilayer using NMI
MNI_aggregated_sum <- MNI_aggregated_mean <- NMI_layers <- c()
for (l in 1:6){
  cat(l);cat('\t')
  # time slice in the multilayer:
  multilayer_affiliation <- membershipRun[,2+l]
  multilayer_affiliation <- multilayer_affiliation[multilayer_affiliation!=0]
  # Module affiliation in the isolated network's layer:
  Q.layer <- unname(data.matrix(read.table(paste(output_folder,'/Q_obs_zero_layer_',l,'.csv',sep=''),he=F,sep=',')))
  maxrun <- which.max(Q.layer)
  monolayer_affiliation <- unname(data.matrix(read.table(paste(output_folder,'/S_obs_zero_layer_',l,'.csv',sep=''),he=F,sep=',')))
  monolayer_affiliation <- monolayer_affiliation[,maxrun]
  #length(monolayer_affiliation)
  #length(speciesPresence[,l])
  monolayer_affiliation <- monolayer_affiliation*speciesPresence[,l]
  monolayer_affiliation <- monolayer_affiliation[monolayer_affiliation!=0]
  cat(ifelse(length(monolayer_affiliation)!=length(multilayer_affiliation),'Warning: lengths of species lists do not match',""));cat('\t')
  cat(length(monolayer_affiliation));cat('\t')
  NMI_layers <- c(NMI_layers, MI(xtabs(~monolayer_affiliation+multilayer_affiliation)))
  cat(NMI_layers[length(NMI_layers)]);cat('\t')
  # Module affiliation in the aggregated networks
  Q.agg_sum <- unname(data.matrix(read.table(paste(output_folder,'/Q_obs_aggregated_sum.csv',sep=''),he=F,sep=',')))
  maxrun <- which.max(Q.agg_sum)
  aggregated_sum_affiliation <- unname(data.matrix(read.table(paste(output_folder,'/S_obs_aggregated_sum.csv',sep=''),he=F,sep=',')))
  aggregated_sum_affiliation <- aggregated_sum_affiliation[,maxrun]
  aggregated_sum_affiliation <- aggregated_sum_affiliation[which(membershipRun[,2+l]!=0)] # the aggregated networks contains all the species so only take thse species that appear in the layer
  cat(ifelse(length(aggregated_sum_affiliation)!=length(multilayer_affiliation),'Warning: lengths of species lists do not match',""));cat('\t')
  MNI_aggregated_sum <- c(MNI_aggregated_sum, MI(xtabs(~aggregated_sum_affiliation+multilayer_affiliation)))
  cat(MNI_aggregated_sum[length(MNI_aggregated_sum)]);cat('\t')
  Q.agg_mean <- unname(data.matrix(read.table(paste(output_folder,'/Q_obs_aggregated_mean.csv',sep=''),he=F,sep=',')))
  maxrun <- which.max(Q.agg_mean)
  aggregated_mean_affiliation <- unname(data.matrix(read.table(paste(output_folder,'/S_obs_aggregated_mean.csv',sep=''),he=F,sep=',')))
  aggregated_mean_affiliation <- aggregated_mean_affiliation[,maxrun]
  aggregated_mean_affiliation <- aggregated_mean_affiliation[which(membershipRun[,2+l]!=0)] # the aggregated networks contains all the species so only take thse species that appear in the layer
  cat(ifelse(length(aggregated_mean_affiliation)!=length(multilayer_affiliation),'Warning: lengths of species lists do not match',""));cat('\t')
  MNI_aggregated_mean <- c(MNI_aggregated_mean, MI(xtabs(~aggregated_mean_affiliation+multilayer_affiliation)))
  cat(MNI_aggregated_mean[length(MNI_aggregated_mean)]);cat('\n')
}
paste(mean(NMI_layers),sd(NMI_layers))
paste(mean(MNI_aggregated_sum),sd(MNI_aggregated_sum))
paste(mean(MNI_aggregated_mean),sd(MNI_aggregated_mean))
range(MNI_aggregated_mean)

#################################################################################################
# Reducibility
#################################################################################################
library(d3heatmap)
library(rCharts)
library(gplots)

# Create an edge list in the general format: node layer node layer weight
# Only intralayer edges are necessary
edges <- data.frame(node_i=as.numeric(), layer_i=as.numeric(), node_j=as.numeric(), layer_j=as.numeric(), w=as.numeric())
for (l in 1:6){
  m = data.list[[l]]
  NZ <- which(m!=0, arr.ind = T) # non-zero elements
  COO <- data.frame(i=rep(NA,nrow(NZ)),j=rep(NA,nrow(NZ)),w=rep(NA,nrow(NZ)))
  COO$i <- rownames(m)[NZ[,1]]
  COO$j <- colnames(m)[NZ[,2]]
  for (n in 1:nrow(COO)){
    COO[n,'w'] <- m[as.character(COO[n,'i']),as.character(COO[n,'j'])]
  }
  COO[,1] <- nodes.df$nodeID[match(COO[,1], nodes.df$nodeLabel)]
  COO[,2] <- nodes.df$nodeID[match(COO[,2], nodes.df$nodeLabel)]
  for (z in 1:nrow(COO)){
    row=nrow(edges)+1
    edges[row,1] <- COO[z,1]
    edges[row,2] <- l
    edges[row,3] <- COO[z,2]
    edges[row,4] <- l
    edges[row,5] <- COO[z,3]
  }
}

#write edgelist for reducibility analysis
currfolder <- getwd()
setwd(paste(currfolder,'/Reducibility',sep=''))
write.table(edges,'host_parasite.edges', sep = ' ', row.names = F, col.names = F, quote = F)
# prepare the configuration file of the reducibilty analysis
f <- 'muxOctaveConfig.m'
projectname <- 'host_parasite_temporal'
write(paste("AnalysisName = \"",projectname,"\";",sep=""), f, append=F)
write('isExtendedEdgesList = 1;', f, append=T)
write('isExtendedEdgesList = 1;', f, append=T)
write(paste0("MultiLayerEdgesListFile = \"",normalizePath('host_parasite.edges',winslash = "/"),"\";"),file=f,append=T)
write(paste("Flags = \"",'DW',"\";",sep=""), f, append=T)
write(paste("MaxNodes = ",nrow(nodes.df),";",sep=""), f, append=T)
write(paste("FirstNodeLabel = ",min(nodes.df$nodeID),";",sep=""), f, append=T)
write("MultisliceType = \"ordered\";", f, append=T)
# Run reducibility
system('octave muxMultisliceReducibility.m')

# Look at the results
resultFile <- paste(projectname,'_reducibility_jsd.txt',sep='')
distanceMatrix <- matrix(scan(resultFile, n = length(data.list)^2), ncol=length(data.list), nrow=length(data.list), byrow = TRUE)
colnames(distanceMatrix) <- names(data.list)
rownames(distanceMatrix) <- names(data.list)
distanceMatrix.df <- as.data.frame(distanceMatrix)


# For plotting:
# svg('reducibility_heatmap.svg')
# my_palette <- brewer.pal(100,"Blues")
# heatmap.2(x=data.matrix(distanceMatrix.df), hclustfun = function(x) hclust(x,method='ward.D2'),
#                 density.info="none",  # turns off density plot inside color legend
#                 trace="none",         # turns off trace lines inside the heat map
#                 margins =c(12,9),     # widens margins around plot
#                 col=my_palette,       # use on color palette defined earlier
#                 dendrogram="col"     # only draw a col dendrogram
#                 )
# dev.off()

d3heatmap(
  distanceMatrix.df,
  color = 'Blues',
  labRow=names(data.list),
  labCol=names(data.list),
  cexRow=1,
  cexCol=1,
  hclustfun=function(x) hclust(x,method='ward.D2'),
  symm=F,
  dendrogram="both"
)

# svg('reducibility_dendrogram.svg')
plot(hclust(as.dist(distanceMatrix),
                  method='ward.D2'),
           col = "#1F77B4", col.main = "#1F77B4", col.lab = "#E08400", 
           col.axis = "#E08400", lwd = 2, 
           labels=names(data.list),
           cex=1,
           main="Reducibility Dendrogram",
           sub="", 
           xlab="")
# dev.off()

resultFile <- paste(projectname,'_reducibility_quality.txt',sep='')
d <- read.table(resultFile, header=T, sep=" ")[,1:2]
colnames(d) <- c("Step", "Q")
linechart <- nPlot(Q ~ Step, data = d, type = 'lineChart')
linechart$addParams(width = 600, height = 400, title="Quality function") 
linechart$xAxis(axisLabel="Step")
linechart$yAxis(axisLabel="Q")
linechart$chart(forceY=c(floor(min(d$Q)),floor(max(d$Q))+1), 
                forceX=c(floor(min(d$Step)),floor(max(d$Step))+1))
linechart
ggplot(d, aes(x=step, y=Q))+geom_point()
svg('reducibility_q.svg')
plot(d$Q~d$Step, type='b', xlab='Step', ylab='q',pch=19, cex=2)
dev.off()
if(file.exists("hclust_merge.txt")) file.remove("hclust_merge.txt")
if(file.exists("jsd_distance_matrix.txt")) file.remove("jsd_distance_matrix.txt")
if(file.exists(resultFile)) file.remove(resultFile)
resultFile <- paste(projectname,'_reducibility_jsd.txt',sep='')
if(file.exists(resultFile)) file.remove(resultFile)
setwd(currfolder)

#################################################################################################
# Analyses of modularity for the SI
#################################################################################################

## Plot distribution of interlayer edges
d <- data.frame(interval=rep(c('layer 1-->layer 2','layer 2-->layer 3','layer 3-->layer 4','layer 4-->layer 5','layer 5-->layer 6'),each=78),value=ile)
d <- d[d$value!=0,]
p=ggplot(d, aes(x=value))+geom_histogram(fill='royalblue3', bins=50)+facet_wrap(~interval)
ggplotToBrowser(p)
svg('interlayer_edges_distribution.svg', width = 15, height = 10)
p
dev.off()


## Analysis of modularity for each layer separately

for (l in 1:6){
  Q.layer <- unname(data.matrix(read.table(paste(output_folder,'/Q_obs_zero_layer_',l,'.csv',sep=''),he=F,sep=',')))
  n.layer <- unname(data.matrix(read.table(paste(output_folder,'/n_obs_zero_layer_',l,'.csv',sep=''),he=F,sep=',')))
  cat(paste('Layer:',l));cat('\t')
  cat(round(mean(Q.layer),2));cat('\t\t')
  cat(round(mean(n.layer),2));cat('\n')
}

## Interlayer edges are inifinity
Q.multilayerlObs_inf <- unname(data.matrix(read.table(paste(output_folder,'/Q_obs_infinity.csv',sep=''),he=F,sep=',')))
S.multilayerlObs_inf <- unname(data.matrix(read.table(paste(output_folder,'/S_obs_infinity.csv',sep=''),he=F,sep=',')))
# Value of the maximum modularity across 100 runs:
mean(Q.multilayerlObs_inf)
# Mean number of modules across runs
mean(apply(S.multilayerlObs_inf, 2, FUN = max))





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
extinctionProcess=matrix(unlist(extantSpeciesMonolayer(monolayerNetwork)),1,s1+s3)
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


## Scenario 3 -- Multilayer (remove flower visitors + uniformly random plant removal and follow parasitoids’ tertiary extinction)
flowerExtinctionProbability <- 0.3 # This parameter is to set our low (0.3) and high (0.8) scenarios
extinctionProcessList=list()
# because removal of plants is random we do this a 100 times and take average at the end.
for (run in 1:20){
  print(run)
  multilayerNetwork <- multilayerNetworkOrig
  extinctionProcess=matrix(unlist(extantSpecies(multilayerNetwork,s1,s2,s3,F)),1,length(colnames(multilayerNetwork))) 
  colnames(extinctionProcess) <- colnames(multilayerNetwork)
  extinctionOrder <- extantSpecies(multilayerNetwork,s1,s2,s3,F)$extantSet2
  extinctionOrder <- names(sort(extinctionOrder, decreasing = T))
  for (x in 1:length(extinctionOrder)){
    #print(paste(run,x,sep='--'))
    #1. Remove a node from Set2 in layer 1
    multilayerNetwork <- removeSpeciesMultilayer(multilayerNetwork,extinctionOrder[x]) 
    #2. see if any nodes from Set1 in layer 1 have a degree of 0
    extinctSet1_1 <- which(extantSpecies(multilayerNetwork,s1,s2,s3,F)$extantSet1_1==0)
    if (length(extinctSet1_1)>0){
      #3. If so, remove these nodes from layer 2 as well
      for (E in names(extinctSet1_1)){multilayerNetwork <- removeSpeciesMultilayer(multilayerNetwork,E)}
    }
    
    # 4. Randomly select a node from Set1 in layer 1 with probability 0.2, and remove it
    if (runif(1)<flowerExtinctionProbability){
      randNode <- names(sample(extantSpecies(multilayerNetwork,s1,s2,s3,F)$extantSet1_1,1))
      multilayerNetwork <- removeSpeciesMultilayer(multilayerNetwork,randNode)
    }
    
    #5. Remove nodes from Set3 in layer 2 which are disconnected
    extinctSet3 <- which(extantSpecies(multilayerNetwork,s1,s2,s3,F)$extantSet3==0)
    if (length(extinctSet3)>0){
      for (E in names(extinctSet3)){multilayerNetwork <- removeSpeciesMultilayer(multilayerNetwork,E)} 
    }
    #6. Record
    extinctionProcess <- rbind(extinctionProcess,unlist(extantSpecies(multilayerNetwork,s1,s2,s3,F)))
  }
  rownames(extinctionProcess) <- c('Initial',extinctionOrder)
  extinctionProcessList[[run]] <- extinctionProcess
}
sapply(extinctionProcessList,dim)
extinctionProcessArray <- list2Array(extinctionProcessList)
extinctionProcess <- apply(extinctionProcessArray,1:2,mean)

# Plot
#extinctionProcess <- extinctionProcess[1:which(rowSums(extinctionProcess)==0)[1],]
Set1Extinction <- extinctionProcess[,1:s1]
Set1Extinction <- apply(Set1Extinction, 1, function(x) sum(x!=0))/s1
Set3Extinction <- extinctionProcess[,(2*s1+s2+1):N]
Set3Extinction <- apply(Set3Extinction, 1, function(x) sum(x!=0))/s3
dScenario3 <- data.frame(order = rep(1:nrow(extinctionProcess),2),
                         extant = c(Set1Extinction,Set3Extinction),
                         Extinction=c(rep('Plants secondary extinction',nrow(extinctionProcess)),rep('Parasitoids tertiary extinction',nrow(extinctionProcess))),
                         Guild=c(rep('Plants',nrow(extinctionProcess)),rep('Parasitoids',nrow(extinctionProcess))))
dScenario3 <- dScenario3[with(dScenario3, order(order,Extinction)), ]
dScenario3$extant <- as.numeric(dScenario3$extant)
dScenario3$order <- as.numeric(dScenario3$order)
dScenario3$order <- dScenario3$order/max(dScenario3$order)
plt=ggplot(dScenario3, aes(x=order, y=extant, group=Extinction, color=Extinction))+geom_point(size=4)+geom_line()+facet_wrap(~Extinction,scales="free_y")+
  labs(title='Flower visitors removed by degree',x='Primary extinction', y='Number of species remaining')+ scale_colour_manual(values=c('orange','blue'))
# ggplotToBrowser(plt, w = 50,20)

# Now repeat the exact same code but with a different extinction probability of the flowers.
extinctionProcessList=list()
flowerExtinctionProbability <- 0.8 # This parameter is to set our low (0.3) and high (0.8) scenarios
for (run in 1:20){
  print(run)
  multilayerNetwork <- multilayerNetworkOrig
  extinctionProcess=matrix(unlist(extantSpecies(multilayerNetwork,s1,s2,s3,F)),1,length(colnames(multilayerNetwork))) 
  colnames(extinctionProcess) <- colnames(multilayerNetwork)
  extinctionOrder <- extantSpecies(multilayerNetwork,s1,s2,s3,F)$extantSet2
  extinctionOrder <- names(sort(extinctionOrder, decreasing = T))
  for (x in 1:length(extinctionOrder)){
    #print(paste(run,x,sep='--'))
    #1. Remove a node from Set2 in layer 1
    multilayerNetwork <- removeSpeciesMultilayer(multilayerNetwork,extinctionOrder[x]) 
    #2. see if any nodes from Set1 in layer 1 have a degree of 0
    extinctSet1_1 <- which(extantSpecies(multilayerNetwork,s1,s2,s3,F)$extantSet1_1==0)
    if (length(extinctSet1_1)>0){
      #3. If so, remove these nodes from layer 2 as well
      for (E in names(extinctSet1_1)){multilayerNetwork <- removeSpeciesMultilayer(multilayerNetwork,E)}
    }
    
    # 4. Randomly select a node from Set1 in layer 1 with probability 0.2, and remove it
    if (runif(1)<flowerExtinctionProbability){
      randNode <- names(sample(extantSpecies(multilayerNetwork,s1,s2,s3,F)$extantSet1_1,1))
      multilayerNetwork <- removeSpeciesMultilayer(multilayerNetwork,randNode)
    }
    
    #5. Remove nodes from Set3 in layer 2 which are disconnected
    extinctSet3 <- which(extantSpecies(multilayerNetwork,s1,s2,s3,F)$extantSet3==0)
    if (length(extinctSet3)>0){
      for (E in names(extinctSet3)){multilayerNetwork <- removeSpeciesMultilayer(multilayerNetwork,E)} 
    }
    #6. Record
    extinctionProcess <- rbind(extinctionProcess,unlist(extantSpecies(multilayerNetwork,s1,s2,s3,F)))
  }
  rownames(extinctionProcess) <- c('Initial',extinctionOrder)
  extinctionProcessList[[run]] <- extinctionProcess
}
sapply(extinctionProcessList,dim)
extinctionProcessArray <- list2Array(extinctionProcessList)
extinctionProcess <- apply(extinctionProcessArray,1:2,mean)

# Plot
#extinctionProcess <- extinctionProcess[1:which(rowSums(extinctionProcess)==0)[1],]
Set1Extinction <- extinctionProcess[,1:s1]
Set1Extinction <- apply(Set1Extinction, 1, function(x) sum(x!=0))/s1
Set3Extinction <- extinctionProcess[,(2*s1+s2+1):N]
Set3Extinction <- apply(Set3Extinction, 1, function(x) sum(x!=0))/s3
dScenario3a <- data.frame(order = rep(1:nrow(extinctionProcess),2),
                          extant = c(Set1Extinction,Set3Extinction),
                          Extinction=c(rep('Plants secondary extinction',nrow(extinctionProcess)),rep('Parasitoids tertiary extinction',nrow(extinctionProcess))),
                          Guild=c(rep('Plants',nrow(extinctionProcess)),rep('Parasitoids',nrow(extinctionProcess))))
dScenario3a <- dScenario3a[with(dScenario3a, order(order,Extinction)), ]
dScenario3a$extant <- as.numeric(dScenario3a$extant)
dScenario3a$order <- as.numeric(dScenario3a$order)
dScenario3a$order <- dScenario3a$order/max(dScenario3a$order)

dScenario1$Scenario <- 'Scenario1'
dScenario2$Scenario <- 'Scenario2'
dScenario3$Scenario <- 'Scenario3_0.3'
dScenario3a$Scenario <- 'Scenario3_0.8'
d_Plot <- rbind(dScenario1,dScenario2,dScenario3,dScenario3a)

plt=ggplot(d_Plot, aes(x=order, y=extant, group=Scenario, color=Scenario))+geom_point(size=3)+geom_line()+facet_wrap(~Guild)+
  labs(title='Removal by degree',x='Proportion of species removed', y='Proportion of leafminer parasitoids remaining')+
  scale_colour_manual(values=c('dark orange','#914200','mediumpurple1','forestgreen'))+theme_Publication()
# ggplotToBrowser(plt, w = 40,20)


#################################################################################################
# Analyses of robustness -- Supplementary Information
# The algorithms and rational for the removal scenarios are explained in detail in the paper
#################################################################################################

# Create the empirical multilayer network
supraAdjMatrix <- createSupraAdjMatrix(visitorsLayer,parasitoidsLayer)
s1 <- nrow(visitorsLayer);s2 <- ncol(visitorsLayer);s3 <- ncol(parasitoidsLayer)
N <- nrow(supraAdjMatrix)
multilayerNetwork <- supraAdjMatrix  
multilayerNetworkOrig <- supraAdjMatrix

# Run extinctions
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
Set1Extinction <- extinctionProcess[,1:s1] # The monolayer scenario - only plants
Set1Extinction <- apply(Set1Extinction, 1, function(x) sum(x!=0))/s1
Set1_3Extinction <- extinctionProcess[,c(1:s1,(2*s1+s2+1):N)] # The multilayer scenario - plants + parasitoids
Set1_3Extinction <- apply(Set1_3Extinction, 1, function(x) sum(x!=0))/(s1+s3)
d <- data.frame(order = rep(1:nrow(extinctionProcess),2),
                extant = c(Set1Extinction,Set1_3Extinction),
                Extinction=c(rep('Plants secondary extinction',nrow(extinctionProcess)),rep('Plants scnd + Parasitoids tertiary extinction',nrow(extinctionProcess))),
                Guild=c(rep('Plants',nrow(extinctionProcess)),rep('Plants+Parasitoids',nrow(extinctionProcess))))
d <- d[with(d, order(order,Extinction)), ]
d$extant <- as.numeric(d$extant)
d$order <- as.numeric(d$order)
d$order <- d$order/max(d$order)
plt=ggplot(d, aes(x=order, y=extant, group=Extinction, color=Extinction))+geom_point(size=4)+geom_line()+
  labs(title='Removal by degree (lowest first)',x='Proportion of flower visitors removed', y='Proportion of species remaining')+
  scale_colour_manual(values=c('orange','blue'))
ggplotToBrowser(plt, w = 50,20)
