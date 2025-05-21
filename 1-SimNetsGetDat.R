### real script for getting trees: 
#### MAKING A LOOP TO EXTRACT DATA fROM TREES: 
# packages: 
library(SiPhyNetwork)
library(ape)
library(tidyverse)
### run your simulation: 
#### we're going to simulate networks under three different hybridization regimes: 
## we're keeping extinction, speciation, hybridization rates, and inherritance probs THE SAME ACROSS SIMULATIONS: 

rm(list=ls())
#### ######
### ##
####
###  SET YOUR CONSTANT THINGS: 
## extinction rate: 
mu<-.2
## speciation rate: 
lambda<-1

## Hybridization rate: 
nu<-.25

## inherritance probabilities: 
inheritance.fxn <- make.beta.draw(3,8)



######## 1: LINEAGE NEUTRAL ONLY: 

hybrid_proportions <-c(0,  ##Lineage Generative
                       0, ##Lineage Degenerative
                       1) ##Lineage Neutral


## if 50 takes way too long, try to just go down to like 10 or something for the test.
LG_nets <-sim.bdh.taxa.gsa(m=30,n=10,numbsim=600,
                           lambda=lambda,mu=mu,
                           nu=nu, hybprops = hybrid_proportions,
                           hyb.inher.fxn = inheritance.fxn)

TreeDat<- data.frame(NetNum=1:length(LG_nets), nTaxa=NA, nSur=NA,termMax=NA,termMin=NA,termMean=NA,termVar=NA,
                     intMax=NA,intMin=NA,intMean=NA,intVar=NA,termRat=NA,nHyb=NA,
                     PropTermH=NA,PropExtH=NA,HybType=NA,ExtRate=NA,HybRate=NA,
                     spRate=NA,inherMean=NA,inherMax=NA,inherVar=NA,
                     isTB=NA,isTC=NA,isFU=NA,level=NA)

for(i in 1:length(LG_nets)){
  # get your full network:
  tree<- LG_nets[[i]]
  xi<-length(tree)
  if(xi>1){
    # get your pruned network:
    treeP<-reconstructedNetwork(tree)
    
    
    ## get number of taxa: 
    TreeDat$nTaxa[i]<-length(tree$tip.label)
    # get the names of the taxa and associated ttiplables: 
    labs_F<- data.frame(names=tree$tip.label, nodes=1:length(tree$tip.label))
    # get a df of extinct taxa: 
    Exctinct<- labs_F %>% filter(!(names %in% treeP$tip.label))
    
    TreeDat$nSur[i]<- paste(nrow(labs_F)-nrow(Exctinct))
    ## get average branch lengths and shit for your network: 
    terminal<-setNames(tree$edge.length[sapply(1:length(tree$tip.label),function(x,y) which (y==x),y=tree$edge[,2])],tree$tip.label)
    TreeDat$termMax[i]<-max(terminal)
    TreeDat$termMin[i]<-min(terminal)
    TreeDat$termMean[i]<-mean(terminal)
    TreeDat$termVar[i]<-var(terminal)
    
    ## get the length of internal branches: 
    distances <- cbind(tree$edge, tree$edge.length)
    colnames(distances) <- c("from", "to", "distance")
    distances<-as.data.frame(distances)
    distances<- distances %>% filter(!(to %in% labs_F$nodes))
    TreeDat$intMax[i]<- max(distances$distance)
    TreeDat$intMin[i]<-min(distances$distance)
    TreeDat$intMean[i]<- mean(distances$distance)
    TreeDat$intVar[i]<-mean(distances$distance)
    
    ## Compare terminal to internal branch lengths: 
    TreeDat$termRat[i]<-mean(terminal)/mean(distances$distance)
    
    
    #### Information about hybridization events: 
    ret<- tree[["reticulation"]]
    ret<- as.data.frame(ret)
    
    ### number of reticulation events: 
    TreeDat$nHyb[i]<-nrow(ret)
    
    dat<-nodepath(tree)
    NODEZ<- data.frame()
    for( l in 1:length(dat)){
      node<- data.frame(node=dat[[l]][length(dat[[l]])-1])
      NODEZ<- rbind(NODEZ, node)
    }
    
    NODEZ$tip<- rownames(NODEZ)
    
    
    
    TreeDat$PropTermH[i]<- paste((length(which(c(ret$from,ret$to) %in% NODEZ$node)))/(nrow(ret)*2))
    
    check<- ret %>% slice(which(to|from %in% NODEZ$node))
    
    
    ### find the number of nodes that occur on extinct branches (or involve extinct branches IDK.): 
    
    extDat<-dat[Exctinct$nodes]
    if(length(extDat)>0){
      EN<- data.frame()
      for(b in 1:length(extDat)){
        Enodes<- data.frame(node=extDat[[b]])
        EN<- rbind(Enodes, EN)
      }
      TreeDat$PropExtH[i]<- paste((length(which(c(ret$from,ret$to) %in% EN$node)))/(nrow(ret)*2)) 
    }
    else{
      TreeDat$PropExtH[i]<-paste(0)
    }
    
    
    
    ## look at inheritance: 
    TreeDat$inherMean[i]<- mean(tree$inheritance)
    TreeDat$inherMax[i]<-max(tree$inheritance)
    if(length(tree$inheritance)>0){
      TreeDat$inherVar[i]<-var(tree$inheritance)
    }
    else{}
    
    ### get summary information for your network: 
    TreeDat$isTB[i]<-isTreeBased(tree)
    TreeDat$isTC[i]<-isTreeChild(tree)
    TreeDat$isFU[i]<-isFUstable(tree)
    TreeDat$level[i]<-getNetworkLevel(tree)
    
    
    ## add the information from your simulations: 
    ## hybrid type: 
    TreeDat$HybType<-paste("LN")
    
    
    ## extinction rate: 
    TreeDat$ExtRate[i]<-paste(mu)
    TreeDat$HybRate[i]<-paste(nu)
    TreeDat$spRate[i]<-paste(lambda)
  }
  
  else{
    
  }
}

LN_dat<- TreeDat %>% filter(isTB==T & nHyb >=1 & nHyb <=3 & nSur ==10 & PropExtH<.9)
LN_dat$finalNetNum<- 1:nrow(LN_dat)
i<-1
setwd("BLN/")
for(i in 1:100){
  num<- LN_dat$NetNum[i]
  net<- LG_nets[[num]]
  net2<- reconstructedNetwork(net)
  write.net(net2,file = paste("BLN_", i, sep=""))
  write.net(net, file= paste("BLNF_", i, sep=""))
}

setwd("../")
write.csv(LN_dat, "Data_BLN.csv", row.names=F)

rm(list=ls())
#### ######
### ##
####
###  SET YOUR CONSTANT THINGS: 
## extinction rate: 
mu<-.2
## speciation rate: 
lambda<-1

## Hybridization rate: 
nu<-.25

## inherritance probabilities: 
inheritance.fxn <- make.beta.draw(3,8)



######## 2: 75% LINEAGE GENERATIVE  

hybrid_proportions <-c(.25,  ##Lineage Generative
                       0, ##Lineage Degenerative
                       .75) ##Lineage Neutral


## if 50 takes way too long, try to just go down to like 10 or something for the test.
LN75_nets <-sim.bdh.taxa.gsa(m=30,n=10,numbsim=500,
                             lambda=lambda,mu=mu,
                             nu=nu, hybprops = hybrid_proportions,
                             hyb.inher.fxn = inheritance.fxn)

TreeDat<- data.frame(NetNum=1:length(LN75_nets), nTaxa=NA, nSur=NA,termMax=NA,termMin=NA,termMean=NA,termVar=NA,
                     intMax=NA,intMin=NA,intMean=NA,intVar=NA,termRat=NA,nHyb=NA,
                     PropTermH=NA,PropExtH=NA,HybType=NA,ExtRate=NA,HybRate=NA,
                     spRate=NA,inherMean=NA,inherMax=NA,inherVar=NA,
                     isTB=NA,isTC=NA,isFU=NA,level=NA)

for(i in 1:length(LN75_nets)){
  # get your full network:
  tree<- LN75_nets[[i]]
  xi<-length(tree)
  if(xi>1){
    # get your pruned network:
    treeP<-reconstructedNetwork(tree)
    
    
    ## get number of taxa: 
    TreeDat$nTaxa[i]<-length(tree$tip.label)
    # get the names of the taxa and associated ttiplables: 
    labs_F<- data.frame(names=tree$tip.label, nodes=1:length(tree$tip.label))
    # get a df of extinct taxa: 
    Exctinct<- labs_F %>% filter(!(names %in% treeP$tip.label))
    
    TreeDat$nSur[i]<- paste(nrow(labs_F)-nrow(Exctinct))
    ## get average branch lengths and shit for your network: 
    terminal<-setNames(tree$edge.length[sapply(1:length(tree$tip.label),function(x,y) which (y==x),y=tree$edge[,2])],tree$tip.label)
    TreeDat$termMax[i]<-max(terminal)
    TreeDat$termMin[i]<-min(terminal)
    TreeDat$termMean[i]<-mean(terminal)
    TreeDat$termVar[i]<-var(terminal)
    
    ## get the length of internal branches: 
    distances <- cbind(tree$edge, tree$edge.length)
    colnames(distances) <- c("from", "to", "distance")
    distances<-as.data.frame(distances)
    distances<- distances %>% filter(!(to %in% labs_F$nodes))
    TreeDat$intMax[i]<- max(distances$distance)
    TreeDat$intMin[i]<-min(distances$distance)
    TreeDat$intMean[i]<- mean(distances$distance)
    TreeDat$intVar[i]<-mean(distances$distance)
    
    ## Compare terminal to internal branch lengths: 
    TreeDat$termRat[i]<-mean(terminal)/mean(distances$distance)
    
    
    #### Information about hybridization events: 
    ret<- tree[["reticulation"]]
    ret<- as.data.frame(ret)
    
    ### number of reticulation events: 
    TreeDat$nHyb[i]<-nrow(ret)
    
    dat<-nodepath(tree)
    NODEZ<- data.frame()
    for( l in 1:length(dat)){
      node<- data.frame(node=dat[[l]][length(dat[[l]])-1])
      NODEZ<- rbind(NODEZ, node)
    }
    
    NODEZ$tip<- rownames(NODEZ)
    
    
    
    TreeDat$PropTermH[i]<- paste((length(which(c(ret$from,ret$to) %in% NODEZ$node)))/(nrow(ret)*2))
    
    check<- ret %>% slice(which(to|from %in% NODEZ$node))
    
    
    ### find the number of nodes that occur on extinct branches (or involve extinct branches IDK.): 
    
    extDat<-dat[Exctinct$nodes]
    if(length(extDat)>0){
      EN<- data.frame()
      for(b in 1:length(extDat)){
        Enodes<- data.frame(node=extDat[[b]])
        EN<- rbind(Enodes, EN)
      }
      TreeDat$PropExtH[i]<- paste((length(which(c(ret$from,ret$to) %in% EN$node)))/(nrow(ret)*2)) 
    }
    else{
      TreeDat$PropExtH[i]<-paste(0)
    }
    
    
    
    ## look at inheritance: 
    TreeDat$inherMean[i]<- mean(tree$inheritance)
    TreeDat$inherMax[i]<-max(tree$inheritance)
    if(length(tree$inheritance)>0){
      TreeDat$inherVar[i]<-var(tree$inheritance)
    }
    else{}
    
    ### get summary information for your network: 
    TreeDat$isTB[i]<-isTreeBased(tree)
    TreeDat$isTC[i]<-isTreeChild(tree)
    TreeDat$isFU[i]<-isFUstable(tree)
    TreeDat$level[i]<-getNetworkLevel(tree)
    
    
    ## add the information from your simulations: 
    ## hybrid type: 
    TreeDat$HybType<- paste("LN75")
    
    
    ## extinction rate: 
    TreeDat$ExtRate[i]<-paste(mu)
    TreeDat$HybRate[i]<-paste(nu)
    TreeDat$spRate[i]<-paste(lambda)
  }
  
  else{
    
  }
}

LN75_dat<- TreeDat %>% filter(isTB==T & nHyb >=1 & nHyb <=3 & nSur ==10 & PropExtH<.9)

LN75_dat$finalNetNum<-paste(1:nrow(LN75_dat))



setwd("BLN75/")

for(i in 1:100){
  num<- LN75_dat$NetNum[i]
  net<- LN75_nets[[num]]
  net2<- reconstructedNetwork(net)
  write.net(net2,file = paste("BLN75_", i, sep=""))
  write.net(net, file= paste("BLN75F_", i, sep=""))
  
}
setwd("../")

write.csv(LN75_dat, "Data_BLN75.csv", row.names=F)

rm(list=ls())
#### ######
### ##
####
###  SET YOUR CONSTANT THINGS: 
## extinction rate: 
mu<-.2
## speciation rate: 
lambda<-1

## Hybridization rate: 
nu<-.25

## inherritance probabilities: 
inheritance.fxn <- make.beta.draw(3,8)



######## 3: 50% LINEAGE GENERATIVE 50% LINEAGE NEUTRAL  

hybrid_proportions <-c(.5,  ##Lineage Generative
                       0, ##Lineage Degenerative
                       .5) ##Lineage Neutral


## if 50 takes way too long, try to just go down to like 10 or something for the test.
LN50_nets <-sim.bdh.taxa.gsa(m=30,n=10,numbsim=500,
                             lambda=lambda,mu=mu,
                             nu=nu, hybprops = hybrid_proportions,
                             hyb.inher.fxn = inheritance.fxn)

TreeDat<- data.frame(NetNum=1:length(LN50_nets), nTaxa=NA, nSur=NA,termMax=NA,termMin=NA,termMean=NA,termVar=NA,
                     intMax=NA,intMin=NA,intMean=NA,intVar=NA,termRat=NA,nHyb=NA,
                     PropTermH=NA,PropExtH=NA,HybType=NA,ExtRate=NA,HybRate=NA,
                     spRate=NA,inherMean=NA,inherMax=NA,inherVar=NA,
                     isTB=NA,isTC=NA,isFU=NA,level=NA)

for(i in 1:length(LN50_nets)){
  # get your full network:
  tree<- LN50_nets[[i]]
  xi<-length(tree)
  if(xi>1){
    # get your pruned network:
    treeP<-reconstructedNetwork(tree)
    
    
    ## get number of taxa: 
    TreeDat$nTaxa[i]<-length(tree$tip.label)
    # get the names of the taxa and associated ttiplables: 
    labs_F<- data.frame(names=tree$tip.label, nodes=1:length(tree$tip.label))
    # get a df of extinct taxa: 
    Exctinct<- labs_F %>% filter(!(names %in% treeP$tip.label))
    
    TreeDat$nSur[i]<- paste(nrow(labs_F)-nrow(Exctinct))
    ## get average branch lengths and shit for your network: 
    terminal<-setNames(tree$edge.length[sapply(1:length(tree$tip.label),function(x,y) which (y==x),y=tree$edge[,2])],tree$tip.label)
    TreeDat$termMax[i]<-max(terminal)
    TreeDat$termMin[i]<-min(terminal)
    TreeDat$termMean[i]<-mean(terminal)
    TreeDat$termVar[i]<-var(terminal)
    
    ## get the length of internal branches: 
    distances <- cbind(tree$edge, tree$edge.length)
    colnames(distances) <- c("from", "to", "distance")
    distances<-as.data.frame(distances)
    distances<- distances %>% filter(!(to %in% labs_F$nodes))
    TreeDat$intMax[i]<- max(distances$distance)
    TreeDat$intMin[i]<-min(distances$distance)
    TreeDat$intMean[i]<- mean(distances$distance)
    TreeDat$intVar[i]<-mean(distances$distance)
    
    ## Compare terminal to internal branch lengths: 
    TreeDat$termRat[i]<-mean(terminal)/mean(distances$distance)
    
    
    #### Information about hybridization events: 
    ret<- tree[["reticulation"]]
    ret<- as.data.frame(ret)
    
    ### number of reticulation events: 
    TreeDat$nHyb[i]<-nrow(ret)
    
    dat<-nodepath(tree)
    NODEZ<- data.frame()
    for( l in 1:length(dat)){
      node<- data.frame(node=dat[[l]][length(dat[[l]])-1])
      NODEZ<- rbind(NODEZ, node)
    }
    
    NODEZ$tip<- rownames(NODEZ)
    
    
    
    TreeDat$PropTermH[i]<- paste((length(which(c(ret$from,ret$to) %in% NODEZ$node)))/(nrow(ret)*2))
    
    check<- ret %>% slice(which(to|from %in% NODEZ$node))
    
    
    ### find the number of nodes that occur on extinct branches (or involve extinct branches IDK.): 
    
    extDat<-dat[Exctinct$nodes]
    if(length(extDat)>0){
      EN<- data.frame()
      for(b in 1:length(extDat)){
        Enodes<- data.frame(node=extDat[[b]])
        EN<- rbind(Enodes, EN)
      }
      TreeDat$PropExtH[i]<- paste((length(which(c(ret$from,ret$to) %in% EN$node)))/(nrow(ret)*2)) 
    }
    else{
      TreeDat$PropExtH[i]<-paste(0)
    }
    
    
    
    ## look at inheritance: 
    TreeDat$inherMean[i]<- mean(tree$inheritance)
    TreeDat$inherMax[i]<-max(tree$inheritance)
    if(length(tree$inheritance)>0){
      TreeDat$inherVar[i]<-var(tree$inheritance)
    }
    else{}
    
    ### get summary information for your network: 
    TreeDat$isTB[i]<-isTreeBased(tree)
    TreeDat$isTC[i]<-isTreeChild(tree)
    TreeDat$isFU[i]<-isFUstable(tree)
    TreeDat$level[i]<-getNetworkLevel(tree)
    
    
    ## add the information from your simulations: 
    ## hybrid type: 
    TreeDat$HybType<- paste("LN50")
    
    
    ## extinction rate: 
    TreeDat$ExtRate[i]<-paste(mu)
    TreeDat$HybRate[i]<-paste(nu)
    TreeDat$spRate[i]<-paste(lambda)
  }
  
  else{
    
  }
}

LN50_dat<- TreeDat %>% filter(isTB==T & nHyb >=1 & nHyb <=3 & nSur ==10 & PropExtH<.9)

LN50_dat$finalNetNum<-paste(1:nrow(LN50_dat))



setwd("BLN50/")

for(i in 1:100){
  num<- LN50_dat$NetNum[i]
  net<- LN50_nets[[num]]
  net2<- reconstructedNetwork(net)
  write.net(net2,file = paste("BLN50_", i, sep=""))
  write.net(net, file= paste("BLN50F_", i, sep=""))
  
}
setwd("../")

write.csv(LN50_dat, "Data_BLN50.csv", row.names=F)




