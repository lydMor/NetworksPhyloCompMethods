# load your libraries: 
rm(list=ls())
library(randomForest)
library(lmer)
library(phylotools)
library(SiPhyNetwork)
library(ape)
library(tidyverse)

####### 1: read in your list of topologies, with relevant data: 

TreeBM<- read.csv("TreeTopologiesForSim.csv")
#### set up your simulations: 
library(foreach)
library(doParallel)

RATES=c(.1,.5,1,1.5,2)
RATES<- as.list(RATES)
rm(list=setdiff(ls(), c("RATES", "TreeBM")))
    ##make your cluster: 
cluster <- makeCluster(4)
registerDoParallel(cluster)

Tree_SimBM22<- foreach(i=1:length(RATES), .combine=rbind, .inorder=TRUE,.packages=c("tidyverse", "phytools", "phylotools", 
                                                                                   "SiPhyNetwork", "geiger", "ape", "foreach")) %dopar% {
  rate<- RATES[[i]]
  
  DF<- foreach::foreach(j=1:nrow(TreeBM), .combine=rbind, .inorder=TRUE) %do% {
    Tree<- read.tree(text=TreeBM$TreeTopology[j])
    #do a fast BM and get Chars: 
    x <- fastBM(Tree, sig2 = rate, internal = TRUE)
    Chars<- as.data.frame(x)
    names(Chars)[1]<-"True_val"
    Chars$node<- rownames(Chars)
    #get tip states in a named vector: 
    tips<- Tree$tip.label
    tipStates<- Chars %>% filter(node %in% tips)
    
    Chars3<- tipStates$True_val
    names(Chars3)=tipStates$node
    ##4: Estimate ancestral characters: 
    
    ## first test for the best mode of evolution: 
    fitBM<-fitContinuous(Tree,Chars3)
    fitOU<-fitContinuous(Tree,Chars3,model="OU")
    aic<-setNames(c(fitBM$opt$aicc,fitOU$opt$aicc),
                  c("BM","OU"))
    aic<- as.data.frame(aic)
    best<-aic %>% filter(aic== min(aic))
    model<- rownames(best)
    
    ### now do your ACE and get values: 
    test<- anc.ML(Tree, Chars3, model=model, CI=T)
    
    ### now get it in a dataframe: 
    Estimated<- data.frame(Node=Tree$node.label, Est=test$ace, LogLik=test$logLik)
    
    True2<- Chars %>% filter(!(node %in% tipStates$node))
    
    ### get all of it in a dataframe: 
    df<- cbind(True2, Estimated)
    df$error<- df$True_val-df$Est
    #make a dataframe of info: 
    Things<- data.frame(mean_error=mean(df$error), sig2_Est=test$sig2, sig2_True=rate, 
                        ModelChosen=model, TreeTopology=TreeBM$TreeTopology[j], rateError=test$sig2-rate,
                        BM_aic=aic$aic[1], OU_aic=aic$aic[2], logLik=test$logLik)
    
    return(Things)
  }
  return(DF)
}

write.csv(Tree_SimBM22, "TREE_CharSims_10K.csv", row.names=F)