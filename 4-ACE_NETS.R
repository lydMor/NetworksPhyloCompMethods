### SECOND TRY- PARALLELIZATION: 
library(phytools)
library(vegan)
library(phylotools)
library(tidyverse)
library(reshape2)
library(geiger)
library(foreach)
library(doParallel)

## List your files: 
list<- list.files(pattern=".csv")

# Make the "cluster"
cluster <- makeCluster(10)
registerDoParallel(cluster)
print("CLUSTEREDD")

## make the loop: 
ACE_TESTS2<- foreach(i=1:length(list), .combine=rbind, .inorder=TRUE, .packages=c("foreach", "phytools", "vegan", "phylotools", 
                                                                                     "tidyverse", "reshape2", "geiger", "doParallel")) %dopar%{
  Nets<- read.csv(list[[i]])
  
  DF<-foreach::foreach(j=1:nrow(Nets), .combine=rbind, .inorder = TRUE) %do% {
    Tree<- read.tree(text=Nets$MajTreeTopology[j])
    rep<- Nets$replicate[j]
    ##3: Get a vector with leaves and character values: 
    tips<- Tree$tip.label
    Chars<- Nets %>% slice(as.numeric(j))
    Chars<- Chars %>% select(which(colnames(Nets) %in% tips))
    
    Chars2<- melt(Chars)
    Chars3<- Chars2$value
    names(Chars3)=Chars2$variable
    ##4: Estimate ancestral characters: 
    
    ## first test for the best mode of evolution: 
    fitBM<-fitContinuous(Tree,Chars3)
    fitOU<-fitContinuous(Tree,Chars3,model="OU")
    aic<-setNames(c(fitBM$opt$aicc,fitOU$opt$aicc),
                  c("BM","OU"))
    aic<- as.data.frame(aic)
    best<-aic %>% filter(aic== min(aic))
    model<- rownames(best)
    
    #### Run your ACE with the chosen model: 
    test<- anc.ML(Tree, Chars3, model=model, CI=T)
    
    
    ### now get it in a dataframe: 
    Estimated<- data.frame(Node=Tree$node.label, Est=test$ace, LogLik=test$logLik)
    
    ## now get your true values: 
    True<- Nets[j,]
    True2<- True %>% select(which(colnames(True) %in% Estimated$Node))
    True2<- melt(True2)
    names(True2)<- c("Node", "TrueVal")
    True2$Replicate<- True$replicate
    True2$Netname<- True$filename
    True2$TreeTopology<- True$MajTreeTopology
    True2$NetTopology<-True$NetTopology
    True2$sig2_True<- True$variance
    True2$ModelUsed<- model
    True2$numHyb<- True$numHybrids
    True2$ancMean<- True$ancestralmean
    True2$numTrans<- True$numTransgressive
    
    ### get all of it in a dataframe: 
    df<- full_join(True2, Estimated)
    
    ## 5: Estimate rate of evolution: 
    df$sig2<- test$sig2
    
    
    ## 6: Calculate Error at each node, and for the tree overall: 
    df$error<- df$Est-df$TrueVal
    df$rate_error<- df$sig2-df$sig2_True
    df$meanError<- mean(df$error)
    
    ### display the model of evolution fit statistics: 
    fitBM<-fitContinuous(Tree,Chars3)
    fitOU<-fitContinuous(Tree,Chars3,model="OU")
    fitEB<-fitContinuous(Tree,Chars3, model="EB")
    
    aic.vals<- data.frame(BM_aic=fitBM$opt$aicc, OU_aic=fitOU$opt$aicc, EB_aic=fitEB$opt$aicc, Netname=True2$Netname[1], Replicate=True2$Replicate[1])
    
    df2<- full_join(df, aic.vals)
    return(df2)
    
  }
  return(DF)
  save.image("ACE_TESTS3.RData")
}

## okay, now combine this with the old ACE_TESTS: 
rm(list=setdiff(ls(), "ACE_TESTS2"))
gc()
## now save your workspace: 
save.image("ACE_TESTS_FINAL.RData")

