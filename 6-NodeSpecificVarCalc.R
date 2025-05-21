#### get node specific variables:
load("ACE_TESTS_JOINED_FINAL.RData")

## select just the nodes: 
rm(ACE_TESTS_both)
gc()

#load your libraries: 
library(ape)
library(SiPhyNetwork)
library(foreach)
library(doParallel)
library(tidyverse)

## set up your clusters: 
## see how many you have available: 
detectCores()
cl <- makeCluster(2)
registerDoParallel(cl)


### do the loop: 
list<- split(ACE_TESTS_nodes, f=ACE_TESTS_nodes$Netname)
#make your selection function: 
select.tip.or.node <- function(element, net) {
  ifelse(element < Ntip(net)+1,
         net$tip.label[element],
         net$node.label[element-Ntip(net)])
}
# create your dataframe: 

Node_Vars<- foreach( i = 1:length(list), .combine=rbind, .inorder=TRUE, .packages= c("SiPhyNetwork", "ape", "tidyverse", "foreach")) %dopar%{
  temp<- list[[i]]
  net<- read.net(text=temp$NetTopology[1])
  ## get branch lengths: 
  distances <- cbind(net$edge, net$edge.length)
  colnames(distances) <- c("from", "to", "distance")
  distances<-as.data.frame(distances)
  ## create your edge reference table: 
  edge_table <- data.frame(
    "node" = net$edge[,1],
    "node_label" = sapply(net$edge[,1],
                          select.tip.or.node,
                          net = net)
  )
  edge_table<- edge_table %>% distinct()
  if(nrow(edge_table>0)){
    hyb_edge<- filter(edge_table, grepl('H', node_label))
    ## get node ages: 
    depths<-data.frame(depth=node.depth.edgelength(net))
    depths$node<- rownames(depths)
    Hyb_depths<- depths %>% filter(node %in% hyb_edge$node)
    ## now get your shit: 
    foreach(j=1:nrow(temp), .combine=rbind, .inorder=TRUE) %do% {
      NodeLab<- temp$Node[j]
      NodeNum<- edge_table %>% filter(node_label==NodeLab) 
      NodeNum<- NodeNum$node
      #branch length leading to and from the node of interest: 
      BrLenTo<- distances %>% filter(to==NodeNum) %>% select(distance)
      BrLenFrom<- distances %>% filter(from==NodeNum) %>% summarize(from=mean(distance))
      ## age of node: 
      age<- depths %>% filter(node==NodeNum) %>% select(depth)
      Hyb_Dist<- Hyb_depths
      Hyb_Dist$dist<- age$depth-Hyb_Dist$depth
      Hyb_Dist$absDist<- abs(Hyb_Dist$dist)
      NODE<- Hyb_Dist %>% dplyr::filter(absDist== min(absDist))
      Hyb_node<- hyb_edge %>% filter(node==NODE$node)
      
      ###OK NOW PUT ALL OF YOUR THINGS TOGETHER: 
      #node ages and distances: 
      temp$NodeAge[j]<- paste(age$depth[1])
      temp$MinDistH[j]<- paste(NODE$dist[1])
      temp$MinPosDistH[j]<- paste(min(Hyb_Dist$dist[Hyb_Dist$dist>0]))
      temp$MinNegDistH[j]<-paste(min(Hyb_Dist$dist[Hyb_Dist$dist<0]))
      temp$ClosestHyb[j]<-paste(Hyb_node$node_label[1])
      #branch lengths before and after the node: 
      temp$brLenTo[j]<- paste(BrLenTo$distance[1])
      temp$brLenFrom[j]<- paste(BrLenFrom$from[1])
    }
  }
  else{
    temp$NodeAge<- paste("FIX")
    temp$MinDistH<- paste("FIX")
    temp$MinPosDistH<- paste("FIX")
    temp$MinNegDistH<-paste("FIX")
    #branch lengths before and after the node: 
    temp$brLenTo<- paste("FIX")
    temp$brLenFrom<- paste("FIX")
    temp$ClosestHyb<- paste("FIX")
  }
  
  return(temp)
}

Node_Vars[Node_Vars=="Inf"|Node_Vars=="NA"]<-NA

save.image("NODE_SPECIFIC_WITHCLOSESTHYBNODE.RData")
## look at the ones you need to fix: 
Fix_Nodes<- Node_Vars %>% filter(brLenTo=="FIX") %>% distinct(Netname)




######### NOW GET THE INFO ABOUT TRANSGRESSIVE HYBRIDIZATION: 
ACE_TESTS_NODES<- ACE_TESTS_both %>% select(Node_index, error, sig2_True, numHyb, numTrans, TrueVal, Tree_index)

## get something for transgressive events ONLY: (DEAL WITH THESE LATERRR)

Trans_Events<- ACE_TESTS_NODES %>% filter(numTrans>0)

rm(list=setdiff(ls(), "Node_Vars"))

load("ACE_TESTS_JOINED_NO.ZEROHybEvents.RData")


ACE_TESTS_NODES<- ACE_TESTS_both %>% select(Node_index, error, sig2_True, numHyb, numTrans, TrueVal, Tree_index)
Node_ind<- Node_Vars %>% filter(!(NodeAge=="FIX"))
Node_ind<- Node_ind %>% select(-Node, -Netname)


##### combine for analysis: 
Nodes_all<- full_join(ACE_TESTS_NODES, Node_ind)

#### okay, now save this For later: 
rm(list=setdiff(ls(), "Nodes_all"))
gc()
save.image("NodeSpecificIndices_REAL_NoTransMagnitude.Rdata")


### Calculate Magnitude of Transgressive trait evolution:

### okay, now make a list of files that contain any transgressive events and loop over your trait sim dataframes: 
ACE_TESTS_TRANS<- ACE_TESTS_both %>% filter(numTrans>0) %>% distinct(Netname)
rm(list=setdiff(ls(), "ACE_TESTS_TRANS"))

#this should be the same directory as your list of files for running ancestral character estimation

TRANSGRESSIVE<- data.frame()
for(i in 1:nrow(ACE_TESTS_TRANS)){
  print(i)
  temp<- read.csv(paste(ACE_TESTS_TRANS$Netname[i], ".nwk.txt.csv", sep=""))
  ### get your tree index and shit: 
  TEST<- temp %>% filter(numTransgressive>0)
  TEST$Netname<- paste(word(TEST$filename, 1, sep="\\."))
  TEST$Tree_index<- paste(TEST$Netname, TEST$replicate, TEST$variance, sep="_")
  ## get just the transgressive events and just what you want: 
  Trans<- select(TEST, Tree_index, contains("brH"))
  Trans<- Trans %>% distinct(Tree_index, .keep_all=T)
  Trans2<- Trans %>% select(-Tree_index) %>% as.matrix()
  rownames(Trans2)<- Trans$Tree_index
  Trans3<- melt(Trans2)
  TRANSGRESSIVE<- rbind(TRANSGRESSIVE, Trans3)
}

write.csv(TRANSGRESSIVE, "TRANSGRESSIVE_MAGNITUDE.csv")
rm(list=ls())
load("NodeSpecificIndices_REAL_NoTransMagnitude.Rdata")
TRANS<- read.csv("TRANSGRESSIVE_MAGNITUDE.csv")

#### get just transgressive nodes: 
names(TRANS)<- c("x", "Tree_index", "ClosestHyb", "TransMagnitude")
TRANS<- TRANS %>% select(-x)
TRANS$ClosestHyb<- str_replace_all(TRANS$ClosestHyb, "br", "#")

Nodes_trans<- Nodes_all %>% filter(Tree_index %in% TRANS$Tree_index)

### make like a new column for joining: 
## first get your ones that don't do transgressive evolution: 
Nodes_noTrans<- anti_join(Nodes_all, Nodes_trans)
Nodes_noTrans$TransMagnitude<- 0

### NOW make your little column: 
TRANS$Mag_Tree<- paste(TRANS$Tree_index, TRANS$ClosestHyb, sep="_")
Nodes_trans$Mag_Tree<- paste(Nodes_trans$Tree_index, Nodes_trans$ClosestHyb, sep="_")
TRANS<- TRANS %>% select(Mag_Tree, TransMagnitude)

Nodes_Tjoin<- full_join(Nodes_trans, TRANS)

library(janitor)
test<- Nodes_Tjoin %>% filter(!(Mag_Tree %in% TRANS$Mag_Tree))
test2<- TRANS %>%filter(!(Mag_Tree %in% Nodes_trans$Mag_Tree))
### OKAY, SO I THINK THAT THERE ARE SOME HYBRID EVENTS THAT ARE NOT CLOSEST TO A NODE, SO THEY AREN't IN THE NODES THING: 
TRANS_join<- TRANS %>% filter(Mag_Tree %in% Nodes_trans$Mag_Tree)
Nodes_Tjoin<- full_join(Nodes_trans, TRANS_join)


### okay, now join back all the stuff: 
Nodes_Tjoin<- Nodes_Tjoin %>% select(-Mag_Tree)

Nodes_allT<- rbind(Nodes_Tjoin, Nodes_noTrans)
### OKAY! Now save this shit: 
rm(list=setdiff(ls(), "Nodes_allT"))
save.image("NodeSpecificVars_WITHTRANSMAG.RData")
