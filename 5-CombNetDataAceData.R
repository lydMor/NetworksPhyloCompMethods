### testing some things with your data: 
### see if you have everything: 
## read in your dataframe: 
load("ACE_TESTS_FINAL.RData")

## check it out: 
TEST<- ACE_TESTS_all %>% distinct(Replicate, Netname, sig2_True)
TEST$sig2_True<-as.factor(TEST$sig2_True)
summary(TEST)
TEST$type<- word(TEST$Netname, 1, sep="_")
TEST$type<- as.factor(TEST$type)
summary(TEST)

## check these: 
Types<- TEST %>% group_by(type) %>% summarize(count=n())
## check on what the net types say: 
check<- unique(TEST$type)
check


### okay, make a summary that you can add to your other things: (ACE summaries)
SAMP<- ACE_TESTS_all %>% slice(1:5)
names<- data.frame(names=colnames(ACE_TESTS_all))
## the F needs to be fixed: (the full thing vs the not full thing) (the way I collected the tree)
## with and without knowledge of extinction. 


NetSummaries<- ACE_TESTS_all %>% select(Netname, NetTopology, sig2_True, Replicate, ModelUsed, rate_error, meanError, numTrans, Replicate)

NetSummaries$Netname<- paste(word(NetSummaries$Netname, 1, sep="\\."))
Names<- data.frame(Netname=unique(NetSummaries$Netname))
NetSummaries<-NetSummaries %>% distinct(Netname, sig2_True, Replicate, .keep_all=T)




## read in your basic network stats: 
Data_AllNets <- read.csv("NetsData_variableAndstableInher.csv")


### get rid of things you don't want and combine these: 
Data_all_Comb<- Data_AllNets %>% select(-NetNum, -finalNetNum)
check1<- Data_all_Comb %>% distinct(Data_all_Comb$Netname)
check2<- NetSummaries %>% distinct(Netname)
check3<- Data_all_Comb %>% filter(!(Netname %in% check2$Netname)) %>% select(Netname) %>% distinct()

Data_comb_test<- merge(Data_all_Comb, NetSummaries, all=T)
## I think combining is weird bc you don't have all the matches, and so theyre just making crazy duplicates all over.
## okay, see what happens here:  (get rid of any data that aren't in your ACEs so far)

## this should be the same length as the number of networks in your ACE summaries(NEtsummaries)

## add F to the netNames: 
Data_Comb3<- Data_all_Comb
Data_Comb3$Netname2<- paste(word(Data_Comb3$Netname, 1, sep="_"), "F", sep="")
Data_Comb3$Netname2<- paste(Data_Comb3$Netname2, word(Data_Comb3$Netname, 2, sep="_"), sep="_")
Data_Comb3$Netname2[1]
Data_Comb3$Netname<- Data_Comb3$Netname2
Data_Comb3<- Data_Comb3 %>% select(-Netname2)

Data_Comb4<- rbind(Data_all_Comb, Data_Comb3)
Data_Comb5<- Data_Comb4 %>% filter(Netname %in% NetSummaries$Netname)
ACE_summaries_all<- full_join(NetSummaries, Data_Comb5)

### check on these bc wtf (zeros and NAs??)
Zeroes<- ACE_summaries_all %>% filter(is.na(nHyb))
#### you should also have nrow(names)*500 entries in your summary sheet(AKA the number of unique networks * 500 replicates for each network,)
## ok so everything is good. 
## BUT let's see if there's even a difference in topology between the full and not full networks that Kevin used: 
Nets_test<- ACE_summaries_all %>% filter(Netname=="BLN_13" | Netname=="BLNF_13")
Nets_test<- Nets_test %>% distinct(Netname, .keep_all=T)
library(SiPhyNetwork)
library(ape)

Net1<- read.net(text=Nets_test$NetTopology[1])
plot(Net1)
Net2<- read.net(text=Nets_test$NetTopology[2])
plot(Net2)
## ok, so the exported network is different. 

### save this: 
write.csv(ACE_summaries_all, "netData_AceData_AllInherFULL.csv", row.names=F)

