##### BIG DATA ANALYSIS SCRIPT FOR HPRC:::: 
# load your libraries: 
rm(list=ls())
library(randomForest)
library(lmer)
library(phylotools)
library(SiPhyNetwork)
library(ape)
library(tidyverse)


##### 1: Compare networks and Trees: 
############################################################
Tree_SimBM22<- read.csv("TREE_CharSims_10K.csv")
###### ok, now look at some stuff: 
summary(as.factor(Tree_SimBM22$ModelChosen))
### model chosen is almost entirely BM: only 7% of the time is the lowest AIC model OU!! 


##### rate error: (plotted by true rate:)
Tree_SimBM2$lrate_error<- Tree_SimBM2$rateError + (1- min(Tree_SimBM2$rateError))
#log transform: 
Tree_SimBM2$lrate_error<- log(Tree_SimBM2$lrate_error)


Tree_rErrorPlot<-ggplot(Tree_SimBM2, aes(x=as.factor(sig2_True), y=lrate_error, color=as.factor(sig2_True))) +
  geom_boxplot()+ 
  xlab("True rate")+
  ylab("standardized rate error")+
  theme(legend.position="none")+
  scale_color_viridis_d()+
  ggtitle("Tree-Based Evolution")


###### mean tree-wide estimation error: 

Tree_tErrorPlot<-ggplot(Tree_SimBM2, aes(x=as.factor(sig2_True), y=mean_error, color=as.factor(sig2_True))) +
  geom_boxplot()+ 
  xlab("True rate")+
  ylab("Tree-wide trait estimation error")+
  theme(legend.position="none")+
  scale_color_viridis_d()+
  ggtitle("Tree-Based Evolution")

###### BM VS OU: 
prop<- data.frame(value=c(7, 93), group=c("OU", "BM"))
pie(prop$value, labels=c("OU-7%", "BM-93%"), col=c("#440154FF","#22A884FF"), border = "white", radius=1.25)
## do it in ggplot:
prop<- Tree_SimBM2 %>% group_by(ModelChosen) %>% summarize(value=n())

prop <- prop %>% 
  arrange(desc(ModelChosen)) %>%
  mutate(prop = value / sum(prop$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

ggplot(prop, aes(x="", y=prop, fill=ModelChosen)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("#440154FF","#22A884FF"))+
  theme_void()+
  theme(legend.position="none")+ 
  geom_text(aes(label = paste(ModelChosen, paste(prop,"%",sep=""), sep="-"), y=ypos), color = "white", size=5, nudge_y = 2)

##### look at them together: ############## 
save.image("TREE_ANALYSES_10K_ONLYTREEPLOTS.RData")



##### load in ALL of your things (your ACE_summaries_ALL from step 5): 
ACE_summaries_all<- read.csv("netData_AceData_AllInherFULL.csv")

bm_rates<- ACE_summaries_all%>% filter(ModelUsed=="BM")


TreeSim<- read.csv("TREE_CharSims_10K.csv")
bm_tree<- TreeSim %>% filter(ModelChosen=="BM")


#### first, get your rates errors: 
nets_rates<- bm_rates %>% select(sig2_True, rate_error)
nets_rates$group<- "Nets"
trees_rates<- bm_tree %>% select(sig2_True, rateError)
trees_rates$group<- "Trees"
names(trees_rates)[2]<-"rate_error"

all_rates<- rbind(trees_rates, nets_rates)
all_rates$lrate_error<- all_rates$rate_error +(1-all_rates$rate_error)
all_rates$lrate_error<- log(all_rates$lrate_error)

ggplot(all_rates, aes(x=as.factor(sig2_True), y=rate_error, color=as.factor(group))) +
  geom_boxplot()+ 
  xlab("True rate")+
  ylab("rate estimation error")+
  stat_compare_means()+
  scale_color_manual(values=c("#440154FF","#22A884FF"), name="Evolutionary Scenario")


#### Now get your char estimation errors: 
nets_traits<- ACE_summaries_all %>% select(meanError, sig2_True)
nets_traits$group<- "Nets"
trees_traits<- TreeSim %>% select(mean_error, sig2_True)
names(trees_traits)[1]<-"meanError"
trees_traits$group<- "Trees"

all_traits<- rbind(nets_traits, trees_traits)

ggplot(all_traits, aes(x=as.factor(sig2_True), y=meanError, color=as.factor(group))) +
  geom_boxplot()+ 
  xlab("True rate")+
  ylab("Tree-wide trait estimation error")+
  stat_compare_means()+
  scale_color_manual(values=c("#440154FF","#22A884FF"), name="Evolutionary Scenario")
### ok now do a levennes test to see about similarity in variance: 
##### traits: 
##HOMOSCEDASTICITY: 
library(car)
result = leveneTest(meanError ~ as.factor(group), all_traits) ## variable ~ group
result
### variance in MEANS: 
t.test(trees_rates$rate_error, nets_rates$rate_error, alternative = "less", var.equal = FALSE)
###### rates: 
##homoescadasticity: 
result_r= leveneTest(rate_error~as.factor(group), all_rates)
result_r
## variance in MEANS: 
t.test(trees_traits$meanError, nets_traits$meanError, alternative = "two.sided", var.equal = FALSE)




###### DO some T-Tests to look at differences in proportions of OU choice: 
Model_Net<- ACE_summaries_all %>% group_by(ModelUsed) %>% summarize(count=n())
Model_Tree<- TreeSim %>% group_by(ModelChosen) %>% summarize(count=n())
### 2-tailed t-test: Are proportions significantly DIFFERENT?
prop.test(x=c(Model_Net$count[2], Model_Tree$count[2]), n=c(sum(Model_Net$count), sum(Model_Tree$count)), p = NULL, alternative = "two.sided",
          correct = TRUE)
prop.test(x=c(Model_Net$count[2], Model_Tree$count[2]), n=c(sum(Model_Net$count), sum(Model_Tree$count)), p = NULL, alternative = "greater",
          correct = TRUE)
#### 1-tailed t-test: Are proportions significantly GREATER in Network group?: 
## save this shit: 
save.image("TreesVsNets.RData")


##### 2: TREE-WIDE ERROR: ##################### ################ ########### ######## 
######################### ##### ####### ############ ############## #####################
#####
## read in your dataframes: 
### read in the df or remove the things you don't need:
rm(list=setdiff(ls(), c("ACE_summaries_all", "bm_rates")))
ACE_summaries_all<- read.csv("netData_AceData_AllInherFULL.csv")


PREDICTORS<- ACE_summaries_all %>% select(sig2_True, numTrans, termMax, termMin, termMean, termVar,
                                          termRat, intMax, intMin, intMean, nHyb, intVar, PropTermH,
                                          PropExtH, HybType, inherMax, inherMean, inherVar)
PRED_vect<- colnames(PREDICTORS)
rm(PREDICTORS)





#### A: MODEL CHOICE: ########## ########

MODEL<- ACE_summaries_all %>% group_by(ModelUsed) %>% summarize(count=n())
MODEL$count[2]/MODEL$count[1]+MODEL$count[2]
#do a pie chart: 
prop<- c(20, 80)
pie(prop, labels=c("OU-20%", "BM-80%"), col=c("#440154FF","#22A884FF"), border = "white")


#### run some analyses about model choice: 
library(randomForest)
Choice_RF<- ACE_summaries_all %>% select(c(PRED_vect), ModelUsed)
#check for NAs: 
Choice_NA<- Choice_RF[!complete.cases(Choice_RF), ]
##inheritance variance is naturally 0 when there's only one hybridization event: 
Choice_RF$inherVar[is.na(Choice_RF$inherVar)]<-0
#check to see if you have any NAs now: 
anyNA(Choice_RF)

rf_choice1<- randomForest(as.factor(ModelUsed)~., data=Choice_RF, importance=TRUE)
#for permutation importance: 
varImpPlot(rf_choice1, type=1)
rf_choice1

### do it with a smaller set of variables: 
##mt
## this model is mod2 in MCLUST: (Mod1 is ALL variables.)
Choice_RF2<- Choice_RF %>% select(-termMax, -termMin, -termMean, -intMax, -intMean, -intMin, -intVar, -termVar)



## MCLUST: 
library(mclust)
class<- Choice_RF2$ModelUsed
X<- Choice_RF2 %>% select(-ModelUsed, -HybType)
mod2 <- MclustDA(X, class, modelType = "EDDA")

summary.MclustDA(mod2)
plot(mod2, what="classification")
##### calculate permutation importance for MCLUST: 
library(foreach)
library(doParallel)
cluster<- makeCluster(2)
registerDoParallel(cluster)
list<- as.list(colnames(X))

PERMS<- foreach(i=1:length(list), .combine=rbind, .inorder=T, .packages=c('mclust', 'foreach', 'tidyverse')) %dopar% {
  variable=list[[i]]
  class<- Choice_RF2$ModelUsed
  df<- Choice_RF2 %>% select(-ModelUsed, -HybType)
  
  LIK<- foreach(j=1:100, .combine=rbind, .inorder=T ) %do%{
    df[,i] <- sample(df[,i],length(df[,i]), replace=F)
    Mod_T<-MclustDA(df, class, modelType = "EDDA")
    likelihood<-data.frame(lik=Mod_T$loglik)
    likelihood$var<- paste(variable)
    likelihood$rep<- paste(j)
    return(likelihood)
  }
  return(LIK)
}

PERM_means<- PERMS %>% group_by(var) %>% summarize(meanLik=mean(lik))
PERM_means$DecreaseLogLik<- mod2$loglik-PERM_means$meanLik
PERM_means$scDec<- PERM_means$DecreaseLogLik/max(PERM_means$DecreaseLogLik)

### do a boxplot of variable importance: 
ggplot(data=PERM_means, aes(x=var, y=scDec, fill=var))+
  geom_bar(stat="identity")+
  scale_fill_viridis_d()+
  xlab("variable")+
  ylab("permutation importance")+
  theme(legend.position="none")


## get out the mean values for each IMPORTANT parameter: 
OU<- mod2[["models"]][["OU"]][["parameters"]][["mean"]]
OU<- as.data.frame(OU)
names(OU)[1]<-"Mean"
OU$Var<- rownames(OU)
OU$Group<- "OU"

BM<- mod2[["models"]][["BM"]][["parameters"]][["mean"]]
BM<- as.data.frame(BM)
names(BM)[1]<-"Mean"
BM$Var<- rownames(BM)
BM$Group<- "BM"
### bind these: 
ModMeans<- rbind(BM, OU)

## simple plot: 
ggplot(ModMeans, aes(x=Var, y=Mean, group=factor(Group))) + 
  geom_line(aes(colour=factor(Group)))+
  scale_color_manual(values=c("#22A884FF","#440154FF"), name="Model Choice")



### Load your Tree-Net comparisons: 
rm(list=ls())
load("TreesVsNets_FULLDATASET_inREAL_ANALYSIS_upTo237.RData")



####### ### ######## ## RATE ERROR #### ##### ###### ###### 

#### look at categorical variation: 
### true rate: 
#### log transform your thing with an added constant: 
## translate so that your minimum value= 1: 
## in this dataset: 2.997578 = ERROR OF 0 (anything below that = UNDERESTIMATE RATES; above = OVERSTIMATE RATES)
# translate
bm_rates$lrate_error<- bm_rates$rate_error + 1-min(bm_rates$rate_error)
#log transform: 
bm_rates$lrate_error<- log(bm_rates$lrate_error)

### do a thing where we look at scaled rate error: 
bm_rates$absRate_error<- abs(bm_rates$rate_error)
#### DO SOME BOXPLOTS: 
## error by true rate: 
netRerror_TrueRate<-ggplot(bm_rates, aes(x=as.factor(sig2_True), y=rate_error, color=as.factor(sig2_True))) +
  geom_boxplot()+ 
  xlab("True rate")+
  ylab("rate error")+
  scale_color_viridis_d()+
  theme(legend.position="none")

## error by hyb type: 
netRerror_HybType<- ggplot(bm_rates, aes(x=as.factor(HybType), y=lrate_error, color=as.factor(HybType))) +
  geom_boxplot()+ 
  xlab("network category")+
  ylab("rate error")+
  scale_color_viridis_d()+
  theme(legend.position="none")

## error by number of hybridization events: 
netRerror_NumHyb<- ggplot(bm_rates, aes(x=as.factor(nHyb), y=lrate_error, color=as.factor(nHyb))) +
  geom_boxplot()+ 
  xlab("number of hybridization events")+
  ylab("rate error")+
  scale_color_viridis_d()+
  theme(legend.position="none")

## error rate by number of transgressive events: 
netRerror_NumTrans<-ggplot(bm_rates, aes(x=as.factor(numTrans), y=lrate_error, color=as.factor(numTrans))) +
  geom_boxplot()+ 
  xlab("number of transgressive events")+
  ylab("rate error")+
  scale_color_viridis_d()+
  theme(legend.position="none")
### plot number of hybridizations by number of transgressive events: 
netRerror_nhybBynTrans<- ggplot(bm_rates, aes(x=as.factor(nHyb), y=lrate_error, color=as.factor(numTrans))) +
  geom_boxplot()+ 
  xlab("number of hybridization events")+
  ylab("rate estimation error")+
  stat_compare_means()+
  scale_color_viridis_d( name="number of transgressive events")



### do a 3-way ANOVA: 
#perform three-way ANOVA
Dat<- bm_rates %>% filter(numTrans<2)
model <- aov(lrate_error ~ as.factor(numTrans)*as.factor(nHyb), data=Dat)
## nested anova: 
model2<- aov(lrate_error ~factor(nHyb), data=Dat)
summary(model2)

TUKEY<- TukeyHSD(model2)
## ok the tukey test is huge, and basically they're all significant?? I'll have to look closer at that. 
PREDICTORS<- ACE_summaries_all %>% select(sig2_True, numTrans, termMax, termMin, termMean, termVar,
                                          termRat, intMax, intMin, intMean, nHyb, intVar, PropTermH,
                                          PropExtH, HybType, inherMax, inherMean, inherVar)
PRED_vect<- colnames(PREDICTORS)
rm(PREDICTORS)
#### Do a GLM: 
Rate_errorDF<- bm_rates %>% select(c(PRED_vect), lrate_error, absRate_error)
Rate_errorDF$nHyb<- as.factor(Rate_errorDF$nHyb)
Rate_errorDF$sig2_True<- as.factor(Rate_errorDF$sig2_True)
Rate_errorDF$HybType<- as.factor(Rate_errorDF$HybType)
Rate_errorDF$numTrans<- as.factor(Rate_errorDF$numTrans)
Rate_errorDF[is.na(Rate_errorDF)]<-0
## run the glm: 
test1<- glm(lrate_error~.-absRate_error, data=Rate_errorDF)
##run this with absolute value of rate error: 

## get your summaries: 
summary(test1)
#AIC: 1175139
#Pseudo R squared: 1-(212594/328851) (.354)

test2<- glm(absRate_error~.-lrate_error, data=Rate_errorDF) 
summary(test2)
#AIC: 6214744
#Pseudo R squared: 1-(198863195/228588586) (.1300)

### let's partial out the effects of hyb type: 
hyb_vect<- glm(lrate_error~nHyb, data=Rate_errorDF)
Rate_errorDF$resid<-hyb_vect$residuals

test2<- glm(resid~. -lrate_error-nHyb, data=Rate_errorDF)
summary(test2)
AIC: 1178327
#### okau, see if you do fewer variables: 
test3<- glm(lrate_error~nHyb+numTrans, data=Rate_errorDF)
summary(test3)

### okay, now do another glm, but only with the significant variables: 
sigs<- as.data.frame(summary(test1)$coefficients)
sigs<- sigs %>% filter(`Pr(>|t|)` <.0001)

### sigVars vector: 
RateError_sig<- rownames(sigs)
NotSig<- Rate_errorDF %>% select(which(!colnames(Rate_errorDF)%in% RateError_sig)) %>% select(where(is.numeric))

rateError_2<- Rate_errorDF %>% select(-colnames(NotSig))
rateError_2$lrate_error<- Rate_errorDF$lrate_error

test4<- glm(lrate_error~., data=rateError_2)
summary(test4)
AIC: 1175314

#### now get just the variables you care about: 
RateError_3<- Rate_errorDF %>% select(nHyb, numTrans, HybType, 
                                      inherMax, inherMean, sig2_True,
                                      intMean, termMean, PropTermH, PropExtH, lrate_error)


test5<- glm(lrate_error~., data=RateError_3)
summary(test5)
#AIC: 1175465
#Pseudo r squared: 1-(212692/328851) (.353)
#P val: 1-pchisq(328851-212692, df=(736677-736658)) (0- model is better than null)
### basically same amount of variance explained (MCFADDENS PSEUDO R SQUARED)


### RANDOM FOREST: 
library(ranger)


# use range
rf_RangerRate<- ranger(lrate_error~.-rate_error, data=Rate_errorDF, write.forest =F, importance="permutation", scale.permutation.importance = T)
rf_RangerRate$variable.importance



### let's do a plot of important variables: 
ImpVar<- data.frame(rf_RangerRate$variable.importance)
names(ImpVar)[1]<-"Perm_Importance"
ImpVar$Variable<- rownames(ImpVar)
ImpVar$Group<- NA
ImpVar$Group[1:2]<-"Trait history"
ImpVar$Group[3:nrow(ImpVar)]<- "network structure"

ImpVarPlot<-ggplot(ImpVar, aes(x=reorder(Variable,-Perm_Importance), y=Perm_Importance)) +
  geom_segment( aes(x=reorder(Variable,Perm_Importance), xend=Variable, y=0, yend=Perm_Importance, color=Group)) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) + 
  xlab("Variable")+
  ylab("Permutation Importance")+ 
  scale_color_manual(values=c("#22A884FF","#440154FF"))

### now do intMax: 
RateError_intMax<-ggplot(Rate_errorDF, aes(x = intMax, y = lrate_error), ) + #the basic mapping (if you want to do a scale color, you need to add color here)
  geom_point(aes(color="#414487FF"), size=.8) +#the geometric funtion (like is point or historgram) (you can add color here, but it won't plot a legend with the scale mapping?)
  labs(x= paste("maximum internal branch length"), y=paste("rate estimation error")) + #adds title and axis labels
  theme(plot.title = element_text(size = 17), legend.position="none")+
  scale_color_manual(values="#414487FF")


### now make your panel: 
plotlist<-list(ImpVarPlot, netRerror_NumTrans,netRerror_TrueRate,netRerror_HybType, RateError_MeanInheritance, RateError_intMax)
ggarrange(plotlist=plotlist, ncol=3, nrow=2)
## now save this: 
save.image("RATES_rf_top5Vars_PanelPlot.RData")



### now see if your top important variables are significantly interacting: 
rate_impGLM<- glm(lrate_error~numTrans+sig2_True+(numTrans*sig2_True), data=Rate_errorDF)

Coeffs_Comb<- data.frame(summary(rate_impGLM)$coefficients)
names(Coeffs_Comb)[4]<- "Pval"
Coeffs_CombI<- Coeffs_Comb %>% filter(Pval < .005)


### MAKE A BOXPLOT: 
Rate_TansRateInt<-ggplot(Rate_errorDF, aes(x=as.factor(sig2_True), y=lrate_error, color=as.factor(numTrans))) +
  geom_boxplot()+ 
  xlab("True rate")+
  ylab("rate estimation error")+
  scale_color_viridis_d(name="number of Transgressive Events")+
  ggtitle("Transgressive-rate interaction")+
  geom_hline(yintercept=2.99)


Rate_TansRateInt
## save this: 
save.image("Rate_Error_interactions_RandomForest_Plots.Rdata")









######## ### ########### ### ## TREE-WIDE TRAIT ESTIMATION ERROR: ######### ##### ###### ###### ##### 
#### get a look at the distribution: 
hist(ACE_summaries_all$meanError)
min(ACE_summaries_all$meanError)
max(ACE_summaries_all$meanError)
mean(ACE_summaries_all$meanError)
#### get a log transformed mean error: 

1-min(ACE_summaries_all$meanError)
ACE_summaries_all$lmeanError<- ACE_summaries_all$meanError + (1-min(ACE_summaries_all$meanError))
## minimum value should be 1
min(ACE_summaries_all$lmeanError)
## now log transform: 
ACE_summaries_all$lmeanError<- log(ACE_summaries_all$lmeanError)
hist(ACE_summaries_all$lmeanError)

##### do some boxplots: 
## error by true rate: 
ggplot(ACE_summaries_all, aes(x=as.factor(sig2_True), y=meanError, color=as.factor(sig2_True))) +
  geom_boxplot()+ 
  xlab("True rate")+
  ylab("mean tree-wide trait estimation error")+
  scale_color_viridis_d()+
  theme(legend.position="none")


## error by hyb type: 
TreeTrait_HybType<- ggplot(ACE_summaries_all, aes(x=as.factor(HybType), y=meanError, color=as.factor(HybType))) +
  geom_boxplot()+ 
  xlab("net structure")+
  ylab("Tre-wide trait estimation error")+
  scale_color_viridis_d()+
  theme(legend.position="none")

## error by number of hybridization events: 
TreeTrait_nHyb<- ggplot(ACE_summaries_all, aes(x=as.factor(nHyb), y=meanError, color=as.factor(nHyb))) +
  geom_boxplot()+ 
  xlab("number of hybridization events")+
  ylab("Tree-wide trait estimation error")+
  scale_color_viridis_d()+
  theme(legend.position="none")

## error rate by number of transgressive events: 
ggplot(ACE_summaries_all, aes(x=as.factor(numTrans), y=meanError, color=as.factor(numTrans))) +
  geom_boxplot()+ 
  xlab("number of transgressive events")+
  ylab("Tree-wide trait estimation error")+
  scale_color_viridis_d()+
  theme(legend.position="none")




## the means are basically the same, so tree-wide there isn't much error: 

#### do a glm: 

Char_errorDF<- ACE_summaries_all %>% select(c(PRED_vect), meanError)
Char_errorDF$nHyb<- as.factor(Char_errorDF$nHyb)
Char_errorDF$sig2_True<- as.factor(Char_errorDF$sig2_True)
Char_errorDF$HybType<- as.factor(Char_errorDF$HybType)
Char_errorDF$numTrans<- as.factor(Char_errorDF$numTrans)
Char_errorDF[is.na(Char_errorDF)]<-0

## run the glm: 

#3 run it with absolute error (general distance from correct answer) ** AKIN TO MSE
Char_errorDF$scError<- abs(Char_errorDF$meanError)
test2C<- glm(scError~.-meanError, data=Char_errorDF)
summary(test2C)
#AIC: -4651850
#Pseudo R squared: 1-(397879/555862) (.284)
#P val = 1-pchisq(534.99 -382.94 , df=(936999- 936973)) (0 = model is significantly better than null)


### RANDOM FOREST: 


rm(list=setdiff(ls(), "Char_errorDF"))
### try to do it in ranger: 
library(ranger)

rf_TreeTraitError<- ranger(scError~.-meanError, data=Char_errorDF, write.forest =F, importance="permutation", scale.permutation.importance = T)
rf_TreeTraitError$variable.importance
rf_TreeTraitError$prediction.error # .0003
rf_TreeTraitError$r.squared #50% of variance explained

### OKAY: this one is better! good to know. 
save.image("Absolute_TreeWideTraitError_RF_RangerResults.RData")

### do some plots of important variables: 
ImpVar<- data.frame(rf_TreeTraitError$variable.importance)
names(ImpVar)[1]<-"Perm_Importance"
ImpVar$Variable<- rownames(ImpVar)

ImpVar$Group<- NA
ImpVar$Group[3:15]<-"Network Structure"
ImpVar$Group[1:2]<-"Trait History"
ImpVar$Group[16:18]<-"Trait History"

ggplot(ImpVar, aes(x=reorder(Variable,-Perm_Importance), y=Perm_Importance)) +
  geom_segment( aes(x=reorder(Variable,Perm_Importance), xend=Variable, y=0, yend=Perm_Importance, color=Group)) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) + 
  xlab("Variable")+
  ylab("Coefficients")+ 
  scale_color_manual(values=c("#22A884FF","#440154FF"))

## okay, now I think we need to look at combinations of hyb type, num Trans, and sig2_True: 
comb_glm_Trait<- glm(scError~(numTrans*sig2_True)+numTrans+sig2_True, data=Char_errorDF)

Coeffs_Comb<- data.frame(summary(comb_glm_Trait)$coefficients)
names(Coeffs_Comb)[4]<- "Pval"
Coeffs_CombI<- Coeffs_Comb %>% filter(Pval < .005)

### also, run a random forest with the top two variables: 
rf_Top2_Char<- ranger(scError~numTrans+sig2_True, data=Char_errorDF, write.forest =F, importance="permutation", scale.permutation.importance = T)
rf_Top2_Char$r.squared ##34 % of variance explained! 
## okay, so these two variables alone explain 34% of the variance! So let's do a by-by plot and look at them: 

### MAKE A BOXPLOT: 
Char_TansRateInt<-ggplot(Char_errorDF, aes(x=as.factor(sig2_True), y=scError, color=as.factor(numTrans))) +
  geom_boxplot()+ 
  xlab("True rate")+
  ylab("standardized trait estimation error")+
  scale_color_viridis_d(name="number of Transgressive Events")+
  ggtitle("Transgressive-rate interaction")

Char_TansRateInt

### save this: 
save.image("TreeWide_Error_RF_Interactions_Boxplots.RData")

##### MAKE TREE-WIDE PANEL PLOTS: ###########
### do some plots of important variables: 
ImpVar<- data.frame(rf_TreeTraitError$variable.importance)
names(ImpVar)[1]<-"Perm_Importance"
ImpVar$Variable<- rownames(ImpVar)

ImpVar$Group<- NA
ImpVar$Group[3:15]<-"Network Structure"
ImpVar$Group[1:2]<-"Trait History"
ImpVar$Group[16:18]<-"Trait History"

VarImp<-ggplot(ImpVar, aes(x=reorder(Variable,-Perm_Importance), y=Perm_Importance)) +
  geom_segment( aes(x=reorder(Variable,Perm_Importance), xend=Variable, y=0, yend=Perm_Importance, color=Group)) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="none",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) + 
  xlab("Variable")+
  ylab("Permuatation Importance")+ 
  scale_color_manual(values=c("#22A884FF","#440154FF"))

NumTrans<- ggplot(Char_errorDF, aes(x=as.factor(numTrans), y=scError, color=as.factor(numTrans))) +
  geom_boxplot()+ 
  xlab("number of transgressive events")+
  ylab("tree-wide trait error")+
  scale_color_viridis_d()+
  theme(legend.position="none")

TrueRate<- ggplot(Char_errorDF, aes(x=as.factor(sig2_True), y=scError, color=as.factor(sig2_True))) +
  geom_boxplot()+ 
  xlab("True rate")+
  ylab("tree-wide trait error")+
  scale_color_viridis_d()+
  theme(legend.position="none")

HybType<- ggplot(Char_errorDF, aes(x=as.factor(HybType), y=scError, color=as.factor(HybType))) +
  geom_boxplot()+ 
  xlab("Network structure")+
  ylab("tree-wide trait error")+
  scale_color_viridis_d()+
  theme(legend.position="none")

TermMin<- ggplot(Char_errorDF, aes(x = termMin, y = scError) ) + #the basic mapping (if you want to do a scale color, you need to add color here)
  geom_point(size=.8, color="#2A788EFF") +#the geometric funtion (like is point or historgram) (you can add color here, but it won't plot a legend with the scale mapping?)
  labs(x= paste("minimum terminal branch length"), y=paste("trait estimation error")) + #adds title and axis labels
  theme(plot.title = element_text(size = 17), legend.position="none")+
  scale_color_manual(values="#2A788EFF")

IntMax<- ggplot(Char_errorDF, aes(x = intMax, y = scError) ) + #the basic mapping (if you want to do a scale color, you need to add color here)
  geom_point(size=.8, color="#414487FF") +#the geometric funtion (like is point or historgram) (you can add color here, but it won't plot a legend with the scale mapping?)
  labs(x= paste("maximum internal branch length"), y=paste("trait estimation error")) + #adds title and axis labels
  theme(plot.title = element_text(size = 17), legend.position="none")+
  scale_color_manual(values="#414487FF")

library(ggpubr)

plotlist<- list(VarImp, NumTrans, TrueRate, HybType, TermMin, IntMax)
ggarrange(plotlist=plotlist, ncol=3, nrow=2)
save.image("tree-wide.PanelPlots.RData")





######## NODE SPECIFIC TRAIT ESTIMATION ERROR: ######### ######## ####### ####### 
######
##
rm(list=ls())
gc()


Nodes_all<- read.csv("Node_specificVariables_WithTransMagnitude.csv")
Nodes_Anal<- Nodes_all %>% select(-Node_index, -Tree_index, -TrueVal,-NetTopology, 
                                  -ClosestHyb, -numHyb, -sig2_True, -numTrans, -MinNegDistH, -MinPosDistH)

## we want absolute value of error:
Nodes_Anal$scError<- abs(Nodes_Anal$error)
rm(Nodes_all)
gc()

## format your cells so they're all numeric 
Nodes_Anal<- as.data.frame(lapply(Nodes_Anal, as.numeric))


### NOW do a glm: (first one: )
test1<- glm(scError~.-error, data=Nodes_Anal)
summary(test1)
AIC: 29766537
pseudo R squared: 1-(16010346/16033729) (.0015)

test2<- glm(error~.-scError, data=Nodes_Anal)
summary(test2)
AIC: 30414139
Pseudo R squared: 1-(17261529/18742249) (.08)
## do one with a gamma distribution: 
Nodes_Gamma<- na.omit(Nodes_Anal)
rm(Nodes_Anal)
gc()

#### set up a list for making starting values: (you have to specify starting coefficients)
G1<- glm(scError~NodeAge, family=Gamma, data=Nodes_Gamma)

test3<- glm(scError~NodeAge+MinDistH+brLenTo+brLenFrom+TransMagnitude, data=Nodes_Gamma, family=Gamma(link="identity"), start=c(coef(G1), 0,0,0,0))
summary(test3)

### test if my model is better thanthe null: (IT LOOKS LIKE GAMMA FAMILY IS BAD- OR AT LEAST IDENTITY LINK)
1-pchisq(16033729 -16010346, df=(8606552- 8606547))
## pseudo R squared: 
1 - (16010346/16033729)
## .001458 ## BAD MODEL THO. I SHOULD CHECK MY OTHER ONES. 
## SO: My model performs significantly better than the null model (THE GAUSSIAN ONE NOT THE OTHER ONE.)
#### plot this shit: 
Coeffs_2<- data.frame(summary(test2)$coefficients)
Coeffs_2$Var.Names<- rownames(Coeffs_2)
Coeffs_2<- as.data.frame(Coeffs_2)
Coeffs_2<- Coeffs_2 %>% slice(-1)
names(Coeffs_2)[4]<-"Pval"

ggplot(Coeffs_2, aes(x=Var.Names, y=Estimate)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=Estimate, color=Pval)) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) + 
  xlab("Variable")+
  ylab("Coefficients")+ 
  scale_color_viridis_c()
### let's save this: 
save.image("GLM_NodeSpecific_WithTransMagnitude.RData")





### RANDOM FOREST: ############# #### 
library(randomForest)



Nodes_RF<- Nodes_Anal %>% drop_na()
rm(Nodes_Anal)
nodes_rf<- randomForest(scError~.-error, data=Nodes_RF)

### with ranger: 
rm(list=ls())
load("~/Documents/NetworkSIMtraits/NodeSpecificVars_For.RF.RData")

RF_NodeError<- ranger(scError~.-error, data=Nodes_RF, write.forest =F, importance="permutation", scale.permutation.importance = T)
RF_NodeError$variable.importance
RF_NodeError$r.squared ### .422 so 42% of variance explained
RF_NodeError$prediction.error ## error= 1.07

save.image("RF_NodeSpecificError_Finished.RData")

### do RF with top variables: 
RF_nodeError_Imp<- ranger(scError~TransMagnitude+MinDistH, data=Nodes_RF, write.forest =F, importance="permutation", scale.permutation.importance = T)
RF_nodeError_Imp$r.squared ### 25% of variance explained

### let's do a plot of important variables: 
ImpVar<- data.frame(RF_NodeError$variable.importance)
names(ImpVar)[1]<-"Perm_Importance"
ImpVar$Variable<- rownames(ImpVar)
ImpVar$Group<- NA

ggplot(ImpVar, aes(x=reorder(Variable,-Perm_Importance), y=Perm_Importance)) +
  geom_segment( aes(x=reorder(Variable,Perm_Importance), xend=Variable, y=0, yend=Perm_Importance, color=Perm_Importance)) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) + 
  xlab("Variable")+
  ylab("Permutation Importance")+ 
  scale_color_viridis_c()

### now do a glm to see if there is significant interaction between the variables: 
test<- glm(scError~TransMagnitude+MinDistH+(TransMagnitude*MinDistH), data=Nodes_RF)
summary(test)
1-(16030085/16033729) ## really bad pseudo r squared (.00002)

## save this: 
save.image("NodeSpecific_RandomForestResults.RData")
rm(list=ls())
load("NodeSpecific_RandomForestResults.RData")



#### PANEL PLOT FOR NODE SPECIFIC: 
### plot relationships in a big panel: (the only one you still need is for inherMean): 
#f) Rate and other shit: 
Node_TransMag<-ggplot(Nodes_RF, aes(x = abs(TransMagnitude), y = scError) ) + #the basic mapping (if you want to do a scale color, you need to add color here)
  geom_point(size=.8, color="#22A884FF") +#the geometric funtion (like is point or historgram) (you can add color here, but it won't plot a legend with the scale mapping?)
  labs(x= paste("magnitude of nearest transgressive event"), y=paste("standardized trait estimation error")) + #adds title and axis labels
  theme(plot.title = element_text(size = 17), legend.position="none")+
  scale_color_manual(values="#22A884FF")

Node_NearestHyb<- ggplot(Nodes_RF, aes(x = MinDistH, y = scError) ) + #the basic mapping (if you want to do a scale color, you need to add color here)
  geom_point(size=.8, color="#440154FF") +#the geometric funtion (like is point or historgram) (you can add color here, but it won't plot a legend with the scale mapping?)
  labs(x= paste("distance to nearest hybridization event"), y=paste("standardized trait estimation error")) + #adds title and axis labels
  theme(plot.title = element_text(size = 17), legend.position="none")+
  scale_color_manual(values="#440154FF")

Node_Age<-  ggplot(Nodes_RF, aes(x = NodeAge, y = scError) ) + #the basic mapping (if you want to do a scale color, you need to add color here)
  geom_point(size=.8, color="#7AD151FF") +#the geometric funtion (like is point or historgram) (you can add color here, but it won't plot a legend with the scale mapping?)
  labs(x= paste("node age"), y=paste("standardized trait estimation error")) + #adds title and axis labels
  theme(plot.title = element_text(size = 17), legend.position="none")+
  scale_color_manual(values="#7AD151FF")

Node_brTo<- ggplot(Nodes_RF, aes(x = brLenTo, y = scError) ) + #the basic mapping (if you want to do a scale color, you need to add color here)
  geom_point(size=.8, color="#414487FF") +#the geometric funtion (like is point or historgram) (you can add color here, but it won't plot a legend with the scale mapping?)
  labs(x= paste("leading branch length"), y=paste("standardized trait estimation error")) + #adds title and axis labels
  theme(plot.title = element_text(size = 17), legend.position="none")+
  scale_color_manual(values="#414487FF")

Node_brFrom<- ggplot(Nodes_RF, aes(x = brLenFrom, y = scError) ) + #the basic mapping (if you want to do a scale color, you need to add color here)
  geom_point(size=.8, color="#2A788EFF") +#the geometric funtion (like is point or historgram) (you can add color here, but it won't plot a legend with the scale mapping?)
  labs(x= paste("trailing branch length"), y=paste("standardized trait estimation error")) + #adds title and axis labels
  theme(plot.title = element_text(size = 17), legend.position="none")+
  scale_color_manual(values="#2A788EFF")

Node_Imp<- ggplot(ImpVar, aes(x=reorder(Variable,-Perm_Importance), y=Perm_Importance, color=Variable)) +
  geom_segment( aes(x=reorder(Variable,Perm_Importance), xend=Variable, y=0, 
                    yend=Perm_Importance, color=Variable)) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) + 
  xlab("Variable")+
  ylab("Permutation Importance")+ 
  scale_color_manual(values=c("#7AD151FF","#2A788EFF","#440154FF","#414487FF","#22A884FF"))
## plotting order: 
# 1: brLenFrom "#7AD151FF"
# 2: brLenTo   "#2A788EFF"
# 3: minDistH "#440154FF"
# 4: nodeAge "#414487FF"
# 5: TransMag "#22A884FF"





### PANEL THEM
plotlist<-list(Node_Imp, Node_TransMag, Node_NearestHyb, Node_brFrom, Node_brTo, Node_Age)

ggarrange(plotlist=plotlist, ncol=3, nrow=2)



