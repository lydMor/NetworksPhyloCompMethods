library(cowplot)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)
library(stringr)
library(forcats)
library(gridExtra)
library(grid)
#=====================================================================Figure 2
load("TreesVsNets_Julia.RData")

all_traits$group[all_traits$group=="Nets"]  <- "Network-based"
all_traits$group[all_traits$group=="Trees"] <- "Tree-based"
all_traits$group <- factor(all_traits$group,levels = c("Tree-based", "Network-based"))

TN_traits <- ggplot(all_traits,aes(x = as.factor(sig2_True), y = meanError, fill = group)) +
  geom_violin(position = position_dodge(width = 0.85),
              color    = "black",
              size     = 0.5,
              alpha    = 1,
              trim     = TRUE) +
  labs(x = "True evolutionary rate",
       y = "Mean tree-wide trait estimation error",
    fill = "Trait simulation scenario") +
  #stat_compare_means()+
  #scale_y_log10() +
  scale_fill_manual(values = c("white", "gray35")) +
  theme_cowplot(font_size = 13) +
  theme(axis.title = element_text(size = 15),
        axis.text  = element_text(size = 13),
        legend.position = c(0.01, 0.99),
        legend.justification = c(0, 1),
        legend.background = element_rect(
                     fill = alpha("white", 0.9),
                    color = "black",
                linewidth = 0.4),
        legend.key = element_rect(fill = NA),
        legend.margin = margin(6, 8, 6, 8))
fig2<-ggdraw(TN_traits)
fig2

#=====================================================================Figure 3
load("TreeWide_Error_RF_Interactions_Boxplots.RData")

Char_errorDF$HybType <- factor(Char_errorDF$HybType, levels = c("LG", "LG75", "LN50", "LN75", "LN"))
#### important variables
VarImpPlot <- ggplot(ImpVar, aes(x = reorder(Variable, -Perm_Importance), y = Perm_Importance)) +
  geom_segment(aes(x = reorder(Variable, Perm_Importance),
                xend = Variable,
                   y = 0,
                yend = Perm_Importance,
               color = Group),
              linewidth = 1) +
  theme_cowplot() +
  coord_flip() +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 13, family = "mono")
  ) +
  xlab("Variable") +
  ylab("Permutation importance") +
  scale_color_grey(start = 0.2, end = 0.2) +
  labs(tag = "A")

#### number of transgressive events
NumTrans <- ggplot(Char_errorDF,
                   aes(x = as.factor(numTrans),
                       y = scError,
                       fill = as.factor(numTrans))) +
  geom_violin(fill = "gray35", color = "black") +
  theme_cowplot() +
  xlab("Number of transgressive events") +
  ylab("Error") +
  theme(legend.position = "none") +
  labs(tag = "B")

#### true rate
TrueRate <- ggplot(Char_errorDF,
                   aes(x = as.factor(sig2_True),
                       y = scError,
                       fill = as.factor(sig2_True))) +
  geom_violin(fill = "gray35", color = "black") +
  theme_cowplot() +
  xlab("True rate") +
  ylab("Error") +
  theme(legend.position = "none") +
  labs(tag = "C")

#### network category
HybType <- ggplot(Char_errorDF,
                  aes(x = as.factor(HybType),
                      y = scError,
                      fill = as.factor(HybType))) +
  geom_violin(fill = "gray35", color = "black") +
  theme_cowplot() +
  xlab("Network category") +
  ylab("Error") +
  theme(legend.position = "none") +
  labs(tag = "D")

#### minimum terminal branch length
TermMin <- ggplot(Char_errorDF,
                  aes(x = termMin, y = scError)) +
  geom_point(color = "gray35", alpha = 0.35, size = 0.8) +
  theme_cowplot() +
  labs(x = "Minimum terminal branch length",
       y = "Error",
       tag = "E") +
  theme(legend.position = "none")

#### maximum internal branch length
IntMax <- ggplot(Char_errorDF,
                 aes(x = intMax, y = scError)) +
  geom_point(color = "gray35", alpha = 0.35, size = 0.8) +
  theme_cowplot() +
  labs(x = "Maximum internal branch length",
       y = "Error",
       tag = "F") +
  theme(legend.position = "none")

#### combine plots
right_top <- ggarrange(NumTrans,TrueRate,HybType, ncol = 1)
top <- ggarrange(VarImpPlot, right_top, ncol = 2, widths = c(1, 1))
bottom <- ggarrange(TermMin, IntMax, ncol = 2)
fig3 <- ggarrange(top, bottom, ncol = 1, heights = c(2, 1))
fig3


#=====================================================================Figure 4
fig4 <- ggplot(Char_errorDF,
            aes(x = as.factor(sig2_True),
                y = scError,
                fill = as.factor(numTrans))) +
  geom_boxplot(width = 0.7,
               outlier.size = 0.8,
               color = "black") +
  scale_fill_grey(start = 0.85, end = 0.3,
                  name = "Number of transgressive events") +
  scale_y_log10() +
  labs(
    x = "True rate",
    y = "Mean tree-wide trait estimation error (log-transformed)"
  ) +
  theme_cowplot(font_size = 14) +
  theme(
    legend.position = c(0.01, 0.96),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.margin = margin(2,2,2,2),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

fig4

#=====================================================================Figure 5
load("NodeSpecific_RandomForestResults.RData")

#### PANEL PLOT FOR NODE SPECIFIC: 
Node_Imp<- ggplot(ImpVar, aes(x=reorder(Variable,-Perm_Importance), y=Perm_Importance, color=Variable)) +
  geom_segment( aes(x=reorder(Variable,Perm_Importance), xend=Variable, y=0, 
                    yend=Perm_Importance, color=Variable)) +
  theme_cowplot()+
  coord_flip() +
  theme(
    legend.position="none",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) + 
  xlab("Variable")+
  ylab("Permutation importance")+ 
  scale_color_grey(start = 0.2, end = 0.2) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 13, family = "mono")  # only text inside plot
  ) +
  labs(tag="A")

Node_TransMag<-ggplot(Nodes_RF, aes(x = abs(TransMagnitude), y = scError) ) + #the basic mapping (if you want to do a scale color, you need to add color here)
  geom_point(color = "gray35", alpha = 0.35, size = 0.8) +
  labs(x= paste("Magnitude of nearest transgressive event"), y=paste("Error")) + #adds title and axis labels
  theme_cowplot()+
  theme(plot.title = element_text(size = 17), legend.position="none")+
  labs(tag="B")

Node_NearestHyb<- ggplot(Nodes_RF, aes(x = MinDistH, y = scError) ) + #the basic mapping (if you want to do a scale color, you need to add color here)
  geom_point(color = "gray35", alpha = 0.35, size = 0.8) +
  labs(x= paste("Distance to nearest hybridization event"), y=paste("Error")) + #adds title and axis labels
  theme_cowplot()+
  theme(plot.title = element_text(size = 17), legend.position="none")+
  labs(tag="C")

Node_brFrom<- ggplot(Nodes_RF, aes(x = brLenFrom, y = scError) ) + #the basic mapping (if you want to do a scale color, you need to add color here)
  geom_point(color = "gray35", alpha = 0.35, size = 0.8) +
  labs(x= paste("Trailing branch length"), y=paste("Error")) + #adds title and axis labels
  theme_cowplot()+
  theme(plot.title = element_text(size = 17), legend.position="none")+
  labs(tag="D")

Node_brTo<- ggplot(Nodes_RF, aes(x = brLenTo, y = scError) ) + #the basic mapping (if you want to do a scale color, you need to add color here)
  geom_point(color = "gray35", alpha = 0.35, size = 0.8) +
  labs(x= paste("Leading branch length"), y=paste("Error")) + #adds title and axis labels
  theme_cowplot()+
  theme(plot.title = element_text(size = 17), legend.position="none")+
  labs(tag="E")

Node_Age<-  ggplot(Nodes_RF, aes(x = NodeAge, y = scError) ) + #the basic mapping (if you want to do a scale color, you need to add color here)
  geom_point(color = "gray35", alpha = 0.35, size = 0.8) +
  labs(x= paste("Node depth"), y=paste("Error")) + #adds title and axis labels
  theme_cowplot()+
  theme(plot.title = element_text(size = 17), legend.position="none")+
  labs(tag="F")

### PANEL THEM
plotlist=list(Node_Imp, Node_TransMag, Node_NearestHyb, Node_brFrom, Node_brTo, Node_Age)

fig5<-ggarrange(plotlist=plotlist, ncol=3, nrow=2)
fig5

#=====================================================================Figure 6
#load("TreesVsNets_Julia.RData")  #data file loaded for fig1
all_rates$group[all_rates$group=="Nets"]<-"Network-based"
all_rates$group[all_rates$group=="Trees"]<-"Tree-based"
all_rates$group <- factor(all_rates$group,levels = c("Tree-based", "Network-based"))

TN_rates <- ggplot(all_rates, aes(x = as.factor(rate),y = rate_error,fill = as.factor(group))) +
  geom_violin() + 
  theme_light() +
  xlab("True rate") +
  ylab("Rate estimation error") +
  theme(legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  scale_fill_manual(
    values = c("white", "gray35"),
    name = "Trait simulation scenario") +
  #stat_compare_means() +
  #coord_cartesian(ylim = c(-5, 5))+
  #scale_y_log10()+ 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.7) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10))+
  theme_cowplot(font_size = 13) +
  theme(
    axis.title = element_text(size = 15),
    axis.text  = element_text(size = 13),
    
    legend.position = c(0.01, 0.99),
    legend.justification = c(0, 1),
    
    legend.background = element_rect(
      fill = alpha("white", 0.9),
      color = "black",
      linewidth = 0.4
    ),
    
    legend.key = element_rect(fill = NA),
    
    legend.margin = margin(6, 8, 6, 8)   # top, right, bottom, left
  )

fig6<-ggdraw(TN_rates)
fig6

#=====================================================================Figure 7
load("Rate_Error_interactions_RandomForest_Plots.Rdata")

bm_rates$HybType<- factor(bm_rates$HybType, levels=c("LG", "LG75", "LN50", "LN75", "LN"))
#### rate factors: 
ImpVarPlot <- ggplot(ImpVar, aes(x = reorder(Variable, -Perm_Importance),y = Perm_Importance)) +
  geom_segment(aes(x = reorder(Variable, Perm_Importance),
                xend = Variable,
                   y = 0,
                yend = Perm_Importance,
               color = Group),
              linewidth = 1) +
  theme_cowplot() +
  coord_flip() +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 13, family = "mono")  # only text inside plot
  ) +
  xlab("Variable") +
  ylab("Permutation importance") +
  scale_color_grey(start = 0.2, end = 0.2) +
  scale_linetype_manual(values = c("solid", "dashed"))+
  labs(tag = "A")

## error rate by number of transgressive events: 
netRerror_NumTrans<-ggplot(bm_rates, aes(x=as.factor(numTrans), y=lrate_error, fill=as.factor(numTrans))) +
  geom_violin(fill = "gray35", color = "black")+ 
  theme_cowplot()+
  xlab("Number of transgressive events")+
  ylab("Error")+
  theme(legend.position="none")+
  labs(tag="B")

## error by true rate: 
netRerror_TrueRate<-ggplot(bm_rates, aes(x=as.factor(sig2_True), y=lrate_error, fill=as.factor(sig2_True))) +
  geom_violin(fill = "gray35", color = "black")+ 
  theme_cowplot()+
  xlab("True rate")+
  ylab("Error")+
  theme(legend.position="none")+
  labs(tag="C")

## error by hyb type: 
netRerror_HybType<- ggplot(bm_rates, aes(x=as.factor(HybType), y=lrate_error, fill=as.factor(HybType))) +
  geom_violin(fill = "gray35", color = "black")+ 
  theme_cowplot()+
  xlab("Network category")+
  ylab("Error")+
  theme(legend.position="none")+
  labs(tag="D")

### mean inheritance: 
RateError_MeanInheritance<-ggplot(bm_rates, aes(x = inherMean, y = lrate_error)) + #the basic mapping (if you want to do a scale color, you need to add color here)
  geom_point(color = "gray35", alpha = 0.35, size = 0.8) +
  theme_cowplot()+
  labs(x= paste("Mean parental inheritance"), 
       y=paste("Error")) + #adds title and axis labels
  scale_color_manual(values="gray35")+
  labs(tag="E")

### INT MAX: 
RateError_intMax <- ggplot(bm_rates, aes(x = intMax, y = lrate_error)) +
  geom_point(color = "gray35", alpha = 0.35, size = 0.8) +
  theme_cowplot() +
  labs(x = "Maximum internal branch length",
       y = "Error",
       tag = "F") +
  theme(plot.title = element_text(size = 17),
        legend.position = "none")

right_top <- ggarrange(netRerror_NumTrans, netRerror_TrueRate, netRerror_HybType, ncol = 1)
top <- ggarrange(ImpVarPlot, right_top, ncol = 2, widths = c(1,1))
bottom <- ggarrange(RateError_MeanInheritance, RateError_intMax, ncol = 2)
fig7 <- ggarrange(top, bottom, ncol = 1, heights = c(2,1))
fig7

#=====================================================================Figure 8
#load("Rate_Error_interactions_RandomForest_Plots.Rdata")  #data file loaded for fig7
Rate_TansRateInt <- ggplot(Rate_errorDF,
                           aes(x = as.factor(sig2_True),
                               y = lrate_error,
                               fill = as.factor(numTrans))) +
  geom_boxplot() +
  theme_light() +
  xlab("True rate") +
  ylab("Rate estimation error") +
  scale_fill_grey(start = 0.4, end = 0.9,
                  name = "Number of transgressive events") +
  geom_hline(yintercept = 2.99, linetype = "dashed") +
  theme_cowplot() +
  theme(
    legend.position = c(0.01, 0.99),
    legend.justification = c(0, 1),
    legend.direction = "horizontal",
    legend.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.4
    ),
    legend.margin = margin(2,2,2,2)  # increase legend padding
  )

fig8<-Rate_TansRateInt
fig8

#=====================================================================Figure 9
load("MCLUST_results.RData")

# make a line plot to go with the rest of your permuatation importance variables: 
PERM_means$scDec<- PERM_means$DecreaseLogLik/max(PERM_means$DecreaseLogLik)
PERM_USE<- PERM_means %>% filter(!(var=="inherMean" |var=="inherVar"))

VarImp_Mod <- ggplot(PERM_USE, aes(x = reorder(var, -scDec), y = scDec)) +
  geom_segment(
    aes(x = reorder(var, scDec),
        xend = var,
        y = 0,
        yend = scDec),
    color = "gray35",
    linewidth = 1
  ) +
  theme_cowplot() +
  coord_flip() +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(family = "mono", size = 12)   # texttt style
  ) +
  xlab("Variable") +
  ylab("Permutation importance") +
  labs(tag = "A")

ModMeans<- ACE_summaries_all %>% select(ModelUsed, inherMax, nHyb, numTrans, PropExtH, PropTermH)
ModMeans$num<- paste(1:nrow(ModMeans))
ModMeans$Model<- paste(ModMeans$ModelUsed, ModMeans$num, sep="_")
rownames(ModMeans)<- ModMeans$Model
ModMeans<- ModMeans %>% select(-num, -Model, -ModelUsed)
test<- melt(as.matrix(ModMeans))
test$ModelUsed<- paste(word(test$Var1, 1, sep="_"))

means<- test %>% group_by(ModelUsed, Var2) %>% summarize(mean=mean(value), sd=sd(value))
means$normMean<- means$mean/means$sd

wilcox.test(
  value ~ ModelUsed,
  data = test %>% filter(Var2 == "inherMax"),
  alternative = "two.sided"
)

t.test(
  value ~ ModelUsed,
  data = test %>% filter(Var2 == "numTrans"),
  alternative = "two.sided"
)

#wilcox.test(BM$inherMax, OU$inherMax, alternative = "two.sided")
#t.test(BM$numTrans, OU$numTrans, alternative= "two.sided")

### show differences in means: 
BM_m<- means %>% filter(ModelUsed=="BM")
BM_m<- BM_m %>% ungroup() %>% select(Var2, normMean)
names(BM_m)[2]<- "mean_BM"

OU_m<- means %>% filter(ModelUsed=="OU")
OU_m<- OU_m %>% ungroup()%>% select(Var2, normMean)
names(OU_m)[2]<- "mean_OU"

means_diff<- full_join(BM_m, OU_m)
means_diff$diff<- means_diff$mean_OU - means_diff$mean_BM

# plot this: 
SMD_model <- ggplot(means_diff, aes(x = Var2, y = diff)) +
  geom_bar(stat = "identity", fill = "gray35", color = "black") +
  xlab("Variable") +
  ylab("Standardized mean difference") +
  theme_cowplot() +
  theme(
    axis.title = element_text(size = 14),
    axis.text.x = element_text(family = "mono", size = 12),
    legend.position = "none"
  ) +
  labs(tag = "B")

plotlist<- list(VarImp_Mod, SMD_model)
fig9<-ggarrange(plotlist=plotlist, ncol=2, nrow=1)
fig9

#=====================================================================Figure 10
load("TreesVsNets_FULLDATASET_inREAL_ANALYSIS_upTo237.RData")

##### let's get it in a numeric situation: 
Mods_tests<- ACE_summaries_all %>% select(ModelUsed, HybType, sig2_True, inherMean, numTrans, nHyb, inherMax)
Mods_tests$Mod_Fact<- as.integer(as.factor(Mods_tests$ModelUsed))
Mods_tests<- Mods_tests %>% mutate(inher_brk = fct_rev(cut(inherMax, seq(0,ceiling(max(inherMax)),.2))))
#inheritance probs: 
inher_props <- ggplot(data = Mods_tests, aes(x = inher_brk, fill = ModelUsed)) +
  geom_bar(position = "fill") +
  theme_cowplot() +
  xlab("Maximum inheritance probabilities") +
  ylab("Proportion") +
  scale_fill_manual(values = c("gray70", "gray30")) +
  labs(fill="Evolutionary model used", tag="A")+
  theme(legend.background = element_rect(fill = alpha("white", 0.9), color = "black", linewidth = 0.4),
    legend.margin = margin(5,5,5,5))

#numHyb: 
numHyb_props<- ggplot(data=Mods_tests, aes(x=nHyb, fill=ModelUsed))+
  geom_bar(position="fill")+
  theme_cowplot() +
  xlab("Number of hybridization events")+
  ylab("Proportion")+
  scale_fill_manual(values = c("gray70", "gray30")) +
  labs(tag="B")

#numTrans: 
numTrans_props<- ggplot(data=Mods_tests, aes(x=numTrans, fill=ModelUsed))+
  geom_bar(position="fill")+
  theme_cowplot() +
  xlab("Number of transgressive events")+
  ylab("Proportion")+
  scale_fill_manual(values = c("gray70", "gray30")) +
  labs(tag="C")

#HybType: 
HybType_props<- ggplot(data=Mods_tests, aes(x=HybType, fill=ModelUsed))+
  geom_bar(position="fill")+
  theme_cowplot() +
  xlab("Network category")+
  ylab("Proportion")+
  scale_fill_manual(values = c("gray70", "gray30")) +
  labs(tag="D")

###### function for sharing legend: 

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "top", "right", "left")
) {
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     
                     "bottom" = arrangeGrob(
                       do.call(arrangeGrob, gl),
                       legend,
                       ncol = 1,
                       heights = unit.c(unit(1, "npc") - lheight, lheight)
                     ),
                     
                     "top" = arrangeGrob(
                       legend,
                       do.call(arrangeGrob, gl),
                       ncol = 1,
                       heights = unit.c(lheight, unit(1, "npc") - lheight)
                     ),
                     
                     "right" = arrangeGrob(
                       do.call(arrangeGrob, gl),
                       legend,
                       ncol = 2,
                       widths = unit.c(unit(1, "npc") - lwidth, lwidth)
                     ),
                     
                     "left" = arrangeGrob(
                       legend,
                       do.call(arrangeGrob, gl),
                       ncol = 2,
                       widths = unit.c(lwidth, unit(1, "npc") - lwidth)
                     )
  )
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}

fig10 <- grid_arrange_shared_legend(inher_props,
                                 numHyb_props,
                                 numTrans_props,
                                 HybType_props,
                                 ncol=2,nrow=2,
                                 position="top")



fig10
