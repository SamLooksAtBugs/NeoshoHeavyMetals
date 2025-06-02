# Me First ####
setwd("C:/Users/Sam/OneDrive/Elements24")
library(vegan)
library(tidyverse)
library(readxl)
library(biotools)
library(glmmTMB)
library(ggeffects)
library(PMCMRplus)
library(rrcov)

# These should be available in the Github repository. 
boogs <- read_excel("ICPoutput.xlsx", sheet= "bugs")
wtr <- read_excel("ICPoutput.xlsx", sheet= "water")
info <- read_excel("elementalbugs.xlsx", sheet= "Sheet1")
waterelements <- wtr[,3:ncol(wtr)]
wtr[,3:ncol(wtr)] <- waterelements/1000
data <- rbind(boogs,wtr)

data <- left_join(info,data, by="Sample")


head(data)
View(data)

#PCA 1: lower bugs vs upper bugs  ####
bug.sub <- data %>% filter(Taxa %in% c("Hexagenia","Naididae","Orthocladiinae")& 
                             !Site %in% c("Elk creek","Hickory creek","Honey creek"))

set.seed(123)
rpca <- PcaHubert(bug.sub[,7:ncol(bug.sub)],scale = T,k=4)

pca <- data.frame(Taxa=bug.sub$Taxa,rpca@scores,Site=bug.sub$Site)
loadings <- data.frame(rpca@loadings)
loadings$Element <- row.names(loadings)
rownames(loadings) <- NULL

ggplot()+stat_ellipse(data=pca,aes(x=PC1,y=PC4,group=interaction(Taxa,Site),color=Taxa,linetype = Site), level=0.68)+
  geom_point(data=pca,aes(x=PC1,y=PC4,color=Taxa,shape=Site),size=3,alpha=0.8)+
  theme(legend.title = element_blank(),legend.key.size = unit(0.5,units = "cm"))+
  geom_segment(data=loadings,aes(x = 0, y = 0, xend = (7*PC1), yend = (7*PC4), label=Element),
               arrow = arrow(length = unit(0.2, "cm")),  # Arrowhead size
               linewidth = 0.8,  # Adjust thickness here
               color = "black")  +# Optional: Customize color+
  ggrepel::geom_label_repel(data=loadings,aes(x = (7*PC1), y= (7*PC4), label=Element ),size=3, 
                            max.overlaps = 100, min.segment.length = 4, label.size = 0.3,alpha=0.8)+
  
  scale_color_brewer(palette = "Dark2")+scale_fill_brewer(palette = "Dark2")+theme_classic()+
  facet_wrap(~Taxa)+
  theme(legend.title = element_blank(),legend.key.size = unit(0.5,units = "cm"), text=element_text(size=16),legend.position = "none")+
  labs(color = "Taxa", shape = "Site", linetype="Site", x= "PC1 (56.2%)", y="PC4 (6.4%)")


#PCA 2: Hexagenia ####
hex.sub <- data %>% filter(Taxa=="Hexagenia")
hex.sub <- hex.sub %>% mutate(Site=case_when(Site %in% c("Elk creek","Hickory creek","Honey creek")~"Grand Lake",
                                             !Site %in% c("Elk creek","Hickory creek","Honey creek")~Site))
# hex.sub <- hex.sub %>% mutate(Site= case_when(Sample %in% c("Hex1p"  ,  "Hex2p"  ,  "Hex3p"  ,  "Hex4p"  ,  "Hex5p" )~"Lower Neosho (past)",
#                                               !Sample %in% c("Hex1p"  ,  "Hex2p"  ,  "Hex3p"  ,  "Hex4p"  ,  "Hex5p" )~Site))

rpca <- PcaHubert( hex.sub[,7:ncol( hex.sub )],scale = T)

pca <- data.frame(Taxa= hex.sub$Taxa,rpca@scores,Site= hex.sub$Site)
loadings <- data.frame(rpca@loadings)
loadings$Element <- row.names(loadings)
rownames(loadings) <- NULL

ggplot()+stat_ellipse(data=pca,aes(x=PC1,y=PC2,color=Site), level=0.68)+
  geom_point(data=pca,aes(x=PC1,y=PC2,color=Site, shape=Site),size=3,alpha=0.8)+
  theme_classic()+
  theme(legend.title = element_blank(),legend.key.size = unit(0.5,units = "cm"))+
  geom_segment(data=loadings,aes(x = 0, y = 0, xend = (5*PC1), yend = (5*PC2), label=Element),
               arrow = arrow(length = unit(0.2, "cm")),  # Arrowhead size
               linewidth = 0.8,  # Adjust thickness here
               color = "black")  +# Optional: Customize color+
  ggrepel::geom_label_repel(data=loadings,aes(x = (5*PC1), y= (5*PC2), label=Element ),size=3, 
                            max.overlaps = 100, min.segment.length = 4, label.size = 0.3,alpha=0.8)+
  
  scale_color_manual(values=c("gray60","forestgreen","red4"))+scale_fill_brewer(palette = "Dark2")+
  theme(legend.title = element_blank(),legend.key.size = unit(0.5,units = "cm"), text=element_text(size=16))+
  labs(color = "Taxa", shape = "Taxa", fill="Taxa", x= "PC1 (65.6%)", y="PC2 (12.4%)")

#Boxplot 1: Upper vs Lower Neosho ####

bug.sub <- data %>% filter(Taxa %in% c("Hexagenia","Naididae","Orthocladiinae","Sediment","Unfiltered water","Filtered water")& 
                             !Site %in% c("Elk creek","Hickory creek","Honey creek")&
                             !Sample %in% c("Hex1p","Hex2p","Hex3p","Hex4p","Hex5p" ))
bug.sub$Taxa <- factor(bug.sub$Taxa, levels=c("Hexagenia","Naididae","Orthocladiinae","Filtered water","Unfiltered water","Sediment"))

#separate plots with boxes
zinc <- ggplot(data=bug.sub,aes(x=Taxa,y=Zn,group=interaction(Taxa,Site),color=Site))+geom_boxplot()+geom_vline(aes(xintercept = 3.5), linetype="dashed")+theme_bw()+xlab("")+theme(legend.position = "none",axis.text.x = element_text(angle=45, hjust=1)) + scale_color_manual(values = c("forestgreen","red4"))
copper <- ggplot(data=bug.sub,aes(x=Taxa,y=Cu,group=interaction(Taxa,Site),color=Site))+geom_boxplot()+geom_vline(aes(xintercept = 3.5), linetype="dashed")+theme_bw()+xlab("")+theme(legend.position = "none",axis.text.x = element_text(angle=45, hjust=1))+ scale_color_manual(values = c("forestgreen","red4"))
chromium <- ggplot(data=bug.sub,aes(x=Taxa,y=Cr,group=interaction(Taxa,Site),color=Site))+geom_boxplot()+geom_vline(aes(xintercept = 3.5), linetype="dashed")+theme_bw()+xlab("")+theme(legend.position = "none",axis.text.x = element_text(angle=45, hjust=1))+ scale_color_manual(values = c("forestgreen","red4"))
manganese <- ggplot(data=bug.sub,aes(x=Taxa,y=Mn,group=interaction(Taxa,Site),color=Site))+geom_boxplot()+geom_vline(aes(xintercept = 3.5), linetype="dashed")+theme_bw()+xlab("")+theme(legend.position = "none",axis.text.x = element_text(angle=45, hjust=1))+ scale_color_manual(values = c("forestgreen","red4"))
nickel <- ggplot(data=bug.sub,aes(x=Taxa,y=Ni,group=interaction(Taxa,Site),color=Site))+geom_boxplot()+geom_vline(aes(xintercept = 3.5), linetype="dashed")+theme_bw()+ylim(0,0.025)+xlab("")+theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))+ scale_color_manual(values = c("forestgreen","red4"))
stront <- ggplot(data=bug.sub,aes(x=Taxa,y=Sr,group=interaction(Taxa,Site),color=Site))+geom_boxplot()+geom_vline(aes(xintercept = 3.5), linetype="dashed")+theme_bw()+xlab("")+theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))+ scale_color_manual(values = c("forestgreen","red4"))

gridExtra::grid.arrange(zinc, nickel,copper,manganese,nrow=2)                     


#faceted plot w/ brackets
bug.long <- bug.sub %>%
  pivot_longer(cols = c(Cu, Ni, Mn, Zn), names_to = "Metal", values_to = "Concentration")

#change concentration to mg/kg
bug.long$Concentration <- bug.long$Concentration *1000
# Step 2: Compute summary stats for each group
bug.summary <- bug.long %>%
  group_by(Taxa, Site, Metal) %>%
  summarise(
    ymin = min(Concentration, na.rm = TRUE),
    ymax = max(Concentration, na.rm = TRUE),
    ymen = mean(Concentration, na.rm = TRUE),
    ymed = median(Concentration, na.rm = TRUE),
    .groups = "drop"
  )

bug.summary$Metal <- factor(bug.summary$Metal, levels=c("Zn", "Ni","Cu","Mn"))

# Step 3: Plot with range bars, whiskers, and median points, faceted by metal
ggplot(bug.summary, aes(x = Taxa, group = interaction(Taxa, Site), color = Site)) +
  # Range bar (min to max)
  # geom_linerange(aes(ymin = ymin, ymax = ymax), position = position_dodge(width = 0.6), size = 1) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), 
                width = 0.5, 
                position = position_dodge(width = 0.6), 
                size = 1)+# Median point
  geom_point(aes(y = ymed), position = position_dodge(width = 0.6), size = 2) +
  # Vertical line (same place in every facet)
  #geom_hline(aes(yintercept = 320), linetype = "dashed") +
  facet_wrap(~Metal, scales = "free_y") +
  theme_classic() +
  xlab("") +
  ylab("Concentration (mg/kg)")+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_color_manual(values = c("forestgreen", "red4"))
#Boxplot 3: Hexagenia ####
hex.sub <- data %>% filter(Taxa=="Hexagenia")
hex.sub <- hex.sub %>% mutate(Site=case_when(Site %in% c("Elk creek","Hickory creek","Honey creek")~"Grand Lake",
                                             !Site %in% c("Elk creek","Hickory creek","Honey creek")~Site))
hex.sub <- hex.sub %>% mutate(Site=case_when(Sample %in% c("Hex1p","Hex2p","Hex3p","Hex4p","Hex5p" )~"Lower Neosho '22",
                                             !Site %in% c("Hex1p","Hex2p","Hex3p","Hex4p","Hex5p" )~Site))
hex.sub$Site <- factor(hex.sub$Site, levels = c("Grand Lake", "Lower Neosho '22","Lower Neosho","Upper Neosho"))


zinc <- ggplot(data=hex.sub)+geom_boxplot(aes(x=Site,y=Zn,group=Site,color=Site))+theme_bw()+xlab("")+theme(legend.position = "none",axis.text.x = element_text(angle=45, hjust=1,size=10))+scale_color_manual(values = c("gray20","lawngreen","forestgreen","red4"))
copper <- ggplot(data=hex.sub)+geom_boxplot(aes(x=Site,y=Cu,group=Site,color=Site))+theme_bw()+xlab("")+theme(legend.position = "none",axis.text.x = element_text(angle=45, hjust=1,size=10))+scale_color_manual(values = c("gray20","lawngreen","forestgreen","red4"))
chromium <- ggplot(data=hex.sub)+geom_boxplot(aes(x=Site,y=Cr,group=Site,color=Site))+theme_bw()+xlab("")+theme(legend.position = "none",axis.text.x = element_text(angle=45, hjust=1,size=10))+scale_color_manual(values = c("gray20","lawngreen","forestgreen","red4"))
manganese <- ggplot(data=hex.sub)+geom_boxplot(aes(x=Site,y=Mn,group=Site,color=Site))+theme_bw()+xlab("")+theme(legend.position = "none",axis.text.x = element_text(angle=45, hjust=1,size=10))+scale_color_manual(values = c("gray20","lawngreen","forestgreen","red4"))
nickel <- ggplot(data=hex.sub)+geom_boxplot(aes(x=Site,y=Ni,group=Site,color=Site))+theme_bw()+xlab("")+theme(legend.position = "none",axis.text.x = element_text(angle=45, hjust=1,size=10))+scale_color_manual(values = c("gray20","lawngreen","forestgreen","red4"))
iron <- ggplot(data=hex.sub)+geom_boxplot(aes(x=Site,y=Fe,group=Site,color=Site))+theme_bw()+xlab("")+theme(legend.position = "none",axis.text.x = element_text(angle=45, hjust=1,size=10))+scale_color_manual(values = c("gray20","lawngreen","forestgreen","red4"))

gridExtra::grid.arrange(zinc, nickel, copper, manganese) 


#faceted plot w/ brackets
hex.long <- hex.sub %>%
  pivot_longer(cols = c(Cu, Ni, Mn, Zn), names_to = "Metal", values_to = "Concentration")

#change concentration to mg/kg
hex.long$Concentration <- hex.long$Concentration *1000

hex.long$Metal <- factor(hex.long$Metal, levels=c("Zn", "Ni","Cu","Mn"))
hex.long$Site <- factor(hex.long$Site, levels=c("Grand Lake", "Lower Neosho '22", "Lower Neosho","Upper Neosho"))
ggplot(data=hex.long)+
  geom_boxplot(aes(x=Site,y=Concentration,group=Site,color=Site))+
  theme_classic()+
  facet_wrap(~Metal, scales = "free_y")+
  xlab("")+ylab("Concentration (mg/kg)")+
  theme(legend.position = "none",axis.text.x = element_text(angle=45, hjust=1,size=10))+scale_color_manual(values = c("gray20","lawngreen","forestgreen","red4"))

k<- kruskal.test(Zn~Site, data=hex.sub)
dun<- dunn.test::dunn.test(hex.sub$Zn,hex.sub$Site)

#Modeling Upper vs. Lower ##### 
bug.sub <- data %>% filter(Taxa %in% c("Hexagenia","Naididae","Orthocladiinae","Tanypodinae","Ceratopoginae")& 
                             !Site %in% c("Elk creek","Hickory creek","Honey creek")&
                             Date > "2023-01-01 12:00:00" )

#remove nickel outliers
bug.sub <- bug.sub %>% filter(Sample != "Tany2u"&Sample != "Naid8u")

bugp <- list()
bugt <- list()
bugse <- list()
bugc <- list()
models <- list()
for (metal in c("Cu","Mn","Ni","Zn")){
  
  #concentration to mg/kg
  bug.sub[[metal]] <- bug.sub[[metal]]*1000
  model <- glmmTMB(bug.sub[[metal]]~1+Site+
                      (1|Taxa),ziformula = ~0,
                    control=glmmTMBControl(optimizer=optim,
                                           optCtrl=list(maxit= 1000),
                                           optArgs=list(method="BFGS")),
                    family=Gamma(link="log"),data=bug.sub)
  
  modelsum <- summary(model)
  bugp[[metal]] <- modelsum$coefficients$cond[2,4]
  bugt[[metal]] <- modelsum$coefficients$cond[2,3]
  bugse[[metal]] <- modelsum$coefficients$cond[2,2]
  bugc[[metal]] <- modelsum$coefficients$cond[2,1]
  
  modeldf <- ggpredict(model,terms=c("Site","Taxa"),type = "random")
  models[[metal]] <- data.frame(Element=metal,Site=modeldf$x,Estimate=modeldf$predicted,
                              se = ifelse(is.null(modeldf$std.error),0,modeldf$std.error),Taxa=modeldf$group,
                              p=modelsum$coefficients$cond[2,4])  # upper confidence interval)
  
  print(metal)
  print(modelsum)
  }

ps <- do.call(rbind,bugp)
ts <- do.call(rbind, bugt)
ses <- do.call(rbind,bugse)
cs <- do.call(rbind, bugc)

modelsummary <- data.frame(Element= rownames(ts),Coefficient = cs, Std.Error=ses,
                           Zstat=ts,Pvalue=ps)
rownames(modelsummary) <- NULL

modelplot <- do.call(rbind, models)

bug.sub <- data %>% filter(Taxa %in% c("Hexagenia","Naididae","Orthocladiinae","Ceratopoginae","Tanypodinae")& 
                             !Site %in% c("Elk creek","Hickory creek","Honey creek")&
                             Date > "2023-01-01 12:00:00" )
#remove nickel outliers
bug.sub <- bug.sub %>% filter(Sample != "Tany2u"& Sample != "Naid8u")

low2 <- bug.sub[,c("Taxa","Site","Cu","Cr","Ni","Zn","Mn","Fe")]
lowlong <- pivot_longer(low2, cols=3:ncol(low2),names_to = "Element",values_to = "Estimate")

lowlong$Estimate <- lowlong$Estimate *1000

plotted <- list()
for (metal in c("Cr", "Cu","Mn","Ni","Zn","Fe")){
  df1 <- lowlong %>% filter(Element==metal)
  df2 <- modelplot %>% filter(Element==metal)

  plotted[[metal]]<- ggplot() +
  # Raw data points
  geom_jitter(data = df1, aes(x = Site, y = Estimate, group =Site, color = Taxa),
             size=2, width = 0.1, height = 0,alpha=0.7) +
  # Model output as lines
  geom_line(data = df2, aes(x = Site, y = Estimate, group = Taxa, color = Taxa), size=0.8) +
  # A clean theme
  theme_classic()+theme(axis.text.x = element_text(hjust=1,angle=45),text=element_text(size=16),legend.position = "none")+labs(y=metal, x="",color="Taxa",group="Taxa")+scale_color_brewer(palette = "Dark2")
    
}

gridExtra::grid.arrange(plotted$Zn,plotted$Ni,plotted$Cu,plotted$Mn)


#modeling enviro ####
env.sub <- data %>% filter(Taxa %in% c("Sediment", "Filtered water", "Unfiltered water")& 
                             !Site %in% c("Elk creek","Hickory creek","Honey creek")&
                             !Sample %in% c("Hex1p","Hex2p","Hex3p","Hex4p","Hex5p" ))

env.sub$Site <- factor(env.sub$Site, levels=c("Lower Neosho","Upper Neosho"))

envp <- list()
envt <- list()
envse <- list()
envc <- list()
models <- list()
for (metal in c("Cr", "Cu","Mn","Ni","Zn","Fe")){
  
  #concentration to mg/kg
  env.sub[[metal]] <- env.sub[[metal]]*1000
  
  model <- glmmTMB(env.sub[[metal]]~Site+
                     (1|Taxa),ziformula = ~0,
                   control=glmmTMBControl(optimizer=optim,
                                          optCtrl=list(maxit= 1000),
                                          optArgs=list(method="BFGS")),
                   family=Gamma(link="log"),data=env.sub)
  
  modelsum <- summary(model)
  envp[[metal]] <- modelsum$coefficients$cond[2,4]
  envt[[metal]] <- modelsum$coefficients$cond[2,3]
  envse[[metal]] <- modelsum$coefficients$cond[2,2]
  envc[[metal]] <- modelsum$coefficients$cond[2,1]
  
  modeldf <- ggpredict(model,terms=c("Site","Taxa"),type = "random")
  models[[metal]] <- data.frame(Element=metal,Site=modeldf$x,Estimate=modeldf$predicted,
                                se = ifelse(is.null(modeldf$std.error),0,modeldf$std.error),Taxa=modeldf$group,
                                p=modelsum$coefficients$cond[2,4])  # upper confidence interval)
  
}

ps <- do.call(rbind,envp)
ts <- do.call(rbind, envt)
ses <- do.call(rbind,envse)
cs <- do.call(rbind, envc)

modelsummary <- data.frame(Element= rownames(ts),Coefficient = cs, Std.Error=ses,
                           Zstat=ts,Pvalue=ps)
rownames(modelsummary) <- NULL

modelplot <- do.call(rbind, models)

env.sub <- data %>% filter(Taxa  %in% c("Sediment", "Filtered water", "Unfiltered water")& 
                             !Site %in% c("Elk creek","Hickory creek","Honey creek")&
                             Date > "2023-01-01 12:00:00" )
#remove nickel outliers
env.sub <- env.sub %>% filter(Sample != "Tany2u"& Sample != "Naid8u")

env.sub$Site <- factor(env.sub$Site, levels=c("Lower Neosho","Upper Neosho"))

low2 <- env.sub[,c("Taxa","Site","Cu","Cr","Ni","Zn","Mn","Fe")]
lowlong <- pivot_longer(low2, cols=3:ncol(low2),names_to = "Element",values_to = "Estimate")

lowlong$Estimate <- lowlong$Estimate *1000

plotted <- list()
for (metal in c("Cr", "Cu","Mn","Ni","Zn","Fe")){
  df1 <- lowlong %>% filter(Element==metal)
  df2 <- modelplot %>% filter(Element==metal)
  
  plotted[[metal]]<- ggplot() +
    # Raw data points
    geom_jitter(data = df1, aes(x = Site, y = Estimate, group =Site, color = Taxa),
                size=2, width = 0.1, height = 0,alpha=0.7) +
    # Model output as lines
    geom_line(data = df2, aes(x = Site, y = Estimate, group = Taxa, color = Taxa), size=0.8) +
    # A clean theme
    theme_classic()+theme(axis.text.x = element_text(hjust=1,angle=45),
                          text=element_text(size=16),
                          legend.position = "none")+
    labs(y=metal, x="",color="Taxa",group="Taxa")+
    scale_color_brewer(palette = "Dark2")
  
}

gridExtra::grid.arrange(plotted$Zn,plotted$Ni,plotted$Cu,plotted$Mn)
#Games Howell randomizations with Hexagenia ####
hex.sub <- data %>% filter(Taxa=="Hexagenia")
hex.sub <- hex.sub %>% mutate(Site=case_when(Site %in% c("Elk creek","Hickory creek","Honey creek")~"Grand Lake",
                                             !Site %in% c("Elk creek","Hickory creek","Honey creek")~Site))
hex.sub <- hex.sub %>% mutate(Site=case_when(Sample %in% c("Hex1p","Hex2p","Hex3p","Hex4p","Hex5p" )~"Lower Neosho (past)",
                                             !Site %in% c("Hex1p","Hex2p","Hex3p","Hex4p","Hex5p" )~Site))
hex.sub$Site <- factor(hex.sub$Site, levels = c("Grand Lake", "Lower Neosho (past)","Lower Neosho","Upper Neosho"))


#check normality #everything is normal except for Cr, so probably lets switch to parametric
hex.shaps <- list()
for (x in unique(hex.sub$Site)) {
  df <- hex.sub %>% filter(Site == x)
  metal.stat <- list()
  for (metal in  c("Zn","Cu","Cr","Mn","Ni")) {
    dff <- shapiro.test(df[[metal]])
    metal.stat[[metal]] <- data.frame(Site=x,Metal=metal,pvalue=dff$p.value)
  }
  
  hex.shaps[[x]] <- do.call(rbind, metal.stat)
}

hex.shap <- do.call(rbind, hex.shaps)

pairwises <- c(
  "Grand Lake - Lower Neosho", 
  "Grand Lake - Lower Neosho (past)",  
  "Grand Lake - Upper Neosho", 
  "Lower Neosho - Lower Neosho (past)", 
  "Lower Neosho - Upper Neosho",        
  "Lower Neosho (past) - Upper Neosho"
)

metalperms <- list()
for (metal in c("Zn","Cu","Cr","Mn","Ni")) {
  
  set.seed(123) 
  values <- hex.sub[[metal]]
  group_labels <- rep(c("Grand Lake", "Lower Neosho (past)","Lower Neosho","Upper Neosho"), times = c(3, 5, 12, 9)) 
  
  randomize_groups_with_labels <- function(values, group_labels) {
    shuffled <- sample(values)  # Shuffle values
    data.frame(Value = shuffled, Group = group_labels)  # Assign back group labels
  }
  
  randomizations <- replicate(999, randomize_groups_with_labels(values, group_labels), simplify = FALSE)
  
  game <- list()
  for (n in 1:999) {
    perm <- randomizations[[n]]
    perm$Group <- as.factor(perm$Group)
    dung <- PMCMRplus::gamesHowellTest(Value~Group,data=perm)
    deeta <- as.data.frame(matrix(dung$statistic, nrow = 1))
    deeta <- deeta[,c("V1","V2","V3","V5","V6","V9")]
    colnames(deeta) <- pairwises
    game[[n]] <-  data.frame(deeta)
  }
  metalperms[[metal]] <- do.call(rbind,game)
}

metal.intervals <- list()
for (metal in c("Zn","Cu","Cr","Mn","Ni")) {
  perms <- metalperms[[metal]]
  
  intervals <- list()
  for (sterts in colnames(perms)) {
    quants <- quantile(perms[[sterts]],probs=c(0.025,0.975))
    intervals[[sterts]] <- data.frame(Element=metal, Stat=sterts,data.frame(t(quants)))
  }
  metal.intervals[[metal]] <- do.call(rbind,intervals)
}
interballs <- do.call(rbind,metal.intervals)

ts2 <- list()
for (metal in c("Zn","Cu","Cr","Mn","Ni")) {
  dung2 <- PMCMRplus::gamesHowellTest(hex.sub[[metal]]~Site,data=hex.sub)
  deeta2 <- as.data.frame(matrix(dung2$statistic, nrow = 1))
  deeta2 <- deeta2[,c("V1","V2","V3","V5","V6","V9")]
  colnames(deeta2) <- pairwises
  ts2[[metal]] <- data.frame(deeta2)
}
tss <- do.call(rbind,ts2)
tss$Element <- row.names(tss)

ttests2 <- pivot_longer(tss,cols = 1:6,names_to = "Stat")

finaldata <- left_join(interballs,ttests2,by=c("Element","Stat"))

finaldata <- finaldata %>%  mutate(Stat=case_when(
                                      Stat=="Grand.Lake...Lower.Neosho"~"GL-LN",
                                      Stat=="Grand.Lake...Lower.Neosho..past."~"GL-LNp",
                                      Stat=="Grand.Lake...Upper.Neosho"~"GL-UN",
                                      Stat=="Lower.Neosho...Lower.Neosho..past."~"LN-LNp",
                                      Stat=="Lower.Neosho...Upper.Neosho"~"LN-UN",
                                      Stat=="Lower.Neosho..past....Upper.Neosho"~"LNp-UN"))

finaldata <- finaldata %>% filter(Element != "Cr")
highlight <- finaldata %>% filter(value<X2.5.|value>X97.5.)
ggplot(data=finaldata, aes(x=Stat))+geom_errorbar(aes(ymin=X2.5.,ymax=X97.5.,color=Stat))+
  geom_point(aes(y=value),size=1)+
  geom_point(data=highlight, aes(x=Stat,y=value),size=3,shape=16,color="red2")+
  facet_wrap(~Element)+theme_classic()+
  labs(x="",y="Statistic from Games-Howell pairwise test")+
  scale_color_brewer(palette = "Paired")+
  theme(text=element_text(size=14), axis.text.x = element_text(angle=45,hjust=1),
        legend.position = "none")

write.csv(finaldata, "randomreassignoutput.csv")

# T-test/Mann Whitney Upper vs. Lower ####
bug.sub <- data %>% filter(Taxa %in% c("Hexagenia","Naididae","Orthocladiinae")& 
                             !Site %in% c("Elk creek","Hickory creek","Honey creek")&
                             !Sample %in% c("Hex1p","Hex2p","Hex3p","Hex4p","Hex5p" ))

bug.shaps <- list()
for (x in unique(bug.sub$Site)) {
  df <- bug.sub %>% filter(Site == x)
  buggers <- list()
  for (berg in unique(bug.sub$Taxa)) {
    dff <- df %>% filter(Taxa==berg)
  metal.stat <- list()
  for (metal in  c("Zn","Cu","Cr","Mn","Ni")) {
    dfff <- shapiro.test(dff[[metal]])
    metal.stat[[metal]] <- data.frame(Site=x,Metal=metal,Taxa=berg,pvalue=dfff$p.value)
  }
  buggers[[berg]] <- do.call(rbind, metal.stat)
  }
  bug.shaps[[x]] <- do.call(rbind, buggers)
}

bug.shap <- do.call(rbind, bug.shaps)

#nonparametric for Naididae-Ni, Hex-Cr, and Orth-Zn
naid <- bug.sub %>% filter(Taxa=="Naididae")
wilcox.test(Ni~Site,data=naid)

hex <- bug.sub %>% filter(Taxa=="Hexagenia")
wilcox.test(Cr~Site,data=hex)

orth <- bug.sub %>% filter(Taxa=="Orthocladiinae")
wilcox.test(Zn~Site,data=orth)

#parametric for the rest

t.tests <- list()
for (bug in c("Hexagenia","Naididae","Orthocladiinae")){
  bugz <- bug.sub %>% filter(Taxa==bug)

  print(bug) 
  
output <- list()
for (metal in c( "Cu","Mn","Ni","Zn")){
  t <- t.test(bugz[[metal]]~Site, data=bugz, var.equal=FALSE)
  output[[metal]] <- data.frame(Taxa=bug,Element=metal,T=t$statistic,P=t$p.value, conf.low=t$conf.int[1], conf.high=t$conf.int[2])

  print(metal)
  print(t)
  
  }
outputs <- do.call(rbind,output)
t.tests[[bug]] <- outputs
}
ttests <- do.call(rbind,t.tests)

write.csv(ttests, "upvlowstatoutput.csv")

# Water vs Sediment vs Up/Lower bugs accumulation ####
bug.sub <-  data %>% filter(Taxa %in% c("Hexagenia","Naididae","Orthocladiinae","Sediment","Unfiltered water","Filtered water")& 
                             !Site %in% c("Elk creek","Hickory creek","Honey creek")&
                             !Sample %in% c("Hex1p","Hex2p","Hex3p","Hex4p","Hex5p" ))

alldiffs <- list()
for (tax in c("Hexagenia","Naididae","Orthocladiinae","Sediment","Unfiltered water","Filtered water")){
  deta <- bug.sub %>% filter(Taxa==tax)
  
  diffs <- list()
  for (metal in  c("Cr", "Cu","Mn","Ni","Zn","Fe")){
    ups <- deta %>% filter(Site == "Upper Neosho")
    ups <- ups[[metal]]
    lows <- deta %>% filter(Site == "Lower Neosho")
    lows <- lows[[metal]]
    diffs[[metal]] <- data.frame(Metal=metal,Taxa=tax,Differences=c(outer(ups,lows, FUN = "-") ))
  }
  
  alldiffs[[tax]] <- do.call(rbind,diffs)
}

alldiffs <- do.call(rbind,alldiffs)

###Stats


###PLOT 
alldiffs$Taxa <- factor(alldiffs$Taxa,levels=c("Hexagenia","Naididae","Orthocladiinae","Filtered water","Unfiltered water","Sediment"))

#Convert to mg/kg
alldiffs$Differences2 <- alldiffs$Differences *1000

#remove outliers
alldiffs <- alldiffs %>% mutate(Differences2 = case_when(Metal == "Ni" & Differences2 > 200 ~ NA_real_,
                                                         TRUE ~ Differences2))

plotz <- list()
for (metal in  c("Cr", "Cu","Mn","Ni","Zn","Fe")){
  df <- alldiffs %>% filter(Metal==metal)
  plotz[[metal]]  <- ggstatsplot::ggbetweenstats(df,Taxa,Differences2, 
                            violin.args = list(fill = NA, color = "transparent"),
                            plot.type = "box",results.subtitle = FALSE,
                            centrality.plotting = F, ggsignif.args = list(textsize=0),
                            caption.args = list(title = NULL),
                            
                            ylab="Difference in concentration (mg/kg)", xlab="",title = metal,
                            ggplot.component = list(theme_classic(),
                                                    theme(text = element_text(size=16),
                                                          axis.text.x = element_text(angle=45, hjust=1), 
                                                          plot.title  = element_text(size = 24),
                                                          axis.text.y.right = element_blank(), 
                                                          axis.ticks.y.right = element_blank(),
                                                          legend.position = "none"),
                                                      annotate(geom = "text", label="")))+
                          theme(axis.title.y.right = element_blank(), 
                         axis.text.y.right = element_blank(), 
                         axis.ticks.y.right = element_blank())


  }

gridExtra::grid.arrange(plotz$Zn,plotz$Ni,plotz$Cu,plotz$Mn)



#without nickel outliers
df <- alldiffs %>% filter(Metal=="Ni")
df <- df %>% filter(Differences < 0.1)
ggstatsplot::ggbetweenstats(df,Taxa,Differences, 
                               violin.args = list(fill = NA, color = "transparent"),
                               plot.type = "box",results.subtitle = FALSE,centrality.plotting = F, ggsignif.args = list(textsize=0),
                               ylab="Difference in concentration", xlab="",title = "Ni",
                               ggplot.component = list(theme(text = element_text(size=20),axis.text.x = element_text(angle=45, hjust=1), plot.title  = element_text(size = 24))))

# t-tests for environmental readings ####
env.sub <- data %>% filter(Taxa %in% c("Sediment", "Filtered water", "Unfiltered water")& 
                             !Site %in% c("Elk creek","Hickory creek","Honey creek")&
                             !Sample %in% c("Hex1p","Hex2p","Hex3p","Hex4p","Hex5p" ))

t.tests <- list()
for (bug in c("Sediment", "Filtered water", "Unfiltered water")){
  bugz <- env.sub %>% filter(Taxa==bug)
  
  print(bug)
  output <- list()
  for (metal in c("Cr", "Cu","Mn","Ni","Zn")){
    t <- t.test(bugz[[metal]]~Site, data=bugz, var.equal=FALSE)
    output[[metal]] <- data.frame(Taxa=bug,Element=metal,T=t$statistic,P=t$p.value, conf.low=t$conf.int[1], conf.high=t$conf.int[2])
  
    print(metal)
    print(t)
    }
  outputs <- do.call(rbind,output)
  t.tests[[bug]] <- outputs
  
}
ttests <- do.call(rbind,t.tests)

write.csv(ttests, "heavymetalenvirottests.csv")
# permutations for environment - taxa comparisons ####

# NOT FOR PUBLICATION Boxplot 2: Lower Neosho ####
#Ignore unless you're interesting in seeing whats in the other bugs
low.sub <- data %>% filter(Taxa %in% c("Dineutus","Hexagenia",
                                       "Naididae","Orthocladiinae",
                                       "Tanypodinae") & Site=="Lower Neosho"&
                             !Sample %in% c("Hex1p","Hex2p","Hex3p","Hex4p","Hex5p" ))

zinc <- ggplot(data=low.sub)+geom_boxplot(aes(x=Taxa,y=Zn,color=Taxa))+theme_bw()+xlab("")+theme(legend.position = "none")
copper <- ggplot(data=low.sub)+geom_boxplot(aes(x=Taxa,y=Cu,color=Taxa))+theme_bw()+xlab("")+theme(legend.position = "none")
chromium <- ggplot(data=low.sub)+geom_boxplot(aes(x=Taxa,y=Cr,color=Taxa))+theme_bw()+xlab("")+theme(legend.position = "none")
manganese <- ggplot(data=low.sub)+geom_boxplot(aes(x=Taxa,y=Mn,color=Taxa))+theme_bw()+xlab("")+theme(legend.position = "none")
nickel <- ggplot(data=low.sub)+geom_boxplot(aes(x=Taxa,y=Ni,color=Taxa))+theme_bw()+xlab("")+theme(legend.position = "none")

gridExtra::grid.arrange(zinc,copper, chromium, manganese, nickel, nrow=5) 

# NOT FOR PUBLICATION Kruskall Dunn for lower bugs ####
low.sub <- data %>% filter(Taxa %in% c("Dineutus","Hexagenia",
                                       "Naididae","Orthocladiinae",
                                       "Tanypodinae") & Site=="Lower Neosho"&
                             !Sample %in% c("Hex1p","Hex2p","Hex3p","Hex4p","Hex5p" ))

kruskstats <- list()
for (metal in c("Cr", "Cu","Mn","Ni","Zn")){
  dung2 <- dunn.test::dunn.test(low.sub[[metal]],low.sub$Taxa)
  deeta2 <- as.data.frame(matrix(dung2$Z, nrow = 1))
  colnames(deeta2) <- dung2$comparisons
  kruskstats[[metal]] <- data.frame(Chi2=dung2$chi2,deeta2)
}
lowkrusk <- do.call(rbind,kruskstats)
lowkrusk$Value <- "Stats"

kruskpvalue<- list()
for (metal in c("Cr", "Cu","Mn","Ni","Zn")){
  dung2 <- dunn.test::dunn.test(low.sub[[metal]],low.sub$Taxa)
  deeta2 <- as.data.frame(matrix(dung2$P.adjusted, nrow = 1))
  colnames(deeta2) <- dung2$comparisons
  kruskpvalue[[metal]] <- data.frame(Chi2=dung2$chi2,deeta2)
}
lowkrusk2 <- do.call(rbind,kruskpvalue)
lowkrusk2$Value <- "Pvalue"

lowkrusks <- rbind(lowkrusk,lowkrusk2)

write.csv(lowkrusks, "lowbugsstats.csv")






# NOT FOR PUBLICATION Hexagenia v. Chaoborus v. Rheumobates ####

bugs22 <- data %>% filter(Taxa %in% c("Hexagenia","Chaoborus","Rheumobates")&
                           Date < "2023-01-01 12:00:00" )
bugs22 <- bugs22 %>% mutate(Site = case_when(Site != "Lower Neosho" ~ "Grand Lake",
                                            Site == "Lower Neosho" ~ Site  ))

#PCA
pca22 <- prcomp(bugs22[,7:ncol(bugs22)],scale. = T, center = T)
p <- ggbiplot::ggbiplot(pca22, choices=c(1:2),scale=0,groups = bugs22$Taxa, 
                        point.size = -1,
                        shape=bugs22$Taxa, ellipse = T, ellipse.alpha=0.1, 
                        ellipse.linewidth = 0.7, var.factor=1.3,
                        varname.adjust = 3.3, varname.size = -1, var.axes = F)
d <- ggbiplot::ggbiplot(pca22, choices=c(1:2),scale=0,groups = bugs22$Taxa, point.size = -1,
                        shape=bugs22$Taxa, ellipse = T, ellipse.alpha=0.1, 
                        ellipse.linewidth = 0.7, var.factor=1.3,varname.adjust = 3.3, varname.size = -1, var.axes = T)
p+ geom_point(aes(shape=bugs22$Site,color=bugs22$Taxa),size=2.5,alpha=0.8)+xlab("PC1 (51.1%)")+ylab("PC2 (28.6%)")+
  geom_segment(aes(x = 0, y = 0, xend = xvar, yend = yvar),
               data = d$layers[[2]]$data,  # Access arrow data
               arrow = arrow(length = unit(0.2, "cm")),  # Arrowhead size
               linewidth = 0.5,  # Adjust thickness here
               color = "black")  +# Optional: Customize color+
  ggrepel::geom_label_repel(data = d$layers[[4]]$data,  # Extract the layer with variable labels
                            aes(x = xvar, y = yvar, label = varname),size=2, 
                            max.overlaps = 100, min.segment.length = 4, label.size = 0.05,alpha=0.8)+
  #scale_color_manual(values = c("forestgreen","red4"))+scale_fill_manual(values = c("forestgreen","red4"))+
  theme_bw()+  theme(legend.title = element_blank(),legend.key.size = unit(0.5,units = "cm"))



#boxplot
bugs22 <- data %>% filter(Taxa %in% c("Hexagenia","Chaoborus","Rheumobates")&
                            Date < "2023-01-01 12:00:00" )
bugs22 <- bugs22 %>% filter(Site == "Lower Neosho")

zinc <- ggplot(data=bugs22 )+geom_boxplot(aes(x=Taxa,y=Zn,color=Taxa))+theme_bw()+xlab("")+theme(legend.position = "none")
copper <- ggplot(data=bugs22 )+geom_boxplot(aes(x=Taxa,y=Cu,color=Taxa))+theme_bw()+xlab("")+theme(legend.position = "none")
chromium <- ggplot(data=bugs22 )+geom_boxplot(aes(x=Taxa,y=Cr,color=Taxa))+theme_bw()+xlab("")+theme(legend.position = "none")
manganese <- ggplot(data=bugs22 )+geom_boxplot(aes(x=Taxa,y=Mn,color=Taxa))+theme_bw()+xlab("")+theme(legend.position = "none")
nickel <- ggplot(data=bugs22 )+geom_boxplot(aes(x=Taxa,y=Ni,color=Taxa))+theme_bw()+xlab("")+theme(legend.position = "none")

gridExtra::grid.arrange(zinc,copper, chromium, manganese, nickel, nrow=5) 




#analyses
bugs22 <- data %>% filter(Taxa %in% c("Hexagenia","Chaoborus","Rheumobates")&
                            Date < "2023-01-01 12:00:00" )
bugs22 <- bugs22 %>% filter(Site == "Lower Neosho")

kruskstats <- list()
for (metal in c("Cr", "Cu","Mn","Ni","Zn")){
  dung2 <- dunn.test::dunn.test(bugs22[[metal]],bugs22$Taxa)
  deeta2 <- as.data.frame(matrix(dung2$Z, nrow = 1))
  colnames(deeta2) <- dung2$comparisons
  kruskstats[[metal]] <- data.frame(Chi2=dung2$chi2,deeta2)
}
lowkrusk <- do.call(rbind,kruskstats)
lowkrusk$Value <- "Stats"

kruskpvalue<- list()
for (metal in c("Cr", "Cu","Mn","Ni","Zn")){
  dung2 <- dunn.test::dunn.test(bugs22[[metal]],bugs22$Taxa)
  deeta2 <- as.data.frame(matrix(dung2$P.adjusted, nrow = 1))
  colnames(deeta2) <- dung2$comparisons
  kruskpvalue[[metal]] <- data.frame(Chi2=dung2$chi2,deeta2)
}
lowkrusk2 <- do.call(rbind,kruskpvalue)
lowkrusk2$Value <- "Pvalue"

lowkrusks <- rbind(lowkrusk,lowkrusk2)

write.csv(lowkrusks, "lowbugs22stats.csv")


# NOT FOR PUBLICATION Old code ####
allpca <- prcomp(data[,7:ncol(data)],scale. = T, center = T)

ggbiplot::ggbiplot(allpca, groups = data$Taxa, ellipse = T)

ggbiplot::ggbiplot(allpca,groups = data$Site, ellipse = T)

data.all <- data  %>%
filter(Taxa %in% c("Hexagenia","Naididae","Orthocladiinae"))


pca.shared <- prcomp(data.all[,7:ncol(data.all)],scale. = T, center = T)



bug.sub <- data %>% filter(Taxa %in% c("Hexagenia","Naididae","Orthocladiinae")& !Site %in% c("Elk creek","Hickory creek","Honey creek"))

zinc <- ggplot(data=bug.sub)+geom_boxplot(aes(x=Taxa,y=Zn,group=interaction(Taxa,Site),color=Site))+theme(legend.position = "none")
copper <- ggplot(data=bug.sub)+geom_boxplot(aes(x=Taxa,y=Cu,group=interaction(Taxa,Site),color=Site))+theme(legend.position = "none")
chromium <- ggplot(data=bug.sub)+geom_boxplot(aes(x=Taxa,y=Cr,group=interaction(Taxa,Site),color=Site))+theme(legend.position = "none")
manganese <- ggplot(data=bug.sub)+geom_boxplot(aes(x=Taxa,y=Mn,group=interaction(Taxa,Site),color=Site))+theme(legend.position = "none")
nickel <- ggplot(data=bug.sub)+geom_boxplot(aes(x=Taxa,y=Ni,group=interaction(Taxa,Site),color=Site))+ylim(0,0.025)

gridExtra::grid.arrange(zinc,copper, chromium, manganese, nickel)                     

t.test(Cu~Site, data=filter(bug.sub, bug.sub$Taxa=="Naididae"))


hex.sub <- data %>% filter(Taxa=="Hexagenia")
hex.sub <- hex.sub %>% mutate(Site=case_when(Site %in% c("Elk creek","Hickory creek","Honey creek")~"Grand Lake",
                                             !Site %in% c("Elk creek","Hickory creek","Honey creek")~Site))


zinc <- ggplot(data=hex.sub)+geom_boxplot(aes(x=Taxa,y=Zn,group=interaction(Taxa,Site),color=Site))+theme(legend.position = "none")
copper <- ggplot(data=hex.sub)+geom_boxplot(aes(x=Taxa,y=Cu,group=interaction(Taxa,Site),color=Site))+theme(legend.position = "none")
chromium <- ggplot(data=hex.sub)+geom_boxplot(aes(x=Taxa,y=Cr,group=interaction(Taxa,Site),color=Site))+theme(legend.position = "none")
manganese <- ggplot(data=hex.sub)+geom_boxplot(aes(x=Taxa,y=Mn,group=interaction(Taxa,Site),color=Site))+theme(legend.position = "none")
nickel <- ggplot(data=hex.sub)+geom_boxplot(aes(x=Taxa,y=Ni,group=interaction(Taxa,Site),color=Site))

gridExtra::grid.arrange(zinc,copper, chromium, manganese, nickel) 



zinc <- ggplot(data=data)+geom_boxplot(aes(x=Taxa,y=Zn,color=Taxa))+theme(legend.position = "none")
copper <- ggplot(data=data)+geom_boxplot(aes(x=Taxa,y=Cu,color=Taxa))+theme(legend.position = "none")
chromium <- ggplot(data=data)+geom_boxplot(aes(x=Taxa,y=Cr,color=Taxa))+theme(legend.position = "none")
manganese <- ggplot(data=data)+geom_boxplot(aes(x=Taxa,y=Mn,color=Taxa))+theme(legend.position = "none")
nickel <- ggplot(data=data)+geom_boxplot(aes(x=Taxa,y=Ni,color=Taxa))

gridExtra::grid.arrange(zinc,copper, chromium, manganese, nickel) 


ggplot(data=bug.sub)+geom_boxplot(aes(x=Taxa,y=Ca,group=interaction(Taxa,Site),color=Site))


low.sub <- data %>% filter(Taxa %in% c("Chaoborus","Dineutus","Hexagenia",
                                       "Naididae","Orthocladiinae","Rheumobates",
                                       "Tanypodinae") & Site=="Lower Neosho"&
                             !Sample %in% c("Hex1p","Hex2p","Hex3p","Hex4p","Hex5p" )))

zinc <- ggplot(data=low.sub)+geom_boxplot(aes(x=Taxa,y=Zn,color=Taxa))+theme(legend.position = "none")
copper <- ggplot(data=low.sub)+geom_boxplot(aes(x=Taxa,y=Cu,color=Taxa))+theme(legend.position = "none")
chromium <- ggplot(data=low.sub)+geom_boxplot(aes(x=Taxa,y=Cr,color=Taxa))+theme(legend.position = "none")
manganese <- ggplot(data=low.sub)+geom_boxplot(aes(x=Taxa,y=Mn,color=Taxa))+theme(legend.position = "none")
nickel <- ggplot(data=low.sub)+geom_boxplot(aes(x=Taxa,y=Ni,color=Taxa))+theme(legend.position = "none")

gridExtra::grid.arrange(zinc,copper, chromium, manganese, nickel) 
