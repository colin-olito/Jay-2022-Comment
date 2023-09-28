library(tidyr)
library(patchwork)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(gghalves)
library(dplyr)
library(colorRamps)
library(cmocean)
library(viridis)
library(ggnewscale)
library(tidyverse)
library(directlabels)
options(scipen=999) #non-scientific notation
ThemeSobr=  theme(
  panel.border = element_blank(),  
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  text = element_text(size=12),
  axis.line = element_line(colour = "grey"),
  legend.spacing.y= unit(0, 'cm'),
  axis.title = element_text(face="bold"),
  plot.title=element_text(size=12, face="bold",hjust=0.5, vjust=2),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10),
  axis.text = element_text(size = 10)
)

### Figure 3C ###
Simul=read.table(paste("~/Paper/ModelSexChrom/V3/CleanDataset/InversionTrajectories_N=1000_Fig3-S13-S15-S19.txt",sep=""), stringsAsFactors = F) #File containing all simulation with N=1000
colnames(Simul)=c("N", "u", "r", "h", "s", "Gen", "DebInv", "FinInv", "Rep", 
                  "MeanMutInv","MinMutInv","MaxMutInv","sdMutInv","FreqMutInv",
                  "MeanMutNoInv","MinMutNoInv","MaxMutNoInv","sdMutNoInv","FreqMutNoInv",
                  "InvFit", "NoInvFit","Freq","Chromosome")
Simul$Position="Y" #Create a position column
Simul[Simul$DebInv>10000000,]$Position="Autosome" #Inversion starting after position 10Mb are on the second chromosome
Simul[Simul$Position=="Y",]$Freq=Simul[Simul$Position=="Y",]$Freq * 4 #frequency of Y inversions in the population of Y chromosome, not the overall frequency
Simul$InvSize=Simul$FinInv - Simul$DebInv # Inversion size

Summary=Simul %>% group_by(N,u,r,h,s,InvSize,Position, Rep) %>% summarise(maxFreq=max(Freq), maxGen=max(Gen), InitMutNumb=min(MeanMutInv), MinSegMut=min(MeanMutNoInv)) # For each simulation,  grep its last generation (when the inversion was lost or fixed) and its maximum frequency
Summary$State="LostEarly" #Define an state defining inversion
Summary[Summary$maxGen>15020,]$State="LostLate" #Inversion lost after 20 generation or more
Summary[Summary$maxGen==24991,]$State="Segregating" #Inversion still segregating at simulation end
Summary[Summary$maxFreq>0.95,]$State="Fixed" #Inversion that reached above 0.95 frequency are considered fixed (for computation purpose, simulation stop when inversion fix,so sometime we do not observe inversion at 1.0)
Summary$StateCode=1 #For estimating the proportion of inversion fixed, note as 1 inversion fixed and 0 otherwise
Summary[Summary$State=="LostEarly",]$StateCode=0 #Define a code for each state
Summary[Summary$State=="LostLate",]$StateCode=0
Summary[Summary$State=="Segregating",]$StateCode=0
SumNoLostEarly=subset(Summary, (Summary$State!="LostEarly" & Summary$InitMutNumb>0)) #Remove inversions that were lost in fewer than 20 generations and those that are mutation free
DataSummary=SumNoLostEarly %>% group_by(N,u,r,h,s,InvSize, Position) %>% summarise(ProbSpread=mean(StateCode), MeanSegMut=mean(MinSegMut)) # For each set of parameter, compute the fraction of mutation fixed (only for not mutation-free inversion)

options(scipen=0) #non-scientific notation
Col=scales::viridis_pal(begin=0.0, end=0.8, option="A")(5) #Define color palette
n=2000000 #Focus on 2Mb inversions
DataSummary$u=factor(DataSummary$u, labels = c("mu==1 %*% 10^{-09}","mu==1 %*% 10^{-08}"), )

Base=ggplot(DataSummary[(DataSummary$InvSize == n & DataSummary$s<0),], aes( y=ProbSpread)) #Plot the result
PlotPropInvSpread=Base+geom_point(aes(x=h,color=as.factor(s)), size=1)+
  geom_line(aes(x=h, color=as.factor(s)), size=1)+
  scale_color_manual("Selection coefficient (s)", values=Col)+
  scale_size_manual("Inversion size", labels=c("0.5Mb", "1Mb", "2Mb", "5Mb"), values=c(1,2,3,4))+
  ylab("Fraction of inversions fixed")+
  xlab("h")+
  ThemeSobr+
  theme(panel.border = element_blank(), 
        legend.position = c(0.5,1.13),
        legend.direction = "horizontal",
        legend.text = element_text(face="bold"),
        legend.key = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        legend.background = element_blank(),
        strip.background.x  = element_blank(),
        axis.title = element_text(face="bold"),
        strip.text.x = element_text(face="bold"),
        plot.title=element_text(size=12, face="bold",hjust=0.5, vjust=4),
        plot.margin = margin(22, 3, 3, 3, "pt"))+
  facet_grid(u~Position, scales = "free_y", labeller = label_parsed)+
  labs(title = "Fraction of inversions fixed after 10,000 generations")

### Figure S15 ### Same as figure 3C but with different inversion size
U="mu==1 %*% 10^{-08}"
DataSummary$InvSize=as.factor(paste0("Inversion size=",DataSummary$InvSize))
DataSummary$InvSize=relevel(DataSummary$InvSize,"Inversion size=500000")
Base=ggplot(DataSummary[(DataSummary$u == U & DataSummary$s<0),], aes( y=ProbSpread))
PlotPropInvSpread_AllInv=Base+geom_point(aes(x=h,color=as.factor(s)), size=4)+
  geom_line(aes(x=h, color=as.factor(s)), size=2)+
  scale_color_manual("Selection coefficient (s)", values=Col)+
  scale_size_manual("Inversion size", labels=c("0.5Mb", "1Mb", "2Mb", "5Mb"), values=c(1,2,3,4))+
  ylab("Fraction of inversions fixed after 10,000 generations")+
  xlab("h")+
  ThemeSobr+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))+
  facet_grid(InvSize ~Position)
PlotPropInvSpread_AllInv
save_plot(paste0("~/Paper/ModelSexChrom/V3/Plot/FigS15.png"),PlotPropInvSpread_AllInv, ncol=2, nrow=3) #Fig. S?
save_plot(paste0("~/Paper/ModelSexChrom/V3/Plot/FigS15.pdf"),PlotPropInvSpread_AllInv, ncol=2, nrow=3) #Fig. S?

### Figure 3A-B ###
#Define a set or parameter to focus on #
U=1e-08
H=0.1
S=-0.01
INVSize=2000000
Simul$Code=paste(Simul$N,Simul$u,Simul$r,Simul$h,Simul$s,Simul$InvSize,Simul$Position, Simul$Rep, sep="_")
write.table(Simul, "~/Paper/ModelSexChrom/V3/CleanDataset/InversionTrajectories_N=1000_CodeB.txt",  quote=F, row.names = F) #To be used for figure S13

SimulSubEx=subset(Simul, (Simul$u==U & Simul$h==H & Simul$s==S & Simul$InvSize==INVSize)) #Extract the simulation corresponding to a certain set of parameter
summarySub=SimulSubEx %>% group_by(Code) %>% summarise(maxFreq=max(Freq), maxGen=max(Gen), InitMutNumb=min(MeanMutInv)) # For each simulation, grep its max frequency and it end generation (which can be different from 25000 in case the inversion is lost or fixed)
GoodSimulSub=SimulSubEx #Just in case, copy the file...
FixedSimul=summarySub[summarySub$maxFreq>0.95,]$Code #Grep the code of the inversion that have fixed
LostSimul=summarySub[summarySub$maxFreq<0.95,]$Code #Grep the code of the inversion that have not fixed
LineToAdd=SimulSubEx[0,] #For computation reason, we stop record population state after inversion fixation or lost. For plotting purpose, recreate ending states.
for (i in LostSimul){ #For each inversion that have not fixed 
  LastGen=summarySub[summarySub$Code==i,]$maxGen 
  FalseEndGoodSimulSub=SimulSubEx[(SimulSubEx$Code==i & SimulSubEx$Gen==LastGen),] #Grep the last generation this inversion was recorded
  FalseEndGoodSimulSub$Gen=24991 #Change it to 24991 (last recorded generation)
  FalseEndGoodSimulSub$Freq=0 #Set its frequency to 0
  FalseEndGoodSimulSub$MeanMutInv=0 # Set its mutation number to 0
  FalseMiddleGoodSimulSub=FalseEndGoodSimulSub #Same thing, but just for 10 generation after the last generation recorded
  FalseMiddleGoodSimulSub$Gen=LastGen+10
  LineToAdd=rbind(LineToAdd, FalseEndGoodSimulSub, FalseMiddleGoodSimulSub)
}
LineToAdd2=SimulSubEx[0,]
for (i in FixedSimul){ #Do the same thing for fixed inversions
  LastGen=summarySub[summarySub$Code==i,]$maxGen
  FalseEndGoodSimulSub=SimulSubEx[(SimulSubEx$Code==i & SimulSubEx$Gen==LastGen),]
  FalseEndGoodSimulSub$Gen=24991
  FalseEndGoodSimulSub$Freq=1 #Except that in this case, the inversion go to freq 1.0
  FalseEndGoodSimulSub$MeanMutInv=SimulSubEx[(SimulSubEx$Code==i & SimulSubEx$Gen==LastGen),]$MeanMutInv
  FalseEndGoodSimulSub$MinMutInv=SimulSubEx[(SimulSubEx$Code==i & SimulSubEx$Gen==LastGen),]$MinMutInv
  FalseMiddleGoodSimulSub=FalseEndGoodSimulSub
  FalseMiddleGoodSimulSub$Gen=LastGen+10
  LineToAdd2=rbind(LineToAdd2, FalseEndGoodSimulSub, FalseMiddleGoodSimulSub)
}

GoodSimulSubComplete=rbind(GoodSimulSub,LineToAdd,LineToAdd2) #Add these false simulation end to the simulation data.frame

Col=scales::viridis_pal(begin=0.0, end=0.6, option="A")(2) #Color Palette
base=ggplot(GoodSimulSubComplete) #Plot the data
PlotA=base+geom_line(aes(x=Gen, y=Freq, group=Rep, color=as.factor(Position)), size=0.5, alpha=0.3)+ #Inversion frequency
  geom_hline(yintercept = 0, linetype=2, size=0.1)+
  scale_color_manual("", values=Col, label=c("Autosomal","Y-linked"))+
  xlab("Generation")+ylab("Inversion frequency")+
  ThemeSobr+
  labs(title = "Change in inversion frequency in a finite population")+
  theme(legend.position = c(0.50,0.98),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        axis.title = element_text(face="bold"),
        plot.title=element_text(size=12, face="bold",hjust=0.5, vjust=2),
        strip.placement = "outside",
        legend.key=element_blank(),
        panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 1)))+
  scale_y_continuous(breaks = c(0.0,0.25,0.5,0.75,1.0), limits=c(-0.05,1.1))+
  facet_grid( ~ Position)

PlotB=base+geom_line(data=GoodSimulSubComplete[GoodSimulSubComplete$MeanMutInv>0,], aes(x=Gen, y=MeanMutInv, group=Rep, color=as.factor(Position)), size=0.5, alpha=0.3)+ #Mutation number 
  geom_line(data=GoodSimulSubComplete[GoodSimulSubComplete$Gen<24991,],
            aes(x=Gen, y=MeanMutInv, group=Rep, color=as.factor(Position)), size=0.5, alpha=0.3, linetype="dotted")+
  geom_hline(yintercept = 0, linetype=2, size=0.1)+
  scale_color_manual("", values=Col, label=c("Autosomal","Y-linked"))+
  xlab("Generation")+ylab("Number of mutations in inversions")+
  labs(title = "Change in inversion mutation load in a finite population")+
  ThemeSobr+
  theme(legend.position = c(0.50,0.98),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(),
        axis.title = element_text(face="bold"),
        panel.spacing = unit(1, "lines"),
        plot.title=element_text(size=12, face="bold",hjust=0.5, vjust=2),
        axis.line = element_line(colour = "grey"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  facet_grid( ~ Position)

MergedPlot=plot_grid(PlotA, PlotB, PlotPropInvSpread, ncol=1, labels = c('a', 'b', 'c')) #Merge the plot to create the figure 3.
save_plot(paste0("~/Paper/ModelSexChrom/V3/Plot/Fig3.png"),MergedPlot, nrow=3, ncol=2, base_aspect_ratio = 1) #Fig. S?
save_plot(paste0("~/Paper/ModelSexChrom/V3/Plot/Fig3.pdf"),MergedPlot, nrow=3, ncol=2, base_aspect_ratio = 1) #Fig. S?

### Figure S16 # Same as before but for N=10000 ###
Simul10k=read.table(paste("~/Paper/ModelSexChrom/V3/CleanDataset/InversionTrajectories_N=10000_FigS11-S12-S16.txt",sep=""), stringsAsFactors = F)
colnames(Simul10k)=c("N", "u", "r", "h", "s", "Gen", "DebInv", "FinInv", "Rep", 
                  "MeanMutInv","MinMutInv","MaxMutInv","sdMutInv","FreqMutInv",
                  "MeanMutNoInv","MinMutNoInv","MaxMutNoInv","sdMutNoInv","FreqMutNoInv",
                  "InvFit", "NoInvFit","Freq","Chromosome")
Simul10k$Position="Y" #Define position 
Simul10k[Simul10k$DebInv>10000000,]$Position="Autosome" 
Simul10k[Simul10k$Position=="Y",]$Freq=Simul10k[Simul10k$Position=="Y",]$Freq * 4 # Define Y inversion frequency
Simul10k$InvSize=Simul10k$FinInv - Simul10k$DebInv #Inversion size

Summary10k=Simul10k %>% group_by(N,u,r,h,s,InvSize,Position, Rep) %>% summarise(maxFreq=max(Freq), maxGen=max(Gen), InitMutNumb=min(MeanMutInv), MinSegMut=min(MeanMutNoInv))
Summary10k$State="LostEarly"
Summary10k[Summary10k$maxGen>15020,]$State="LostLate"
Summary10k[Summary10k$maxGen==24991,]$State="Segregating"
Summary10k[Summary10k$maxFreq>0.95,]$State="Fixed"
Summary10k$StateCode=1
Summary10k[Summary10k$State=="LostEarly",]$StateCode=0
Summary10k[Summary10k$State=="LostLate",]$StateCode=0
Summary10k[Summary10k$State=="Segregating",]$StateCode=0
SumNoLostEarly10k=subset(Summary10k, (Summary10k$State!="LostEarly" & Summary10k$InitMutNumb>0))
DataSummary= SumNoLostEarly10k %>% group_by(N,u,r,h,s,InvSize, Position) %>% summarise(ProbSpread=mean(StateCode), MeanSegMut=mean(MinSegMut))
Col=scales::viridis_pal(begin=0.0, end=0.8, option="A")(5)
options(scipen=0)
DataSummary$InvSize=as.factor(paste0("Inversion size=",DataSummary$InvSize))
DataSummary$InvSize=relevel(DataSummary$InvSize,"Inversion size=500000")

Base=ggplot(DataSummary[(DataSummary$s<0),], aes( y=ProbSpread))
PlotPropInvSpread_AllInv=Base+geom_point(aes(x=h,color=as.factor(s)), size=4)+
  geom_line(aes(x=h, color=as.factor(s)), size=2)+
  scale_color_manual("Selection coefficient (s)", values=Col)+
  scale_size_manual("Inversion size", labels=c("0.5Mb", "1Mb", "2Mb", "5Mb"), values=c(1,2,3,4))+
  ylab("Fraction of inversions fixed after 10,000 generations")+
  xlab("h")+
  ThemeSobr+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))+
  facet_grid(InvSize ~Position)

save_plot(paste0("~/Paper/ModelSexChrom/V3/Plot/FigS16.png"),PlotPropInvSpread_AllInv, ncol=2, nrow=3) #Fig. S?
save_plot(paste0("~/Paper/ModelSexChrom/V3/Plot/FigS16.pdf"),PlotPropInvSpread_AllInv, ncol=2, nrow=3) #Fig. S?

### Figure S11-12 ### Same as before (Figure 3A-B) but for N=10000
U=1e-08
H=0.1
S=-0.001 #For S12
#S=-0.01 #For S11
INVSize=2000000
SimulSubEx10k=subset(Simul10k, (Simul10k$u==U & Simul10k$h==H & Simul10k$s==S & Simul10k$InvSize==INVSize))
SimulSubEx10k$Code=paste(SimulSubEx10k$N,SimulSubEx10k$u,SimulSubEx10k$r,SimulSubEx10k$h,SimulSubEx10k$s,SimulSubEx10k$InvSize,SimulSubEx10k$Position, SimulSubEx10k$Rep, sep="_")

summarySub=SimulSubEx10k %>% group_by(Code) %>% summarise(maxFreq=max(Freq), maxGen=max(Gen), InitMutNumb=min(MeanMutInv))
LostSimulSub=summarySub[summarySub$maxGen<24991,]$Code # Inversion that have been lost or fixed
GoodSimulSub=SimulSubEx10k #Keep only non-lost Inversion
FixedSimul=summarySub[summarySub$maxFreq>0.95,]$Code #Grep the code of the inversion that have fixed
LostSimul=summarySub[summarySub$maxFreq<0.95,]$Code
LineToAdd=SimulSubEx10k[0,] #For computation reason, we stop record population state after inversion fixation or lost. For plotting purpose, recreate endind states.
for (i in LostSimul){ #For each inversion that have not fixed 
  LastGen=summarySub[summarySub$Code==i,]$maxGen 
  FalseEndGoodSimulSub=SimulSubEx10k[(SimulSubEx10k$Code==i & SimulSubEx10k$Gen==LastGen),] #Grep the last generation this inversion was recorded
  FalseEndGoodSimulSub$Gen=24991 #Change it to 24991 (last recorded generation)
  FalseEndGoodSimulSub$Freq=0 #Set its frequency to 0
  FalseEndGoodSimulSub$MeanMutInv=0 # Set its mutation number to 0
  FalseMiddleGoodSimulSub=FalseEndGoodSimulSub #Same thing, but just for 10 generation after the last generation recorded
  FalseMiddleGoodSimulSub$Gen=LastGen+10
  LineToAdd=rbind(LineToAdd, FalseEndGoodSimulSub, FalseMiddleGoodSimulSub)
}
LineToAdd2=SimulSubEx10k[0,]
for (i in FixedSimul){ #Do the same thing for fixed inversions
  LastGen=summarySub[summarySub$Code==i,]$maxGen
  FalseEndGoodSimulSub=SimulSubEx10k[(SimulSubEx10k$Code==i & SimulSubEx10k$Gen==LastGen),]
  FalseEndGoodSimulSub$Gen=24991
  FalseEndGoodSimulSub$Freq=1 #Except that in this case, the inversion go to freq 1.0
  FalseEndGoodSimulSub$MeanMutInv=SimulSubEx10k[(SimulSubEx10k$Code==i & SimulSubEx10k$Gen==LastGen),]$MeanMutInv
  FalseEndGoodSimulSub$MinMutInv=SimulSubEx10k[(SimulSubEx10k$Code==i & SimulSubEx10k$Gen==LastGen),]$MinMutInv
  FalseMiddleGoodSimulSub=FalseEndGoodSimulSub
  FalseMiddleGoodSimulSub$Gen=LastGen+10
  LineToAdd2=rbind(LineToAdd2, FalseEndGoodSimulSub, FalseMiddleGoodSimulSub)
}
GoodSimulSubComplete=rbind(GoodSimulSub,LineToAdd,LineToAdd2) #Add these false simulation end to the simulation data.frame
Col=scales::viridis_pal(begin=0.0, end=0.6, option="A")(2)#Color Palette#

base=ggplot(GoodSimulSubComplete) #Plot
PlotA=base+geom_line(aes(x=Gen, y=Freq, group=Rep, color=as.factor(Position)), size=0.5, alpha=0.3)+
  geom_hline(yintercept = 0, linetype=2, size=0.1)+
  scale_color_manual("", values=Col, label=c("Autosomal","Y-Linked"))+
  xlab("Generation")+ylab("Inversion frequency")+
  ThemeSobr+
  labs(title = "Change in inversion frequency in a finite population")+
  theme(legend.position = c(0.50,0.98),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        axis.title = element_text(face="bold"),
        plot.title=element_text(size=12, face="bold",hjust=0.5, vjust=2),
        strip.placement = "outside",
        legend.key=element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 1)))+
  scale_y_continuous(breaks = c(0.0,0.25,0.5,0.75,1.0), limits=c(-0.05,1.1))+
  facet_grid( ~ Position)

PlotB=base+geom_line(data=GoodSimulSubComplete[GoodSimulSubComplete$MeanMutInv>0,], aes(x=Gen, y=MeanMutInv, group=Rep, color=as.factor(Position)), size=0.5, alpha=0.3)+
  geom_line(data=GoodSimulSubComplete,
            aes(x=Gen, y=MeanMutInv, group=Rep, color=as.factor(Position)), size=0.5, alpha=0.3, linetype="dotted")+
  geom_hline(yintercept = 0, linetype=2, size=0.1)+
  scale_color_manual("", values=Col, label=c("Autosomal","Y-Linked"))+
  xlab("Generation")+ylab("Number of mutations in inversions")+
  labs(title = "Change in inversion mutation load in a finite population")+
  ThemeSobr+
  theme(legend.position = c(0.50,0.98),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(),
        axis.title = element_text(face="bold"),
        plot.title=element_text(size=12, face="bold",hjust=0.5, vjust=2),
        axis.line = element_line(colour = "grey"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  facet_grid( ~ Position)

MergedPlot=plot_grid(PlotA, PlotB, ncol=1, labels = c('a', 'b'))
save_plot(paste0("~/Paper/ModelSexChrom/V3/Plot/FigS12.png"),MergedPlot, nrow=2) 
save_plot(paste0("~/Paper/ModelSexChrom/V3/Plot/FigS12.pdf"),MergedPlot, nrow=2) 

### Figure 2 (deterministic model) ###
library(directlabels)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(viridis)
options(scipen = 0)

ThemeSobr=  theme(
  panel.border = element_blank(),  
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  text = element_text(size=12),
  axis.line = element_line(colour = "grey"),
  legend.spacing.y= unit(0, 'cm'),
  panel.spacing = unit(0.8, "lines"),
  legend.title = element_text(size = 11),
  legend.text = element_text(size = 9)
)

ResultHeatLessLoad=data.frame(h=double(), s=double(), u=double(), n=double(), q=double(), ### Table containing the result to display. 
                              nq=double(), ExpectedFreqY=double(), ExpectedFreqAutosome=double(),
                              FracFixedY=double(), FracFixedAuto=double(), TotProb=double())

for (u in c(1e-08, 1e-09)) ## Here studying inversion under various u, n, h, and s.
{
  for (n in c(2000000))
  {
    for (h in seq(0.001,0.25,0.0025))
    {
      for (s in  seq(0.001,0.1,0.001))
      {
        q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2)))) #Mutation frequencies
        nq=n*q #Average number of mutation of segment of size n.
        #B1=(q*n*h)/(1-h) #Threshold under which inversions fix on the autosome. (WNI > WNN)
        CumulProbY=0 ## Expected equilibrium Frequency
        CumulProbAuto=0 ## Expected equilibrium Frequency
        FracFixedY=0 #Fraction of fixed inversion
        FracFixedAuto=0 #Fraction of fixed autosomal inversion
        TotProb=0 #Sum of inversion occurring probability. This allows to determine the fraction of mutations that are less loaded (if we use (m in seq(0,nq))). it must equal 1 at the end if we consider (m in seq(0,n))
        for (m in seq(0,nq)) #here, consider less-loaded inversion. Replace by "m in seq(0,n)" to consider all inversions (Figure S9)
        {
          Pm=dbinom(m,n,q) #Probabilities of occuring inversion with m mutations
          if (m<nq) #If the inversion is less loaded... (which is always true here because "m in seq(0,nq-1)")
          {
            FequilY=1 #Inversion on chromosome Y fix. 
            FracFixedY=FracFixedY+Pm
            WNI=(q*(1-s) + (1-q)*(1-h*s))^m * (q*(1-h*s) + 1-q)^(n-m) #Fitness of individual heterozygous for the inversion.
            WNN=(1-2*q*(1-q)*h*s - q*q*s)^n #Fitness of individual homozygous for the absence of inversion.
            WII=(1-s)^m #Fitness of individual homozygous for the inversion.
            if (WNI>WII)
            {   
                s1=WNI-WII
                s2=WNI-WNN
                FequilAuto=s2/(s1+s2) #Autosomal inversion equlibrium frequency
            }
            else #m<B3
            {
              FequilAuto=1 #Inversion fix
              FracFixedAuto=FracFixedAuto+Pm
            }
          }
          else ## Should never occurs here !
          { #If m> nq, the inversion never rise in frequency 
             FequilY=0 
             FequilAuto=0
          }
          ProbEquilY=Pm*FequilY 
          ProbEquilAuto=Pm*FequilAuto
          CumulProbY=CumulProbY+ProbEquilY
          CumulProbAuto=CumulProbAuto+ProbEquilAuto
          TotProb=TotProb + Pm #Probabilities of less loaded inversion
        }
        ResultHeatLessLoad[nrow(ResultHeatLessLoad)+1,]=c(h,s,u, n, q, n*q, CumulProbY,CumulProbAuto,FracFixedY, FracFixedAuto, TotProb) #Fill the table
      }
    }
  }
}

ResultHeatLessLoad$Y=ResultHeatLessLoad$ExpectedFreqY/ ResultHeatLessLoad$TotProb #Because we do not considere all inversion (i.e. TotProb != 1), to get the expected frequencies of less-loaded inversion, devide the summed equilibrium frequency by the probabilities of less-loaded inversions
ResultHeatLessLoad$Autosome=ResultHeatLessLoad$ExpectedFreqAutosome/ ResultHeatLessLoad$TotProb #Same for autosomal inversion
LongTable=ResultHeatLessLoad %>% pivot_longer(cols=c(Y, Autosome), names_to="Position", values_to="ProbEquil") # Pivot the table for plotting
LongTable[LongTable$Position=="Y",]$Position="Y-linked" #Change labels
LongTable[LongTable$Position=="Autosome",]$Position="Autosomal"

# To display the set of parameter that produde nq=1, nq=2, nq=... mutations (white and red lines on the figure), create a table
nqLimit=data.frame(h=double(),u=double(), s_nq1=double(),
                   s_nq2=double())
n=2000000
for (h in seq(0.001,0.30,0.0025)) # Depend only on h and mu, and not on s
{
  for( u in c(1e-08, 1e-09))
  {
    s_nq1=(n*u)/(h*1)
    s_nq2=(n*u)/(h*2)
    nqLimit[nrow(nqLimit)+1,]=c(h,u,s_nq1,s_nq2)
  }
}

colnames(nqLimit)[3]="nq=1" #Change label for plotting purpose
colnames(nqLimit)[4]="nq=2"
nqLimitLong=nqLimit %>% pivot_longer(cols=c("nq=1","nq=2"), names_to="n", values_to="Values")
LongTable$u=factor(LongTable$u, labels = c("mu==1 %*% 10^{-09}","mu==1 %*% 10^{-08}"), )
nqLimitLong$u=factor(nqLimitLong$u, labels = c("mu==1 %*% 10^{-09}","mu==1 %*% 10^{-08}"))

base=ggplot(LongTable[LongTable$n==n,])
# Figure 2d #
InversionFrequency=base+  geom_tile(aes(x=h, y=s, fill=ProbEquil))+
  geom_step(data=nqLimitLong, aes(x=h, y=Values, color=n), direction='mid', alpha=0.5)+
  geom_dl(data=nqLimitLong, aes(x=h, y=Values, label=n, color=n),
          method=list('last.bumpup', cex = 0.8, vjust=-2.8, hjust = 3.1))+
  scale_color_manual(values=c("red","white"), guide=F)+
  coord_cartesian(ylim=c(0,0.1),xlim=c(0,0.25), expand = F)+
  facet_grid(u~Position, labeller = label_parsed)+
  ylab("Selection coefficient (s)")+
  xlab("Dominance coefficient (h)")+
  scale_fill_viridis( option="D","", direction = -1)+
  labs(title="Expected equilibrium frequency of less-loaded inversions")+
  ThemeSobr+
  theme(panel.border = element_blank(), 
        legend.position = c(0.5,1.13),
        legend.direction = "horizontal",
        legend.text = element_text(face="bold"),
        strip.text = element_text(face="bold"),
        #legend.key.width = unit(1.2,"cm"),
        legend.background = element_blank(),
        strip.background.x  = element_blank(),
        panel.spacing.x = unit(1.0, "lines"),
        strip.text.x = element_text(face="bold"),
        plot.title=element_text(size=11, face="bold",hjust=0.5, vjust=8),
        plot.margin = margin(22, 3, 3, 5, "pt"))+
  guides(fill = guide_colorbar( title = NULL,
                                # hjust = 0.5 centres the title horizontally
                                label.position = "top",
                                label.vjust =1
  ))
InversionFrequency

# Figure 2 b #
LessLoadProb=base+geom_tile(aes(x=h, y=s, fill=TotProb))+
  geom_step(data=nqLimitLong, aes(x=h, y=Values, color=n), direction='mid', alpha=0.5)+
  geom_dl(data=nqLimitLong, aes(x=h, y=Values, label=n, color=n),
          method=list('last.bumpup', cex = 0.8, vjust=-2.8, hjust =5))+
  scale_color_manual(values=c("red","white"), guide=F)+
  ylab("Selection coefficient (s)")+
  xlab("Dominance coefficient (h)")+
  coord_cartesian(ylim=c(0,0.1),xlim=c(0,0.25), expand = F)+
  facet_grid(u~., labeller = label_parsed)+
  scale_fill_viridis(option = "B", "")+
  ThemeSobr+
  theme(panel.border = element_blank(),  
        legend.position = c(0.5,1.08),
        legend.direction = "horizontal",
        legend.text = element_text(face="bold"),
        strip.text = element_text(face="bold"),
        legend.background = element_blank(),
        plot.title=element_text(size=11, face="bold",hjust=0.5, vjust=12),
        plot.margin = margin(29, 3, 3, 5, "pt"))+
  labs(title="Probability of occurence of less-loaded inversions")+
  guides(fill = guide_colorbar( title = NULL,
                                label.position = "top",
                                label.vjust =1
  ))
LessLoadProb

### Figure 2a ###
my_breaks = c(0.1, 1, 10, 100)
NMut=base+geom_tile(aes(x=h, y=s, fill=nq))+
  geom_step(data=nqLimitLong, aes(x=h, y=Values, color=n), direction='mid', alpha=0.5)+
  geom_dl(data=nqLimitLong, aes(x=h, y=Values, label=n, color=n),
          method=list('last.bumpup', cex = 0.8, vjust=-2.8, hjust = 5))+
  scale_color_manual(values=c("red","white"), guide=F)+
  coord_cartesian(ylim=c(0,0.1),xlim=c(0,0.25), expand = F)+
  ylab("Selection coefficient (s)")+
  xlab("Dominance coefficient (h)")+
  facet_grid(u~., labeller = label_parsed)+
  scale_fill_viridis(option = "A", "", trans="log10", breaks = my_breaks, labels = my_breaks)+
  ThemeSobr+
  theme(panel.border = element_blank(),  
        legend.position = c(0.5,1.08),
        legend.direction = "horizontal",
        legend.text = element_text(face="bold"),
        strip.text = element_text(face="bold"),
        #legend.key.width = unit(1.2,"cm"),
        legend.background = element_blank(),
        plot.title=element_text(size=11, face="bold",hjust=0.5, vjust=12),
        plot.margin = margin(29, 5, 3, 3, "pt"))+
  labs(title="Expected number of segregating mutations (nq)")+
  guides(fill = guide_colorbar( title = NULL,
                                # hjust = 0.5 centres the title horizontally
                                label.position = "top",
                                label.vjust =1
  ))
NMut
## Figure 2c #
library(data.table)
Line=paste("time","h","s","u","n","r","Recomb","P",
           "FXN","FXI","FXNm","FXIm","FXNf","FXIf","FXm","FXf",
           "FYN","FYI","FY","Wm","Wf","D","q", sep=" ")
u=1e-08
n=2000000
File=paste0("~/Paper/ModelSexChrom/V3/CleanDataset/TimeSimul_h_s_u=",u, "n=", n, "_XYsyst_NoMutAccumul_Example.txt")
fwrite(list(Line), File) #Directly write the result in a file, for computing purpose
for (Recomb in c(0.0,0.5)) #Linked or unlinked
{
  for (h in c(0.1))
  {
    for (s in c(0.001,0.01,0.1))
    {
      q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2))))
      P=0.95 #Inversion with 5% less mutation than average
      m=floor(P*n*q)
      WNI=(q*(1-s) + (1-q)*(1-h*s))^m * (q*(1-h*s) + 1-q)^(n-m) #Fitness of individual heterozygous for the inversion.
      WNN=(1-2*q*(1-q)*h*s - q*q*s)^n #Fitness of individual homozygous for the absence of inversion.
      WII=(1-s)^m #Fitness of individual homozygous for the inversion.
      FXNm=1.00 #Frequency of X chromosomes in males with the non-inverted segment
      FXIm=0.00  #Frequency of X chromosomes in males with the inverted segment
      FXNf=1.00  #Frequency of X chromosomes in females with the non-inverted segment
      FXIf=0.00 #Frequency of X chromosomes in females with the inverted segment
      FYN=0.99 #Frequency of Y chromosome with the non-inverted segment
      FYI=0.01 #Frequency of Y chromosome with the non-inverted segment #Introduce the inversion in 1% of Y chromosomes
      FY=FYN+FYI #Frequency of Y chromosome (must equal 1, used for check)
      FXm=FXNm+FXIm  #Frequency of X chromosome in males(must equal 1, used for check)
      FXf=FXNf+FXIf #Frequency of X chromosome in females (must equal 1, used for check)
      FXI=(2/3)*FXIf + (1/3)*FXIm #Two third of the X chromosome are in females and one third in males
      FXN=(2/3)*FXNf + (1/3)*FXNm
      Wm=FXNf*FYN*WNN + FXNf*FYI*WNI + FXIf*FYN*WNI + FXIf*FYI*WII # Mean fitness of the males
      Wf=FXNf*FXNm*WNN + FXNf*FXIm*WNI + FXNm*FXIf*WNI + FXIf*FXIm*WII  #mean fitness of the females
      D=FXI*FYN - FXN*FYI #Linkage disequilibrium
      time=0 #Initial time
      Line=paste(time,h,s,u,n,m,Recomb,P,FXN,FXI,FXNm,FXIm,FXNf,FXIf,FXm,FXf,FYN,FYI,FY,Wm,Wf,D,q,sep=" ") #Table row
      fwrite(list(Line), File, append = T)
      for (time in seq(2,10000,1)) #During 10000 generation
      {  #Recalculate the frequency of each inversion in each chromosome (see appendix for detail)
        FYI_t=(FYI*FXIf*WII +
                 FYI*FXNf*WNI*(1-Recomb) + 
                 FYN*FXIf*WNI*Recomb)/Wm
        FYN_t=(FYN*FXNf*WNN +
                 FYN*FXIf*WNI*(1-Recomb) + 
                 FYI*FXNf*WNI*Recomb)/Wm
        FXIm_t=(FXIf*FYI*WII +
                  FXIf*FYN*WNI*(1-Recomb) + 
                  FXNf*FYI*WNI*Recomb)/Wm
        FXNm_t=(FXNf*FYN*WNN +
                  FXNf*FYI*WNI*(1-Recomb) + 
                  FXIf*FYN*WNI*Recomb)/Wm
        FXIf_t=(FXIf*FXIm*WII +
                  (1/2)*(FXIm*FXNf + FXIf*FXNm)*WNI)/Wf
        FXNf_t=(FXNf*FXNm*WNN +
                  (1/2)*(FXIm*FXNf + FXIf*FXNm)*WNI)/Wf
        FYI=FYI_t #For the next generation, define the new inversion frequency
        FYN=FYN_t
        FXIm=FXIm_t
        FXNm=FXNm_t
        FXIf=FXIf_t
        FXNf=FXNf_t
        FXI=(2/3)*FXIf + (1/3)*FXIm
        FXN=(2/3)*FXNf + (1/3)*FXNm
        FY=FYN+FYI
        FX=FXN+FXI
        D=FXI*FYN - FXN*FYI
        Wm=FXNf*FYN*WNN + FXNf*FYI*WNI + FXIf*FYN*WNI + FXIf*FYI*WII
        Wf=FXNf*FXNm*WNN + FXNf*FXIm*WNI + FXNm*FXIf*WNI + FXIf*FXIm*WII
        Line=paste(time,h,s,u,n,m,Recomb,P,FXN,FXI,FXNm,FXIm,FXNf,FXIf,FXm,FXf,FYN,FYI,FY,Wm,Wf,D,q,sep=" ")
        fwrite(list(Line), File, append = T)
      }
    }
  }
}

Table=read.table(File, sep=" ", header=T, stringsAsFactors = F) #Read the dataset
Table$Linkage="Autosomal" #Add labels
Table[Table$Recomb=="0",]$Linkage="Y-linked"
Table$FI=0.25*Table$FYI + 0.75*Table$FXI #Frequency of the inversion when not linked to the Y chromosome
Table[Table$Linkage=="Y-Linked",]$FI=Table[Table$Linkage=="Y-Linked",]$FYI #Frequency of the inversion when linked to the Y chromosome (frequency among the Y chromosome)
Col=scales::viridis_pal(begin=0.0, end=0.6, option="A")(2) #Color palette
base=ggplot(Table)
PlotA=base+
  geom_hline(yintercept = 0, linetype=2, size=0.1)+
  geom_line(aes(x=time, y=FYI, linetype=as.factor(s), color=Linkage), size=1, alpha=1.0)+
  ylab("Inversion frequency")+
  xlab("Generation")+
  scale_linetype_manual("s=",values=c("solid", "dashed", "dotted"))+
  scale_color_manual("Linkage", values=Col, guide=F)+
  labs(title = "Deterministic trajectory of less-loaded inversions")+
  ThemeSobr+
  theme(panel.border = element_blank(),  
        legend.position = c(0.5,1.03),
        legend.direction = "horizontal",
        legend.key.width = unit(1.2,"cm"),
        legend.background = element_blank(),
        legend.key=element_blank(),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title=element_text(size=11, face="bold",hjust=0.5, vjust=8),
        strip.placement = "outside",
        plot.margin = margin(22, 3, 3, 3, "pt"),
        axis.line = element_line(colour = "grey"))+
  guides(linetype = guide_legend(title.position = "left", 
                                 label.position = "top",
                                 label.vjust = -2,
  ))+
  scale_y_continuous(breaks=c(0.0, 0.25, 0.5, 0.75, 1.0), limits=c(-0.05,1.05))+
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000), limits=c(0,12000))+
  geom_dl(aes(x=time, y=FYI, label=Linkage, color=Linkage),method=list('last.bumpup', cex =0.8, hjust = 0.00, vjust=0.1))
PlotA
##

# Align the four panel
plotsAB <- align_plots(NMut,LessLoadProb,align = 'hv', axis = 'rltp') #Align plots
MergedPlot1=plot_grid(plotsAB[[1]],plotsAB[[2]],  ncol=2, labels = c('a', 'b'), rel_widths = c(1,1,2))
MergedPlot2=plot_grid(PlotA, InversionFrequency,  ncol=2, labels = c('c', 'd'), rel_widths = c(1.2,2))
MergedPlot3=plot_grid(MergedPlot1, MergedPlot2,  ncol=1, labels = c('', ''))
save_plot(paste0("~/Paper/ModelSexChrom/V3/Plot/Fig2.png"),MergedPlot3, nrow=2, ncol=2) #Fig. S?
save_plot(paste0("~/Paper/ModelSexChrom/V3/Plot/Fig2.pdf"),MergedPlot3, nrow=2, ncol=2) #Fig. S?

### Figure 4 and S24###

Dir="~/Paper/ModelSexChrom/V3/CleanDataset/"
ID="1750881723915" #Seed Id of the simulation to plot (figure 4) #Used to easily plot multiple simulation
#ID="1953648752606" #Seed Id of the simulation to plot (uncomment to create figure S24) 

Precision=10000 #Plot precision
Freq=intersect(intersect(list.files(Dir, pattern = ID, full.names = T),  #file containing inversion frequency
                         list.files(Dir, pattern = "Freq", full.names = T)),
               list.files(Dir, pattern = "parsed", full.names = T))
Recom=intersect(intersect(list.files(Dir, pattern = ID, full.names = T), #file containing recombination rate
                          list.files(Dir, pattern = "Recom", full.names = T)),
                list.files(Dir, pattern = "parsed", full.names = T)) 

Data=read.table(Recom, stringsAsFactors = F, header = F) #Read the recombination file (panel B)
Pos=seq(1:(length(colnames(Data))-1))
colnames(Data)=c("Generation", Pos) #Change col names with position in Mb
Data$"101"=Data$"101"-0.5 #The bin 101 contain the position that split the two chromosome in SLiM. This position recombine with a probability 0.5 each generation (as in nature). So remove these recombination event from the recombination rate
DataLong=gather(Data, Position, Value, "1":"200") #Pivot the plot
DataLong$Position=as.numeric(as.character(DataLong$Position)) #Define position as numeric for plotting
DataLong$Chrom="Chromosome 1 (X/Y)" #Change labels
DataLong[DataLong$Position>100,]$Chrom="Chromosome 2"
DataLong$NormVal=DataLong$Value/mean(DataLong[DataLong$Generation<10000,]$Value) #Define recombination rate relative to the recombination rate during the burn-in phase
DataLong[DataLong$Chrom=="Chromosome 2",]$Position=DataLong[DataLong$Chrom=="Chromosome 2",]$Position - 100 #Change the position of the autosome (from 101-200 to 1-100)
DataSub=subset(DataLong, (DataLong$Chrom=="Chromosome 1 (X/Y)" & DataLong$Generation < 101000)) #Only sex chrom data (ugly, just used for plotting the sex locus position)
base=ggplot(DataLong[DataLong$Generation < 101000,]) #Limit plot length
PlotRecomb=base+geom_tile(aes(x=Position, y=Generation, fill=NormVal),colour="white",size=0.02)+
  scale_y_reverse(expand=c(0.005,0),breaks=pretty(DataLong[DataLong$Generation < 101000,]$Generation, 10))+
  geom_vline(data=DataSub, aes(xintercept = 50), color="red2")+
  scale_color_manual("", values=c("red"))+
  scale_x_continuous(expand=c(0.01,0.01))+
  scale_fill_viridis("Relative \n recombination rate", option = "inferno", direction = -1, limits=c(0,max(DataLong$NormVal)))+
  facet_wrap(.~fct_relevel(as.factor(Chrom), "Chromosome 1 (X/Y)", "Chromosome 2"), scale="free_x")+
  xlab("Position (Mb)")+
  theme(panel.border = element_blank(),
        plot.margin = margin(0,0,0,0, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.title = element_text(face="bold"),
        axis.title = element_text(face="bold", size=14),
        axis.text = element_text(size=10),
        strip.text = element_blank())


DataInv=read.table(Freq, stringsAsFactors = F, header = F)  #Read the inversion frequency file (panel A)
colnames(DataInv)=c("Generation", "InvStart", "InvEnd", "InitMut", "p__", "pI_", "pII", "Yfreq", "Xfreq","MeanMutInv", "MeanMutNoInv")
DataInv$InvStart=DataInv$InvStart/1000000 #Define the inversion position on the new scale (in Mb instead of in bp)
DataInv$InvEnd=DataInv$InvEnd/1000000
DataInv$InvFreq=(DataInv$Xfreq + DataInv$Yfreq)/2000 #Overall inversion frequency (note that the Yfreq and Xfreq column contain the number of chromosome with inversion )
DataInv$Yfreq=DataInv$Yfreq/500 #Define inversion frequency on the Y and X chromosome (note their is n=2000 haploid genome, so 500 Y chromosomes)
DataInv$Xfreq=DataInv$Xfreq/1500
DataInv$Chrom="Chromosome 1 (X/Y)" #Define label
DataInv[DataInv$InvStart>100,]$Chrom="Chromosome 2"
DataInv[DataInv$Chrom=="Chromosome 2",]$InvStart=DataInv[DataInv$Chrom=="Chromosome 2",]$InvStart - 100 #Change the position of the autosome (from 101-200 to 1-100)
DataInv[DataInv$Chrom=="Chromosome 2",]$InvEnd=DataInv[DataInv$Chrom=="Chromosome 2",]$InvEnd - 100

Col=scales::viridis_pal(begin=0.00, end=0.80, option="cividis")(3)
DataInvSub=subset(DataInv, (DataInv$Chrom=="Chromosome 1 (X/Y)" & DataInv$Generation %% Precision == 0 & DataInv$Generation <101000))
base2=ggplot(DataInv[(DataInv$Generation %% Precision == 0 & DataInv$Generation <101000), ])
PlotFreq=base2+
  geom_rect(data=DataInv[(DataInv$Generation %% Precision == 0 & DataInv$Generation <101000 & DataInv$Chrom=="Chromosome 1 (X/Y)"), ],
            aes(xmin=InvStart, xmax=InvEnd, ymin=0, ymax=Yfreq, fill="Y-linked"),color="black", alpha=0.6, size=0.2)+
  geom_rect(data=DataInv[(DataInv$Generation %% Precision == 0 & DataInv$Generation <101000 & DataInv$Chrom=="Chromosome 1 (X/Y)"), ], 
            aes(xmin=InvStart, xmax=InvEnd, ymin=0, ymax=Xfreq, fill="X-linked"),color="black", alpha=0.6, size=0.2)+
  geom_rect(data=DataInv[(DataInv$Generation %% Precision == 0 & DataInv$Generation <101000 & DataInv$Chrom=="Chromosome 2"), ],
            aes(xmin=InvStart, xmax=InvEnd, ymin=0, ymax=InvFreq, fill="Autosomal"),color="black", alpha=0.6, size=0.2)+
  scale_fill_manual("Inversions", values = Col)+
  ylab("Inversion Frequency")+
  scale_color_manual("", values=c("red"))+
  geom_vline(data=DataInvSub, aes(xintercept = 50), color="red2")+
  scale_x_continuous(expand=c(0.01,0.01), limits = c(0,100))+
  scale_y_continuous(expand=c(0,0),breaks = c(0.0,0.5,1.0), position = "right")+
  facet_grid(Generation~fct_relevel(as.factor(Chrom), "Chromosome 1 (X/Y)", "Chromosome 2"), scale="free_x", drop=F, switch="y")+
  theme(
    strip.text.y.left = element_text(angle = 0, size=10),
    strip.text.x = element_text(face="bold", size=14),
    panel.border = element_blank(),  
    plot.margin = margin(0,0,0,0,"pt"),
    panel.grid.major = element_line(color="grey", size=0.2),
    panel.grid.minor = element_line(color="grey", size=0.2),
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.text = element_text(face="bold", size=12),
    legend.title = element_text(face="bold", size=12),
    axis.title.y = element_text(face="bold", size=14),
    axis.text.x = element_blank(),
    axis.text.y.right = element_text(size=8))

plots <- align_plots(PlotFreq, PlotRecomb, align = 'vh', axis = 'lrtb')

MergedPlot=plot_grid(plots[[1]],plots[[2]],  ncol=1, labels = c('a', 'b'))
Sub=sub(".*//", "",sub("_InvFreq.*", "", Freq))
 save_plot(paste("~/Paper/ModelSexChrom/V3/Plot/Fig4.png"), MergedPlot, ncol = 2, nrow=2)
 save_plot(paste("~/Paper/ModelSexChrom/V3/Plot/Fig4.pdf"), MergedPlot, ncol = 2, nrow=2)
 #save_plot(paste("~/Paper/ModelSexChrom/V3/Plot/FigS24.png"), MergedPlot, ncol = 2, nrow=2)
 #save_plot(paste("~/Paper/ModelSexChrom/V3/Plot/FigS24.pdf"), MergedPlot, ncol = 2, nrow=2)

### Figure 5 ###
# To get the fraction of the sex chromosome not recombining,
# in the directory containing the result of all simulations, do (in a bash terminal):
# for i in *NRecomb_IndivSimulation_OnlyXY_Optim.txt ; do base=${i%%.txt} ; ../../Util/ParseRecombinationOutput_SepSex.pl -i $i -o $base.parsed.txt ; done
# for i in *_NRecomb_IndivSimulation_OnlyXY_Optim.parsed.txt ; do N_BP=`echo $i | sed 's/.*-BP=//' | sed 's/_.*//'`; echo $N_BP ;  REP=`echo $i | sed 's/.*_Rep_//' | sed 's/_.*//'`; echo $REP; echo -ne "$N_BP\t$REP\t" >> LastGenRecomb.txt ; tail -n 1 $i >> LastGenRecomb.txt ; done
Data=read.table("~/Paper/ModelSexChrom/V3/CleanDataset/LastGenRecomb_Fig5.txt", stringsAsFactors = F, header=F)
Pos=seq(1:(length(colnames(Data))-3))
colnames(Data)=c("N_BP","Rep","Generation", Pos) # Set the position along the genome in Mb
Data$NumberNonRecomb=0 #Number of Mb with no recombination
for (i in 1:nrow(Data)){ #For all simulation
  for (o in 4:103){ #For all position on the sex-chromosome
    if (Data[i,o]==0.0){ #If no recombination is recorded
      Data$NumberNonRecomb[i]=Data$NumberNonRecomb[i]+1 #Increment the number of non-recombining region by 1
    }
  }
}

Data$FracNonRecom=Data$NumberNonRecomb/100 #Define here the fraction of the sex-chromosme not recombining (the sex chromosome is 100Mb long)
Col4=scales::viridis_pal(begin=0.3, end=0.8, option="B")(4) #Color Palette
# Panel A #
baseNRecomb=ggplot(Data)
PlotNRecomb=baseNRecomb + geom_boxplot(aes(x=as.factor(N_BP), y=FracNonRecom, fill=as.factor(N_BP)), outlier.shape=NA) +
  geom_jitter(aes(x=as.factor(N_BP), y=FracNonRecom, fill=as.factor(N_BP)), size=1, alpha=0.3, shape=21, width = 0.1)+
  scale_fill_manual(values = Col4, guide=F)+
  ThemeSobr+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(3, 3, 3, 3, "pt"))+
   ylab("Fraction of the sex chromosome \n not recombining")

#Panel B#
# To get a summary of the inversion and reversion that have appeared, 
# in the directory containing the result of all simulations, do (in a bash terminal):
# for i in *_NewInv_Optim.txt ; do base=${i%%_NewInv_Optim.txt}; NInv=` cat $i | wc -l`; RevFile=$base"_Reversion_IndivSimulation_OnlyXY_NbMut_Optim.txt"; if [ -f "$RevFile" ]; then NRev=`cat $RevFile | wc -l` ; else NRev=0; fi; BlockRevFile=$base"_BlockedReversion_IndivSimulation_OnlyXY_NbMut_Optim.txt"; if [ -f "$BlockRevFile" ]; then NBlockRev=`cat $BlockRevFile | wc -l` ; else NBlockRev=0; fi;  N_BP=`echo $i | sed 's/.*-BP=//' | sed 's/_.*//'`; echo $N_BP ;  REP=`echo $i | sed 's/.*_Rep_//' | sed 's/_.*//'`; echo -ne "$N_BP\t$REP\t$NRev\t$NInv\t$NBlockRev\n" >> Nreversion_N=1000.txt ; done
DataRev=read.table("~/Paper/ModelSexChrom/V3/CleanDataset/Nreversion_N=1000_Fig5.txt", stringsAsFactors = F, header=F)
colnames(DataRev)=c("N_BP","Rep","N_Revers", "N_Inv", "NBlockRev")
baseNRev=ggplot(DataRev)
PlotNRev=baseNRev + geom_boxplot(aes(x=as.factor(N_BP), y=N_Revers, fill=as.factor(N_BP)), outlier.shape=NA) +
  geom_jitter(aes(x=as.factor(N_BP), y=N_Revers, fill=as.factor(N_BP)), size=1, alpha=0.3, shape=21, width = 0.1)+
  scale_fill_manual(values = Col4, guide=F)+
  ThemeSobr+xlab("Number of potential inversion breakpoints (k)")+
  theme(plot.margin = margin(3, 30, 3, 3, "pt"))+
  ylab("Number of reversions")

plotsAB <- align_plots(PlotNRecomb, PlotNRev,align = 'hv', axis = 'rltp') #Align plots
MergedPlot=plot_grid(plotsAB[[1]], plotsAB[[2]], ncol=1, labels = c('a', 'b'))
save_plot("~/Paper/ModelSexChrom/V3/Plot/Fig5.png", MergedPlot, nrow=2, base_aspect_ratio = 1.1)
save_plot("~/Paper/ModelSexChrom/V3/Plot/Fig5.pdf", MergedPlot, nrow=2, base_aspect_ratio = 1)

### Figure S9 ###
#Similar to figure 2d for for all inversion (m in seq(1, n))
ResultHeat=data.frame(h=double(), s=double(), u=double(), n=double(), q=double(), nq=double(), Y=double(), Autosome=double(), TotProb=double())
for (u in c(1e-08, 1e-09))
{
  for (n in c(2000000))
  {
    for (h in seq(0.001,0.5,0.005))
    {
      for (s in  seq(0.001,0.5,0.005))
      {
        q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2))))
        nq=n*q
        CumulProbY=0
        CumulProbAuto=0
        TotProb=0 #Sum of inversion occurring probability: must equal 1  at the end
        for (m in seq(0,nq*5)) #For computing purpose, do not consider all possible inversions, but inversion with m in 0,1,2,...,5nq (qualitatively similar to all inversions (m in seq(1 n) but hundred time quicker)
        {
          Pm=dbinom(m,n,q) #Proba of occuring inversion with m mutations
          if (m<nq)
          {
            FequilY=1 
            if(m>0)
            {
              WNI=(q*(1-s) + (1-q)*(1-h*s))^m * (q*(1-h*s) + 1-q)^(n-m) #Fitness of individual heterozygous for the inversion.
              WNN=(1-2*q*(1-q)*h*s - q*q*s)^n #Fitness of individual homozygous for the absence of inversion.
              WII=(1-s)^m #Fitness of individual homozygous for the inversion.
              if (WII > WNI)
              {
                FequilAuto=1
              }
              else
              {
                s1=WNI-WII
                s2=WNI-WNN
                FequilAuto=s2/(s1+s2)
              }
            }
            else #m=0
            {
              FequilAuto=1
            }
          }
          else { 
            FequilY=0 
            FequilAuto=0
          }
          ProbEquilY=Pm*FequilY
          ProbEquilAuto=Pm*FequilAuto
          CumulProbY=CumulProbY+ProbEquilY
          CumulProbAuto=CumulProbAuto+ProbEquilAuto
          TotProb=TotProb + Pm
        }
        ResultHeat[nrow(ResultHeat)+1,]=c(h,s,u, n, q, n*q, CumulProbY,CumulProbAuto, TotProb)
      }
    }
  }
}

ResultHeat$Y=ResultHeat$Y/ ResultHeat$TotProb
ResultHeat$Autosome=ResultHeat$Autosome/ ResultHeat$TotProb
LongTable=ResultHeat %>% pivot_longer(cols=c(Y, Autosome), names_to="Position", values_to="ProbEquil")
n=2000000 #Focus on n=2000000 
LongTable[LongTable$Position=="Y",]$Position="Y-linked" #Change label
LongTable[LongTable$Position=="Autosome",]$Position="Autosomal"
nqLimit=data.frame(h=double(),u=double(), s_nq1=double(),
                   s_nq2=double())
for (h in seq(0.001,0.5,0.005))
{
  for( u in c(1e-08, 1e-09))
  {
    s_nq1=(n*u)/(h*1)
    s_nq2=(n*u)/(h*2)
    nqLimit[nrow(nqLimit)+1,]=c(h,u,s_nq1,s_nq2)
  }
}

colnames(nqLimit)[3]="nq=1" #Label of the line 
colnames(nqLimit)[4]="nq=2"
nqLimitLong=nqLimit %>% pivot_longer(cols=c("nq=1","nq=2"), names_to="n", values_to="Values")
base=ggplot(LongTable[LongTable$n==n,])

InversionFrequency=base+  geom_tile(aes(x=h, y=s, fill=ProbEquil))+
  geom_step(data=nqLimitLong, aes(x=h, y=Values, color=n), direction='mid', alpha=0.5)+
  geom_dl(data=nqLimitLong, aes(x=h, y=Values, label=n, color=n),
          method=list('last.bumpup', cex = 0.8,vjust=-1.1, hjust = 1))+
  scale_color_manual(values=c("red","white"), guide=F)+
  coord_cartesian(ylim=c(0,0.5),xlim=c(0,0.5))+
  facet_grid(paste0("u=",u)~Position)+
  ylab(expression(paste("Selection coefficient (", italic("s"),")")))+
  xlab(expression(paste("Dominance coefficient (", italic("h"),")")))+
  scale_fill_viridis( option="B","", direction = -1)+
  labs(title="Expected equilibrium frequency of all inversions")+
  ThemeSobr+
  theme(panel.border = element_blank(), 
        legend.position = c(0.5,1.13),
        legend.direction = "horizontal",
        legend.text = element_text(face="bold"),
        legend.background = element_blank(),
        strip.background.x  = element_blank(),
        strip.text.x = element_text(face="bold"),
        plot.title=element_text(size=10, face="bold",hjust=0.5, vjust=10),
        plot.margin = margin(22, 3, 3, 3, "pt"))+
  guides(fill = guide_colorbar( title = NULL,
                                # hjust = 0.5 centres the title horizontally
                                label.position = "top",
                                label.vjust =1
  ))

save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS9.png", InversionFrequency)
save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS9.pdf", InversionFrequency)

### Figure S10 ###
Simul=read.table(paste("~/Paper/ModelSexChrom/V3/CleanDataset/X_Y_Autosomal_Inversion_Trajectories_N=1000_FigS10.txt",sep=""), stringsAsFactors = F)#Simulation runned with inversion on the X or Y chromosome, or on the autosome. We recorded simulation state every 10 generations.
colnames(Simul)=c("N", "u", "r", "h", "s", "Gen", "DebInv", "FinInv","Chrom", "Rep", "MutInv", "FreqMutInv", "InvFit", "MutNoInv","FreqMutNoInv","NoInvFit","Freq")
Simul=Simul[!(Simul$DebInv>10000000 & Simul$Chrom=="Y"),] # 20,000 inversion are runned in autosome (10,000 in Y-bearing genome, 10,000 in X-bearing genomes). Removing here all the simulations that were perform with autosomal inversion on the Y-bearing slim genome. (Same result are obtain by removed autosomal inverison on X-bearing genome)
Simul[Simul$DebInv>10000000,]$Chrom="Autosome" #Label the autosomal inversions
Simul$InvSize=Simul$FinInv - Simul$DebInv #Inversion size
SimulSub=Simul #Copyt, just in case...
SimulSub=SimulSub[!(SimulSub$Freq==0.0),] #Remove generation with lost inversion.
SimulSub$Code=paste(SimulSub$u,SimulSub$r,SimulSub$h,SimulSub$s,SimulSub$DebInv,SimulSub$FinInv,SimulSub$Chrom, SimulSub$Rep, sep="_") #Define Code
FocS=c(-0.001, -0.01,-0.1) # s value to focus on
SimulSub=SimulSub[SimulSub$s %in% FocS,] 
SimulSub20G=unique(SimulSub[SimulSub$Gen>15020,]$Code)# Keep only inversion that have not been lost the first 20 generations
SimulSub=SimulSub[SimulSub$Code %in% SimulSub20G,]
summarySub=SimulSub %>% group_by(Code) %>%summarise(LastFreq=last(Freq), maxGen=max(Gen)) # For each simulation, get its last recorded frequency and and generation
LostSimulSub=summarySub[summarySub$maxGen<25000,]$Code # Inversion that have been lost or fixed # As before, the simulation end when the inversion is lost or fixed. We grep this simulations by looking at simulation that did not end at generation 25000
GoodSimulSub=SimulSub #Full dataset to be completed with the end state of each simulation
FalseEndGoodSimulSub=SimulSub[(SimulSub$Code %in% LostSimulSub & SimulSub$Gen==15010),]#For lost or fixed inversion, grep their initial state
FalseEndGoodSimulSub$Gen=25000 #Modify their initial state to create a false end state
FalseEndGoodSimulSub$Freq=0.0 ## Defined their end frequency as 0 (they have been lost) or
if (length(summarySub[(summarySub$LastFreq>0.95 & summarySub$maxGen != 25000) ,]$Code)>0)  # Simulation for which the last recorded frequency was above 0.95 were considered to have fix (as they it is very very unlikely that they may have been lost in less than 10 generations)
{
  FixedSimul=summarySub[summarySub$LastFreq>0.95,]$Code #Grep the code of the inversions that have fixed
  FalseEndGoodSimulSub[FalseEndGoodSimulSub$Code %in% FixedSimul,]$Freq=1.0 #For all inversion that have fixed, define their end frequency as 1.0
}
GoodSimulSub=rbind(GoodSimulSub,FalseEndGoodSimulSub) #Concatenate the dataset
DataAll=GoodSimulSub #Copy the dataset, just in case
DataAll$s=-DataAll$s #Change s value, to be homogenous with the deterministic analyses
DataAllEnd=subset(DataAll, DataAll$Gen==25000) #Grep the end generation of each simulation
DataAllEnd$state="Fixed" # As before, label the end state of each simulation
DataAllEnd[DataAllEnd$Freq==0.0,]$state="Lost"
DataAll$Gen=DataAll$Gen-15000 #Simulation now start at generation 0 (burn-in generation removed for plotting)
SumEnd=DataAllEnd %>% count(h,s,Chrom,state, InvSize, sort = TRUE)  #Summarize the date for plotting summary
SumEnd$Pos=0.0
SumEnd[SumEnd$state=="Fixed",]$Pos=1.0 #Define position for plotting the text
SumEnd[(SumEnd$state=="Fixed" & SumEnd$Chrom=="Y"),]$Pos=0.25
SumEnd[(SumEnd$state=="Fixed" & SumEnd$Chrom=="X"),]$Pos=0.75
Col=scales::viridis_pal(begin=0, end=0.6, option="A")(3)

base=ggplot(DataAll[(DataAll$h==0.01 & DataAll$InvSize==2000000),])
Plot2Mb0.01=base+geom_line(aes(x=Gen, y=Freq, group=Rep, color=as.factor(Chrom)), size=0.5, alpha=0.3)+
  geom_hline(yintercept = 0, linetype=2, size=0.1)+
  scale_color_manual("", values=Col, labels=c( "Autosome","X-linked", "Y-linked"))+
  xlab("Generation")+ylab("Inversion frequency")+
  geom_text(data=SumEnd[(SumEnd$h==0.01 & SumEnd$InvSize==2000000),], aes(x=8000, y=Pos, label=paste0("c=",n), color=as.factor(Chrom)), vjust = -0.5, hjust = 0, size=6, show.legend = FALSE)+
  theme(legend.position = c(0.50,1.01),
        legend.direction = "horizontal",
        panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(size=25),
        legend.text = element_text(size=25),
        axis.line = element_line(colour = "grey"),
        strip.placement = "outside",
        strip.text.x = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 2, alpha=1.0)))+
  scale_y_continuous(breaks = c(0.0,0.25,0.5,0.75,1.0), limits=c(-0.05,1.1))+
  facet_grid(paste0("s=",s)  ~ Chrom)+ggtitle("2Mb inversions, h=0.01")

base=ggplot(DataAll[(DataAll$h==0.1  & DataAll$InvSize==2000000),])
Plot2Mb0.1=base+geom_line(aes(x=Gen, y=Freq, group=Rep, color=as.factor(Chrom)), size=0.5, alpha=0.3)+
  geom_hline(yintercept = 0, linetype=2, size=0.1)+
  scale_color_manual("", values=Col, labels=c( "Autosome","X-linked", "Y-linked"))+
  xlab("Generation")+ylab("Inversion frequency")+
  geom_text(data=SumEnd[(SumEnd$h==0.1 & SumEnd$InvSize==2000000),], aes(x=8000, y=Pos, label=paste0("c=",n), color=as.factor(Chrom)), vjust = -0.5, hjust = 0, size=6, show.legend = FALSE)+
  theme(legend.position = c(0.50,1.01),
        legend.direction = "horizontal",
        panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(size=25),
        legend.text = element_text(size=25),
        axis.line = element_line(colour = "grey"),
        strip.placement = "outside",
        strip.text.x = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 2, alpha=1.0)))+
  scale_y_continuous(breaks = c(0.0,0.25,0.5,0.75,1.0), limits=c(-0.05,1.1))+
  facet_grid(paste0("s=",s) ~ Chrom)+ggtitle("2Mb inversions, h=0.1")

base=ggplot(DataAll[(DataAll$h==0.001 & DataAll$InvSize==2000000),])
Plot2Mb0.001=base+geom_line(aes(x=Gen, y=Freq, group=Rep, color=as.factor(Chrom)), size=0.5, alpha=0.3)+
  geom_hline(yintercept = 0, linetype=2, size=0.1)+
  scale_color_manual("", values=Col, labels=c( "Autosome","X-linked", "Y-linked"))+
  xlab("Generation")+ylab("Inversion frequency")+
  geom_text(data=SumEnd[(SumEnd$h==0.001 & SumEnd$InvSize==2000000),], aes(x=8000, y=Pos, label=paste0("c=",n), color=as.factor(Chrom)), vjust = -0.5, hjust = 0, size=6, show.legend = FALSE)+
  theme(legend.position = c(0.50,1.01),
        legend.direction = "horizontal",
        panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(size=25),
        axis.line = element_line(colour = "grey"),
        legend.text = element_text(size=25),
        strip.placement = "outside",
        strip.text.x = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 2, alpha=1.0)))+
  scale_y_continuous(breaks = c(0.0,0.25,0.5,0.75,1.0), limits=c(-0.05,1.1))+
  facet_grid(paste0("s=",s) ~ Chrom)+ggtitle("2Mb inversions, h=0.001")

Plot2Mb=plot_grid(Plot2Mb0.001,Plot2Mb0.01,Plot2Mb0.1, nrow=3, labels=c("a", "b","c"), label_size = 30)
save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS10.png", Plot2Mb, ncol=3, nrow=9)
save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS10.pdf", Plot2Mb, ncol=3, nrow=9)

### Figure S13 ###
SimulSub=read.table("~/Paper/ModelSexChrom/V3/CleanDataset/InversionTrajectories_N=1000_CodeB.txt",  header=T, stringsAsFactors = F)
SimulSubEx=subset(SimulSub, (SimulSub$h==0.1 & SimulSub$s==-0.01 & SimulSub$InvSize==2000000)) #Value to focus on
summarySub=SimulSubEx %>% group_by(Code) %>% summarise(maxFreq=max(Freq), maxGen=max(Gen), InitMutNumb=min(MeanMutInv))
LostSimulSub=summarySub[summarySub$maxGen<24991,]$Code # Inversion that have been lost or fixed
GoodSimulSub=SimulSubEx #Copy, just in case
FixedSimul=summarySub[(summarySub$maxFreq>0.95 & summarySub$maxGen<24991),]$Code #Grep the code of the inversion that have fixed
LostSimul=summarySub[summarySub$maxFreq<0.95,]$Code #Mutation that have been lost
LineToAdd=SimulSubEx[0,]
for (i in LostSimul){ # As done before (Figure 3A), recreate false end state for plotting
  LastGen=summarySub[summarySub$Code==i,]$maxGen
  FalseEndGoodSimulSub=SimulSubEx[(SimulSubEx$Code==i & SimulSubEx$Gen==LastGen),]
  FalseEndGoodSimulSub$Gen=24991
  FalseEndGoodSimulSub$Freq=0
  FalseEndGoodSimulSub$MeanMutInv=0
  FalseEndGoodSimulSub$MinMutInv=0
  FalseMiddleGoodSimulSub=FalseEndGoodSimulSub
  FalseMiddleGoodSimulSub$Gen=LastGen+10
  LineToAdd=rbind(LineToAdd, FalseEndGoodSimulSub, FalseMiddleGoodSimulSub)
}

LineToAdd2=SimulSubEx[0,]
for (i in FixedSimul){
  LastGen=summarySub[summarySub$Code==i,]$maxGen
  FalseEndGoodSimulSub=SimulSubEx[(SimulSubEx$Code==i & SimulSubEx$Gen==LastGen),]
  FalseEndGoodSimulSub$Gen=24991
  FalseEndGoodSimulSub$Freq=1
  FalseEndGoodSimulSub$MeanMutInv=SimulSubEx[(SimulSubEx$Code==i & SimulSubEx$Gen==LastGen),]$MeanMutInv
  FalseEndGoodSimulSub$MinMutInv=SimulSubEx[(SimulSubEx$Code==i & SimulSubEx$Gen==LastGen),]$MinMutInv
  FalseMiddleGoodSimulSub=FalseEndGoodSimulSub
  FalseMiddleGoodSimulSub$Gen=LastGen+10
  LineToAdd2=rbind(LineToAdd2, FalseEndGoodSimulSub, FalseMiddleGoodSimulSub)
}


GoodSimulSubComplete=rbind(GoodSimulSub,LineToAdd,LineToAdd2) #Add the false end to the dataset
GoodSimulSubComplete=GoodSimulSub
GoodSimulSubComplete$Gen = GoodSimulSubComplete$Gen - 15000
options(scipen=0) #Non scientific notation
Col=scales::viridis_pal(begin=0.0, end=0.6, option="A")(2) #Color Palette

# Panel a
base=ggplot(GoodSimulSubComplete[(GoodSimulSubComplete$Position=="Autosome"),])
FreqAuto=base+geom_line(aes(x=Gen, y=Freq, group=Rep), color=Col[1],  size=0.5, alpha=0.8)+
  facet_wrap(.~paste0("u=",u), scales = "free")+
  xlab("Generation")+
  ylab("Autosomal inversion frequency")+
  labs(title = "Inversion frequency")+
  theme(strip.text.x = element_text(size=12, face="bold"))+
  ThemeSobr

# Panel c
MinMutAuto=base+geom_line(data=GoodSimulSubComplete[(GoodSimulSubComplete$Position=="Autosome"),],
                          aes(x=Gen, y=MinMutInv, group=Rep), color=Col[1],  size=0.5, alpha=0.8)+
  facet_wrap(.~paste0("u=",u), scales = "free")+
  xlab("Generation")+
  ylab("Number of mutations")+
  theme(strip.text.x = element_text(size=12, face="bold"))+
  labs(title = "Minimum number of mutations in autosomal inversions")+
  ThemeSobr

# Panel e
MeanMutAuto=base+geom_line(data=GoodSimulSubComplete[(GoodSimulSubComplete$Position=="Autosome"),],
                           aes(x=Gen, y=MeanMutInv, group=Rep), color=Col[1],  size=0.5, alpha=0.8)+
  facet_wrap(.~paste0("u=",u), scales = "free")+
  xlab("Generation")+
  ylab("Number of mutations")+
  theme(strip.text.x = element_text(size=12, face="bold"))+
  labs(title = "Mean number of mutations in autosomal inversions")+
  ThemeSobr

# Panel b
base=ggplot(GoodSimulSubComplete[(GoodSimulSubComplete$Position=="Y"),])
FreqY=base+geom_line(aes(x=Gen, y=Freq, group=Rep), color=Col[2],  size=0.5, alpha=0.8)+
  facet_wrap(.~paste0("u=",u), scales = "free")+
  xlab("Generation")+
  theme(strip.text.x = element_text(size=12, face="bold"))+
  ylab("Y-linked inversion frequency")+
  labs(title = "Inversion frequency")+
  ThemeSobr
# Panel d
MinMutY=base+geom_line(data=GoodSimulSubComplete[(GoodSimulSubComplete$Position=="Y" ),],
                       aes(x=Gen, y=MinMutInv, group=Rep),  size=0.5, alpha=0.8, color=Col[2])+
  facet_wrap(.~paste0("u=",u), scales = "free")+
  xlab("Generation")+
  theme(strip.text.x = element_text(size=12, face="bold"))+
  ylab("Number of mutations")+
  labs(title = "Minimum number of mutations in Y-linked inversions")+
  ThemeSobr

# Panel f
MeanMutY=base+geom_line(data=GoodSimulSubComplete[(GoodSimulSubComplete$Position=="Y"),],
                        aes(x=Gen, y=MeanMutInv, group=Rep), color=Col[2],  size=0.5, alpha=0.8)+
  facet_wrap(.~paste0("u=",u), scales = "free")+
  theme(strip.text.x = element_text(size=12, face="bold"))+
  xlab("Generation")+
  ylab("Number of mutations")+
  labs(title = "Mean number of mutations in Y-linked inversions")+
  ThemeSobr

# Panel g
u=1e-09
n=2000000
s=0.01
h=0.1
q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2))))
NmutAUT=n*q
t=seq(1,15000,1)
NMut=n*(u/(s*h))*(1-exp(-s*h*t))
DataTableU1=data.frame(t,u,NMut)
u=1e-08
q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2))))
NmutAUT=n*q
t=seq(1,15000,1)
NMut=n*(u/(s*h))*(1-exp(-s*h*t))
DataTableU2=data.frame(t,u,NMut)
DataTable=rbind(DataTableU1,DataTableU2)
NMutDeterministic=ggplot(DataTable)+geom_line(aes(x=t, y=NMut))+  facet_grid(.~Position)+
  ThemeSobr+
  xlab("Generation")+
  ylab("Number of mutations")+
  theme(strip.text.x = element_text(size=12, face="bold"))+
  facet_wrap(.~paste0("u=",u), scales = "free")+labs(title = "Deterministic change in mutation number")

#Merged plot
MergedPlot=plot_grid(FreqAuto, FreqY ,MinMutAuto, MinMutY ,MeanMutAuto, MeanMutY, ncol=2, labels = c('a', 'b', 'c','d','e','f'))
MergedPlot2=plot_grid(MergedPlot, NMutDeterministic, ncol=1, labels=c('','g'),rel_heights = c(2,1),rel_widths = c(2,1))
save_plot(paste0("~/Paper/ModelSexChrom/V3/Plot/FigS13.png"),MergedPlot2, nrow=4, ncol=2)
save_plot(paste0("~/Paper/ModelSexChrom/V3/Plot/FigS13.pdf"),MergedPlot2, nrow=4, ncol=2)

### Figure S17 ### Same as Figure S15 but considering that mutation have their selection coefficient drawn from a distribution
SimulLamb=read.table(paste("~/Paper/ModelSexChrom/V3/CleanDataset/MutationLambdaDistribution_InversionTrajectories_FigS17.txt",sep=""), stringsAsFactors = F)
colnames(SimulLamb)=c("N", "u", "r", "s", "Gen", "DebInv", "FinInv", "Rep", 
                      "MeanMutInv","MinMutInv","MaxMutInv","sdMutInv","FreqMutInv",
                      "MeanMutNoInv","MinMutNoInv","MaxMutNoInv","sdMutNoInv","FreqMutNoInv",
                      "InvFit", "NoInvFit","Freq","Chromosome")

SimulLamb$Position="Y"
SimulLamb[SimulLamb$DebInv>10000000,]$Position="Autosome"
SimulLamb[SimulLamb$Position=="Y",]$Freq=SimulLamb[SimulLamb$Position=="Y",]$Freq * 4
SimulLamb$InvSize=SimulLamb$FinInv - SimulLamb$DebInv
SimulLambSub=SimulLamb
summarySub_BLamb=SimulLambSub %>% group_by(N,u,r,s,InvSize,Position, Rep) %>% summarise(maxFreq=max(Freq), maxGen=max(Gen), InitMutNumb=min(MeanMutInv), MinSegMut=min(MeanMutNoInv))
summarySub_BLamb$State="LostEarly"
summarySub_BLamb[summarySub_BLamb$maxGen>15020,]$State="LostLate"
summarySub_BLamb[summarySub_BLamb$maxGen==24991,]$State="Segregating"
summarySub_BLamb[summarySub_BLamb$maxFreq>0.95,]$State="Fixed"
summarySub_BLamb$StateCode=1
summarySub_BLamb[summarySub_BLamb$State=="LostEarly",]$StateCode=0
summarySub_BLamb[summarySub_BLamb$State=="LostLate",]$StateCode=0
summarySub_BLamb[summarySub_BLamb$State=="Segregating",]$StateCode=0
SumNoLostEarly=subset(summarySub_BLamb, (summarySub_BLamb$State!="LostEarly" & summarySub_BLamb$InitMutNumb>0))
DataSummary=SumNoLostEarly %>% group_by(N,u,r,s,InvSize, Position) %>% summarise(ProbSpread=mean(StateCode), MeanSegMut=mean(MinSegMut))

options(scipen=0)
Col=scales::viridis_pal(begin=0.2, end=0.6, option="A")(2)
options(scipen=0)
DataSummary$InvSize=as.factor(paste0("Inversion size=",DataSummary$InvSize))
DataSummary$InvSize=relevel(DataSummary$InvSize,"Inversion size=500000")
DataSub=DataSummary[(DataSummary$s<0),]
DataSub$s=-DataSub$s

Base=ggplot(DataSub, aes( y=ProbSpread))
PlotPropInvSpread_AllInv=Base+geom_point(aes(x=s, color=as.factor(u)), size=4)+
  geom_line(aes(x=s, color=as.factor(u)), size=2)+
  scale_color_manual("mutation rate (u)", values=Col)+
  scale_x_log10()+
  ylab("Fraction of inversions fixed after 10,000 generations")+
  xlab("mean(s)")+
  ThemeSobr+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))+
  facet_grid(InvSize ~ Position)

save_plot(paste0("~/Paper/ModelSexChrom/V3/Plot/FigS17.png"),PlotPropInvSpread_AllInv, ncol=2, nrow=3) #Fig. S?
save_plot(paste0("~/Paper/ModelSexChrom/V3/Plot/FigS17.pdf"),PlotPropInvSpread_AllInv, ncol=2, nrow=3) #Fig. S?

### Figure S14 ###
Col=scales::viridis_pal(begin=0, end=0.6, option="A")(2)
Simul=read.table(paste("~/Paper/ModelSexChrom/V3/CleanDataset/RecombinationSuppressorEvolution_n=2Mb_FigS14.txt",sep=""), stringsAsFactors = F) #Evolution of recombination suppressors (suppress recombination in heterozygote and in homozygote, unlike inversions)
colnames(Simul)=c("N", "u", "r", "h", "s", "Gen", "DebInv", "FinInv", "Rep", "MutInv", "FreqMutInv", "InvFit", "MutNoInv","FreqMutNoInv","NoInvFit","Freq")
Simul=subset(Simul, Simul$Gen %% 10 == 0) #Keep only every 10 generation for computing purpose
Simul$InvSize=Simul$FinInv - Simul$DebInv #Size of the region affected (n)
SimulSub=Simul #Copy, just in case
SimulSub$Link="Linked" #Define linkage
SimulSub[SimulSub$DebInv > 10000000,]$Link="Unlinked" 
SimulSub$Code=paste(SimulSub$u,SimulSub$r,SimulSub$h,SimulSub$s,SimulSub$DebInv,SimulSub$FinInv, SimulSub$Rep, sep="_") #Define Code
FocS=c(-0.01,-0.05,-0.1) #Value of s to focus on 
SimulSub=SimulSub[(SimulSub$s %in% FocS),] 
SimulSub20G=unique(SimulSub[SimulSub$Gen>15020,]$Code) #Keep only inversion not lost during the first 20 generation
SimulSub=SimulSub[SimulSub$Code %in% SimulSub20G,]
summarySub=SimulSub %>% group_by(Code) %>% summarise(LastFreq=last(Freq), maxGen=max(Gen)) # As before, recreate end state
LostSimulSub=summarySub[summarySub$maxGen<25000,]$Code # Inversion that have been lost or fixed
GoodSimulSub=SimulSub #Keep only non-lost Inversion
FalseEndGoodSimulSub=SimulSub[(SimulSub$Code %in% LostSimulSub & SimulSub$Gen==15010),]#For lost inversion, grep their initial state
FalseEndGoodSimulSub$Gen=25000 # Modify their initial state
FalseEndGoodSimulSub$Freq=0.0 ## Defined their end frequency as 0 (they have been lost) or
if (length(summarySub[summarySub$LastFreq>0.95,]$Code)>0)
{
  FixedSimul=summarySub[summarySub$LastFreq>0.95,]$Code #Grep the code of the inversion that have fixed
  FalseEndGoodSimulSub[FalseEndGoodSimulSub$Code %in% FixedSimul,]$Freq=1.0 #For all inversion that have fixed, define their end frequency as 1.0
}
GoodSimulSub=rbind(GoodSimulSub,FalseEndGoodSimulSub) #Concatenate the dataset
DataAll=GoodSimulSub
DataAll$s=-DataAll$s
DataAllEnd=subset(DataAll, DataAll$Gen==25000)
DataAllEnd$state="Lost"
DataAllEnd[DataAllEnd$Freq==0.5,]$state="Segregate"
DataAllEnd[DataAllEnd$Freq==1.0,]$state="Fixed"
SumEnd=DataAllEnd %>% count(h,s,Link, InvSize, state, sort = TRUE) 
SumEnd$Pos=0.0
SumEnd[SumEnd$state=="Segregate",]$Pos=0.5
SumEnd[SumEnd$state=="Fixed",]$Pos=1.0

DataAll$Gen=DataAll$Gen-15000 #Suppress the burn-in generation.
themeInvFreq=theme(legend.position = c(0.50,0.98), #Change a bit the theme.
                   legend.direction = "horizontal",
                   panel.border = element_blank(),  
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   text = element_text(size=25),
                   axis.line = element_line(colour = "grey"),
                   strip.text.x  = element_blank(),
                   strip.text.y  = element_text(size=25),
                   strip.background.x = element_blank() )

base=ggplot(DataAll[(DataAll$h==0.1 & DataAll$InvSize==2000000),])
Plot0.1_2Mb=base+geom_line(aes(x=Gen, y=Freq, group=Rep, color=fct_relevel(as.factor(Link), "Unlinked", "Linked")), size=0.5, alpha=0.3)+
  geom_hline(yintercept = 0, linetype=2, size=0.1)+
  scale_color_manual("", values=Col)+
  xlab("Generation")+ylab("Inversion frequency")+
  themeInvFreq+
  geom_text(data=SumEnd[(SumEnd$h==0.1 & SumEnd$InvSize==2000000),], 
            aes(x=8000, y=Pos, label=paste0("c=",n), color=fct_relevel(as.factor(Link), "Unlinked", "Linked")), 
            vjust = -0.5, hjust = 0, size=6, show.legend = FALSE)+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  scale_y_continuous(breaks = c(0.0,0.25,0.5,0.75,1.0), limits=c(-0.05,1.1))+
  facet_grid(paste0("s=",s) ~ fct_relevel(as.factor(Link), "Unlinked", "Linked"))+ggtitle("2Mb recombination suppressors, h=0.1")
Plot0.1_2Mb

base=ggplot(DataAll[(DataAll$h==0.01 & DataAll$InvSize==2000000),])
Plot0.01_2Mb=base+geom_line(aes(x=Gen, y=Freq, group=Rep, color=fct_relevel(as.factor(Link), "Unlinked", "Linked")), size=0.5, alpha=0.3)+
  geom_hline(yintercept = 0, linetype=2, size=0.1)+
  scale_color_manual("", values=Col)+
  xlab("Generation")+ylab("Inversion frequency")+
  themeInvFreq+
  geom_text(data=SumEnd[(SumEnd$h==0.01 & SumEnd$InvSize==2000000),], 
            aes(x=8000, y=Pos, label=paste0("c=",n), color=fct_relevel(as.factor(Link), "Unlinked", "Linked")), 
            vjust = -0.5, hjust = 0, size=6, show.legend = FALSE)+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  scale_y_continuous(breaks = c(0.0,0.25,0.5,0.75,1.0), limits=c(-0.05,1.1))+
  facet_grid(paste0("s=",s) ~ fct_relevel(as.factor(Link), "Unlinked", "Linked"))+ggtitle("2Mb recombination suppressors, h=0.01")

base=ggplot(DataAll[(DataAll$h==0.001 & DataAll$InvSize==2000000),])
Plot0.001_2Mb=base+geom_line(aes(x=Gen, y=Freq, group=Rep, color=fct_relevel(as.factor(Link), "Unlinked", "Linked")), size=0.5, alpha=0.3)+
  geom_hline(yintercept = 0, linetype=2, size=0.1)+
  scale_color_manual("", values=Col)+
  xlab("Generation")+ylab("Inversion frequency")+
  themeInvFreq+
  geom_text(data=SumEnd[(SumEnd$h==0.001 & SumEnd$InvSize==2000000),], 
            aes(x=8000, y=Pos, label=paste0("c=",n), color=fct_relevel(as.factor(Link), "Unlinked", "Linked")), 
            vjust = -0.5, hjust = 0, size=6, show.legend = FALSE)+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  scale_y_continuous(breaks = c(0.0,0.25,0.5,0.75,1.0), limits=c(-0.05,1.1))+
  facet_grid(paste0("s=",s) ~ fct_relevel(as.factor(Link), "Unlinked", "Linked"))+ggtitle("2Mb recombination suppressors, h=0.001")

Plot2Mb=plot_grid(Plot0.001_2Mb,Plot0.01_2Mb,Plot0.1_2Mb, nrow=2, labels=c("a", "b","c"), label_size = 30)
# save_plot("~/Analysis/DelSheltering/V2/Output/Clean/FigSRecombSuppress_2Mb_17-05.png", Plot2Mb, ncol=4, nrow=6)
save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS14.png", Plot2Mb, ncol=4, nrow=6)
save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS14.pdf", Plot2Mb, ncol=4, nrow=6)

### Figure S25 ###  
data=read.table("~/Paper/ModelSexChrom/V3/CleanDataset/BurnIn_Stat_FigS25.txt", stringsAsFactors = F) #Stat computed during the burn in of each simulations
colnames(data)=c("N", "mu","r","h","s","generation", "meanNbMut", "Nmut1", "MeanFreq1", "NbMutXY", "FreqMutXY")
Col=scales::viridis_pal(begin=0.0, end=0.8, option="A")(7)#Color Palette 
options(scipen=0) #Non-scientific notation
base=ggplot(data) #Plot number of mutation
nMut=base+geom_line(aes(x=generation, y=Nmut1, color=as.factor(h)))+ 
  facet_grid(paste0("u=",mu)~paste0("s=",s), scale="free")+
  scale_color_manual("h=", values=Col)+
  ylab("Number of mutations genome-wide")+
  xlab("Generation (burn-in)")+
  ThemeSobr

FreqMut=base+geom_line(aes(x=generation, y=MeanFreq1, color=as.factor(h)))+
  facet_grid(paste0("u=",mu)~paste0("s=",s), scale="free")+
  scale_color_manual("h=", values=Col)+
  ylab("Mean frequency of mutations")+
  xlab("Generation (burn-in)")+
  ThemeSobr

PlotMerg=plot_grid(nMut, FreqMut,  ncol=1, labels = c('a', 'b'))
save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS25.png", nrow = 2, ncol=2,PlotMerg)
save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS25.pdf", nrow = 2, ncol=2,PlotMerg)
### Figure S20-21 ###

ResultHeatLessLoad=data.frame(h=double(), s=double(), u=double(), n=double(), q=double(),
                              nq=double(), HC=double(), ExpectedFreqY=double(), ExpectedFreqAutosome=double(),
                              FracFixedY=double(), FracFixedAuto=double(), TotProb=double())
for (u in c(1e-08))
{
  for (n in c(1000000,2000000,5000000))
  {
    for (h in seq(0.001,0.25,0.0025))
    {
      for (s in  seq(0.001,0.1,0.001))
      {
        for (HC in c(0, 0.001,0.005,0.01,0.1))
        {
          q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2))))
          nq=n*q
          CumulProbY=0
          CumulProbAuto=0
          FracFixedY=0
          FracFixedAuto=0
          TotProb=0 #Sum of inversion occurring probability: must equal 1 at the end if all inversion are considered (m in seq(0,n))
          for (m in seq(0,nq)) # Only consider less-loaded inversions
          {
            Pm=dbinom(m,n,q) #Proba of occuring inversion with m mutations
            WNI=(exp(m*q*(2*h*s-s)-h*s*(n*q+m)))*(1-HC)
            WNN=exp(n*q*q*(2*h*s-s)-h*s*(n*q*2))
            WII=exp(-s*m)
            if (WNI > WNN)
            {
              FequilY=1 
              FracFixedY=FracFixedY+Pm
              if(WNI>WII)
              {
                s1=WNI-WII
                s2=WNI-WNN
                FequilAuto=s2/(s1+s2)
              }
              else #WII>WNI
              {
                FequilAuto=1
                FracFixedAuto=FracFixedAuto+Pm
              }
            }
            else #WNN>WNI
            { 
              FequilY=0 
              FequilAuto=0
            }
            ProbEquilY=Pm*FequilY
            ProbEquilAuto=Pm*FequilAuto
            CumulProbY=CumulProbY+ProbEquilY
            CumulProbAuto=CumulProbAuto+ProbEquilAuto
            TotProb=TotProb + Pm
          }
          ResultHeatLessLoad[nrow(ResultHeatLessLoad)+1,]=c(h,s,u, n, q, n*q,HC, CumulProbY,CumulProbAuto,FracFixedY, FracFixedAuto, TotProb)
        }
      }
    }
  }
}


ResultHeatLessLoad$FracFixedY=ResultHeatLessLoad$FracFixedY/ ResultHeatLessLoad$TotProb
ResultHeatLessLoad$FracFixedAuto=ResultHeatLessLoad$FracFixedAuto/ ResultHeatLessLoad$TotProb
LongTableFracFix=ResultHeatLessLoad %>% pivot_longer(cols=c(FracFixedY, FracFixedAuto), names_to="Position", values_to="FracFixed")
LongTableFracFix[LongTableFracFix$Position=="FracFixedY",]$Position="Y-linked"
LongTableFracFix[LongTableFracFix$Position=="FracFixedAuto",]$Position="Autosomal"

nqLimit=data.frame(h=double(),u=double(),n=double(),HC=double(), s_nq1=double(),
                   s_nq2=double(), s_nq3=double(), s_nq4=double(), s_nq5=double(), s_nq6=double())

for (h in seq(0.001,0.30,0.0025)) # Depend only on h and mu, and not on s
{
  for (n in c(1000000,2000000,5000000))
  {
    for (HC in c(0, 0.001,0.005,0.01))
    {
    s_nq1=(n*u)/(h*1)
    s_nq2=(n*u)/(h*2)
    s_nq3=(n*u)/(h*3)
    s_nq4=(n*u)/(h*4)
    s_nq5=(n*u)/(h*5)
    s_nq6=(n*u)/(h*6)
    nqLimit[nrow(nqLimit)+1,]=c(h,u,n,HC,s_nq1,s_nq2,s_nq3,s_nq4,s_nq5,s_nq6)
    }
  }
}

colnames(nqLimit)[5]="nq=1" #Change label for plotting purpose
colnames(nqLimit)[6]="nq=2"
colnames(nqLimit)[7]="nq=3"
colnames(nqLimit)[8]="nq=4"
colnames(nqLimit)[9]="nq=5"
colnames(nqLimit)[10]="nq=6"
nqLimitLong=nqLimit %>% pivot_longer(cols=c("nq=1","nq=2","nq=3","nq=4","nq=5","nq=6"), names_to="nq", values_to="Values")

base=ggplot(LongTableFracFix[(LongTableFracFix$HC %in% c(0, 0.001,0.005, 0.01) & LongTableFracFix$n==1000000),])
InversionFrac=base+  geom_tile(aes(x=h, y=s, fill=FracFixed))+
  geom_step(data=nqLimitLong[nqLimitLong$n==1000000,], aes(x=h, y=Values, color=nq), direction='mid', alpha=0.5)+
  scale_color_manual(values=rep("red",6), guide=F)+
  coord_cartesian(ylim=c(0,0.1),xlim=c(0,0.25), expand = F)+
  facet_grid(HC~Position)+
  ylab("Selection coefficient (s)")+
  xlab("Dominance coefficient (h)")+
  scale_fill_viridis( option="D","", direction = -1)+
  labs(title="Fraction of less-loaded inversions that went to fixation")+
  ThemeSobr+
  theme(panel.border = element_blank(), 
        legend.position = c(0.5,1.04),
        legend.key.size = unit(1.2, 'cm'),
        legend.direction = "horizontal",
        legend.text = element_text(face="bold", size=16),
        strip.text = element_text(face="bold", size=18),
        legend.background = element_blank(),
        strip.background.x  = element_blank(),
        axis.title = element_text(face="bold", size=16),
        strip.text.x = element_text(face="bold"),
        plot.title=element_text(size=18, face="bold",hjust=0.5, vjust=9),
        plot.margin = margin(35, 3, 3, 3, "pt"))+
  guides(fill = guide_colorbar( title = NULL,
                                label.position = "top",
                                label.vjust =1
  ))
InversionFrac
save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS20.pdf", InversionFrac, nrow=4, ncol = 2)
save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS20.png", InversionFrac, nrow=4, ncol = 2)

base=ggplot(LongTableFracFix[(LongTableFracFix$HC %in% c(0, 0.001,0.005, 0.01) & LongTableFracFix$n==5000000),])
InversionFrac=base+ geom_tile(aes(x=h, y=s, fill=FracFixed))+
  geom_step(data=nqLimitLong[nqLimitLong$n==5000000,], aes(x=h, y=Values, color=nq), direction='mid', alpha=0.5)+
  scale_color_manual(values=rep("red",6), guide=F)+
  coord_cartesian(ylim=c(0,0.1),xlim=c(0,0.25), expand = F)+
  facet_grid(HC~Position)+
  ylab("Selection coefficient (s)")+
  xlab("Dominance coefficient (h)")+
  scale_fill_viridis( option="D","", direction = -1)+
  labs(title="Fraction of less-loaded inversions that went to fixation")+
  ThemeSobr+
  theme(panel.border = element_blank(), 
        legend.position = c(0.5,1.04),
        legend.key.size = unit(1.2, 'cm'),
        legend.direction = "horizontal",
        legend.text = element_text(face="bold", size=16),
        strip.text = element_text(face="bold", size=18),
        legend.background = element_blank(),
        strip.background.x  = element_blank(),
        axis.title = element_text(face="bold", size=16),
        strip.text.x = element_text(face="bold"),
        plot.title=element_text(size=18, face="bold",hjust=0.5, vjust=9),
        plot.margin = margin(35, 3, 3, 3, "pt"))+
  guides(fill = guide_colorbar( title = NULL,
                                label.position = "top",
                                label.vjust =1
  ))
InversionFrac
save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS21.pdf", InversionFrac, nrow=4, ncol = 2)
save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS21.png", InversionFrac, nrow=4, ncol = 2)

### Figure S22 ###
Simul=read.table(paste("~/Paper/ModelSexChrom/V3/CleanDataset/HaploDiplo_InvTrajectories_FigS22.txt",sep=""), stringsAsFactors = F) #File containing all simulation with N=1000

colnames(Simul)=c("N", "u", "r", "h", "s", "Gen", "DebInv", "FinInv", "FreqHaplo", "Rep", 
                  "MeanMutInv","MinMutInv","MaxMutInv","sdMutInv","FreqMutInv",
                  "MeanMutNoInv","MinMutNoInv","MaxMutNoInv","sdMutNoInv","FreqMutNoInv",
                  "InvFit", "NoInvFit","Freq")
Simul$InvSize=Simul$FinInv - Simul$DebInv # Inversion size
Simul$Position="Sex-linked inversions"
Simul[Simul$DebInv==500000,]$Position="Autosomal inversions"
Simul$Code=paste(Simul$N,Simul$u,Simul$r,Simul$h,Simul$s,Simul$InvSize,Simul$FreqHaplo,Simul$Position, Simul$Rep, sep="_")
Summary=Simul %>% group_by(N,u,r,h,s,InvSize,FreqHaplo, Position, Rep) %>% summarise(maxFreq=max(Freq), maxGen=max(Gen)) # For each simulation,  grep its last generation (when the inversion was lost or fixed) and its maximum frequency
Summary$State="LostEarly" #Define an state defining inversion
Summary[Summary$maxGen>15020,]$State="LostLate" #Inversion lost after 20 generation or more
Summary[Summary$maxGen==24991,]$State="Segregating" #Inversion still segregating at simulation end
Summary[Summary$maxFreq==0.5,]$State="Fixed" #Inversion that reached above 0.95 frequency are considered fixed (for computation purpose, simulation stop when inversion fix,so sometime we do not observe inversion at 1.0)
Summary$StateCode=1 #For estimating the proportion of inversion fixed, note as 1 inversion fixed and 0 otherwise
Summary[Summary$State=="LostEarly",]$StateCode=0
Summary[Summary$State=="LostLate",]$StateCode=0
Summary[Summary$State=="Segregating",]$StateCode=0

SumNoLostEarly=subset(Summary, Summary$State!="LostEarly") #Remove inversion that were lost in fewer than 20 generations
DataSummary=SumNoLostEarly %>% group_by(N,u,r,h,s,InvSize, FreqHaplo, Position) %>% summarise(ProbSpread=mean(StateCode)) # For each set of parameter, compute the fraction of mutation fixed (only for not mutation-free inversion)
options(scipen=0) #Scientific notation
Col=scales::viridis_pal(begin=0.0, end=0.8, option="A")(4) #Define color

Base=ggplot(DataSummary, aes( y=ProbSpread)) #Plot the result
PlotPropInvSpread=Base+geom_point(aes(x=as.factor(h), color=as.factor(FreqHaplo)), size=3, alpha=0.7)+
  scale_color_manual("Life cycle", values=Col, labels=c("Haplodiplontic (1/2 haploid)","Haplodiplontic (1/3 haploid)","Haplodiplontic (1/10 haploid)", "Fully diploid"))+
  geom_hline(yintercept = 0, linetype="dotted")+
  ylab("Fraction of inversions fixed \n after 10,000 generations")+
  xlab("h")+
  ThemeSobr+
  facet_grid(.~Position)+
  theme(panel.border = element_blank(), 
        legend.text = element_text(face="bold"),
        legend.key = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        legend.background = element_blank(),
        axis.title = element_text(face="bold"),
        strip.text.x = element_text(face="bold"),
        plot.title=element_text(size=12, face="bold",hjust=0.5, vjust=4))
PlotPropInvSpread

save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS22.png", PlotPropInvSpread, ncol=2, base_aspect_ratio = 1)
save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS22.pdf", PlotPropInvSpread, ncol=2, base_aspect_ratio = 1)

### Figure S23 ###
data=read.table("~/Paper/ModelSexChrom/V3/CleanDataset/HaploDiplo_InitStat_FigS23.txt", stringsAsFactors = F)
colnames(data)=c("N", "mu","r","h","s","FreqHaplo", "generation", "meanNbMut", "Nmut1", "MeanFreq1")
Col=scales::viridis_pal(begin=0.0, end=0.8, option="A")(4)
options(scipen=0)
base=ggplot(data)
nMut=base+geom_line(aes(x=generation, y=Nmut1, color=as.factor(FreqHaplo)))+
  facet_grid(paste0("h=",h)~paste0("s=",s), scale="free")+
  scale_color_manual("Frequency of haploid phases \n (every x generation)", values=Col)+
  ylab("Number of mutations \n segregating genome-wide")+
  xlab("Generation (burn-in)")+
  ThemeSobr

Mean_nMut=base+geom_line(aes(x=generation, y=meanNbMut, color=as.factor(FreqHaplo)))+
  facet_grid(paste0("h=",h)~paste0("s=",s), scale="free")+
  scale_color_manual("Frequency of haploid phases \n (every x generation)", values=Col)+
  ylab("Average number of mutations \n per genome")+
  xlab("Generation (burn-in)")+
  ThemeSobr

FreqMut=base+geom_line(aes(x=generation, y=MeanFreq1, color=as.factor(FreqHaplo)))+
  facet_grid(paste0("h=",h)~paste0("s=",s), scale="free")+
  scale_color_manual("Frequency of haploid phases \n (every x generation)", values=Col)+
  ylab("Mean frequency of mutations")+
  xlab("Generation (burn-in)")+
  ThemeSobr

PlotMerg=plot_grid(nMut, Mean_nMut, FreqMut,  ncol=1, labels = c('a', 'b', 'c'))
save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS23.pdf", nrow = 3, ncol=2,PlotMerg)
save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS23.png", nrow = 3, ncol=2,PlotMerg)

### Figure S1 ###
## set parameter to plot ##
u=1e-08
h=0.1
s=0.01
q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2))))
n=2000000
BinomProba=data.frame(m=double(),u=double(), Proba=double()) #Dataset 
for (m in seq(0,n*q*8)) #Only consider m varying from 0 to 8nq to save time
{
  Prob=dbinom(m,n,q) #Calculate binomial probabilities
  BinomProba[nrow(BinomProba)+1,]=c(m,u,Prob) #Fill the table (ugly way of doing that...)
}

dfProb_areaStep <- bind_rows(old = BinomProba, 
                             new = BinomProba %>% mutate(Proba = lag(Proba)),
                             .id = "source") %>% arrange(m, source)

dfProb_areaStep$Proba[1]=0
B3=n*q
B1=(q*h*n)/(1-h)
dfProb_areaStep$m=dfProb_areaStep$m-0.5 # to center values
Col=scales::viridis_pal(begin=0.2, end=0.8, option="A")(3)
base=ggplot(BinomProba, aes(x=m, y=Proba))
PlotProbaStat=base+geom_ribbon(data = dfProb_areaStep[dfProb_areaStep$m<B1+1,],aes(xmin = -1, xmax=B1, ymin = 0, ymax = Proba, fill="WII>WNI (inversion fixation)"))+
  geom_ribbon(data = dfProb_areaStep[(dfProb_areaStep$m>B1 & dfProb_areaStep$m<B3+1),],aes(xmin = -1, xmax=B3+2, ymin = 0, ymax = Proba, fill="WNI>WII (inversion maintained at intermediate frequency)"))+
  geom_ribbon(data = dfProb_areaStep[dfProb_areaStep$m>=B3,],aes(xmin = -1, xmax=max(dfProb_areaStep$m), ymin = 0, ymax = Proba, fill="WNI<WNN (inversion lost)"))+
  xlim(-0.5,n*q*3)+
  scale_fill_manual("Inversion state", values=Col)+
  annotate("segment", x=B1, xend=B1,y=0, yend=max(BinomProba$Proba)+0.05*max(BinomProba$Proba), color="black")+
  annotate("segment", x=B3, xend=B3,y=0, yend=max(BinomProba$Proba)+0.05*max(BinomProba$Proba), color="black")+
  xlab("m (number of mutations in the inversion)")+
  ylab("Probability of occurence")+
  annotate("text", x=B1, y=max(BinomProba$Proba)+0.1*max(BinomProba$Proba), label="m=qhn/(1-h)")+
  annotate("text", x=B3, y=max(BinomProba$Proba)+0.1*max(BinomProba$Proba), label="m=nq")+
  ThemeSobr+
  theme(legend.text = element_text(size=12, face="bold"),
        legend.position = c(0.75, 0.5),
        legend.title = element_text(size=12, face="bold"))+
  ggtitle(paste0("n=",n,", u=", u, ", s=", s, ", h=", h))

save_plot(paste0("~/Paper/ModelSexChrom/V3/Plot/FigS1.png"),PlotProbaStat, nrow=1, ncol=1, base_aspect_ratio = 3) #Fig. S?
save_plot(paste0("~/Paper/ModelSexChrom/V3/Plot/FigS1.pdf"),PlotProbaStat, nrow=1, ncol=1, base_aspect_ratio = 3) #Fig. S?

### Fig S2-3 ###
Data=read.table("~/Paper/ModelSexChrom/V3/CleanDataset/Linked_2MbInv_Gen1000_FigS2.txt", stringsAsFactors = F) #FigS3
# Data=read.table("~/Paper/ModelSexChrom/V3/CleanDataset/Unlinked_2MbInv_Gen1000_FigS3.txt", stringsAsFactors = F) #FigS2
colnames(Data)=c("N", "u", "r", "h", "s", "Gen", "DebInv", "FinInv", "Rep", "MutInv", "FreqMutInv", "InvFit", "MutNoInv","FreqMutNoInv","NoInvFit","Freq")
Data=Data[Data$MutInv != -1, ]
Data$Code=paste(Data$u,Data$r,Data$h,Data$s,Data$DebInv,Data$FinInv, Data$Rep, sep="_")
Data$MutInv=Data$MutInv ### To remove from count the inversion mutation
DataBegin=Data[Data$Gen==15002,]
DataEnd=Data[Data$Gen==16002,]
DataBeginLost=DataBegin[! DataBegin$Code %in% DataEnd$Code, ]
FalseDataEnd=DataBeginLost
FalseDataEnd$Gen=16002
FalseDataEnd$Freq=0.0
DataAllEnd=rbind(DataEnd, FalseDataEnd)
DataAllEnd$State="Lost"
DataAllEnd[DataAllEnd$Freq>0,]$State="Segregate"
DataAllEnd$RelativeMutNumber=DataAllEnd$MutInv/DataAllEnd$MutNoInv
DataAllEnd$RelativeMutFreq=DataAllEnd$FreqMutInv/DataAllEnd$FreqMutNoInv


FocS=c(-0.01,-0.05,-0.1)
FocH=c(0.001,0.01,0.1)
DataAllEnd=DataAllEnd[(DataAllEnd$s %in% FocS & DataAllEnd$h %in% FocH),]
Col=scales::viridis_pal(begin=0.4, end=0.6, option="A")(2)
Plot=DataAllEnd %>% 
  group_split( s) %>% 
  map( ~ggplot(., aes(x=State, y=RelativeMutNumber)) +
         geom_hline(yintercept = 1.0,linetype = 2)+
         geom_half_point_panel(aes(fill=RelativeMutFreq), shape=21,size=2, range_scale = 1.0)+
         scale_fill_gradient2("Relative frequency of\nmutations in inversions",
                              low = "navyblue", 
                              mid = "white", 
                              high = "firebrick", 
                              midpoint = 1,
         ) +
         new_scale("fill") +
         geom_half_boxplot(aes(fill=State),side = "l",  errorbar.draw = FALSE, outlier.shape = NA, alpha=1)+
         scale_fill_manual("", values=Col,guide=FALSE)+
         xlab("Inversion state after 1000 generations")+
         ylab("Relative number of mutations")+
         facet_grid(paste0("s=",s)~paste0("h=", h), scales = "free_y", labeller = function(x) label_value(x, multi_line = FALSE))+
         ggtitle("Dominance")+
         theme(panel.border = element_blank(),  
               legend.background = element_blank(),
               panel.background = element_blank(),
               text = element_text(size=14),
               axis.line = element_line(colour = "black"),
               axis.title.y = element_text(face="bold"),
               axis.text.x=element_text(face="bold"),
               legend.text = element_text(size = 10),
               legend.title = element_text(size = 12),
               legend.key.size = unit(15, 'pt'),
               plot.title = element_text(hjust = 0.5, size=12, vjust=-1, face="bold")
         )
  ) %>% plot_grid(plotlist = ., align = 'hv', ncol = 1)

#save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS3.png", Plot, ncol=3,nrow=3, base_aspect_ratio = 1)
#save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS3.pdf", Plot, ncol=3,nrow=3, base_aspect_ratio = 1)
 save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS2.png", Plot, ncol=3,nrow=3, base_aspect_ratio = 1)
 save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS2.pdf", Plot, ncol=3,nrow=3, base_aspect_ratio = 1)

### Fig S18 ###
Simul=read.table(paste("~/Paper/ModelSexChrom/V3/CleanDataset/N=1000_2MbChromosomeFusion_Trajectory_FigS18.txt",sep=""), stringsAsFactors = F) #File containing all simulation with N=1000
colnames(Simul)=c("N", "u", "r", "h", "s", "Gen", "DebInv", "FinInv", "Rep", 
                  "MeanMutInv","MinMutInv","MaxMutInv","sdMutInv","FreqMutInv",
                  "MeanMutNoInv","MinMutNoInv","MaxMutNoInv","sdMutNoInv","FreqMutNoInv",
                  "InvFit", "NoInvFit","Freq","Chromosome")

Simul[Simul$Chromosome=="Y",]$Freq=Simul[Simul$Chromosome=="Y",]$Freq * 4 #frequency of Y inversions in the population of Y chromosome, not the overall frequency
Simul[Simul$Chromosome=="X",]$Freq=Simul[Simul$Chromosome=="X",]$Freq * (1/0.75) #frequency of Y inversions in the population of Y chromosome, not the overall frequency

Simul$InvSize=Simul$FinInv - Simul$DebInv # Inversion size
Simul$Gen = Simul$Gen - 15000
Col=scales::viridis_pal(begin=0.2, end=0.8, option="A")(2)
base=ggplot(Simul) #Plot the data
PlotA=base+geom_line(aes(x=Gen, y=Freq, group=Rep, color=as.factor(Chromosome)), size=0.5, alpha=0.3)+ #Inversion frequency
  geom_hline(yintercept = 0, linetype=2, size=0.1)+
  scale_color_manual("", values=Col, label=c("X-autosome fusion","Y-autosome fusion"))+
  xlab("Generation")+ylab("Fused chromosome frequency")+
  ThemeSobr+
  theme(legend.position = c(0.50,0.98),
        legend.direction = "horizontal",
        legend.text = element_text(face="bold", size=16),
        legend.background = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        axis.title = element_text(face="bold", size=16),
        plot.title=element_text(size=12, face="bold",hjust=0.5, vjust=2),
        strip.placement = "outside",
        legend.key=element_blank(),
        panel.spacing = unit(1, "lines"),
        strip.background.x = element_blank(),
        strip.text.x = element_blank()
  )+
  guides(color = guide_legend(override.aes = list(size = 1)))+
  scale_y_continuous(breaks = c(0.0,0.25,0.5,0.75,1.0), limits=c(-0.05,1.1))+
  facet_grid(paste0("h=", h) ~ Chromosome)

save_plot(paste0("~/Paper/ModelSexChrom/V3/Plot/FigS18.pdf"),PlotA, nrow=3, ncol=2) #Fig. S18
save_plot(paste0("~/Paper/ModelSexChrom/V3/Plot/FigS18.png"),PlotA, nrow=3, ncol=2) #Fig. S18

### Figure S19 ###
Simul=read.table(paste("~/Paper/ModelSexChrom/V3/CleanDataset/InversionTrajectories_N=1000_Fig3-S13-S15-S19.txt",sep=""), stringsAsFactors = F) #File containing all simulation with N=1000
colnames(Simul)=c("N", "u", "r", "h", "s", "Gen", "DebInv", "FinInv", "Rep", 
                  "MeanMutInv","MinMutInv","MaxMutInv","sdMutInv","FreqMutInv",
                  "MeanMutNoInv","MinMutNoInv","MaxMutNoInv","sdMutNoInv","FreqMutNoInv",
                  "InvFit", "NoInvFit","Freq","Chromosome")

Simul=subset(Simul, Simul$Gen %% 10 == 1) #Keep only every 10 generation Simul=Simul[!(Simul$DebInv>10000000 & Simul$Chrom=="Y"),] # 20000 inversion are runned in autosome. Removing here all the ionversion that were the Y-bearing slim genome to get 10000 inversion
Simul=Simul[!Simul$DebInv>10000000,]
Simul=Simul[!Simul$Chromosome=="X",]
Simul$InvSize=Simul$FinInv - Simul$DebInv
unique(Simul$InvSize)
SimulSub=Simul
SimulSub$cM=SimulSub$r * (SimulSub$DebInv-1) * 100 # Change the unity
SimulSub$Code=paste(SimulSub$u,SimulSub$r,SimulSub$h,SimulSub$s,SimulSub$DebInv,SimulSub$FinInv,SimulSub$Chrom, SimulSub$Rep, sep="_") #Define Code
summarySub=SimulSub %>% group_by(Code) %>% summarise(maxFreq=max(Freq), maxGen=max(Gen)) # For each simulation, grep its max frequency and it end generation (which can be differet from 25000 in case the inversion is lost or fixed)
LostSimulSub=summarySub[summarySub$maxGen<24991,]$Code # Inversion that have been lost or fixed
NonLostSimulSub=summarySub[summarySub$maxGen==24991,]$Code #Inversion still segregating at the end
GoodSimulSub=SimulSub[SimulSub$Code %in% NonLostSimulSub,] #Keep only non-lost Inversion
FalseEndGoodSimulSub=SimulSub[(SimulSub$Code %in% LostSimulSub & SimulSub$Gen==15001),]#For lost inversion, grep their initial state
FalseEndGoodSimulSub$Gen=24991 #Modify their initial state
FalseEndGoodSimulSub$Freq=0.0 ## Defined their end frequency as 0 (they have been lost) or
GoodSimulSubEnd=rbind(GoodSimulSub[GoodSimulSub$Gen==24991,],FalseEndGoodSimulSub) #Concatenate the dataset
DataAllEnd=GoodSimulSubEnd
DataAllEnd$State=0
DataAllEnd[DataAllEnd$Freq>0,]$State=1
DataSummary=DataAllEnd %>% group_by(u,r,h,s,DebInv,FinInv,Chromosome) %>% summarise(ProbSpread=mean(State), MeanSegMut=mean(MeanMutNoInv))
Col=scales::viridis_pal(begin=0.1, end=0.9, option="A")(4)
FocS=c(-0.01,-0.1, -0.001)
FocH=c(0.00,0.01,0.10,0.5)
Base=ggplot(DataSummary[DataSummary$s %in% FocS,], aes( y=ProbSpread))
Base=ggplot(DataSummary[(DataSummary$s %in% FocS & DataSummary$h %in% FocH),], aes( y=ProbSpread))
PlotSegMut=Base+geom_point(aes(x=MeanSegMut,color=as.factor(h), shape=as.factor(s), size=as.factor(FinInv - DebInv)))+
  scale_color_manual("h", values=Col)+
  scale_size_manual("Inversion size", labels=c("0.5Mb", "1Mb", "2Mb", "5Mb"), values=c(1,2,3,4))+
  scale_x_log10(breaks=c(0,1,10,100))+
  scale_shape_manual("s", values=c(15,16,17,18,8))+
  ylab("Fraction of inversions spreading or fixed")+
  xlab("Mean number of mutations segregating at the focal region")+
  theme(
    panel.border = element_blank(),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    text = element_text(size=10),
    axis.line = element_line(colour = "grey"),
    strip.placement = "outside",
    strip.text = element_blank(),
    legend.spacing.y= unit(0, 'cm'),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS19.pdf", PlotSegMut)
save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS19.png", PlotSegMut)

### Figure S4 ###
### Perform the simul ### Same as Figure 3C
TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), r=double(), P=double(),
                      FXN=double(),FXI=double(),FXNm=double(),FXIm=double(),FXNf=double(),FXIf=double(),FXm=double(),FXf=double(),  
                      FYN=double(),FYI=double(), FY=double(), Wm=double(), Wf=double(),D=double(),q=double(), Chrom=character())
write.table(TableLinkR, "~/Paper/ModelSexChrom/V3/CleanDataset/DeterministicInvEvolution_XY_Y.txt", append = F, quote=F, row.names = F)
P=0.80
Chrom="Y-linked"
for (R in c(0.0, 0.001, 0.01, 0.05, 0.5))
{
  for (h in c(0.001, 0.01, 0.1))
  {
    for (s in c(0.01, 0.1))
    {
      TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), r=double(), P=double(),
                            FXN=double(),FXI=double(),FXNm=double(),FXIm=double(),FXNf=double(),FXIf=double(),FXm=double(),FXf=double(),  
                            FYN=double(),FYI=double(), FY=double(), Wm=double(), Wf=double(),D=double(),q=double())
      u=1e-08
      n=2000000
      q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2))))
      m=floor(P*n*q)
      WNI=(q*(1-s) + (1-q)*(1-h*s))^m * (q*(1-h*s) + 1-q)^(n-m) #Fitness of individual heterozygous for the inversion.
      WNN=(1-2*q*(1-q)*h*s - q*q*s)^n #Fitness of individual homozygous for the absence of inversion.
      WII=(1-s)^m #Fitness of individual homozygous for the inversion.
      FXNm=1.00
      FXIm=0.00
      FXNf=1.00
      FXIf=0.00
      FYN=0.99
      FYI=0.01
      FY=FYN+FYI
      FXm=FXNm+FXIm
      FXf=FXNf+FXIf
      FXI=(2/3)*FXIf + (1/3)*FXIm
      FXN=(2/3)*FXNf + (1/3)*FXNm
      Wm=FXNf*FYN*WNN + FXNf*FYI*WNI + FXIf*FYN*WNI + FXIf*FYI*WII
      Wf=FXNf*FXNm*WNN + FXNf*FXIm*WNI + FXNm*FXIf*WNI + FXIf*FXIm*WII
      D=FXI*FYN - FXN*FYI
      FirstGen=c(1,h,s,u,n,R,P,FXN,FXI,FXNm,FXIm,FXNf,FXIf,FXm,FXf,FYN,FYI,FY,Wm,Wf,D,q)
      TableLinkR[1,]=FirstGen
      for (time in seq(2,10000,1))
      {
        FYI=(TableLinkR$FYI[time-1]*TableLinkR$FXIf[time-1]*WII +
               TableLinkR$FYI[time-1]*TableLinkR$FXNf[time-1]*WNI*(1-R) + 
               TableLinkR$FYN[time-1]*TableLinkR$FXIf[time-1]*WNI*R)/TableLinkR$Wm[time-1]
        FYN=(TableLinkR$FYN[time-1]*TableLinkR$FXNf[time-1]*WNN +
               TableLinkR$FYN[time-1]*TableLinkR$FXIf[time-1]*WNI*(1-R) + 
               TableLinkR$FYI[time-1]*TableLinkR$FXNf[time-1]*WNI*R)/TableLinkR$Wm[time-1]
        
        FXIm=(TableLinkR$FXIf[time-1]*TableLinkR$FYI[time-1]*WII +
                TableLinkR$FXIf[time-1]*TableLinkR$FYN[time-1]*WNI*(1-R) + 
                TableLinkR$FXNf[time-1]*TableLinkR$FYI[time-1]*WNI*R)/TableLinkR$Wm[time-1]
        FXNm=(TableLinkR$FXNf[time-1]*TableLinkR$FYN[time-1]*WNN +
                TableLinkR$FXNf[time-1]*TableLinkR$FYI[time-1]*WNI*(1-R) + 
                TableLinkR$FXIf[time-1]*TableLinkR$FYN[time-1]*WNI*R)/TableLinkR$Wm[time-1]
        FXIf=(TableLinkR$FXIf[time-1]*TableLinkR$FXIm[time-1]*WII +
                (1/2)*(TableLinkR$FXIm[time-1]*TableLinkR$FXNf[time-1] + TableLinkR$FXIf[time-1]*TableLinkR$FXNm[time-1])*WNI)/TableLinkR$Wf[time-1]
        FXNf=(TableLinkR$FXNf[time-1]*TableLinkR$FXNm[time-1]*WNN +
                (1/2)*(TableLinkR$FXIm[time-1]*TableLinkR$FXNf[time-1] + TableLinkR$FXIf[time-1]*TableLinkR$FXNm[time-1])*WNI)/TableLinkR$Wf[time-1]
        FXI=(2/3)*FXIf + (1/3)*FXIm
        FXN=(2/3)*FXNf + (1/3)*FXNm
        FY=FYN+FYI
        FX=FXN+FXI
        D=FXI*FYN - FXN*FYI
        Wm=FXNf*FYN*WNN + FXNf*FYI*WNI + FXIf*FYN*WNI + FXIf*FYI*WII
        Wf=FXNf*FXNm*WNN + FXNf*FXIm*WNI + FXNm*FXIf*WNI + FXIf*FXIm*WII
        Gen=c(time,h,s,u,n,R,P,FXN,FXI,FXNm,FXIm,FXNf,FXIf,FXm,FXf,FYN,FYI,FY,Wm,Wf,D,q)
        TableLinkR[nrow(TableLinkR)+1,]=Gen
      }
      TableLinkR$Chrom=Chrom
      write.table(TableLinkR, "~/Paper/ModelSexChrom/V3/CleanDataset/DeterministicInvEvolution_XY_Y.txt", append = T, quote=F, row.names = F, col.names = F)
    }
  }
}

TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), r=double(), P=double(),
                      FXN=double(),FXI=double(),FXNm=double(),FXIm=double(),FXNf=double(),FXIf=double(),FXm=double(),FXf=double(),  
                      FYN=double(),FYI=double(), FY=double(), Wm=double(), Wf=double(),D=double(),q=double(), Chrom=character())
write.table(TableLinkR, "~/Paper/ModelSexChrom/V3/CleanDataset/DeterministicInvEvolution_XY_X.txt", append = F, quote=F, row.names = F)
P=0.80
Chrom="X-linked"
for (R in c(0.0, 0.001, 0.01, 0.05, 0.5))
{
  for (h in c(0.001, 0.01, 0.1))
  {
    for (s in c(0.01, 0.1))
    {
      TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), r=double(), P=double(),
                            FXN=double(),FXI=double(),FXNm=double(),FXIm=double(),FXNf=double(),FXIf=double(),FXm=double(),FXf=double(),  
                            FYN=double(),FYI=double(), FY=double(), Wm=double(), Wf=double(),D=double(),q=double())
      u=1e-08
      n=2000000
      q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2))))
      m=floor(P*n*q)
      WNI=(q*(1-s) + (1-q)*(1-h*s))^m * (q*(1-h*s) + 1-q)^(n-m) #Fitness of individual heterozygous for the inversion.
      WNN=(1-2*q*(1-q)*h*s - q*q*s)^n #Fitness of individual homozygous for the absence of inversion.
      WII=(1-s)^m #Fitness of individual homozygous for the inversion.
      FXNm=0.99
      FXIm=0.01
      FXNf=0.99
      FXIf=0.01
      FYN=1.00
      FYI=0.00
      FY=FYN+FYI
      FXm=FXNm+FXIm
      FXf=FXNf+FXIf
      FXI=(2/3)*FXIf + (1/3)*FXIm
      FXN=(2/3)*FXNf + (1/3)*FXNm
      Wm=FXNf*FYN*WNN + FXNf*FYI*WNI + FXIf*FYN*WNI + FXIf*FYI*WII
      Wf=FXNf*FXNm*WNN + FXNf*FXIm*WNI + FXNm*FXIf*WNI + FXIf*FXIm*WII
      D=FXI*FYN - FXN*FYI
      FirstGen=c(1,h,s,u,n,R,P,FXN,FXI,FXNm,FXIm,FXNf,FXIf,FXm,FXf,FYN,FYI,FY,Wm,Wf,D,q)
      TableLinkR[1,]=FirstGen
      for (time in seq(2,10000,1))
      {
        FYI=(TableLinkR$FYI[time-1]*TableLinkR$FXIf[time-1]*WII +
               TableLinkR$FYI[time-1]*TableLinkR$FXNf[time-1]*WNI*(1-R) + 
               TableLinkR$FYN[time-1]*TableLinkR$FXIf[time-1]*WNI*R)/TableLinkR$Wm[time-1]
        FYN=(TableLinkR$FYN[time-1]*TableLinkR$FXNf[time-1]*WNN +
               TableLinkR$FYN[time-1]*TableLinkR$FXIf[time-1]*WNI*(1-R) + 
               TableLinkR$FYI[time-1]*TableLinkR$FXNf[time-1]*WNI*R)/TableLinkR$Wm[time-1]
        
        FXIm=(TableLinkR$FXIf[time-1]*TableLinkR$FYI[time-1]*WII +
                TableLinkR$FXIf[time-1]*TableLinkR$FYN[time-1]*WNI*(1-R) + 
                TableLinkR$FXNf[time-1]*TableLinkR$FYI[time-1]*WNI*R)/TableLinkR$Wm[time-1]
        FXNm=(TableLinkR$FXNf[time-1]*TableLinkR$FYN[time-1]*WNN +
                TableLinkR$FXNf[time-1]*TableLinkR$FYI[time-1]*WNI*(1-R) + 
                TableLinkR$FXIf[time-1]*TableLinkR$FYN[time-1]*WNI*R)/TableLinkR$Wm[time-1]
        FXIf=(TableLinkR$FXIf[time-1]*TableLinkR$FXIm[time-1]*WII +
                (1/2)*(TableLinkR$FXIm[time-1]*TableLinkR$FXNf[time-1] + TableLinkR$FXIf[time-1]*TableLinkR$FXNm[time-1])*WNI)/TableLinkR$Wf[time-1]
        FXNf=(TableLinkR$FXNf[time-1]*TableLinkR$FXNm[time-1]*WNN +
                (1/2)*(TableLinkR$FXIm[time-1]*TableLinkR$FXNf[time-1] + TableLinkR$FXIf[time-1]*TableLinkR$FXNm[time-1])*WNI)/TableLinkR$Wf[time-1]
        FXI=(2/3)*FXIf + (1/3)*FXIm
        FXN=(2/3)*FXNf + (1/3)*FXNm
        FY=FYN+FYI
        FX=FXN+FXI
        D=FXI*FYN - FXN*FYI
        Wm=FXNf*FYN*WNN + FXNf*FYI*WNI + FXIf*FYN*WNI + FXIf*FYI*WII
        Wf=FXNf*FXNm*WNN + FXNf*FXIm*WNI + FXNm*FXIf*WNI + FXIf*FXIm*WII
        Gen=c(time,h,s,u,n,R,P,FXN,FXI,FXNm,FXIm,FXNf,FXIf,FXm,FXf,FYN,FYI,FY,Wm,Wf,D,q)
        TableLinkR[nrow(TableLinkR)+1,]=Gen
      }
      TableLinkR$Chrom=Chrom
      write.table(TableLinkR, "~/Paper/ModelSexChrom/V3/CleanDataset/DeterministicInvEvolution_XY_X.txt", append = T, quote=F, row.names = F, col.names = F)
    }
  }
}


Data_Y=read.table("~/Paper/ModelSexChrom/V3/CleanDataset/DeterministicInvEvolution_XY_Y.txt", stringsAsFactors = F, header = T)
Data_X=read.table("~/Paper/ModelSexChrom/V3/CleanDataset/DeterministicInvEvolution_XY_X.txt", stringsAsFactors = F, header = T)
FocR=c(0,0.1,1,5,50) # Value of r to focus on
FocS=c(0.01,0.1) # Value of s to focus on
Data_Y=Data_Y[Data_Y$s %in% FocS,]
Data_Y$r=Data_Y$r*100
Data_Y=Data_Y[Data_Y$r %in% FocR,]
Data_Y$Chrom="Y-linked"
Data_X=Data_X[Data_X$s %in% FocS,]
Data_X$r=Data_X$r*100
Data_X=Data_X[Data_X$r %in% FocR,]
Data_X$Chrom="X-linked"
Data=rbind(Data_X,Data_Y)
Data$Freq=Data$FYI+Data$FXI

Col=scales::viridis_pal(begin=0, end=0.8, option="A")(5)
base=ggplot(Data[Data$Chrom=="Y-linked",])
Plot0.8_Y=base+
  geom_hline(yintercept = 0, linetype=2)+
  geom_line(aes(x=time, y=Freq, color=paste0(as.factor(r), " cM")), size=1, alpha=0.8)+
  facet_grid(paste0("s=",s)~paste0("h=",h))+
  ylab("Inversion frequency")+
  xlab("Generation")+
  scale_color_manual("Distance between the inversion and the permanently heterozygous allele:", values=Col)+
  theme(panel.border = element_blank(),  
        legend.position = "top",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.placement = "outside",
        plot.margin = margin(10, 50, 10, 10, "pt"),
        text = element_text(size=24),
        axis.line = element_line(colour = "grey"))+
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1.0), limits=c(-0.08,1.1))+
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000))

base=ggplot(Data[Data$Chrom=="X-linked",])
Plot0.8_X=base+
  geom_hline(yintercept = 0, linetype=2)+
  geom_line(aes(x=time, y=Freq, color=paste0(as.factor(r), " cM")), size=1, alpha=0.8)+
  facet_grid(paste0("s=",s)~paste0("h=",h))+
  ylab("Inversion frequency")+
  xlab("Generation")+
  scale_color_manual("Distance between the inversion and the permanently heterozygous allele:", values=Col)+
  theme(panel.border = element_blank(),  
        legend.position = "top",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.placement = "outside",
        plot.margin = margin(10, 50, 10, 10, "pt"),
        text = element_text(size=24),
        axis.line = element_line(colour = "grey"))+
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1.0), limits=c(-0.08,1.1))+
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000))

Plot=plot_grid(Plot0.8_Y, Plot0.8_X, nrow=2, labels=c('a','b'), label_size = 25)
save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS4.png", Plot, ncol = 3, nrow=4)
save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS4.pdf", Plot, ncol = 3, nrow=4)

### Fig S5 ###
TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), r=double(), P=double(),
                      FXN=double(),FXI=double(),FX=double(),
                      FYN=double(),FYI=double(), FY=double(), 
                      W=double(),D=double(), q=double())
write.table(TableLinkR, "~/Paper/ModelSexChrom/V3/CleanDataset/DeterministicInvEvolution_2MT.txt", append = F, quote=F, row.names = F)

for (P in c(0.8, 0.9, 0.99))
{
  for (R in c(0.0, 0.001, 0.01, 0.05, 0.5))
  {
    for (h in c(0.001, 0.01, 0.1))
    {
      for (s in c(0.01, 0.1))
      {
        TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), r=double(), P=double(),
                              FXN=double(),FXI=double(),FX=double(),
                              FYN=double(),FYI=double(), FY=double(), 
                              W=double(),D=double(),q=double())
        u=1e-08
        n=2000000
        q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2))))
        m=floor(P*n*q)
        WNI=(q*(1-s) + (1-q)*(1-h*s))^m * (q*(1-h*s) + 1-q)^(n-m) #Fitness of individual heterozygous for the inversion.
        WNN=(1-2*q*(1-q)*h*s - q*q*s)^n #Fitness of individual homozygous for the absence of inversion.
        WII=(1-s)^m #Fitness of individual homozygous for the inversion.
        FXN=1.0
        FXI=0.00
        FYN=0.99
        FYI=0.01
        FY=FYN+FYI
        FX=FXN+FXI
        W=FXN*FYN*WNN + FXN*FYI*WNI + FXI*FYN*WNI + FXI*FYI*WII
        D=FXI*FYN - FXN-FYI
        FirstGen=c(1,h,s,u,n,R,P,FXN,FXI,FX,FYN,FYI,FY,W,D,q)
        TableLinkR[nrow(TableLinkR)+1,]=FirstGen
        for (time in seq(2,10000,1))
        {
          FXI2=(FXI*FYI*WII + FXI*FYN*WNI*(1-R) + FXN*FYI*WNI*R)/W
          FXN2=(FXN*FYN*WNN + FXN*FYI*WNI*(1-R) + FXI*FYN*WNI*R)/W
          FYI2=(FYI*FXI*WII + FYI*FXN*WNI*(1-R) + FYN*FXI*WNI*R)/W
          FYN2=(FYN*FXN*WNN + FYN*FXI*WNI*(1-R) + FYI*FXN*WNI*R)/W
          FXI=FXI2
          FXN=FXN2
          FYI=FYI2
          FYN=FYN2
          FY=FYN+FYI
          FX=FXN+FXI
          D=FXI*FYN - FXN-FYI
          W=FXN*FYN*WNN + FXN*FYI*WNI + FXI*FYN*WNI + FXI*FYI*WII
          Gen=c(time,h,s,u,n,R,P,FXN,FXI,FX,FYN,FYI,FY,W,D,q)
          TableLinkR[nrow(TableLinkR)+1,]=Gen
        }
        write.table(TableLinkR, "~/Paper/ModelSexChrom/V3/CleanDataset/DeterministicInvEvolution_2MT.txt", append = T, quote=F, row.names = F, col.names = F)
      }
    }
  }
}

Col=scales::viridis_pal(begin=0, end=0.8, option="A")(5)
Data=read.table("~/Paper/ModelSexChrom/V3/CleanDataset/DeterministicInvEvolution_2MT.txt", stringsAsFactors = F, header = T)
FocS=c(0.01,0.1)
Data=Data[Data$s %in% FocS,]
Data$r=Data$r*100
FoccM=c(0, 0.1, 1, 5, 50)
Data=Data[Data$r %in% FoccM,]
base=ggplot(Data[Data$P==0.8,])
Plot0.8=base+
  geom_hline(yintercept = 0, linetype=2)+
  geom_line(aes(x=time, y=(FYI+FXI)/2, color=paste0(as.factor(r), " cM")), size=1, alpha=0.8)+
  facet_grid(paste0("s=",s)~paste0("h=",h))+
  ylab("Inversion frequency")+
  xlab("Generation")+
  scale_color_manual("Distance between the inversion and the permanently heterozygous allele:", values=Col)+
  theme(panel.border = element_blank(), 
        legend.position = "top",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.placement = "outside",
        plot.margin = margin(10, 50, 10, 10, "pt"),
        text = element_text(size=24),
        axis.line = element_line(colour = "grey"))+
  ylim(-0.08,0.6)+
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000), limits=c(0,10000))


base=ggplot(Data[Data$P==0.9,])
Plot0.9=base+
  geom_hline(yintercept = 0, linetype=2)+
  geom_line(aes(x=time, y=(FYI+FXI)/2, color=paste0(as.factor(r), " cM")), size=1, alpha=0.8)+
  facet_grid(paste0("s=",s)~paste0("h=",h))+
  ylab("Inversion frequency")+
  xlab("Generation")+
  scale_color_manual("Distance between the inversion and the permanently heterozygous allele:", values=Col)+
  theme(panel.border = element_blank(), 
        legend.position = "top",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.placement = "outside",
        plot.margin = margin(10, 50, 10, 10, "pt"),
        text = element_text(size=24),
        axis.line = element_line(colour = "grey"))+
  ylim(-0.08,0.6)+
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000), limits=c(0,10000))

base=ggplot(Data[Data$P==0.99,])
Plot0.99=base+
  geom_hline(yintercept = 0, linetype=2)+
  geom_line(aes(x=time, y=(FYI+FXI)/2, color=paste0(as.factor(r), " cM")), size=1, alpha=0.8)+
  facet_grid(paste0("s=",s)~paste0("h=",h))+
  ylab("Inversion frequency")+
  xlab("Generation")+
  scale_color_manual("Distance between the inversion and the permanently heterozygous allele:", values=Col)+
  theme(panel.border = element_blank(), 
        legend.position = "top",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.placement = "outside",
        plot.margin = margin(10, 50, 10, 10, "pt"),
        text = element_text(size=24),
        axis.line = element_line(colour = "grey"))+
  ylim(-0.08,0.6)+
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000), limits=c(0,10000))

Plot=plot_grid(Plot0.8,Plot0.9,Plot0.99, nrow=3, labels = c('a', 'b', 'c'), label_size = 30)
save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS5.png",Plot, ncol=3, nrow=6)
save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS5.pdf",Plot, ncol=3, nrow=6)

### Figure S6
TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), r=double(), P=double(),
                      FXN=double(),FXI=double(),FX=double(),
                      FYN=double(),FYI=double(), FY=double(), 
                      W=double(),D=double(), q=double(), Des=double(),  InitFreq=double())
write.table(TableLinkR, "~/Paper/ModelSexChrom/V3/CleanDataset/DeterministicInvEvolution_2MT_VarLD.txt", append = F, quote=F, row.names = F)
P=0.80
for( Des in c(1,0.5)) #When Des=1, the inversion is introduced only on Y chromosomes. When Des=0.5, the inversion is introduce in equal proportion in X and Y chromosomes
{
  for ( InitFreq in c(0.01, 0.001, 0.0001)) #Variation in the initial frequency of the inversion
  {
    for (R in c(0.0, 0.001, 0.005))
    {
      for (h in c(0.01))
      {
        for (s in c(0.1, 0.01))
        {
          TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), r=double(), P=double(),
                                FXN=double(),FXI=double(),FX=double(),
                                FYN=double(),FYI=double(), FY=double(), 
                                W=double(),D=double(),q=double(), Des=double(),  InitFreq=double())

          u=1e-08
          n=2000000
          q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2))))
          m=floor(P*n*q)
          WNI=(q*(1-s) + (1-q)*(1-h*s))^m * (q*(1-h*s) + 1-q)^(n-m) #Fitness of individual heterozygous for the inversion.
          WNN=(1-2*q*(1-q)*h*s - q*q*s)^n #Fitness of individual homozygous for the absence of inversion.
          WII=(1-s)^m #Fitness of individual homozygous for the inversion.
          XRelatF=(1-Des)*InitFreq
          YRelatF=Des*InitFreq
          FXN=1.0-XRelatF
          FXI=XRelatF
          FYN=1.0-YRelatF
          FYI=YRelatF
          FY=FYN+FYI
          FX=FXN+FXI
          W=FXN*FYN*WNN + FXN*FYI*WNI + FXI*FYN*WNI + FXI*FYI*WII
          D=FXI*FYN - FXN*FYI
          FirstGen=c(1,h,s,u,n,R,P,FXN,FXI,FX,FYN,FYI,FY,W,D,q, Des, InitFreq)
          TableLinkR[nrow(TableLinkR)+1,]=FirstGen
          for (time in seq(2,10000,1))
          {
            FXI2=(FXI*FYI*WII + FXI*FYN*WNI*(1-R) + FXN*FYI*WNI*R)/W
            FXN2=(FXN*FYN*WNN + FXN*FYI*WNI*(1-R) + FXI*FYN*WNI*R)/W
            FYI2=(FYI*FXI*WII + FYI*FXN*WNI*(1-R) + FYN*FXI*WNI*R)/W
            FYN2=(FYN*FXN*WNN + FYN*FXI*WNI*(1-R) + FYI*FXN*WNI*R)/W
            FXI=FXI2
            FXN=FXN2
            FYI=FYI2
            FYN=FYN2
            FY=FYN+FYI
            FX=FXN+FXI
            D=FXI*FYN - FXN*FYI
            W=FXN*FYN*WNN + FXN*FYI*WNI + FXI*FYN*WNI + FXI*FYI*WII
            Gen=c(time,h,s,u,n,R,P,FXN,FXI,FX,FYN,FYI,FY,W,D,q,Des,InitFreq)
            TableLinkR[nrow(TableLinkR)+1,]=Gen
          }
          write.table(TableLinkR, "~/Paper/ModelSexChrom/V3/CleanDataset/DeterministicInvEvolution_2MT_VarLD.txt", append = T, quote=F, row.names = F, col.names = F)
        }
      }
    }
  }
}

Col=scales::viridis_pal(begin=0, end=0.8, option="A")(3)
Data=read.table("~/Paper/ModelSexChrom/V3/CleanDataset/DeterministicInvEvolution_2MT_VarLD.txt", stringsAsFactors = F, header = T)
Data$r=Data$r*100
Data=Data[Data$Des %in% c(0.5,1.0),]
Data$InitD=0
for (H in unique(Data$h))
{
  for (S in unique(Data$s))
  {
    for (I in unique(Data$InitFreq))
    {
      for (D in unique(Data$Des))
      {
        for (R in unique(Data$r))
        {
          InitDFix=Data[(Data$s==S & Data$h==H & Data$InitFreq==I & Data$Des==D & Data$time==1),]$D
          Data[(Data$s==S & Data$h==H & Data$InitFreq==I & Data$Des==D),]$InitD=InitDFix
        }
      }
    }
  }
}

base=ggplot(Data)
Plot0.8_VarFreq=base+
  geom_hline(yintercept = 0, linetype=2)+
  geom_line(aes(x=time, y=(FYI+FXI)/2, color=paste0(as.factor(r), " cM")), size=1, alpha=0.8)+
  facet_grid(paste0("s=",s)~ fct_relevel(paste0("Initial D=",-InitD), "D=0", "D=0.0001", "D=0.001", "D=0.01"))+
  ylab("Inversion frequency")+
  xlab("Generation")+
  scale_color_manual("Distance between the inversion and the permanently heterozygous allele:", values=Col)+
  theme(panel.border = element_blank(), 
        legend.position = "top",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.placement = "outside",
        plot.margin = margin(10, 50, 10, 10, "pt"),
        text = element_text(size=24),
        axis.line = element_line(colour = "grey"))+
  ylim(-0.08,0.6)+
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000), limits=c(0,10000))

base=ggplot(Data)
Plot0.8_VarLD=base+
  geom_hline(yintercept = 0, linetype=2)+
  geom_line(aes(x=time, y=-D, color=paste0(as.factor(r), " cM")), size=1, alpha=0.8)+
  facet_grid(paste0("s=",s)~ fct_relevel(paste0("Initial D=",-InitD), "D=0", "D=0.0001", "D=0.001", "D=0.01"))+
  ylab("Linkage desequilibrium")+
  xlab("Generation")+
  scale_color_manual("Distance between the inversion and the permanently heterozygous allele:", values=Col)+
  theme(panel.border = element_blank(), 
        legend.position = "top",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.placement = "outside",
        plot.margin = margin(10, 50, 10, 10, "pt"),
        text = element_text(size=24),
        axis.line = element_line(colour = "grey"))+
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000), limits=c(0,10000))

Plot=plot_grid(Plot0.8_VarFreq, Plot0.8_VarLD, nrow=2, labels=c("A", "B"), label_size = 30)
save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS6.png", Plot, ncol=3, nrow=4)

### Fig S7 ###
TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), r=double(),P=double(), 
                      FXN=double(),FXI=double(),FX=double(),
                      FYN=double(),FYI=double(), FY=double(),
                      FZN=double(),FZI=double(), FZ=double(),
                      FXIYI=double(),FXNYI=double(),FXIYN=double(),FXNYN=double(),
                      FYIZI=double(),FYNZI=double(),FYIZN=double(),FYNZN=double(),
                      FXIZI=double(),FXNZI=double(),FXIZN=double(),FXNZN=double(),
                      FT=double(),W=double(),q=double())

write.table(TableLinkR, "~/Paper/ModelSexChrom/V3/CleanDataset/DeterministicInvEvolution_3MT", append = F, quote=F, row.names = F)
P=0.8
for (R in c(0.0, 0.001, 0.01, 0.05, 0.1, 0.5))
{
  for (h in c(0.001, 0.01, 0.1))
  {
    for (s in c(0.01, 0.1))
    {
      TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(),r=double(), P=double(), 
                            FXN=double(),FXI=double(),FX=double(),
                            FYN=double(),FYI=double(), FY=double(),
                            FZN=double(),FZI=double(), FZ=double(),
                            FXIYI=double(),FXNYI=double(),FXIYN=double(),FXNYN=double(),
                            FYIZI=double(),FYNZI=double(),FYIZN=double(),FYNZN=double(),
                            FXIZI=double(),FXNZI=double(),FXIZN=double(),FXNZN=double(),
                            FT=double(),W=double(),q=double())
      u=1e-08
      n=2000000
      q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2))))
      m=floor(P*n*q)
      WNI=(q*(1-s) + (1-q)*(1-h*s))^m * (q*(1-h*s) + 1-q)^(n-m) #Fitness of individual heterozygous for the inversion.
      WNN=(1-2*q*(1-q)*h*s - q*q*s)^n #Fitness of individual homozygous for the absence of inversion.
      WII=(1-s)^m #Fitness of individual homozygous for the inversion.
      F_XN=0.99
      F_XI=0.01
      F_YN=1.00
      F_YI=0.00
      F_ZN=1.00
      F_ZI=0.00
      FY=1/3
      FX=1/3
      FZ=1/3
      FXN=F_XN*FX
      FXI=F_XI*FX
      FYN=F_YN*FY
      FYI=F_YI*FY
      FZN=F_ZN*FZ
      FZI=F_ZI*FZ
      FXIYI=FXI*(FYI/(FY+FZ))+FYI*(FXI/(FX+FZ))
      FXIYN=FXI*(FYN/(FY+FZ))+FYN*(FXI/(FX+FZ))
      FXNYI=FXN*(FYI/(FY+FZ))+FYI*(FXN/(FX+FZ))
      FXNYN=FXN*(FYN/(FY+FZ))+FYN*(FXN/(FX+FZ))
      FYIZI=FZI*(FYI/(FY+FX))+FYI*(FZI/(FX+FZ))
      FYIZN=FZI*(FYN/(FY+FX))+FYN*(FZI/(FX+FZ))
      FYNZI=FZN*(FYI/(FY+FX))+FYI*(FZN/(FX+FZ))
      FYNZN=FZN*(FYN/(FY+FX))+FYN*(FZN/(FX+FZ))
      FXIZI=FZI*(FXI/(FY+FX))+FXI*(FZI/(FY+FZ))
      FXIZN=FZN*(FXI/(FY+FX))+FXI*(FZN/(FY+FZ))
      FXNZI=FZI*(FXN/(FY+FX))+FXN*(FZI/(FY+FZ))
      FXNZN=FZN*(FXN/(FY+FX))+FXN*(FZN/(FY+FZ))
      Sum=FXIYI+FXNYI+FXIYN+FXNYN+
        FYIZI+FYNZI+FYIZN+FYNZN+
        FXIZI+FXNZI+FXIZN+FXNZN
      FT=FX+FY+FZ
      W=FXIYI*WII+FXNYI*WNI+FXIYN*WNI+FXNYN*WNN+
        FYIZI*WII+FYNZI*WNI+FYIZN*WNI+FYNZN*WNN+
        FXIZI*WII+FXNZI*WNI+FXIZN*WNI+FXNZN*WNN
      
      FirstGen=c(1,h,s,u,n,R,P,FXN,FXI,FX,FYN,FYI,FY,FZN,FZI,FZ,
                 FXIYI,FXNYI,FXIYN,FXNYN,FYIZI,FYNZI,FYIZN,FYNZN,FXIZI,FXNZI,FXIZN,FXNZN,
                 FT,W,q)
      TableLinkR[1,]=FirstGen
      for (time in seq(2,10000,1))
      {
        
        FXI=((TableLinkR$FXIYI[time-1]+TableLinkR$FXIZI[time-1])*WII+
               (TableLinkR$FXIYN[time-1]+TableLinkR$FXIZN[time-1])*(1-R)*WNI+
               (TableLinkR$FXNYI[time-1]+TableLinkR$FXNZI[time-1])*R*WNI)/(2*TableLinkR$W[time-1])
        FXN=((TableLinkR$FXNYN[time-1]+TableLinkR$FXNZN[time-1])*WNN+
               (TableLinkR$FXNYI[time-1]+TableLinkR$FXNZI[time-1])*(1-R)*WNI+
               (TableLinkR$FXIYN[time-1]+TableLinkR$FXIZN[time-1])*R*WNI)/(2*TableLinkR$W[time-1])
        FYI=((TableLinkR$FXIYI[time-1]+TableLinkR$FYIZI[time-1])*WII+
               (TableLinkR$FXNYI[time-1]+TableLinkR$FYIZN[time-1])*(1-R)*WNI+
               (TableLinkR$FXIYN[time-1]+TableLinkR$FYNZI[time-1])*R*WNI)/(2*TableLinkR$W[time-1])
        FYN=((TableLinkR$FXNYN[time-1]+TableLinkR$FYNZN[time-1])*WNN+
               (TableLinkR$FXIYN[time-1]+TableLinkR$FYNZI[time-1])*(1-R)*WNI+
               (TableLinkR$FXNYI[time-1]+TableLinkR$FYIZN[time-1])*R*WNI)/(2*TableLinkR$W[time-1])
        FZI=((TableLinkR$FXIZI[time-1]+TableLinkR$FYIZI[time-1])*WII+
               (TableLinkR$FXNZI[time-1]+TableLinkR$FYNZI[time-1])*(1-R)*WNI+
               (TableLinkR$FXIZN[time-1]+TableLinkR$FYIZN[time-1])*R*WNI)/(2*TableLinkR$W[time-1])
        FZN=((TableLinkR$FXNZN[time-1]+TableLinkR$FYNZN[time-1])*WNN+
               (TableLinkR$FXIZN[time-1]+TableLinkR$FYIZN[time-1])*(1-R)*WNI+
               (TableLinkR$FXNZI[time-1]+TableLinkR$FYNZI[time-1])*R*WNI)/(2*TableLinkR$W[time-1])
        FY=FYN+FYI
        FX=FXN+FXI
        FZ=FZN+FZI
        FXIYI=FXI*(FYI/(FY+FZ))+FYI*(FXI/(FX+FZ))
        FXIYN=FXI*(FYN/(FY+FZ))+FYN*(FXI/(FX+FZ))
        FXNYI=FXN*(FYI/(FY+FZ))+FYI*(FXN/(FX+FZ))
        FXNYN=FXN*(FYN/(FY+FZ))+FYN*(FXN/(FX+FZ))
        FYIZI=FZI*(FYI/(FY+FX))+FYI*(FZI/(FX+FZ))
        FYIZN=FZI*(FYN/(FY+FX))+FYN*(FZI/(FX+FZ))
        FYNZI=FZN*(FYI/(FY+FX))+FYI*(FZN/(FX+FZ))
        FYNZN=FZN*(FYN/(FY+FX))+FYN*(FZN/(FX+FZ))
        FXIZI=FZI*(FXI/(FY+FX))+FXI*(FZI/(FY+FZ))
        FXIZN=FZN*(FXI/(FY+FX))+FXI*(FZN/(FY+FZ))
        FXNZI=FZI*(FXN/(FY+FX))+FXN*(FZI/(FY+FZ))
        FXNZN=FZN*(FXN/(FY+FX))+FXN*(FZN/(FY+FZ))
        FT=FX+FY+FZ
        W=FXIYI*WII+FXNYI*WNI+FXIYN*WNI+FXNYN*WNN+
          FYIZI*WII+FYNZI*WNI+FYIZN*WNI+FYNZN*WNN+
          FXIZI*WII+FXNZI*WNI+FXIZN*WNI+FXNZN*WNN
        
        Gen=c(time,h,s,u,n,R,P,FXN,FXI,FX,FYN,FYI,FY,FZN,FZI,FZ,
              FXIYI,FXNYI,FXIYN,FXNYN,FYIZI,FYNZI,FYIZN,FYNZN,FXIZI,FXNZI,FXIZN,FXNZN,
              FT,W,q)
        TableLinkR[nrow(TableLinkR)+1,]=Gen
      }
      write.table(TableLinkR, "~/Paper/ModelSexChrom/V3/CleanDataset/DeterministicInvEvolution_3MT", append = T, quote=F, row.names = F, col.names = F)
    }
  }
}

Col=scales::viridis_pal(begin=0, end=0.8, option="A")(6)
Data=read.table("~/Paper/ModelSexChrom/V3/CleanDataset/DeterministicInvEvolution_3MT", stringsAsFactors = F, header = T)
FocS=c(0.01,0.1)
Data=Data[Data$s %in% FocS,]
FocH=c(0.001,0.01,0.1)
Data=Data[Data$h %in% FocH,]
Data$r=Data$r*100
FoccM=c(0, 0.1, 1, 5, 50)
Data=Data[Data$r %in% FoccM,]
base=ggplot(Data)
Plot0.8=base+
  geom_hline(yintercept = 0, linetype=2)+
  geom_line(aes(x=time, y=FXI, color=paste0(r, " cM")), size=1, alpha=0.8)+
  facet_grid(paste0("s=",s)~paste0("h=",h))+
  ylab("Inversion frequency")+
  xlab("Generation")+
  scale_color_manual("Distance between the inversion and the permanently heterozygous allele:", values=Col)+
  theme(panel.border = element_blank(),  
        legend.position = "top",
        legend.direction = "horizontal",  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.placement = "outside",
        plot.margin = margin(10, 50, 10, 10, "pt"),
        text = element_text(size=25),
        axis.line = element_line(colour = "grey"))+
  scale_y_continuous(breaks=c(0, 0.1,0.2,0.3,0.4), limits=c(0,0.5))+
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000))#+

save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS7.png", Plot0.8, ncol = 3, nrow=2)
save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS7.pdf", Plot0.8, ncol = 3, nrow=2)

### Figure S8 ###
TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), s2=double(), P=double(),
                      FXN=double(),FXI=double(),FX=double(), FYN=double(),FYI=double(), FY=double(),
                      W=double(),D=double(),q=double())

write.table(TableLinkR, "~/Paper/ModelSexChrom/V3/CleanDataset/DeterministicInvEvolution_Over.txt", append = F, quote=F, row.names = F)
P=0.8
R=0.0
for (s2 in c(0.0, 0.5, 0.9, 0.95, 0.99,1))
{
  for (h in c(0.001, 0.01, 0.1))
  {
    for (s in c(0.01, 0.05, 0.1, 0.5))
    {
      TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), s2=double(), P=double(),
                            FXN=double(),FXI=double(),FX=double(), FYN=double(),FYI=double(), FY=double(),
                            W=double(),D=double(),q=double())
      
      u=1e-08
      n=2000000
      q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2))))
      m=floor(P*n*q)
      WNI=(q*(1-s) + (1-q)*(1-h*s))^m * (q*(1-h*s) + 1-q)^(n-m) #Fitness of individual heterozygous for the inversion.
      WNN=(1-2*q*(1-q)*h*s - q*q*s)^n #Fitness of individual homozygous for the absence of inversion.
      WII=(1-s)^m #Fitness of individual homozygous for the inversion.
      FXN=0.49
      FXI=0.01
      FYN=0.50
      FYI=0.0
      FY=FYN+FYI
      FX=FXN+FXI
      W=WII*((1-s2)*(FXI*FXI + FYI*FYI) + 2*FXI*FYI) + 
        WNN*((1-s2)*(FXN*FXN + FYN*FYN) + 2*FXN*FYN) +
        2*WNI*((1-s2)*(FXI*FXN + FYI*FYN) + FXI*FYN + FXN*FYI)
      D=FXI*FYN - FXN*FYI
      FirstGen=c(1,h,s,u,n,s2,P,FXN,FXI,FX,FYN,FYI,FY,W,D,q)
      TableLinkR[1,]=FirstGen
      for (time in seq(2,10000,1))
      {
        FXI=(TableLinkR$FXI[time-1]*(WII*((1-s2)*TableLinkR$FXI[time-1] + TableLinkR$FYI[time-1]) +
                                       WNI*((1-s2)*TableLinkR$FXN[time-1] + TableLinkR$FYN[time-1])) - 
               R*WNI*TableLinkR$D[time-1])/TableLinkR$W[time-1]
        FYI=(TableLinkR$FYI[time-1]*(WII*((1-s2)*TableLinkR$FYI[time-1] + TableLinkR$FXI[time-1]) +
                                       WNI*((1-s2)*TableLinkR$FYN[time-1] + TableLinkR$FXN[time-1])) + 
               R*WNI*TableLinkR$D[time-1])/TableLinkR$W[time-1]
        FXN=(TableLinkR$FXN[time-1]*(WNN*((1-s2)*TableLinkR$FXN[time-1] + TableLinkR$FYN[time-1]) +
                                       WNI*((1-s2)*TableLinkR$FXI[time-1] + TableLinkR$FYI[time-1])) + 
               R*WNI*TableLinkR$D[time-1])/TableLinkR$W[time-1]
        FYN=(TableLinkR$FYN[time-1]*(WNN*((1-s2)*TableLinkR$FYN[time-1] + TableLinkR$FXN[time-1]) +
                                       WNI*((1-s2)*TableLinkR$FYI[time-1] + TableLinkR$FXI[time-1])) - 
               R*WNI*TableLinkR$D[time-1])/TableLinkR$W[time-1]
        FY=FYN+FYI
        FX=FXN+FXI
        W=WII*((1-s2)*(FXI*FXI + FYI*FYI) + 2*FXI*FYI) + 
          WNN*((1-s2)*(FXN*FXN + FYN*FYN) + 2*FXN*FYN) +
          2*WNI*((1-s2)*(FXI*FXN + FYI*FYN) + FXI*FYN + FXN*FYI)
        D=FXI*FYN - FXN*FYI
        Gen=c(time,h,s,u,n,s2,P,FXN,FXI,FX,FYN,FYI,FY,W,D,q)
        TableLinkR[nrow(TableLinkR)+1,]=Gen
      }
      write.table(TableLinkR, "~/Paper/ModelSexChrom/V3/CleanDataset/DeterministicInvEvolution_Over.txt", append = T, quote=F, row.names = F, col.names = F)
    }
  }
}

Col=scales::viridis_pal(begin=0, end=0.8)(6)
Data=read.table("~/Paper/ModelSexChrom/V3/CleanDataset/DeterministicInvEvolution_Over.txt", stringsAsFactors = F, header = T)
FocS=c(0.01,0.1)
Data=Data[Data$s %in% FocS,]
base=ggplot(Data)
Plot0.8=base+
  geom_hline(yintercept = 0, linetype=2)+
  geom_line(aes(x=time, y=FXI, color=fct_rev(fct_infreq(as.factor(s2)))), size=1, alpha=0.8)+
  facet_grid(paste0("s=",s)~paste0("h=",h))+
  ylab("Inversion frequency")+
  xlab("Generation")+
  scale_color_manual("s2=", values=Col)+
  theme(panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.placement = "outside",
        plot.margin = margin(10, 50, 10, 10, "pt"),
        text = element_text(size=25),
        axis.line = element_line(colour = "grey"))+
  scale_y_continuous(breaks=c(0, 0.25, 0.5), limits=c(0,0.6))+
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000))#+

save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS8.png", Plot0.8, ncol = 3, nrow=4)
save_plot("~/Paper/ModelSexChrom/V3/Plot/FigS8.pdf", Plot0.8, ncol = 3, nrow=4)
