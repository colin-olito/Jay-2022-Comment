#' R code to (1) faithfully recreate essential figures
#' from Jay et al. (2022) and also (2) corresponding 
#' figures where inversions that have gone extinct by
#' generation 20 are not excluded prior to plotting. 

#######################################################
#' Code to reproduce figs from Jay et al. is unaltered.
#' Compare with downloaded branch at time of writing
#' in directory ./MutationShelteringTheory-main
#' 
#' Only filenames have been changed to correctly
#' reference within this repo's directory tree.
#' 
#' Original code can be downloaded from
#' GitHub: https://github.com/PaulYannJay/MutationShelteringTheory/tree/main
#' File:   DeleteriousMutationSheltering_Figures_Simulations.R
#' Commit: c2b2f25
#' Data:   https://figshare.com/articles/dataset/Model_of_sex-chromosome_Evolution_-_datasets/19961033

#################
# Dependencies
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
# library(tidyverse) # I get compile error for dependency package `ragg`. However, can run the code here without `tidyverse`.
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


#############################################################
##############################################################
#' Re-create Fig.3c as it appears in Jay et al. (2022)
# Simul=read.table(paste("~/Paper/ModelSexChrom/V3/CleanDataset/InversionTrajectories_N=1000_Fig3-S13-S15-S19.txt",sep=""), stringsAsFactors = F) #File containing all simulation with N=1000
Simul=read.table(paste("./ModelSexChrom/InversionTrajectories_N=1000_Fig3-S13-S15-S19.txt",sep=""), stringsAsFactors = F) #File containing all simulation with N=1000
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
  ylab("Conditional fraction of inversions fixed")+
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
  labs(title = "Fraction of inversions fixed after 10,000 generations (excluding early extinctions)")

PlotPropInvSpread
ggsave(filename="./figures/Jay-etal-2022-Fig3c.pdf", 
	   width = 20, height = 10, units = "cm")


##############################################################
#' Re-create Fig.3c without excluding inversions lost by
#' generation 20

# CREATE DATA SUMMARY WITHOUT EXCLUDING LostEarly INVERSIONS 
# Also calculate poisson error c.i.
DataSummaryAll=Summary %>% group_by(N,u,r,h,s,InvSize, Position) %>% summarise(ProbSpread=mean(StateCode), MeanSegMut=mean(MinSegMut), nFix=sum(StateCode), nSeg=sum(State=="Segregating"), upper=qgamma(0.975, shape=nFix, rate=1)/10000, lower=qgamma(0.025, shape=nFix, rate=1)/10000) # For each set of parameter, compute the fraction of mutation fixed (only for not mutation-free inversion)
DataSummaryAll$u=factor(DataSummaryAll$u, labels = c("mu==1 %*% 10^{-09}","mu==1 %*% 10^{-08}"), )


#' Dummy data to plot benchmarks for 1/Ne
neutralExpect  <-  c(1/(2*DataSummaryAll$N[1]), (2/DataSummaryAll$N[1]))
columns        <-  rep(c("Autosome", "Y"), times=length(levels(DataSummaryAll$u)))
rows           <-  rep(levels(DataSummaryAll$u), each=2)
hLineData      <-  data.frame("Position"=columns, "u"=rows, neutralExpect)
hLineData$u    <- factor(hLineData$u, levels = levels(DataSummaryAll$u))

# Make the plot
BaseAll=ggplot(DataSummaryAll[(DataSummaryAll$InvSize == n & DataSummaryAll$s<0),], aes( y=ProbSpread)) #Plot the result
PlotPropInvSpread_Incl20=BaseAll+geom_point(aes(x=h,color=as.factor(s)), size=1)+
  geom_line(aes(x=h, color=as.factor(s)), size=1)+
  geom_errorbar(aes(x=h, color=as.factor(s), ymin = lower, ymax = upper), width = 0.0125)+
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
  facet_grid(u~Position, scales = "free", labeller = label_parsed)+
  geom_hline(data=hLineData, aes(yintercept= neutralExpect), color="dodgerblue", linewidth=1, linetype='dashed')+
    labs(title = "Fraction of inversions fixed after 10,000 generations")

#' Show plot, save as pdf
PlotPropInvSpread_Incl20
ggsave(filename="./figures/Jay-etal-2022-Fig3c-AllInversions.pdf", 
	   width = 20, height = 11, units = "cm")


#' Figure 3c, including inversions lost by generation 20
#' Reversing grid arrangment so Autosomal and Y inversions
#' are on the same y-axis scale

# Make the plot
BaseAll=ggplot(DataSummaryAll[(DataSummaryAll$InvSize == n & DataSummaryAll$s<0),], aes( y=ProbSpread)) #Plot the result
PlotPropInvSpread_Incl20_newGrid=BaseAll+geom_point(aes(x=h,color=as.factor(s)), size=1)+
  geom_line(aes(x=h, color=as.factor(s)), size=1)+
  geom_errorbar(aes(x=h, color=as.factor(s), ymin = lower, ymax = upper), width = 0.0125)+
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
  facet_grid(Position~u, scales = "free_y", labeller = label_parsed)+
  geom_hline(data=hLineData, aes(yintercept= neutralExpect), color="dodgerblue", linewidth=1, linetype='dashed')+
    labs(title = "Fraction of inversions fixed after 10,000 generations")

#' Show plot, save as pdf
PlotPropInvSpread_Incl20_newGrid
ggsave(filename="./figures/Jay-etal-2022-Fig3c-AllInversions-NewGrid.pdf", 
     width = 20, height = 11, units = "cm")


#' How many of the parameter sets from Fig 3c 
#' had no autosomal inversions go to fixation?
Fig3cData  <-  subset(DataSummaryAll, (Position=="Autosome" & InvSize==2000000))
#' Calculate proportion of these datasets where 
#' autosomal fixations > 0 
1-sum(Fig3cData$nFix==0)/length(Fig3cData$nFix) 
#' Fig 3a presented one of the minority cases where no autosomal 
#' inversions went to fixation. Autosomal fixations occurred in
#'  ~86% of the parameter sets used in 3c.



#############################################################
#############################################################
#' Re-create Fig.S15 as it appears in Jay et al. (2022)
#' 
U="mu==1 %*% 10^{-08}"
DataSummary$InvSize=as.factor(paste0("Inversion size=",DataSummary$InvSize))
DataSummary$InvSize=relevel(DataSummary$InvSize,"Inversion size=500000")
# Make the plot
Base=ggplot(DataSummary[(DataSummary$u == U & DataSummary$s<0),], aes( y=ProbSpread))
PlotPropInvSpread_AllInv=Base+geom_point(aes(x=h,color=as.factor(s)), size=4)+
  geom_line(aes(x=h, color=as.factor(s)), size=2)+
  scale_color_manual("Selection coefficient (s)", values=Col)+
  scale_size_manual("Inversion size", labels=c("0.5Mb", "1Mb", "2Mb", "5Mb"), values=c(1,2,3,4))+
  ylab("Fraction of inversions fixed after 10,000 generations (excl. early extinctions)")+
  xlab("h")+
  ThemeSobr+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))+
  facet_grid(InvSize ~ Position)

#' Show plot, save as pdf
PlotPropInvSpread_AllInv
ggsave(filename="./figures/Jay-etal-2022-FigS15.pdf", 
	   width = 35, height = 30, units = "cm")



########################################################################
#' Re-create Fig.S15 without excluding inversions lost by generation 20
#' Same as figure 3C but with different inversion size


# Make InvSize in complete dataset a factor for plotting
DataSummaryAll$InvSize=as.factor(paste0("Inversion size=",DataSummaryAll$InvSize))
DataSummaryAll$InvSize=relevel(DataSummaryAll$InvSize,"Inversion size=500000")

#' Dummy data to plot benchmarks for 1/Ne
neutralExpect  <-  c(1/(2*DataSummaryAll$N[1]), (2/DataSummaryAll$N[1]))
columns        <-  rep(c("Autosome", "Y"), times=length(levels(DataSummaryAll$InvSize)))
rows           <-  rep(levels(DataSummaryAll$InvSize), each=2)
hLineData      <-  data.frame("Position"=columns, "InvSize"=rows, neutralExpect)
hLineData$InvSize <- factor(hLineData$InvSize, levels = levels(DataSummaryAll$InvSize))

#' Make the plot
Base=ggplot(DataSummaryAll[(DataSummaryAll$u == U & DataSummaryAll$s<0),], aes( y=ProbSpread))
PlotPropInvSpread_AllInv_Incl20=Base+geom_point(aes(x=h,color=as.factor(s)), size=4)+
  geom_line(aes(x=h, color=as.factor(s)), size=2)+
  geom_errorbar(aes(x=h, color=as.factor(s), ymin = lower, ymax = upper), width = 0.025)+
  scale_color_manual("Selection coefficient (s)", values=Col)+
  scale_size_manual("Inversion size", labels=c("0.5Mb", "1Mb", "2Mb", "5Mb"), values=c(1,2,3,4))+
  ylab("Fraction of inversions fixed after 10,000 generations")+
  xlab("h")+
  ThemeSobr+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))+
  facet_grid(Position ~ InvSize, scales = "free_y")+
  geom_hline(data=hLineData, aes(yintercept= neutralExpect), color="dodgerblue", linewidth=1, linetype='dashed')

#' Show plot, save as pdf
PlotPropInvSpread_AllInv_Incl20
ggsave(filename="./figures/Jay-etal-2022-FigS15-AllInversions.pdf", 
	   width = 50, height = 20, units = "cm")





#############################################################
#############################################################
#' Re-create Fig.S16 as it appears in Jay et al. (2022)
#' Same as before but for N=10,000 ###

#Simul10k=read.table(paste("~/Paper/ModelSexChrom/V3/CleanDataset/InversionTrajectories_N=10000_FigS11-S12-S16.txt",sep=""), stringsAsFactors = F)
Simul10k=read.table(paste("./ModelSexChrom/InversionTrajectories_N=10000_FigS11-S12-S16.txt",sep=""), stringsAsFactors = F)
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
DataSummary10k= SumNoLostEarly10k %>% group_by(N,u,r,h,s,InvSize, Position) %>% summarise(ProbSpread=mean(StateCode), MeanSegMut=mean(MinSegMut))
Col=scales::viridis_pal(begin=0.0, end=0.8, option="A")(5)
options(scipen=0)
DataSummary10k$InvSize=as.factor(paste0("Inversion size=",DataSummary10k$InvSize))
DataSummary10k$InvSize=relevel(DataSummary10k$InvSize,"Inversion size=500000")

#' Make Plot
Base=ggplot(DataSummary10k[(DataSummary10k$s<0),], aes( y=ProbSpread))
PlotPropInvSpread_AllInv_10k=Base+geom_point(aes(x=h,color=as.factor(s)), size=4)+
  geom_line(aes(x=h, color=as.factor(s)), size=2)+
  scale_color_manual("Selection coefficient (s)", values=Col)+
  scale_size_manual("Inversion size", labels=c("0.5Mb", "1Mb", "2Mb", "5Mb"), values=c(1,2,3,4))+
  ylab("Fraction of inversions fixed after 10,000 generations (excl. early extinctions)")+
  xlab("h")+
  ThemeSobr+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))+
  facet_grid(InvSize ~Position)

#' Show plot, save as pdf
PlotPropInvSpread_AllInv_10k
ggsave(filename="./figures/Jay-etal-2022-10k-FigS16.pdf", 
	   width = 35, height = 30, units = "cm")


########################################################################
#' Re-create Fig.S16 without excluding inversions lost by generation 20
#' Same as before but for N=10,000 ###

# CREATE DATA SUMMARY WITHOUT EXCLUDING LostEarly INVERSIONS 
DataSummary10kAll= Summary10k %>% group_by(N,u,r,h,s,InvSize, Position) %>% summarise(ProbSpread=mean(StateCode), MeanSegMut=mean(MinSegMut), nFix=sum(StateCode), upper=qgamma(0.975, shape=nFix, rate=1)/10000, lower=qgamma(0.025, shape=nFix, rate=1)/10000)
DataSummary10kAll$InvSize=as.factor(paste0("Inversion size=",DataSummary10kAll$InvSize))
DataSummary10kAll$InvSize=relevel(DataSummary10kAll$InvSize,"Inversion size=500000")

#' Dummy data to plot benchmarks for 1/Ne
neutralExpect  <-  c(1/(2*DataSummary10kAll$N[1]), (2/DataSummary10kAll$N[1]))
columns        <-  rep(c("Autosome", "Y"), times=length(levels(DataSummary10kAll$InvSize)))
rows           <-  rep(levels(DataSummary10kAll$InvSize), each=2)
hLineData      <-  data.frame("Position"=columns, "InvSize"=rows, neutralExpect)
hLineData$InvSize <- factor(hLineData$InvSize, levels = levels(DataSummary10kAll$InvSize))

#' Make Plot
Base=ggplot(DataSummary10kAll[(DataSummary10kAll$s<0),], aes( y=ProbSpread))
PlotPropInvSpread_AllInv_10k_Incl20=Base+geom_point(aes(x=h,color=as.factor(s)), size=4)+
  geom_line(aes(x=h, color=as.factor(s)), size=2)+
  geom_errorbar(aes(x=h, color=as.factor(s), ymin = lower, ymax = upper), width = 0.025)+
  scale_color_manual("Selection coefficient (s)", values=Col)+
  scale_size_manual("Inversion size", labels=c("0.5Mb", "1Mb", "2Mb", "5Mb"), values=c(1,2,3,4))+
  ylab("Fraction of inversions fixed after 10,000 generations")+
  xlab("h")+
  ThemeSobr+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))+
  facet_grid(Position ~ InvSize, scales = "free_y")+
  geom_hline(data=hLineData, aes(yintercept= neutralExpect), color="dodgerblue", linewidth=1, linetype='dashed')

#' Show plot, save as pdf
PlotPropInvSpread_AllInv_10k_Incl20
ggsave(filename="./figures/Jay-etal-2022-FigS16-10k-AllInversions.pdf", 
	   width = 50, height = 20, units = "cm")



#' Modified version of Fig.S16 to zoom in on results for h > 0.01
Base=ggplot(DataSummary10kAll[(DataSummary10kAll$s<0),], aes( y=ProbSpread))
PlotPropInvSpread_AllInv_10k_Incl20_Zoom=Base+geom_point(aes(x=h,color=as.factor(s)), size=4)+
  geom_line(aes(x=h, color=as.factor(s)), size=2)+
  geom_errorbar(aes(x=h, color=as.factor(s), ymin = lower, ymax = upper), width = 0.025)+
  scale_color_manual("Selection coefficient (s)", values=Col)+
  scale_size_manual("Inversion size", labels=c("0.5Mb", "1Mb", "2Mb", "5Mb"), values=c(1,2,3,4))+
  ylab("Fraction of inversions fixed after 10,000 generations")+
  xlab("h")+
  ThemeSobr+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))+
  coord_cartesian(ylim = c(0, (3*(2/DataSummary10kAll$N[1])))) +
  facet_grid(Position ~ InvSize, scales = "free_y")+
  geom_hline(data=hLineData, aes(yintercept= neutralExpect), color="dodgerblue", linewidth=1, linetype='dashed')

#' Show plot, save as pdf
PlotPropInvSpread_AllInv_10k_Incl20_Zoom
ggsave(filename="./figures/Jay-etal-2022-FigS16-10k-AllInversions-Zoom.pdf", 
     width = 50, height = 20, units = "cm")





#############################################################
#############################################################
#' Re-create Fig.S17 as it appears in Jay et al. (2022)
#' 
#' Same as Figure S15 but considering that mutation have 
#' their selection coefficient drawn from a gamma distribution (shape = 0.2)
#' and dominance coefficient drawn from a uniform distribution

#SimulLamb=read.table(paste("~/Paper/ModelSexChrom/V3/CleanDataset/MutationLambdaDistribution_InversionTrajectories_FigS17.txt",sep=""), stringsAsFactors = F)
SimulLamb=read.table(paste("./ModelSexChrom/MutationLambdaDistribution_InversionTrajectories_FigS17.txt",sep=""), stringsAsFactors = F)
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
DataSummaryLamb=SumNoLostEarly %>% group_by(N,u,r,s,InvSize, Position) %>% summarise(ProbSpread=mean(StateCode), MeanSegMut=mean(MinSegMut))

options(scipen=0)
Col=scales::viridis_pal(begin=0.2, end=0.6, option="A")(2)
DataSummaryLamb$InvSize=as.factor(paste0("Inversion size=",DataSummaryLamb$InvSize))
DataSummaryLamb$InvSize=relevel(DataSummaryLamb$InvSize,"Inversion size=500000")
DataSub=DataSummaryLamb[(DataSummaryLamb$s<0),]
DataSub$s=-DataSub$s

#' Make Plot
Base=ggplot(DataSub, aes( y=ProbSpread))
PlotPropInvSpread_AllInv_Lamb=Base+geom_point(aes(x=s, color=as.factor(u)), size=4)+
  geom_line(aes(x=s, color=as.factor(u)), size=2)+
  scale_color_manual("mutation rate (u)", values=Col)+
  scale_x_log10()+
  ylab("Fraction of inversions fixed after 10,000 generations (excl. early extinctions)")+
  xlab("mean(s)")+
  ThemeSobr+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))+
  facet_grid(InvSize ~ Position)

#' Show plot, save as pdf
PlotPropInvSpread_AllInv_Lamb
ggsave(filename="./figures/Jay-etal-2022-FigS17-shDistribution.pdf", 
	   width = 30, height = 40, units = "cm")





########################################################################
#' Re-create Fig.S17 without excluding inversions lost by generation 20
#' 
#' Same as Figure S15 but considering that mutation have 
#' their selection coefficient drawn from a gamma distribution (shape = 0.2)
#' and dominance coefficient drawn from a uniform distribution

# CREATE DATA SUMMARY WITHOUT EXCLUDING LostEarly INVERSIONS 
DataSummaryLambAll=summarySub_BLamb %>% group_by(N,u,r,s,InvSize, Position) %>% summarise(ProbSpread=mean(StateCode), MeanSegMut=mean(MinSegMut), nFix=sum(StateCode), upper=qgamma(0.975, shape=nFix, rate=1)/10000, lower=qgamma(0.025, shape=nFix, rate=1)/10000)
DataSummaryLambAll$InvSize=as.factor(paste0("Inversion size=",DataSummaryLambAll$InvSize))
DataSummaryLambAll$InvSize=relevel(DataSummaryLambAll$InvSize,"Inversion size=500000")
DataSubAll=DataSummaryLambAll[(DataSummaryLambAll$s<0),]
DataSubAll$s=-DataSubAll$s

#' Dummy data to plot benchmarks for 1/Ne
neutralExpect  <-  c(1/(2*DataSubAll$N[1]), (2/DataSubAll$N[1]))
columns        <-  rep(c("Autosome", "Y"), times=length(levels(DataSubAll$InvSize)))
rows           <-  rep(levels(DataSubAll$InvSize), each=2)
hLineData      <-  data.frame("Position"=columns, "InvSize"=rows, neutralExpect)
hLineData$InvSize <- factor(hLineData$InvSize, levels = levels(DataSubAll$InvSize))

# Make Plot
Base=ggplot(DataSubAll, aes( y=ProbSpread))
PlotPropInvSpread_AllInv_Lamb_noExcl=Base+geom_point(aes(x=s, color=as.factor(u)), size=4)+
  geom_line(aes(x=s, color=as.factor(u)), size=2)+
  geom_errorbar(aes(x=s, color=as.factor(u), ymin = lower, ymax = upper), width = 0.1)+
  scale_color_manual("mutation rate (u)", values=Col)+
  scale_x_log10()+
  ylab("Fraction of inversions fixed after 10,000 generations")+
  xlab("mean(s)")+
  ThemeSobr+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))+
  facet_grid(Position ~ InvSize, scales = "free_y")+
  geom_hline(data=hLineData, aes(yintercept= neutralExpect), color="dodgerblue", linewidth=1, linetype='dashed')

#' Show plot, save as .pdf
PlotPropInvSpread_AllInv_Lamb_noExcl
ggsave(filename="./figures/Jay-etal-2022-FigS17-shDistribution-AllInversions.pdf", 
	   width = 50, height = 20, units = "cm")
