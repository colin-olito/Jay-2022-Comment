###############
# DEPENDENCIES
###############
library(extrafont)
library(fontcm)
loadfonts(quiet = TRUE)
library(wesanderson)
library(plotrix)

# Packages from Jay code to correcly load data
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

#######################
# AUXILLIARY FUNCTIONS
#######################

toPdf <- function(expr, filename, ...) {
  toDev(expr, pdf, filename, ...)
}

figPath  <-  function(name) {
  file.path('./figures', name)
}

toDev <- function(expr, dev, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf('Creating %s\n', filename))
    dev(filename, family="Times New Roman", ...)
#    dev(filename, family='Arial', ...)
    on.exit(dev.off())
    eval.parent(substitute(expr))
}



####################
# PLOTTING FUNCTIONS
####################

#' Plot text or points according to relative axis position.
#'
#' @title Plot text or points according to relative axis position
#' @param px Relative x-axis position (in proportion) where character is to be plotted.
#' @param py Relative y-axis position (in proportion) where character is to be plotted.
#' @param lab Plotted text. Works if argument \code{\link[graphics]{text}} is \code{TRUE}.
#' @param adj See argument of same name in R base function \code{\link[graphics]{par}}.
#' @param text Logical. Should text or points be plotted?
#' @param log Used if the original plot uses the argument log, e.g. \code{log='x'}, \code{log='y'} or \code{log='xy'}.
#' @param ... Additional arguments to R base function \code{\link[graphics]{text}}.
#' @export
proportionalLabel <- function(px, py, lab, adj=c(0, 1), text=TRUE, log=FALSE, ...) {
    usr  <-  par('usr')
    x.p  <-  usr[1] + px*(usr[2] - usr[1])
    y.p  <-  usr[3] + py*(usr[4] - usr[3])
    if(log=='x') {
        x.p<-10^(x.p)
    }
    if(log=='y') {
        y.p<-10^(y.p)
    }
    if(log=='xy') {
        x.p<-10^(x.p)
        y.p<-10^(y.p)
    }
    if(text){
        text(x.p, y.p, lab, adj=adj, ...)
    } else {
        points(x.p, y.p, ...)
    }
}



proportionalArrows <- function(px1, py1, px2, py2, adj=c(0, 1), log=FALSE, length=length, ...) {
    usr  <-  par('usr')
    x.p1  <-  usr[1] + px1*(usr[2] - usr[1])
    y.p1  <-  usr[3] + py1*(usr[4] - usr[3])
    x.p2  <-  usr[1] + px2*(usr[2] - usr[1])
    y.p2  <-  usr[3] + py2*(usr[4] - usr[3])
    if(log=='x') {
        x.p1  <-  10^(x.p1)
        x.p2  <-  10^(x.p2)
    }
    if(log=='y') {
        y.p1  <-  10^(y.p1)
        y.p2  <-  10^(y.p2)
    }
    if(log=='xy') {
        x.p1  <-  10^(x.p1)
        y.p1  <-  10^(y.p1)
        x.p2  <-  10^(x.p2)
        y.p2  <-  10^(y.p2)
    }
    arrows(x0=x.p1, y0=y.p1, x1=x.p2, y1=y.p2, length=length,...)
}

#' Draw equally-spaced white lines on plot window.
#'
#' @title Equally-spaced white lines on plot window
#' @param ... Additional arguments to internal function \code{\link{proportionalLabel}}.
#' @author Diego Barneche
#' @export
plotGrid  <-  function(lineCol='white',...) {
    proportionalLabel(rep(0.2, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.4, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.6, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.8, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.2, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.4, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.6, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.8, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
}


#' Internal. Create nice rounded numbers for plotting.
#'
#' @title Rounded numbers for plotting
#' @param value A numeric vector.
#' @param precision Number of rounding digits.
#' @return A character vector.
#' @author Diego Barneche.
rounded  <-  function(value, precision=1) {
  sprintf(paste0('%.', precision, 'f'), round(value, precision))
}


#' Creates transparent colours
#'
#' @title Creates transparent colours
#' @param col Colour.
#' @param opacity equivalent to alpha transparency parameter
#' @export
transparentColor <- function(col, opacity=0.5) {
    if (length(opacity) > 1 && any(is.na(opacity))) {
        n        <-  max(length(col), length(opacity))
        opacity  <-  rep(opacity, length.out=n)
        col      <-  rep(col, length.out=n)
        ok       <-  !is.na(opacity)
        ret      <-  rep(NA, length(col))
        ret[ok]  <-  Recall(col[ok], opacity[ok])
        ret
    } else {
        tmp  <-  col2rgb(col)/255
        rgb(tmp[1,], tmp[2,], tmp[3,], alpha=opacity)
    }
}

fibonacci.scale  <-  function(n) {
    fibs  <-  c(0,1)
    for(i in 2:n) {
        fibs  <-  c(fibs, (fibs[i] + fibs[i-1]))
    }
    (fibs/max(fibs))[-2]
}



####################
# LOAD JAY DATA
####################
load_Jay_Data  <-  function() {

    cat('Loading data from Jay et al. (2022) fig. 3c\n')
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
    DataSummaryAll=Summary %>% group_by(N,u,r,h,s,InvSize, Position) %>% summarise(ProbSpread=mean(StateCode), MeanSegMut=mean(MinSegMut), nFix=sum(StateCode), nSeg=sum(State=="Segregating"), upper=qgamma(0.975, shape=nFix, rate=1)/10000, lower=qgamma(0.025, shape=nFix, rate=1)/10000) # For each set of parameter, compute the fraction of mutation fixed (only for not mutation-free inversion)
    DataSummaryAll$u=factor(DataSummaryAll$u, labels = c("mu==1 %*% 10^{-09}","mu==1 %*% 10^{-08}"), )

    return(DataSummaryAll)
}

JayDat  <-  load_Jay_Data()

#################
# APPROXIMATIONS
#################

wkSel_upper  <-  function(n, theta, h, N, s) {
    rel  <-  1 + (1/144)*n*theta*(1 + 2*h^2)*(1 + n*s*theta*(1 + 2*h))*(N*s)^2
    (2/N)*rel
}
recess_upper  <-  function(n, theta, u, N, s) {
    rel  <-  1 + 0.00923*n*theta * (1 + 2*n*u)*(2*N*s)^(1/2)
    (2/N)*rel
}







#####################
# Approximations Fig

upperBoundApproxFig  <-  function(JayDat = JayDat) {

    # Colors
    COLS  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

    # set plot layout
    layout.mat   <-  matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE)
    layout       <-  layout(layout.mat,respect=TRUE)

    # Subset simulation data for h = 0
    RecDat  <-  subset(JayDat, h == 0 & Position == 'Y' & InvSize == 2000000)
    RecDat$mu  <-  rep(0,length=nrow(RecDat))
    RecDat$mu[RecDat$u == 'mu==1 %*% 10^{-08}']  <-  1e-8
    RecDat$mu[RecDat$u == 'mu==1 %*% 10^{-09}']  <-  1e-9
    RecDat$theta  <-  4*RecDat$N*RecDat$mu
    # Calculate approximate upper bounds
    ss  <-  0:100/200
    RecDat$theta    <-  4*RecDat$N*RecDat$mu
    thetas          <-  unique(RecDat$theta)
    RecUpperBound1  <-  recess_upper(n = RecDat$InvSize[1], 
                                     theta = thetas[1], 
                                     u = unique(RecDat$mu)[1],
                                     N = RecDat$N[1], 
                                     s = ss
                                     )
    RecUpperBound2  <-  recess_upper(n = RecDat$InvSize[1], 
                                     theta = thetas[2], 
                                     u = unique(RecDat$mu)[2],
                                     N = RecDat$N[1], 
                                     s = ss
                                     )
    # Generate Plot
     # Panel (A) Completely Recessive deleterious mutations
    par(omi=c(0.5, 0.5, 0.5, 0.5), mar = c(4,4,4,1), bty='o', xaxt='s', yaxt='s')
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0,0.06), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot Upper Bounds
        lines(RecUpperBound1 ~ ss, lty=1, lwd=3, col=COLS[3])
        points(RecDat$ProbSpread[RecDat$mu==1e-9] ~ abs(RecDat$s[RecDat$mu==1e-9]), pch=21, col=1, bg=transparentColor(COLS[3], opacity=0.6))
        lines(RecUpperBound2 ~ ss, lty=1, lwd=3, col=COLS[2])
        points(RecDat$ProbSpread[RecDat$mu==1e-8] ~ abs(RecDat$s[RecDat$mu==1e-8]), pch=21, col=1, bg=transparentColor(COLS[2], opacity=0.6))
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.5,  1.1,   expression(paste("Recessive mutations (", italic(h)==0,")")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.22, 0.5,  expression(paste("Inversion fixation probability")), cex=1.3, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25, expression(paste("Selection (", italic(s), ")")), cex=1.3, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2]*0.35,
               y       =  usr[4]*0.99,
               legend  =  c(expression(italic(u)==10^-8),
                            expression(italic(u)==10^-9)),
               lty     =  1,
               lwd     =  3,
               col     =  c(COLS[2], COLS[3]),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


    # Subset simulation data for s = 0.001
    wkSelDat     <-  subset(JayDat, s == -0.001 & h !=0 & Position == 'Y' & InvSize == 2000000)
    wkSelDat$mu  <-  rep(0,length=nrow(wkSelDat))
    wkSelDat$mu[wkSelDat$u == 'mu==1 %*% 10^{-08}']  <-  1e-8
    wkSelDat$mu[wkSelDat$u == 'mu==1 %*% 10^{-09}']  <-  1e-9
    wkSelDat$theta  <-  4*wkSelDat$N*wkSelDat$mu
    # Calculate approximate upper bounds
    thetas            <-  unique(wkSelDat$theta)
    hs                <-  1:100/200
    wkSelUpperBound1  <-  wkSel_upper(n = wkSelDat$InvSize[1], 
                                     theta = thetas[1], 
                                     h = hs, 
                                     N = wkSelDat$N[1], 
                                     s = wkSelDat$s[1]
                                     )
    wkSelUpperBound2  <-  wkSel_upper(n = wkSelDat$InvSize[1], 
                                     theta = thetas[2], 
                                     h = hs, 
                                     N = wkSelDat$N[1], 
                                     s = wkSelDat$s[1]
                                     )

     # Panel (B) Weak selection, partially recessive deleterious mutations
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0,0.007), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Plot Upper Bounds
        lines(wkSelUpperBound1 ~ hs, lty=1, lwd=3, col=COLS[3])
        points(wkSelDat$ProbSpread[wkSelDat$mu==1e-9] ~ wkSelDat$h[wkSelDat$mu==1e-9], pch=21, col=1, bg=transparentColor(COLS[3], opacity=0.6))
        lines(wkSelUpperBound2 ~ hs, lty=1, lwd=3, col=COLS[2])
        points(wkSelDat$ProbSpread[wkSelDat$mu==1e-8] ~ wkSelDat$h[wkSelDat$mu==1e-8], pch=21, col=1, bg=transparentColor(COLS[2], opacity=0.6))
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.5,   1.1,  expression(paste("Weak selection (", italic(Ns)==1,")")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5, -0.25, expression(paste("Dominance (", italic(h), ")")), cex=1.3, adj=c(0.5, 0.5), xpd=NA)        

}