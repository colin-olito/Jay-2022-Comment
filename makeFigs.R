################################################################
#'  R code to reproduce figures for:
#' 
#'  Olito, C., and B Charlesworth. 2023. Comment on Jay et al.~(2022). 
#'  PLoS Biology, doi: XXX
#'
#'  Author: Colin Olito
#'
#'  NOTES:  To reproduce the figures for Jay et al. (2022), and
#' 			accompanying figures without excluding inversions 
#'			lost by generation 20, run the following R code,
#'			either interactively, or using source() in terminal.
#'			This will generate a figure directory, and save all
#' 			figures as .pdf files there.
#'		



######################################
#' Create figures directories if they
#' do not already exist
figuresDirectoriesExist  <-  dir.exists("./figures")

if(!figuresDirectoriesExist) {
	dir.create("./figures")
}

source('./R/Recreate-and-Compare-Jay-etal-Figs.R')



######################################
#' Approximation Figures
source('./R/functions-figures.R')

toPdf(upperBoundApproxFig(JayDat = JayDat), 
			figPath(name='upperBoundApproxFig.pdf'), width=10, height=7)
embed_fonts(figPath(name='upperBoundApproxFig.pdf'))

