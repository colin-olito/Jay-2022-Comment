# Comment on Jay et al. (2022): Ignoring inversion extinctions paints a misleading picture of sex chromosome evolution.

## Overview

This is a GitHub repository for R code accompanying a comment article addressing the presentation of simulation results from Jay et al. (2022) "*Sheltering of deleterious mutations explains the stepwise extension of recombination suppression on sex chromosomes and other supergenes*" (doi: [`10.1371/journal.pbio.3001698`](10.1371/journal.pbio.3001698)). In it you will find all the code necessary to reproduce key figures from the original article and also the same figures produced using the full datasets.


## Abstract

In a recent PLoS Biology article, Jay and colleagues (doi: [`10.1371/journal.pbio.3001698`](10.1371/journal.pbio.3001698)) proposed "a novel theory of sex chromosome evolution" based on the sheltering of deleterious mutations by an inversion linked to a sex-determining locus. Their key result was a demonstration that the fraction of fixed inversions among replicate individual-based Wright-Fisher simulations was much higher for inversions expanding the sex-linked region on a proto-Y chromosome compared with the fraction for autosomal inversions. The results were presented and discussed as estimates of the probability of inversion spread and fixation. However, inspection of Jay et al.’s results and their supplementary materials reveals that they presented the fraction of fixed inversions conditioned on those inversions having survived beyond 20 generations. This unusual condition excludes around 90% of the extinctions occuring in the simulations, and means that their values greatly overestimated the fixation probabilities of new inversions. Moreover, the authors simulated completely neutral inversions (i.e., those capturing only neutral mutations), but these data were excluded from all figures and never mentioned in the paper. Correctly plotting the fraction of fixed inversions using their full dataset reveals that they are indistinguishable from neutrality for most parameter conditions studied, with the notable exception of Y-linked inversions that captured many loci segregating for completely recessive (or nearly so) deleterious mutations. Reevaluation of Jay et al.’s results indicate that the sheltering mechanism they propose is unlikely to robustly provide a selective advantage to inversions suppressing recombination between newly evolving X and Y chromosomes.

## Citing information

Citing information will be made available after this manuscript has passed through peer review.

##  Instructions

This repository provides all code necessary to (1) reproduce the original figures from Jay et al. (2022), (2) produce accompanying figures WITHOUT EXCLUDING neutral inversions or those going extinct in the first 20 generations, and (3) export all figures as .pdf's. To do so, please follow these basic steps:

1. Clone the repo using the following: `git clone https://https://github.com/colin-olito/Jay-2022-Comment`. Alternatively, on the project main page on GitHub, click on the green button `clone` or `download` and then click on `Download ZIP`.  
2. Download the original datasets for Jay et al. (2022) from [`here`](https://figshare.com/articles/dataset/Model_of_sex-chromosome_Evolution_-_datasets/19961033) to the main working directory of this repository, and unzip.  
	- If necessary, rename the unzipped directory containing the datasets `./ModelSexChrom` so that it can be correctly referenced by `R`.  
3. Check that you have a recent version of [`R`](https://www.r-project.org/) installed, and that you have installed and loaded the required packages in `./R/Recreate-and-Compare-Jay-etal-Figs.R`. Correctly running this code requires installation of the following packages:  
	- `library(tidyr)`  
	- `library(patchwork)`  
	- `library(cowplot)`  
	- `library(ggplot2)`  
	- `library(RColorBrewer)`  
	- `library(gghalves)`  
	- `library(dplyr)`  
	- `library(colorRamps)`  
	- `library(cmocean)`  
	- `library(viridis)`  
	- `library(ggnewscale)`  
	- `library(tidyverse)` # I get compile error for dependency package `ragg`. However, the code will run without `tidyverse` installed.  
	- `library(directlabels)`  
4. Make sure that the working directory for your `R` session is the root directory of this repo (e.g., `Jay-2022-Comment-master/`).  
5. Run `./makeFigs.R` either interactively in `R` or in terminal.  


## Repository structure and contents 

The directories/files in this repository needed to reproduce the results for this study are as follows:  

- **`R`**   
	- `Recreate-and-Compare-Jay-etal-Figs.R`  
	- `functions-figures.R`  
- **`figures`***  
- **`ModelSexChrom`*** (Download from [Figshare](https://figshare.com/authors/Paul_Jay/12493000))
- **`MutationShelteringTheory-main`**  
	- Contains original code from Jay et al. (2022)
- `makeFigs.R`  
- `LICENSE.txt`   

**Note:** * The `./figures` directory will be created locally the first time `./makeFigs.R` is run (*if needed*). The folder `./ModelSexChrom` is not tracked in this GitHub repository because the file sizes are too large. The data must be downloaded and unzipped locally by the user to match the above directory tree structure.


### File & variable descriptions

Plotting function files
- `Recreate-and-Compare-Jay-etal-Figs.R`: Code to reproduce the figures from Jay et al. (2022), as well as modified versions of the same plots using the full datasets to calculate inversion fixation rates.  
- `functions-figures.R`: Convenience functions for plotting and code to reproduce supplementary figure S1 from Charlesworth & Olito (2024) comparing analytic approximations with the simulation data from Jay et al. (2022).  

Executables
- `makeFigs.R`: executable functions to create .pdf figures using simulation output files.

Archived Code from Jay et al. (2022)
- `./MutationShelteringTheory-main`: Downloaded original code from Jay et al. (2022)'s GitHub (commit #c2b2f25).

Archived Data from Jay et al. (2022)
- `./ModelSexChrom`: Supplementary data from, must be downloaded from [Figshare](https://figshare.com/authors/Paul_Jay/12493000)). Not tracked in this repo because of large file sizes.

License    
- `LICENSE.txt`: MIT license for this repository.  


## Contact & bug reporting

Please report any bugs, problems, or issues by opening an issue on the inversionSize github [issues page](https://github.com/colin-olito/Jay-2022-Comment/issues). If you would like to report a bug/issue, and do not have a github account (and don't want to get one), please send a brief email detailing the problem you encountered to colin.olito at biol dot lu dot se.

## Licence information

This repository is provided by the authors under the MIT License ([MIT](https://opensource.org/licenses/MIT)).