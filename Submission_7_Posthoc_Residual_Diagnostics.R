#########
## Producing the "binarized" surrogate residual diagnostic plots for our
## non-proportional odds fits. We use the enhanced version of the source code for
## "sure" package to allow creating a separate plot (residuals-vs-fitted or QQ)
## for each individual cumulative odds fit of our non-proportional odds model.
#########


library(ordinalNet)
# Loading the enhanced ordinalNet code.
source("Submission_0_ordinalNet_package_source_code_enhanced.R")

library(sure)
# Loading the enhanced sure code.
source("Submission_0p2_sure_package_source_code_enhanced.R")

library(tidyverse)
library(car)
library(fuzzyjoin)
library(ordinal)
library(ggplot2)
library(ggfortify)


# We will be combining the yards for offense and special teams returns.
yds.combine <- TRUE

## Loading all the relevant objects
load("ORDNET_projected.points.intercept.RData")
load("ORDNET_ordnet.fits.RData")


## Looping through the fits for all seasons,
## plotting the "binarized" residuals-vs-fitted and QQ-plot diagnostics
## for each cumulative odds fit separately 
## (via using the modified version of "sure"'s residual generation)

set.seed(101)

for (year in 2014:2020){
  
  print("Year:")
  print(year)
  
  ordnet.obj <- ordnet.fits[[as.character(year)]]$with.complem
  ordnet.paral.obj <- ordnet.fits[[as.character(year)]]$with.complem.paral.only
  
  ###
  # "Binarized" residuals-vs-fitted plots for all 4 cumulative odds fits
  ###
  
  categ.names <- c("\u2265 Safety", "\u2265 No Score", "\u2265 FG", "\u2265 Off. TD")
  plot.list <- list()
  for (categ.num in 1:(ordnet.obj$nLev-1)){
    plot.list[[categ.num]] <- autoplot.ordinalNet(ordnet.obj, what = c("fitted"), which.categ = categ.num, binarized = TRUE) +
      ggtitle(paste0(year, ", Odds(", categ.names[categ.num], ")")) #
  }
  
  print(grid.arrange(plot.list[[1]], plot.list[[2]],
                     plot.list[[3]], plot.list[[4]],
                     ncol=2, nrow=2))
  
  ###
  # "Binarized" QQ plots for all 4 cumulative odds fits
  ###
  
  categ.names <- c("\u2265 Safety", "\u2265 No Score", "\u2265 FG", "\u2265 Off. TD")
  plot.list <- list()
  for (categ.num in 1:(ordnet.obj$nLev-1)){
    set.seed(101)
    plot.list[[categ.num]] <- autoplot.ordinalNet(ordnet.obj, what = c("qq"), which.categ = categ.num, binarized = TRUE) +
      ggtitle(paste0(year, ", Odds(", categ.names[categ.num], ")")) #
    
  }
  
  print(grid.arrange(plot.list[[1]], plot.list[[2]],
                     plot.list[[3]], plot.list[[4]],
                     ncol=2, nrow=2))
  
}

