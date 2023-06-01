#####
## Reproducing the figure from the paper that showed the effect of a non-scoring takeaway
## on the expected points scored per drive (plotted against the baseline of points scored without a takeaway).
#####


library(tidyverse)
library(gridExtra)
library(fuzzyjoin)
library(plotrix)


# We will only focus on the complementary statistics that were consistently picked
# (yards gained and non-scoring takeaways)
all.variables <- FALSE

# We will be combining the yards for offense and special teams returns.
yds.combine <- TRUE

## Loading all the relevant objects
load("ORDNET_projected.points.intercept.RData")
load("ORDNET_ordnet.fits.RData")


######
##  Putting together a full data set of the expected points per drive projected onto a league-average complementary unit
##  for each team across all seasons
######

my.df <- NULL

for (year in 2014:2020){
  
  my.df <- rbind(my.df,
                 data.frame(Year = year, projected.points.intercept[[as.character(year)]]))
  
}

my.df$Year <- as.factor(my.df$Year)


####
# THE MAIN PLOT
####

## Width: 553; Height: 444
print(ggplot(data = my.df %>% mutate(Extra = Offense.Avg.Points.With.Takeaway - Offense.Avg.Points.No.Takeaway),
aes(x=Offense.Avg.Points.No.Takeaway, y=Extra, col=Year)) +
  geom_smooth(linetype=1) + 
  xlab("Points Scored Without A Takeaway (Baseline)") +
  ylab("Extra Points Scored After Takeaway") + 
  ggtitle("Offense"))