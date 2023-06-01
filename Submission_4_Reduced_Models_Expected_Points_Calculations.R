#####
## Fitting the reduced models with consistently selected complementary features,
##  calculating the expected points per drive for various scenarios
##  (e.g. with takeaway vs without, with complementary statistics vs without)
## 
## The resulting expected points projections are saved in "projected.points.intercept.RData" file,
## while the full ordinalNet fitted objects are saved in "ordnet.fits.RData"file.
#####

library(ordinalNet)

# Loading the enhanced ordinalNet code.
source("Submission_0_ordinalNet_package_source_code_enhanced.R")

library(tidyverse)
library(car)
library(fuzzyjoin)
library(sure)
library(ordinal)
library(ggplot2)
library(ggfortify)


# We will only focus on the complementary statistics that were consistently picked
# (yards gained and non-scoring takeaways)
all.variables <- FALSE

# We will be combining the yards for offense and special teams returns.
yds.combine <- TRUE




## alpha trade-off value
alpha.val <- 0.99

## We include the POS_TEAM/DEF_POS_TEAM/HOMEFIELD
incl.pos_team <- TRUE


projected.points <- NULL
projected.points.intercept <- NULL
ordnet.fits <- NULL




for (year in 2014:2020){
  
  print("Year:")
  print(year)
  
  
  load(paste0("pbp_by_drive_",year,".Robj"))
  
  ## Non-scoring FORCED TOs, such as picks, fumbles, on-side kick recovery
  pbp_by_drive$takeaways.nonscor <- pbp_by_drive$drive_result_detailed %in% c("Fumble Recovery (Opponent)",
                                                                              "Interception Return",
                                                                              "Kickoff Team Fumble Recovery",
                                                                              "Punt Team Fumble Lost",
                                                                              "Blocked Punt Team Fumble Lost",
                                                                              "On-Side Kick Lost"
  )
  
  ## Any TOs, including FG misses/blocks, turnover on downs
  pbp_by_drive$turnovers.nonscor <- (pbp_by_drive$takeaways.nonscor | 
                                       pbp_by_drive$drive_result_detailed %in% c("Blocked Field Goal",
                                                                                 "Field Goal Missed",
                                                                                 "Missed Field Goal Return",
                                                                                 "Downs Turnover")
  )
  
  pbp_by_drive$defensive.score <- pbp_by_drive$drive_result_detailed %in% c("Blocked Field Goal Touchdown",
                                                                            "Fumble Return Touchdown",
                                                                            "Interception Return Touchdown",
                                                                            "Punt Team Fumble Recovery Touchdown",
                                                                            "Blocked Punt Team Fumble Recovery Touchdown",
                                                                            "Kickoff Team Fumble Recovery Touchdown"
  )
  
  
  ## Homefield indicator (no NEUTRAL field, evidently)
  pbp_by_drive$pos.homefield <- pbp_by_drive$home == pbp_by_drive$pos_team
  
  ## Punts or safeties (both result in turning the ball over via a 'punting motion', although the latter also has -2 pts)
  pbp_by_drive$punt.safety <- pbp_by_drive$drive_result_detailed %in% c("Punt", "Safety")
  
  ## Number of positive runs
  pbp_by_drive$n.positive.runs <- pbp_by_drive$n.rush - pbp_by_drive$n.stuffed.runs
  
  ## Special team return variables:
  pbp_by_drive$yds_ST_return <- ifelse(!is.na(pbp_by_drive$yds_punt_return),
                                       pbp_by_drive$yds_punt_return,
                                       ifelse(!is.na(pbp_by_drive$yds_kickoff_return),
                                              pbp_by_drive$yds_kickoff_return,
                                              0))
  
  pbp_by_drive$yds_ST_return_net <- ifelse(!is.na(pbp_by_drive$yds_punt_net),
                                           pbp_by_drive$yds_punt_net,
                                           ifelse(!is.na(pbp_by_drive$yds_kickoff_net),
                                                  pbp_by_drive$yds_kickoff_net,
                                                  0))
  
  ####
  ## Combining the yards for offense and special teams.
  ####
  if (yds.combine == TRUE){
    pbp_by_drive$off.yards_gained <- pbp_by_drive$off.yards_gained + pbp_by_drive$yds_ST_return
  }
  
  
  #####
  ## Creating the indicator variable for "complementary drives", as in - anything but
  #     1) Kickoffs at the start of halves, 
  #     2) Drives directly after a defensive touchdown.
  #####
  
  pbp_by_drive$complem.ind <- ifelse(lag(pbp_by_drive$game_id_half) == pbp_by_drive$game_id_half,
                                     ifelse(lag(pbp_by_drive$pos_team) != pbp_by_drive$pos_team,
                                            1,
                                            0),
                                     1)
  pbp_by_drive$complem.ind[1] <- 0   # Cleaning up the "NA" generated by the first lag
  
  
  # Initializing a lagged variable set
  pbp_by_drive_lag_vars <- data.frame(game_id_half = pbp_by_drive$game_id_half,
                                      pos_team = pbp_by_drive$pos_team,
                                      pbp_by_drive[, 8:ncol(pbp_by_drive)])
  
  # For numerical variables, creating their lagged values, 
  for (c.ind in 3:ncol(pbp_by_drive_lag_vars)){
    pbp_by_drive_lag_vars[, c.ind] <- lag(pbp_by_drive_lag_vars[, c.ind])
    pbp_by_drive_lag_vars[1, c.ind] <- 0 # For the initial NA value
  }
  
  # Making the variable names for lagged values explicit (adding "lagged" to them)
  colnames(pbp_by_drive_lag_vars)[3:ncol(pbp_by_drive_lag_vars)] <- paste0("lagged_", colnames(pbp_by_drive_lag_vars)[3:ncol(pbp_by_drive_lag_vars)])
  
  # Merging the lagged values into original data set
  pbp_by_drive <- data.frame(pbp_by_drive, pbp_by_drive_lag_vars[,-1])
  
  ## All the complementary statistics under consideration
  relev.var.names <- c("n.completions", "n.incompletions", "n.stuffed.runs", "n.positive.runs", "n.negative.plays.with.penalties",
                       "off.yards_gained",
                       "score_pts_by_text",
                       "yds_ST_return_net", "yds_ST_return",
                       "n.sacks", "n.fumbles",
                       "first.downs.gained", "third.downs.converted",
                       "takeaways.nonscor", "turnovers.nonscor",
                       "punt.safety"
  )
  
  
  ## If we combine yards, then DROP the "yds_ST_return"
  if (yds.combine) relev.var.names <- relev.var.names[relev.var.names != "yds_ST_return"]
  
  
  
  #######
  ## Complementary statistics to be included in the model for testing:
  #######
  
  
  if (all.variables){
    relev.var.names.lagged <- paste0("lagged_", relev.var.names)
  } else {
    relev.var.names <- c("off.yards_gained", "takeaways.nonscor")
    relev.var.names.lagged <- paste0("lagged_", relev.var.names)
  }
  
  pbp_by_drive$categ.score <- ifelse(pbp_by_drive$score_pts_by_text %in% c(6,7,8),
                                     "Offensive TD",
                                     ifelse(pbp_by_drive$score_pts_by_text == 3,
                                            "Field Goal",
                                            ifelse(pbp_by_drive$score_pts_by_text == -2,
                                                   "Safety",
                                                   ifelse(pbp_by_drive$score_pts_by_text %in% c(-6,-7,-8),
                                                          "Defensive TD",
                                                          "No Score"))))
  
  pbp_by_drive$categ.score <- factor(pbp_by_drive$categ.score,
                                     ordered=T,
                                     levels = c("Defensive TD", "Safety", "No Score", "Field Goal", "Offensive TD"))

  
  
  ########
  ## Cleaning up the team names, including the "Non-Major" category for all the non-FBS teams
  ########
  
  team.names.year <- data.frame(Team=sort(unique(c(pbp_by_drive$pos_team, pbp_by_drive$def_pos_team))))
  CFBSTATS_Team_Names <- data.frame(Team = read.csv("~/Documents/Work/New_College/Research/Play_by_Play_Complementary_Football_Project/CFBSTATS_vs_REFERENCE_Team_Names.csv")$CFBSTATS)
  CFBSTATS_Team_Names_Bar_Last <- data.frame(Team = gsub("\\s\\w+$", "", CFBSTATS_Team_Names$Team))
  
  
  
  # Matching up with standardized team names from "cfbstats.com" website.
  matches.df <- stringdist_join(CFBSTATS_Team_Names_Bar_Last, team.names.year, 
                                by='Team', #match based on team
                                mode='left', #use left join
                                method = "jw", #use jw distance metric
                                max_dist=1, 
                                distance_col='dist') %>%
    group_by(Team.x) %>%
    slice_min(order_by=dist, n=1)
  
  
  ## Fixing up all the identified mismatches, such as:
  ##       "Ole Miss" for "Mississippi"
  ##       "NC State" for "North Carolina State" 
  ##       "North Carolina A&T" matches to "North Carolina" along with "North Carolina Tar"...
  ##        For 2020, "Connecticut" didn't play, gotta drop its match ("Cincinnati")
  ##        etc..
  
  FBS.team.names <- matches.df$Team.y
  FBS.team.names[matches.df$Team.x == "Mississippi"] <- "Ole Miss"
  FBS.team.names[matches.df$Team.x == "North Carolina State"] <- "NC State"
  FBS.team.names[matches.df$Team.x == "Miami (Florida)"] <- "Miami"
  FBS.team.names <- sort(
    FBS.team.names[!(matches.df$Team.x == "Coastal Carolina" & matches.df$Team.y == "East Carolina") & 
                     FBS.team.names != "Charleston Southern" & 
                     matches.df$Team.x != "North Carolina A&T" & 
                     FBS.team.names != "North Carolina A&T" &
                     !(matches.df$Team.x == "UAB" & matches.df$Team.y == "UT San Antonio") & 
                     !(matches.df$Team.x == "Connecticut" & matches.df$Team.y == "Cincinnati") & 
                     !(matches.df$Team.x == "Idaho" & matches.df$Team.y == "Indiana") & 
                     !(matches.df$Team.x == "New Mexico State" & matches.df$Team.y == "New Mexico") & 
                     !(matches.df$Team.x == "Old Dominion" & matches.df$Team.y == "Wyoming")])
  
  # Liberty only became FBS from 2018 onwards:
  if (year <= 2017) FBS.team.names <- FBS.team.names[FBS.team.names != "Liberty"]
  # Idaho stopped being FBS from 2018 onwards:
  if (year >= 2018) FBS.team.names <- FBS.team.names[FBS.team.names != "Idaho"]
  
  
  ## Creating the non-major category
  nonmajor.teams <- team.names.year$Team[!team.names.year$Team %in% FBS.team.names]
  pbp_by_drive$pos_team <- ifelse(pbp_by_drive$pos_team %in% FBS.team.names,
                                  pbp_by_drive$pos_team,
                                  "Non-Major")
  pbp_by_drive$def_pos_team <- ifelse(pbp_by_drive$def_pos_team %in% FBS.team.names,
                                      pbp_by_drive$def_pos_team,
                                      "Non-Major")
  
  
  
  ###########
  ###########
  ### Calculating PROJECTED POINTS for EACH TEAM
  ###########
  ###########
  
  pbp_by_drive$game_id <- as.factor(str_remove(pbp_by_drive$game_id_half, "-1|-2"))
  
  pbp_by_drive.noNA <- na.omit(pbp_by_drive[, c("categ.score", "pos_team", "def_pos_team", relev.var.names.lagged, "game_id", "game_id_half", "pos.homefield", "complem.ind")] %>%
                                 mutate(pos_team = factor(pos_team),
                                        def_pos_team = factor(def_pos_team),
                                        game_id = factor(game_id),
                                        game_id_half = factor(game_id_half)))
  
  
  ## Setting up the contrasts so that the intercept represented a "league-average opponent"
  ## (sum_i alpha_i = sum_j beta_j = 0)
  contrasts(pbp_by_drive.noNA$pos_team) <- contr.sum(nlevels(pbp_by_drive.noNA$pos_team))
  contrasts(pbp_by_drive.noNA$def_pos_team) <- contr.sum(nlevels(pbp_by_drive.noNA$def_pos_team))
  colnames(contrasts(pbp_by_drive.noNA$pos_team)) <- rownames(contrasts(pbp_by_drive.noNA$pos_team))[-nrow(contrasts(pbp_by_drive.noNA$pos_team))]
  colnames(contrasts(pbp_by_drive.noNA$def_pos_team)) <- rownames(contrasts(pbp_by_drive.noNA$def_pos_team))[-nrow(contrasts(pbp_by_drive.noNA$def_pos_team))]
  
  
  
  
  ## First: basic lm fit, just to get the design matrix X & response y
  lm.obj <- lm(formula(paste0("as.numeric(categ.score) ~", paste0(c("pos_team", "def_pos_team", "pos.homefield", "complem.ind", paste0(relev.var.names.lagged, ":complem.ind", sep="")),
                                                                  collapse = " + "), collapse =" ")),
               data=pbp_by_drive.noNA[, c(if(incl.pos_team) c("pos_team", "def_pos_team", "pos.homefield", "complem.ind") else NULL,
                                          relev.var.names.lagged,
                                          "categ.score")])
  X <- model.matrix(lm.obj)[,-1]
  y <- pbp_by_drive.noNA$categ.score

  # Indices of the design matrix X to NOT penalize:
  ind.nonpen <- which(!str_detect(colnames(X), fixed("lagged")))
  penfact.vec <- rep(1, ncol(X))
  penfact.vec[ind.nonpen] <- 0
  

  #####
  ## Fitting the reduced selected model
  #####
  
  ordnet.obj <- ordinalNet.mine(x = as.matrix(X), y = y, alpha=alpha.val, standardize = T,
                                penaltyFactors = penfact.vec,
                                family="cumulative",
                                parallelTerms = T,
                                nonparallelTerms = T,
                                nonparallel.categ = c(FALSE, FALSE, TRUE, FALSE),
                                nonparallel.factor = str_detect(colnames(X), "lagged"),
                                lambdaVals = 0,
                                printIter=TRUE
  )
  
  

  
  ######
  ## NO COMPLEMENTARY (only strength-of-schedule and homefield adjustment)
  ######
  
  
  ## First: basic lm fit, just to get the design matrix X & response y
  lm.no.complem.obj <- lm(as.numeric(categ.score) ~ pos_team + def_pos_team + pos.homefield,
                          data=pbp_by_drive.noNA)
  
  X.no.complem <- model.matrix(lm.no.complem.obj)[,-1]
  y.no.complem <- pbp_by_drive.noNA$categ.score
  
  ind.nonpen <- which(!str_detect(colnames(X.no.complem), fixed("lagged")))
  penfact.vec <- rep(1, ncol(X.no.complem))
  penfact.vec[ind.nonpen] <- 0
  
  
  ordnet.no.complem.obj <- ordinalNet.mine(x = X.no.complem, y = y.no.complem, alpha=alpha.val, standardize = T, 
                                           penaltyFactors = penfact.vec,
                                           family="cumulative",
                                           parallelTerms = T,
                                           nonparallelTerms = F,
                                           lambdaVals = 0,
                                           printIter=TRUE
  )
  
  ordnet.no.complem.obj$coefs
  
  
  
  ###
  # SAVING THE "ordinalNet" FITS
  ###
  ordnet.fits[[as.character(year)]] <- list(with.complem = ordnet.obj,
                                            no.complem = ordnet.no.complem.obj,
                                            with.complem.paral.only = ordnet.paral.obj)
  

  ####
  # Making expected point predictions for each team.
  ####
  
  # Saving the matrix of contrasts for offensive and defensive margin parameters.
  contr.mat.pos_team <- contrasts(pbp_by_drive.noNA$pos_team)
  contr.mat.def_pos_team <- contrasts(pbp_by_drive.noNA$def_pos_team)
  
  # Calculating the mean vector for the values of all complementary statistics,
  # to project the performances onto the league-average complementary unit.
  mean.vec <- as.matrix(t(apply(pbp_by_drive.noNA[, relev.var.names.lagged],
                                2,
                                mean)))
  
  # Making the predictions for each team for a variety of scenarios
  # (e.g. with takeaway vs without, with complementary statistics vs without)
  
  for (i in 1:nrow(contr.mat.pos_team)){
    team.name <- rownames(contr.mat.pos_team)[i]
    
    pos_team.mat <- matrix(contr.mat.pos_team[i,], nrow = nrow(contr.mat.pos_team), ncol = ncol(contr.mat.pos_team), byrow=T)
    def_pos_team.mat <- matrix(contr.mat.def_pos_team[i,], nrow = nrow(contr.mat.def_pos_team), ncol = ncol(contr.mat.def_pos_team), byrow=T)
    
    
    off.vec.points.no_takeaway <- sum(predict.ordinalNet.mine(ordnet.obj, 
                                                              newx = cbind(if (incl.pos_team){
                                                                cbind(as.matrix(t(pos_team.mat[1,])), 
                                                                      as.matrix(t(rep(0, ncol(contr.mat.def_pos_team)))),
                                                                      pos.homefield = mean(pbp_by_drive.noNA$pos.homefield),
                                                                      complem.ind = mean(pbp_by_drive.noNA$complem.ind))} else NULL,
                                                                as.matrix(t(ifelse(relev.var.names.lagged == "lagged_takeaways.nonscor",
                                                                                   0,
                                                                                   mean.vec)))))*c(-7,-2,0,3,7))
    
    
    
    def.vec.points.no_takeaway <- sum(predict.ordinalNet.mine(ordnet.obj, 
                                                              newx = cbind(if (incl.pos_team){
                                                                cbind(as.matrix(t(rep(0, ncol(contr.mat.pos_team)))), 
                                                                      as.matrix(t(def_pos_team.mat[1,])),
                                                                      pos.homefield = mean(pbp_by_drive.noNA$pos.homefield),
                                                                      complem.ind = mean(pbp_by_drive.noNA$complem.ind))} else NULL,
                                                                as.matrix(t(ifelse(relev.var.names.lagged == "lagged_takeaways.nonscor",
                                                                                   0,
                                                                                   mean.vec)))))*c(-7,-2,0,3,7))
    
    
    
    
    off.vec.points <- sum(predict.ordinalNet.mine(ordnet.obj, 
                                                  newx = cbind(if (incl.pos_team){
                                                    cbind(as.matrix(t(pos_team.mat[1,])), 
                                                          as.matrix(t(rep(0, ncol(contr.mat.def_pos_team)))),
                                                          pos.homefield = mean(pbp_by_drive.noNA$pos.homefield),
                                                          complem.ind = mean(pbp_by_drive.noNA$complem.ind))} else NULL,
                                                    as.matrix(t(ifelse(relev.var.names.lagged == "lagged_takeaways.nonscor",
                                                                       mean.vec,
                                                                       mean.vec)))))*c(-7,-2,0,3,7))
    
    
    
    def.vec.points <- sum(predict.ordinalNet.mine(ordnet.obj, 
                                                  newx = cbind(if (incl.pos_team){
                                                    cbind(as.matrix(t(rep(0, ncol(contr.mat.pos_team)))), 
                                                          as.matrix(t(def_pos_team.mat[1,])),
                                                          pos.homefield = mean(pbp_by_drive.noNA$pos.homefield),
                                                          complem.ind = mean(pbp_by_drive.noNA$complem.ind))} else NULL,
                                                    as.matrix(t(ifelse(relev.var.names.lagged == "lagged_takeaways.nonscor",
                                                                       mean.vec,
                                                                       mean.vec)))))*c(-7,-2,0,3,7))
    
    off.vec.points.no.complem <- sum(predict.ordinalNet.mine(ordnet.no.complem.obj, 
                                                             newx = cbind(as.matrix(t(pos_team.mat[1,])), as.matrix(t(rep(0, ncol(contr.mat.def_pos_team)))),
                                                                          pos.homefield = mean(pbp_by_drive.noNA$pos.homefield)))*c(-7,-2,0,3,7))
    
    
    def.vec.points.no.complem <- sum(predict.ordinalNet.mine(ordnet.no.complem.obj, 
                                                             newx = cbind(as.matrix(t(rep(0, ncol(contr.mat.pos_team)))), as.matrix(t(def_pos_team.mat[1,])),
                                                                          pos.homefield = mean(pbp_by_drive.noNA$pos.homefield)))*c(-7,-2,0,3,7))

    off.games.included <- length(unique(pbp_by_drive.noNA[pbp_by_drive.noNA$pos_team == team.name, ]$game_id))
    off.halves.included <- length(unique(pbp_by_drive.noNA[pbp_by_drive.noNA$pos_team == team.name, ]$game_id_half))
    
    def.games.included <- length(unique(pbp_by_drive.noNA[pbp_by_drive.noNA$def_pos_team == team.name, ]$game_id))
    def.halves.included <- length(unique(pbp_by_drive.noNA[pbp_by_drive.noNA$def_pos_team == team.name, ]$game_id_half))
    
    off.vec.points.with_takeaway <- sum(predict.ordinalNet.mine(ordnet.obj, 
                                                                newx = cbind(if (incl.pos_team){
                                                                  cbind(as.matrix(t(pos_team.mat[1,])), 
                                                                        as.matrix(t(rep(0, ncol(contr.mat.def_pos_team)))),
                                                                        pos.homefield = mean(pbp_by_drive.noNA$pos.homefield),
                                                                        complem.ind = mean(pbp_by_drive.noNA$complem.ind))} else NULL,
                                                                  as.matrix(t(ifelse(relev.var.names.lagged == "lagged_takeaways.nonscor",
                                                                                     1,
                                                                                     mean.vec)))))*c(-7,-2,0,3,7))
    
    def.vec.points.with_takeaway <- sum(predict.ordinalNet.mine(ordnet.obj, 
                                                                newx = cbind(if (incl.pos_team){
                                                                  cbind(as.matrix(t(rep(0, ncol(contr.mat.pos_team)))), 
                                                                        as.matrix(t(def_pos_team.mat[1,])),
                                                                        pos.homefield = mean(pbp_by_drive.noNA$pos.homefield),
                                                                        complem.ind = mean(pbp_by_drive.noNA$complem.ind))} else NULL,
                                                                  as.matrix(t(ifelse(relev.var.names.lagged == "lagged_takeaways.nonscor",
                                                                                     1,
                                                                                     mean.vec)))))*c(-7,-2,0,3,7))
    
    ####
    ## SAVING THE PROJECTED POINTS
    ####
    
    projected.points.intercept[[as.character(year)]] <- rbind(projected.points.intercept[[as.character(year)]],
                                                              data.frame(Team=team.name, 
                                                                         Offense.Avg.Points.Avg.Takeaways = mean(off.vec.points),
                                                                         Offense.Avg.Points.No.Complem = mean(off.vec.points.no.complem),
                                                                         Offense.Games.Included = off.games.included,
                                                                         Offense.Halves.Included = off.halves.included,
                                                                         Offense.Avg.Points.No.Takeaway = mean(off.vec.points.no_takeaway),
                                                                         Offense.Avg.Points.With.Takeaway = mean(off.vec.points.with_takeaway),
                                                                         
                                                                         Defense.Avg.Points.Avg.Takeaways = mean(def.vec.points),
                                                                         Defense.Avg.Points.No.Complem = mean(def.vec.points.no.complem),
                                                                         Defense.Games.Included = def.games.included,
                                                                         Defense.Halves.Included = def.halves.included,
                                                                         Defense.Avg.Points.No.Takeaway = mean(def.vec.points.no_takeaway),
                                                                         Defense.Avg.Points.With.Takeaway = mean(def.vec.points.with_takeaway)
                                                              ))
    
    
  }
  
}



save(projected.points.intercept, file= paste0("ORDNET_projected.points.intercept.RData"))
save(ordnet.fits, file=paste0("ORDNET_ordnet.fits.RData"))




