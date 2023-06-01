####
## (NOTE: Takes a LOOOOOOONG time. Best to do year-by-year in parallel.)
##
## For each season (2014-2020), carrying out 3 replicates of 10-fold cross-validation.
## Recording the resulting cross-validation objects into "TUNE_OBJECTS.RData" files,
## and the picked "lambda.1se" estimates into "lambda.1se.est.RData" files.
####


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


# We will be combining the yards for offense and special teams returns.
yds.combine <- TRUE


# Creating the folder to contain the resulting fits/estimates.
dir.create("MULT_REP_NONPARALLEL_ORDNET")


n.lambdas <- 40
n.folds <- 10
n.CV.reps <- 3

## Whether we include POS_TEAM/DEF_POS_TEAM/HOMEFIELD
incl.pos_team <- TRUE

## For the selections
lambda.min.est <- lambda.1se.est <- BIC.est <- list()

## For the # of folds that yielded negative probabilities
bad.folds <- list()

# List of calculated ordnet objects
ordnet.obj.list <- list()


for (year in 2014:2020){
  
  # year <- 2014
  
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
  
  relev.var.names.lagged <- paste0("lagged_", relev.var.names)
  
  ## Defining the categorical ordinal scoring outcome variable
  
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
  ### Preparing the data frame and design matrices for cross-validation
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
  
  
  ########
  ## For CROSS-VALIDATION: Setting the seed, and RE-RUNNING THE K-FOLD CV SEVERAL TIMES
  ########
  
  
  set.seed(1)
  
  for (j in 1:n.CV.reps){
    print("CV REPLICATE #:")
    print(j)
    
    ######
    ## THE MAIN FITTING FUNCTION
    ######
    
    ordnet.obj <- ordinalNetTune.mine(x = as.matrix(X), y = y, alpha=0.99, standardize = T, 
                                      penaltyFactors = penfact.vec, 
                                      family="cumulative",
                                      parallelTerms = T,
                                      nonparallelTerms = T,
                                      nonparallel.factor = str_detect(colnames(X), "lagged"),
                                      nFolds = n.folds,
                                      nLambda= n.lambdas,
                                      printProgress=TRUE,
                                      warn=FALSE
    )
    
    
    ### Saving the full CV object.
    ordnet.obj.list[[as.character(year)]] <- append(ordnet.obj.list[[as.character(year)]], list(ordnet.obj))
    
    ### Saving the lambda.1se estimate
    mean.CVs <- apply(ordnet.obj$loglik, 1, function(x) mean(x[x != -Inf]))
    se.CVs <- apply(ordnet.obj$loglik, 1, function(x) sd(x[x != -Inf]))/sqrt(n.folds)
    
    ## Calculating the lambda.1se estimate
    lambda.min.ind <- which.max(mean.CVs); 
    lambda.min <- ordnet.obj$lambdaVals[lambda.min.ind]
    lambda.1se.ind <- which(ordnet.obj$lambdaVals >= lambda.min & ((mean.CVs + se.CVs) >= mean.CVs[lambda.min.ind]));
    lambda.1se <- max(ordnet.obj$lambdaVals[lambda.1se.ind])
    lambda.1se.ind <- which(ordnet.obj$lambdaVals == lambda.1se)
    
    lambda.1se.est[[as.character(year)]] <- rbind(lambda.1se.est[[as.character(year)]],
                                                  ordnet.obj$fit$coefs[lambda.1se.ind, which(str_detect(colnames(ordnet.obj$fit$coefs), fixed("lagged")))])
    
    
  }
  

  save(ordnet.obj.list, file= paste0("MULT_REP_NONPARALLEL_ORDNET/year=", as.character(year),  "_nlambdas=", n.lambdas, "_nfolds=", n.folds, "_nreps=", n.CV.reps, "_TUNE_OBJECTS.RData"))
  save(lambda.1se.est, file= paste0("MULT_REP_NONPARALLEL_ORDNET/year=", as.character(year), "_nlambdas=", n.lambdas, "_nfolds=", n.folds, "_nreps=", n.CV.reps, "_lambda.1se.est.RData"))

}

