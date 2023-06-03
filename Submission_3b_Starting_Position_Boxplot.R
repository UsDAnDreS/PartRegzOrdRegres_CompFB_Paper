#####
## Reproducing the box plot of starting positions depending on
## whether the defense got a takeaway on the previous drive.
#####


library(tidyverse)
library(car)
library(fuzzyjoin)
library(ggplot2)
library(ggfortify)

# We will only focus on the complementary statistics that were consistently picked
# (yards gained and non-scoring takeaways)
all.variables <- FALSE

## Loading all the seasons

pbp_by_drive_full <- NULL

for (year in 2014:2020){
  print(year)
  load(paste0("pbp_by_drive_",year,".Robj"))
  pbp_by_drive_full <- rbind(pbp_by_drive_full,
                             cbind(Year = year, pbp_by_drive))
}



## Non-scoring FORCED TOs, such as picks, fumbles, on-side kick recovery
pbp_by_drive_full$takeaways.nonscor <- pbp_by_drive_full$drive_result_detailed %in% c("Fumble Recovery (Opponent)",
                                                                                      "Interception Return",
                                                                                      "Kickoff Team Fumble Recovery",
                                                                                      "Punt Team Fumble Lost",
                                                                                      "Blocked Punt Team Fumble Lost",
                                                                                      "On-Side Kick Lost"
)

## Any TOs, including FG misses/blocks, turnover on downs
pbp_by_drive_full$turnovers.nonscor <- (pbp_by_drive_full$takeaways.nonscor | 
                                          pbp_by_drive_full$drive_result_detailed %in% c("Blocked Field Goal",
                                                                                         "Field Goal Missed",
                                                                                         "Missed Field Goal Return",
                                                                                         "Downs Turnover")
)

pbp_by_drive_full$defensive.score <- pbp_by_drive_full$drive_result_detailed %in% c("Blocked Field Goal Touchdown",
                                                                                    "Fumble Return Touchdown",
                                                                                    "Interception Return Touchdown",
                                                                                    "Punt Team Fumble Recovery Touchdown",
                                                                                    "Blocked Punt Team Fumble Recovery Touchdown",
                                                                                    "Kickoff Team Fumble Recovery Touchdown"
)

## Homefield indicator (no NEUTRAL field, evidently)
pbp_by_drive_full$pos.homefield <- pbp_by_drive_full$home == pbp_by_drive_full$pos_team

## Punts or safeties (both result in turning the ball over via a 'punting motion', although the latter also has -2 pts)
pbp_by_drive_full$punt.safety <- pbp_by_drive_full$drive_result_detailed %in% c("Punt", "Safety")

## Number of positive runs
pbp_by_drive_full$n.positive.runs <- pbp_by_drive_full$n.rush - pbp_by_drive_full$n.stuffed.runs

## Special team return variables:
pbp_by_drive_full$yds_ST_return <- ifelse(!is.na(pbp_by_drive_full$yds_punt_return),
                                          pbp_by_drive_full$yds_punt_return,
                                          ifelse(!is.na(pbp_by_drive_full$yds_kickoff_return),
                                                 pbp_by_drive_full$yds_kickoff_return,
                                                 0))

pbp_by_drive_full$yds_ST_return_net <- ifelse(!is.na(pbp_by_drive_full$yds_punt_net),
                                              pbp_by_drive_full$yds_punt_net,
                                              ifelse(!is.na(pbp_by_drive_full$yds_kickoff_net),
                                                     pbp_by_drive_full$yds_kickoff_net,
                                                     0))


## Some "off yards" are > 100 due to FG being blocked and returned back to the team
## (THERE'S ONLY 2 SUCH OBSERVATIONS)
pbp_by_drive_full$off.yards_gained[pbp_by_drive_full$off.yards_gained > 100]

## Way more of the ones where OFF YARDS + YDS_ST_RETURN are >110
sum(pbp_by_drive_full$off.yards_gained + pbp_by_drive_full$yds_ST_return >= 110)



# Initializing a lagged variable set
pbp_by_drive_full_lag_vars <- data.frame(game_id_half = pbp_by_drive_full$game_id_half,
                                         pos_team = pbp_by_drive_full$pos_team,
                                         pbp_by_drive_full[, 9:ncol(pbp_by_drive_full)])

for (c.ind in 3:ncol(pbp_by_drive_full_lag_vars)){
  pbp_by_drive_full_lag_vars[, c.ind] <- lag(pbp_by_drive_full_lag_vars[, c.ind])
}

# Making the variable names for lagged values explicit (adding "lagged" to them)
colnames(pbp_by_drive_full_lag_vars)[3:ncol(pbp_by_drive_full_lag_vars)] <- paste0("lagged_", colnames(pbp_by_drive_full_lag_vars)[3:ncol(pbp_by_drive_full_lag_vars)])

# Merging the lagged values into original data set
pbp_by_drive_full <- data.frame(pbp_by_drive_full, pbp_by_drive_full_lag_vars[,-1])


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



#######
## Complementary statistics to be included in the model for testing:
#######

if (all.variables){
  relev.var.names.lagged <- paste0("lagged_", relev.var.names)
} else {
  relev.var.names <- c("off.yards_gained", "takeaways.nonscor")
  relev.var.names.lagged <- paste0("lagged_", relev.var.names)
}


# All the distinct drive outcomes
all.unique.res <- sort(unique(pbp_by_drive_full$drive_result_detailed))

# Excluding the outcomes that imply a KICK/PUNT to be performed first (which messes with the significance of "stopped position")
all.unique.res.excl.punt.kick.preceeding <- all.unique.res[!str_detect(all.unique.res, "Touchdown") 
                                                           & all.unique.res != "Punt" 
                                                           & !str_detect(all.unique.res, "Safety")
                                                           & !str_detect(all.unique.res, "Field Goal Good")]


########
## Calculating a SPECIAL TEAMS RETURN starting position (so, after having accounted for KICK/PUNT length)
########

pbp_by_drive_full$stop_pos_with.kick.punt.yardage <- pbp_by_drive_full$stopped_position
pbp_by_drive_full$stop_pos_with.kick.punt.yardage[!is.na(lead(pbp_by_drive_full$yds_punted)) | !is.na(lead(pbp_by_drive_full$yds_kickoff))]

pbp_by_drive_full$stop_pos_with.kick.punt.yardage <- ifelse(!is.na(lead(pbp_by_drive_full$yds_punted)) & (lead(pbp_by_drive_full$game_id_half) == pbp_by_drive_full$game_id_half),
                                                            ifelse(str_detect(tolower(lead(pbp_by_drive_full$first.play.text.no.pen)), "touchback"), 
                                                                   80,
                                                                   pbp_by_drive_full$stop_pos_with.kick.punt.yardage + lead(pbp_by_drive_full$yds_punted)),
                                                            ifelse(!is.na(lead(pbp_by_drive_full$yds_kickoff)) & (lead(pbp_by_drive_full$game_id_half) == pbp_by_drive_full$game_id_half),
                                                                   ifelse(str_detect(tolower(lead(pbp_by_drive_full$first.play.text.no.pen)), "touchback"),
                                                                          ifelse(str_detect(tolower(pbp_by_drive_full$drive_result_detailed), "safety"),
                                                                                 80,
                                                                                 75),
                                                                          pbp_by_drive_full$stop_pos_with.kick.punt.yardage + lead(pbp_by_drive_full$yds_kickoff)),
                                                                   pbp_by_drive_full$stop_pos_with.kick.punt.yardage))
pbp_by_drive_full$stop_pos_with.kick.punt.yardage[!is.na(lead(pbp_by_drive_full$yds_kickoff))]


####
# For the TAKEAWAYS IMPACT
####

## width:  height:
ggplot(data = pbp_by_drive_full %>%
         filter(stop_pos_with.kick.punt.yardage <= 110) %>%
         mutate(takeaways.nonscor = ifelse(takeaways.nonscor, "Yes", "No")),
       aes(x=takeaways.nonscor, y=stop_pos_with.kick.punt.yardage, fill = takeaways.nonscor)) +
  geom_boxplot() + 
  xlab("Non-Scoring Takeaway\non Previous Drive") +
  ylab("Yards to Go\n on Current Drive") + 
  theme(legend.position = "none") 