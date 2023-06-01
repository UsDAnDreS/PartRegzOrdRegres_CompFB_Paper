#########
## Generating the table from the paper of ranking/value shifts for 
## offensive and defensive points per drive modifications based on the
## complementary football adjustments. 
#########


library(ordinalNet)
# Loading the enhanced ordinalNet code.
source("Submission_0_ordinalNet_package_source_code_enhanced.R")

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


## Picking the season
year <- 2016

## Picking the side ("Offense" or "Defense")
side <- "Offense"

#####
## Obtaining the values with and without adjustment for complementary football features
#####

adj_opp.rank.table <- projected.points.intercept[[as.character(year)]][, c("Team", paste0(side, ".Avg.Points.No.Complem"))]
adj_opp_and_ballside.rank.table <- projected.points.intercept[[as.character(year)]][, c("Team", paste0(side, c(".Avg.Points.Avg.Takeaways")))]
colnames(adj_opp.rank.table)[1:2] <- colnames(adj_opp_and_ballside.rank.table)[1:2] <- c("Team", "Value")

adj_opp.rank.table$Rank <- if (side == "Defense") rank(adj_opp.rank.table$Value) else rank(-adj_opp.rank.table$Value)
adj_opp_and_ballside.rank.table$Rank <- if (side == "Defense") rank(adj_opp_and_ballside.rank.table$Value) else rank(-adj_opp_and_ballside.rank.table$Value)

opp_vs_ballside_shift <- data.frame(adj_opp.rank.table,
                                    adj_opp_and_ballside.rank.table[,-1],
                                    Rank.Diff = adj_opp.rank.table$Rank - adj_opp_and_ballside.rank.table$Rank,
                                    Value.Diff = ifelse(rep(side, nrow(adj_opp_and_ballside.rank.table)) == "Offense",
                                                        adj_opp_and_ballside.rank.table$Value - adj_opp.rank.table$Value,
                                                        adj_opp.rank.table$Value - adj_opp_and_ballside.rank.table$Value),
                                    Halves = projected.points.intercept[[as.character(year)]][, paste0(side, ".Halves.Included")],
                                    Games = projected.points.intercept[[as.character(year)]][, paste0(side, ".Games.Included")])

opp_vs_ballside_shift <-
  opp_vs_ballside_shift %>%
  mutate(Value=round(Value,2), Value.1=round(Value.1,2), Value.Diff = round(Value.Diff, 2))



#####
# Quick ranking demonstration
#####

### The Top-12 teams
opp_vs_ballside_shift %>% arrange(Rank.1)  %>% head(12)


## Sorting by the strongest (positive or negative) value shifts 
opp_vs_ballside_shift %>% 
  arrange(desc(abs(Value.Diff)))  %>% head(12)



## Prepping for the final table

if (year == 2016){
  top.off <- 10
  top.def <- 10
  if (side == "Offense") which.ind <- c(1:top.off, 15, 59)
  if (side == "Defense") which.ind <- c(1:top.def, 57, 76)
}

## Table sorted by rank
our.table <- opp_vs_ballside_shift %>% arrange(Rank.1) 



######
## Loading all the team statistics from that season season
## to calculate the per-drive complementary statistics values
## (which will be used to demonstrate the underlying mechanism of ranking shifts)
######

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
## !! COMBINING THE YARDS for OFFENSE + SPECIAL TEAMS !!
####
if (yds.combine == TRUE){
  pbp_by_drive$off.yards_gained <- pbp_by_drive$off.yards_gained + pbp_by_drive$yds_ST_return
}


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



nonmajor.teams <- team.names.year$Team[!team.names.year$Team %in% FBS.team.names]
nonmajor.teams


## Creating the non-major category
nonmajor.teams <- team.names.year$Team[!team.names.year$Team %in% FBS.team.names]
pbp_by_drive$pos_team <- ifelse(pbp_by_drive$pos_team %in% FBS.team.names,
                                pbp_by_drive$pos_team,
                                "Non-Major")
pbp_by_drive$def_pos_team <- ifelse(pbp_by_drive$def_pos_team %in% FBS.team.names,
                                    pbp_by_drive$def_pos_team,
                                    "Non-Major")




## Calculating the per-drive average rankings on complementary football stats

if (side == "Offense"){
  our.avg.mat <- pbp_by_drive %>%
    group_by(def_pos_team) %>%
    summarize(n.drives.played = length(unique(game_id_drive)),
              off.yards_gained.per.drive = sum(off.yards_gained)/n.drives.played,
              takeaways.nonscor.per.drive = sum(takeaways.nonscor)/n.drives.played) %>%
    rename(Team = def_pos_team)
} else {
  our.avg.mat <- pbp_by_drive %>%
    group_by(pos_team) %>%
    summarize(n.drives.played = length(unique(game_id_drive)),
              off.yards_gained.per.drive = sum(off.yards_gained)/n.drives.played,
              takeaways.nonscor.per.drive = sum(takeaways.nonscor)/n.drives.played) %>%
    rename(Team = pos_team)
}

our.pred <- colnames(our.avg.mat)[str_detect(colnames(our.avg.mat), "per.drive")]


## Adding the values and rankings of teams based on the complementary football statistics

if (side == "Offense"){
  
  our.rank.mat <- data.frame(Team=our.avg.mat$Team,
                             Yds.Value=format(round(unlist(our.avg.mat[, our.pred[1]]))),
                             Yds.Rank=rank(unlist(our.avg.mat[, our.pred[1]]), tie="first"),
                             TO.Value=format(round(unlist(our.avg.mat[, our.pred[2]]),2), nsmall=2),
                             TO.Rank=rank(-unlist(our.avg.mat[, our.pred[2]]), tie="first"))
  rownames(our.rank.mat) <- NULL

  final.table <- our.table %>% mutate(PPG.Value.Shift = ifelse(Value.Diff > 0, 
                                                               paste0(format(round(Value.1,2), nsmall=2), " ($\\textcolor{darkgreen}{\\Uparrow}$ ", format(round(abs(Value.Diff),2), nsmall=2), ")"),
                                                               ifelse(Value.Diff <0,
                                                                      paste0(format(round(Value.1,2), nsmall=2), " ($\\textcolor{red}{\\Downarrow}$ ", format(round(abs(Value.Diff),2), nsmall=2), ")"),
                                                                      paste0(format(round(Value.1,2), nsmall=2), " ($\\textcolor{blue}{=}", format(round(abs(Value.Diff),2), nsmall=2), "$)"))),
                                      PPG.Rank.Shift = ifelse(Rank.Diff > 0, 
                                                              paste0(Rank.1, " ($\\textcolor{darkgreen}{\\Uparrow}$ ", abs(Rank.Diff), ")"),
                                                              ifelse(Rank.Diff <0,
                                                                     paste0(Rank.1, " ($\\textcolor{red}{\\Downarrow}$ ", abs(Rank.Diff), ")"),
                                                                     paste0(Rank.1, " ($\\textcolor{blue}{=}", abs(Rank.Diff), "$)")))) %>%
    dplyr::select(-Value, -Rank,  -Rank.Diff, -Value.Diff) %>%
    left_join(our.rank.mat)
  
  final.table <- final.table %>% mutate(TO.Rank.Value = paste0(TO.Rank, " (", TO.Value, ")"),
                                        Yds.Rank.Value = paste0(Yds.Rank, " (", Yds.Value, ")"))
  
  
} else {
  
  our.rank.mat <- data.frame(Team=our.avg.mat$Team,
                             Yds.Value=format(round(unlist(our.avg.mat[, our.pred[1]]))),
                             Yds.Rank=rank(-unlist(our.avg.mat[, our.pred[1]]), tie="first"),
                             TO.Value=format(round(unlist(our.avg.mat[, our.pred[2]]),2), nsmall=2),
                             TO.Rank=rank(unlist(our.avg.mat[, our.pred[2]]), tie="first"))
  rownames(our.rank.mat) <- NULL
  
  final.table <- our.table %>% mutate(PPG.Value.Shift = ifelse(Value.Diff > 0, 
                                                               paste0(format(round(Value.1,2), nsmall=2), " ($\\textcolor{darkgreen}{\\Downarrow}$ ", format(round(abs(Value.Diff),2), nsmall=2), ")"),
                                                               ifelse(Value.Diff <0,
                                                                      paste0(format(round(Value.1,2), nsmall=2), " ($\\textcolor{red}{\\Uparrow}$ ", format(round(abs(Value.Diff),2), nsmall=2), ")"),
                                                                      paste0(format(round(Value.1,2), nsmall=2), " ($\\textcolor{blue}{=}", format(round(abs(Value.Diff),2), nsmall=2), "$)"))),
                                      PPG.Rank.Shift = ifelse(Rank.Diff > 0, 
                                                              paste0(Rank.1, " ($\\textcolor{darkgreen}{\\Uparrow}$ ", abs(Rank.Diff), ")"),
                                                              ifelse(Rank.Diff <0,
                                                                     paste0(Rank.1, " ($\\textcolor{red}{\\Downarrow}$ ", abs(Rank.Diff), ")"),
                                                                     paste0(Rank.1, " ($\\textcolor{blue}{=}", abs(Rank.Diff), "$)")))) %>%
    dplyr::select(-Value, -Rank,  -Rank.Diff, -Value.Diff) %>%
    left_join(our.rank.mat)
  
  final.table <- final.table %>% mutate(TO.Rank.Value = paste0(TO.Rank, " (", TO.Value, ")"),
                                        Yds.Rank.Value = paste0(Yds.Rank, " (", Yds.Value, ")"))
  
}


final.table$Value.1 <- format(round(final.table$Value.1, 2), nsmall=2)

## Whether to include ALL 131 the teams in the LATEX-converted table (NOT RECOMMENDED)
include.all <- FALSE



#######
## Generating the Latex code to be copy-pasted into the tabular environment
## to reproduce the tables
#######

  if (side == "Offense"){
    cat(paste0(c("Team \\ ", "Points Scored (Per Drive)", "Defense Statistics (Per Drive)"), collapse = " & "), "\\\\ \\hline")
  }
  
  
  if (side == "Defense"){
    cat(paste0(c("Team \\ ", "Points Allowed (Per Drive)", "Offense Statistics (Per Drive)"), collapse = " & "), "\\\\ \\hline")
  }
  
  
  if (!include.all){
    
    cat(" \n \\begin{tabular}{r} \\\\ \\\\", paste0(final.table$Team[which.ind[1:ifelse(side == "Offense", top.off, top.def)]], collapse = " \\\\ "),  # " \\\\ .... \\\\",
        if (length(which.ind) != top.off) {
          c(" \\\\ .... \\\\", paste0(final.table$Team[which.ind[-c(1:ifelse(side == "Offense", top.off, top.def))]], collapse = " \\\\ .... \\\\ "))} else NULL, " \\\\ .... \\\\",
        "\\end{tabular} &")
    
    ## The Offensive Pts Per Drive shifts
    cat("\n \\begin{tabular}{rr}",
      "\\\\ Val (Shift) \\ \\ & Rk (Shift) \\\\ \\hline ")    
    
    final.string <- NULL
    for (i in which.ind[1:ifelse(side == "Offense", top.off, top.def)]){
      final.string <- paste0(final.string, paste0(final.table[i, c(6,7)],  collapse=" & "), " \\\\ ")
    }
    
    for (i in which.ind[-c(1:ifelse(side == "Offense", top.off, top.def))]){
      final.string <- paste0(final.string, " ... & ... \\\\", paste0(final.table[i, c(6,7)],  collapse=" & "), " \\\\ ")
    }
    
    final.string <- paste0(final.string,  "... & ... \\\\", " \\end{tabular} & ")
    cat(final.string)
    
    
    if (side == "Offense"){
      cat("\n \\begin{tabular}{rr} ",
          "Takeaways & Yards Allowed \\\\ \\hline Rk (Val) & Rk (Val) \\\\ \\hline ")
    } else {
      cat("\n \\begin{tabular}{rr} ",
          "Giveaways & Yards Gained \\\\ \\hline Rk (Val) & Rk (Val) \\\\ \\hline ")
    }
    
    final.string <- NULL
    for (i in which.ind[1:ifelse(side == "Offense", top.off, top.def)]){
      final.string <- paste0(final.string, paste0(c(final.table[i, 12:13]), collapse =" & "), " \\\\ ")
    }
    
    for (i in which.ind[-c(1:ifelse(side == "Offense", top.off, top.def))]){
      final.string <- paste0(final.string, "... & ... \\\\", paste0(c(final.table[i, 12:13]), collapse =" & "), " \\\\ ")
    }
    
    
    final.string <- paste0(final.string, "... & ... \\\\", " \\end{tabular}")
    
    cat(final.string)
    
  } else {
    
    cat(" \n \\begin{tabular}{r} \\\\", paste0(final.table$Team, collapse = " \\\\ "),  "\\end{tabular} &")
    
    ## The Offensive Pts Per Game shifts
    
    cat("\n \\begin{tabular}{rr}",
      "\\\\ Val (Shift) \\ \\ & Rk (Shift) \\\\ \\hline ")
    
    final.string <- NULL
    for (i in 1:nrow(final.table)){
      final.string <- paste0(final.string, paste0(final.table[i, c(6,7)],  collapse=" & "), " \\\\ ")
    }
    
    for (i in which.ind[-c(1:ifelse(side == "Offense", top.off, top.def))]){
      final.string <- paste0(final.string, "... & ... \\\\", paste0(final.table[i, c(6,7)],  collapse=" & "), " \\\\ ")
      
    }
    
    final.string <- paste0(final.string,  "... & ... \\\\", " \\end{tabular} & ")
    
    cat(final.string)
    
    ## The Defensive TDs ranks
    if (side == "Offense"){
      cat("\n \\begin{tabular}{rr} ",
          "Takeaways & Yards Allowed \\\\ \\hline Rk (Val) & Rk (Val) \\\\ \\hline ")
    } else {
      cat("\n \\begin{tabular}{rr} ",
          "Giveaways & Yards Gained \\\\ \\hline Rk (Val) & Rk (Val) \\\\ \\hline ")
    }
    
    final.string <- NULL
    for (i in 1:nrow(final.table)){
      final.string <- paste0(final.string, paste0(final.table[i, 12:13],  collapse=" & "), " \\\\ ")
    }
    
    for (i in which.ind[-c(1:ifelse(side == "Offense", top.off, top.def))]){
      final.string <- paste0(final.string, "... & ...  \\\\", paste0(c(final.table[i, 12:13]), collapse =" & "), " \\\\ ")
    }
    
    
    final.string <- paste0(final.string, "... & ... \\\\", " \\end{tabular}")

    cat(final.string)
  }
