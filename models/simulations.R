# Simulating networks and running extinction simulations
setwd("~/Documents/Uni/M.sc/Master Thesis/Networks/models/")

library(tapnet)
source("tapnet_helper.R")
source("simnetfromtap_aug.R")
source("simulate_tapnet.R")
source("helper_functions.R")
library(bipartite)
library(phytools)
# library(vegan)
library(MPSEM)
library(purrr)
library(dplyr)
library(data.table)

# read species names and make df 
animals <- read.csv("animal-names.csv", colClasses = c("NULL", "character"),
                    stringsAsFactors = F)
plants <- read.csv("flower-and-plant-names.csv", stringsAsFactors = F)
sp_names <- data.frame("plants" = plants,
                       "animals" = animals[1:nrow(plants),],
                       stringsAsFactors = F)

## initial simulation ----
# set no of nets and simulations
n_webs <- 1001
n_nets <- 4
n_sims <- 10
ctrbs <- c("Atl" = 1, "aTl" = 2, "atL" = 3, "ATL" = 4)
rew_names <- c("abund", "trait", "phylo")

# create initial network for community variables
init_sim <- map(seq(n_webs), ~ simulate_tapnet_aug(nlower = 20, nhigher = 30, ntraits_nopem = 2,
                                                   ntraits_pem = 2, abuns = "lognormal", Nobs = 1111,
                                                   names = sp_names, initial_sim = T))


# Create simulated networks with specified contributions ----
sims <- list()

# list w/ ctrb_vecs
ctrb_list <- list("Atl" = c("high", "low", "low"),
                  "aTl" = c("low", "high", "low"),
                  "atL" = c("low", "low", "high"),
                  "ATL" = c("high", "high", "high"))

sims <- map(.x = seq(n_webs), function(x) {
  map(.x = ctrb_list, ~ simnetfromtap_ctrb(ctrb_vec = .x,
                                           traits = pluck(init_sim, x)$traits_all,
                                           abuns = pluck(init_sim, x)$networks[[1]]$abuns,
                                           paramsList = pluck(init_sim, x)$sim_params,
                                           pems = pluck(init_sim, x)$networks[[1]]$pems,
                                           tmatch_type_pem = "normal",
                                           tmatch_type_obs = "normal",
                                           Nobs = 1111))})

## rewiring partner choice ----
# abundances
abunds <- map(seq(n_webs), ~ pluck(init_sim, .x)$networks[[1]]$abuns)

# traits 
traits <- map(seq(n_webs), ~ pluck(init_sim, .x)$networks[[1]]$traits)

# phylogenetic distances
library(ape)

phylos <- map(seq(n_webs), ~ list("low" = cophenetic.phylo(pluck(init_sim, .x)$trees$low), 
                                  "high" = cophenetic.phylo(pluck(init_sim, .x)$trees$high)))

## run extinction models ----
source("rewiring_vizentin-bugoni_2019/Functions/extinction.mod.R")
source("one.second.extinct.mod.R")
source("rewiring_vizentin-bugoni_2019/Functions/calc.mean.one.second.extinct.mod.R")
source("rewiring_vizentin-bugoni_2019/Functions/matrix.p1.R")
source("rewiring_vizentin-bugoni_2019/Functions/IC.R") # calc of 95 percent confidence interval

# only use Atl, aTl, and atL for extc sims
n_nets <- 3
sims_all <- sims # copy to preserve all contribution importances 
sims <- map(seq(n_webs), ~ pluck(sims, .x)[-4]) # delete ATL

# original web lower level
extc_sims_lower_org <- map(seq(n_webs), function(x) {map2(.x = list("abund" = abunds,
                                                                    "trait" = traits,
                                                                    "phylo" = phylos),
                                                          .y = c("abund", "trait", "phylo"),
                                                          ~ replicate(n_sims, simplify = F,
                                                                      one.second.extinct.mod_aug(web = pluck(init_sim, x)$networks[[1]]$web,
                                                                                                 participant = "lower",
                                                                                                 method = "random",
                                                                                                 rewiring = T,
                                                                                                 partner.choice = pluck(.x, x),
                                                                                                 interactions = pluck(init_sim, x)$networks[[1]]$I_mat,
                                                                                                 method.rewiring = .y)))})

# original web higher level
extc_sims_higher_org <- map(seq(n_webs), function(x) {map2(.x = list("abund" = abunds,
                                                                     "trait" = traits,
                                                                     "phylo" = phylos),
                                                           .y = c("abund", "trait", "phylo"),
                                                           ~ replicate(n_sims, simplify = F,
                                                                       one.second.extinct.mod_aug(web = pluck(init_sim, x)$networks[[1]]$web,
                                                                                                  participant = "higher",
                                                                                                  method = "random",
                                                                                                  rewiring = T,
                                                                                                  partner.choice = pluck(.x, x),
                                                                                                  interactions = pluck(init_sim, x)$networks[[1]]$I_mat,
                                                                                                  method.rewiring = .y)))})

# no rewiring, initial extinction on lower level
extc_sims_lower_norew <- map(seq(n_webs), function(x) {
  run_extc(web = pluck(sims, x),
           participant = "lower",
           method = "random",
           rewiring = F,
           partner.choice = NULL, 
           interactions = pluck(sims, x),
           method.rewiring = "NULL",
           n_sims = n_sims,
           multiple.webs = T)})

# no rewiring, initial extinction on higher level
extc_sims_higher_norew <- map(seq(n_webs), function(x) {
  run_extc(web = pluck(sims, x),
           participant = "higher",
           method = "random",
           rewiring = F,
           partner.choice = NULL, 
           interactions = pluck(sims, x),
           method.rewiring = "NULL",
           n_sims = n_sims,
           multiple.webs = T)})

# initial extinction on lower level
extc_sims_lower <- map(seq(n_webs), function(x) {
  map2(.x = list("abund" = abunds,
                 "trait" = traits,
                 "phylo" = phylos),
       .y = c("abund", "trait", "phylo"),
       ~ run_extc(web = pluck(sims, x),
                  participant = "lower",
                  method = "random",
                  rewiring = T,
                  partner.choice = pluck(.x, x), 
                  interactions = pluck(sims, x),
                  method.rewiring = .y,
                  n_sims = n_sims,
                  multiple.webs = T))})

# initial extinction on higher level
extc_sims_higher <- map(seq(n_webs), function(x) {
  map2(.x = list("abund" = abunds,
                 "trait" = traits,
                 "phylo" = phylos),
       .y = c("abund", "trait", "phylo"),
       ~ run_extc(web = pluck(sims, x),
                  participant = "higher",
                  method = "random",
                  rewiring = T,
                  partner.choice = pluck(.x, x), 
                  interactions = pluck(sims, x),
                  method.rewiring = .y,
                  n_sims = n_sims,
                  multiple.webs = T))})

## compute means of all simulations ----
# initial extinction on lower level org
extc_sims_lower_mean_org <- map(seq(n_webs), function(x) {
  list("lower" = map(rew_names,
                     ~ list_mean(pluck(extc_sims_lower_org, x),
                                 y = .x,
                                 original = T)),
       "higher" = map(rew_names,
                      ~ list_mean(pluck(extc_sims_lower_org, x),
                                  y = .x,
                                  lower = F,
                                  original = T)))})

# initial extinction on lower level org
extc_sims_higher_mean_org <- map(seq(n_webs), function(x) {
  list("lower" = map(rew_names,
                     ~ list_mean(pluck(extc_sims_higher_org, x),
                                 y = .x,
                                 original = T)),
       "higher" = map(rew_names,
                      ~ list_mean(pluck(extc_sims_higher_org, x),
                                  y = .x,
                                  lower = F,
                                  original = T)))})

# initial extinction on lower level, no rewiring
extc_sims_lower_mean_norew <- map(seq(n_webs), function(x) {
  list("lower" = map(1:3,
                     ~ list_mean(pluck(extc_sims_lower_norew, x),
                                 y = .x,
                                 original = T)),
       "higher" = map(1:3,
                      ~ list_mean(pluck(extc_sims_lower_norew, x),
                                  y = .x,
                                  lower = F,
                                  original = T)))})

# initial extinction on higher level, no rewiring
extc_sims_higher_mean_norew <- map(seq(n_webs), function(x) {
  list("lower" = map(1:3,
                     ~ list_mean(pluck(extc_sims_higher_norew, x),
                                 y = .x,
                                 original = T)),
       "higher" = map(1:3,
                      ~ list_mean(pluck(extc_sims_higher_norew, x),
                                  y = .x,
                                  lower = F,
                                  original = T)))})

# initial extinction on lower level
extc_sims_lower_mean <- map(seq(n_webs), function(x) {
  list("lower" = map(rew_names,
                     ~ list_mean(pluck(extc_sims_lower, x),
                                 y = .x)),
       "higher" = map(rew_names,
                      ~ list_mean(pluck(extc_sims_lower, x),
                                  y = .x,
                                  lower = F)))})

# initial extiction on higher level
extc_sims_higher_mean <- map(seq(n_webs), function(x) {
  list("lower" = map(rew_names,
                     ~ list_mean(pluck(extc_sims_higher, x),
                                 y = .x)),
       "higher" = map(rew_names,
                      ~ list_mean(pluck(extc_sims_higher, x),
                                  y = .x,
                                  lower = F)))})

## calculate percentages of remaining sp ----
# initial extinction on lower level original web
sp_remain_lower_org <- map(seq(n_webs), function(x) {
  list("lower" = map(1:3, ~ per_surv(pluck(extc_sims_lower_mean_org, x),
                                     y = .x,
                                     original = T)),
       "higher" = map(1:3, ~ per_surv(pluck(extc_sims_lower_mean_org, x),
                                      y = .x,
                                      lower = F,
                                      original = T)))})

# initial extinction on higher level original web
sp_remain_higher_org <- map(seq(n_webs), function(x) {
  list("lower" = map(1:3, ~ per_surv(pluck(extc_sims_higher_mean_org, x),
                                     y = .x,
                                     original = T)),
       "higher" = map(1:3, ~ per_surv(pluck(extc_sims_higher_mean_org, x),
                                      y = .x,
                                      lower = F,
                                      original = T)))})

# initial extinction on lower level, no rewiring
sp_remain_lower_norew <- map(seq(n_webs), function(x) {
  list("lower" = map(1:3, ~ per_surv(pluck(extc_sims_lower_mean_norew, x),
                                     y = .x,
                                     original = T)),
       "higher" = map(1:3, ~ per_surv(pluck(extc_sims_lower_mean_norew, x),
                                      y = .x,
                                      lower = F,
                                      original = T)))})

# initial extinction on higher level, no rewiring
sp_remain_higher_norew <- map(seq(n_webs), function(x) {
  list("lower" = map(1:3, ~ per_surv(pluck(extc_sims_higher_mean_norew, x),
                                     y = .x,
                                     original = T)),
       "higher" = map(1:3, ~ per_surv(pluck(extc_sims_higher_mean_norew, x),
                                      y = .x,
                                      lower = F,
                                      original = T)))})

# initial extinction on lower level
sp_remain_lower <- map(seq(n_webs), function(x) {
  list("lower" = map(1:3, ~ per_surv(pluck(extc_sims_lower_mean, x),
                                     y = .x)),
       "higher" = map(1:3, ~ per_surv(pluck(extc_sims_lower_mean, x),
                                      y = .x,
                                      lower = F)))})

# initial extinction on higher level
sp_remain_higher <- map(seq(n_webs), function(x) {
  list("lower" = map(1:3, ~ per_surv(pluck(extc_sims_higher_mean, x),
                                     y = .x)),
       "higher" = map(1:3, ~ per_surv(pluck(extc_sims_higher_mean, x),
                                      y = .x,
                                      lower = F)))})

## mean over all webs ----
# initial extincion on lower level original web
sp_remain_lower_web_mean_org <- list("lower" =
                                       map(.x = 1:3, function(x) {
                                         map(seq(n_webs), ~ pluck(sp_remain_lower_org, .x, "lower", x))  %>%
                                           as.data.table() %>%
                                           match_lengths() %>%
                                           rowMeans()}),
                                     "higher" = map(.x = 1:3, function(x){
                                       map(seq(n_webs), ~ pluck(sp_remain_lower_org, .x, "higher", x))  %>%
                                         as.data.table() %>%
                                         match_lengths() %>%
                                         rowMeans()}))

# initial extincion on lower level original web
sp_remain_higher_web_mean_org <- list("lower" =
                                        map(.x = 1:3, function(x) {
                                          map(seq(n_webs), ~ pluck(sp_remain_higher_org, .x, "lower", x))  %>%
                                            as.data.table() %>%
                                            match_lengths() %>%
                                            rowMeans()}),
                                      "higher" = map(.x = 1:3, function(x){
                                        map(seq(n_webs), ~ pluck(sp_remain_higher_org, .x, "higher", x))  %>%
                                          as.data.table() %>%
                                          match_lengths() %>%
                                          rowMeans()}))

# initial extincion on lower level, no rewiring
sp_remain_lower_web_mean_norew <-list("lower" =
                                        map(.x = 1:3, function(x) {
                                          map(seq(n_webs), ~ pluck(sp_remain_lower_norew, .x, "lower", x))  %>%
                                            as.data.table() %>%
                                            match_lengths() %>%
                                            rowMeans()}),
                                      "higher" = map(.x = 1:3, function(x){
                                        map(seq(n_webs), ~ pluck(sp_remain_lower_norew, .x, "higher", x))  %>%
                                          as.data.table() %>%
                                          match_lengths() %>%
                                          rowMeans()}))


# initial extincion on higher level
sp_remain_higher_web_mean_norew <- list("lower" =
                                          map(.x = 1:3, function(x) {
                                            map(seq(n_webs), ~ pluck(sp_remain_higher_norew, .x, "lower", x))  %>%
                                              as.data.table() %>%
                                              match_lengths() %>%
                                              rowMeans()}),
                                        "higher" = map(.x = 1:3, function(x){
                                          map(seq(n_webs), ~ pluck(sp_remain_higher_norew, .x, "higher", x))  %>%
                                            as.data.table() %>%
                                            match_lengths() %>%
                                            rowMeans()}))

# initial extincion on lower level
sp_remain_lower_web_mean <- list("lower" =
                                   map(.x = 1:3, function(x) map(.x = 1:3, function(y) {
                                     map(seq(n_webs), ~ pluck(sp_remain_lower, .x, "lower", y, x))  %>%
                                       as.data.table() %>%
                                       match_lengths() %>%
                                       rowMeans()})),
                                 "higher" = map(.x = 1:3, function(x) map(.x = 1:3, function(y) {
                                   map(seq(n_webs), ~ pluck(sp_remain_lower, .x, "higher", y, x))  %>%
                                     as.data.table() %>%
                                     match_lengths() %>%
                                     rowMeans()})))


# initial extincion on higher level
sp_remain_higher_web_mean <- list("lower" =
                                    map(.x = 1:3, function(x) map(.x = 1:3, function(y) {
                                      map(seq(n_webs), ~ pluck(sp_remain_higher, .x, "lower", y, x))  %>%
                                        as.data.table() %>%
                                        match_lengths() %>%
                                        rowMeans()})),
                                  "higher" = map(.x = 1:3, function(x) map(.x = 1:3, function(y) {
                                    map(seq(n_webs), ~ pluck(sp_remain_higher, .x, "higher", y, x))  %>%
                                      as.data.table() %>%
                                      match_lengths() %>%
                                      rowMeans()})))
