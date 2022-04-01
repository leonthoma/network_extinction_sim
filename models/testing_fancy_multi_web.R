# first test using tapnet to simulate network
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

# initial simulation ----
# create initial network for community variables
init_sim <- map(1:10, ~ simulate_tapnet_aug(nlower = 20, nhigher = 30, ntraits_nopem = 2,
                                ntraits_pem = 2, abuns = "lognormal", Nobs = 1111,
                                names = sp_names, initial_sim = T))

# set no of webs, nets and simulations
n_nets <- 4
n_sims <- 10

# Create simulated networks with setting contributions ----
sims <- list()

# list w/ ctrb_vecs
ctrb_list <- list("Atl" = c("high", "low", "low"),
                  "aTl" = c("low", "high", "low"),
                  "atL" = c("low", "low", "high"),
                  "ATL" = c("high", "high", "high"))

sims <- map(.x = 1:10, function(x) {map(.x = ctrb_list, ~ simnetfromtap_ctrb(ctrb_vec = .x,
                                                        traits = pluck(init_sim, x)$traits_all,
                                                        abuns = pluck(init_sim, x)$networks[[1]]$abuns,
                                                        paramsList = pluck(init_sim, x)$sim_params,
                                                        pems = pluck(init_sim, x)$networks[[1]]$pems,
                                                        tmatch_type_pem = "normal",
                                                        tmatch_type_obs = "normal",
                                                        Nobs = 1111))})


# visualize networks ----
# map(map(sims, 1), ~ plotweb(.x))

# basic metrics ----
# # connectance
# connectance <- function(x) {
#   sum(x != 0) / (nrow(x) * ncol(x))
# }
# 
# connectance(sim$networks[[1]]$web)
# connectance(sim1web)
# connectance(sim2web)
# connectance(sim3web)
# connectance(sim4web)
# 
# # nestedness; nested rank
# nestedness(sim$networks[[1]]$web)
# nestedness(sim1web)
# nestedness(sim2web)
# nestedness(sim3web)

# rewiring partner choice ----
  # abundances
  abunds <- map(1:10, ~ pluck(init_sim, .x)$networks[[1]]$abuns)
  
  # traits 
  traits <- map(1:10, ~ pluck(init_sim, .x)$networks[[1]]$traits)
  
  # phylogenetic distances
  library(ape)
  
  phylos <- map(1:10, ~ list("low" = cophenetic.phylo(pluck(init_sim, .x)$trees$low), 
                    "high" = cophenetic.phylo(pluck(init_sim, .x)$trees$high)))
  
# run extinction models ----
source("rewiring_vizentin-bugoni_2019/Functions/extinction.mod.R")
source("one.second.extinct.mod.R")
source("rewiring_vizentin-bugoni_2019/Functions/calc.mean.one.second.extinct.mod.R")
source("rewiring_vizentin-bugoni_2019/Functions/matrix.p1.R")
source("rewiring_vizentin-bugoni_2019/Functions/IC.R") # calc of 95 percent confidence interval

# only use Atl, aTl, and atL for extc sims
n_nets <- 3
sims <- map(1:10, ~ pluck(sims, .x)[-4]) # delete ATL

# # original web lower level
# extc_sims_lower_org <- map(1:10, function(x) {map2(.x = list("abund" = abunds,
#                                       "trait" = traits,
#                                       "phylo" = phylos),
#                             .y = c("abund", "trait", "phylo"),
#                             ~ replicate(n_sims, simplify = F, 
#                                  one.second.extinct.mod_aug(web = pluck(init_sim, x)$networks[[1]]$web, 
#                                                             participant = "lower",
#                                                             method = "random",
#                                                             rewiring = T,
#                                                             partner.choice = pluck(.x, x),
#                                                             interactions = pluck(init_sim, x)$networks[[1]]$I_mat,
#                                                             method.rewiring = .y)))})
# 
# # original web higher level
# extc_sims_higher_org <- map(1:10, function(x) {map2(.x = list("abund" = abunds,
#                                       "trait" = traits,
#                                       "phylo" = phylos),
#                             .y = c("abund", "trait", "phylo"),
# ~ replicate(n_sims, simplify = F,
#            one.second.extinct.mod_aug(web = pluck(init_sim, x)$networks[[1]]$web,
#                                       participant = "higher",
#                                       method = "random",
#                                       rewiring = T,
#                                       partner.choice = pluck(.x, x),
#                                       interactions = pluck(init_sim, x)$networks[[1]]$I_mat,
#                                       method.rewiring = .y)))})
# 

# initial extinction on lower level
extc_sims_lower <- map(1:10, function(x) {
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
extc_sims_higher <- map(1:10, function(x) {
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


# compute means of all simulations
rew_names <- c("abund", "trait", "phylo")

# # initial extinction on lower level org
# extc_sims_lower_org_mean <- map(1:10, function(x) {
#   list("lower" = map(rew_names,
#                      ~ list_mean(pluck(extc_sims_lower_org, x),
#                                  y = .x,
#                                  original = T)),
#        "higher" = map(rew_names,
#                       ~ list_mean(pluck(extc_sims_lower_org, x),
#                                   y = .x,
#                                   lower = F,
#                                   original = T)))})
# 
# # initial extinction on lower level org
# extc_sims_higher_org_mean <- map(1:10, function(x) {
#   list("lower" = map(rew_names,
#                      ~ list_mean(pluck(extc_sims_higher_org, x),
#                                  y = .x,
#                                  original = T)),
#        "higher" = map(rew_names,
#                       ~ list_mean(pluck(extc_sims_higher_org, x),
#                                   y = .x,
#                                   lower = F,
#                                   original = T)))})

# initial extinction on lower level
extc_sims_lower_mean <- map(1:10, function(x) {
  list("lower" = map(rew_names,
                     ~ list_mean(pluck(extc_sims_lower, x),
                                 y = .x)),
       "higher" = map(rew_names,
                      ~ list_mean(pluck(extc_sims_lower, x),
                                  y = .x,
                                  lower = F)))})

# initial extiction on higher level
extc_sims_higher_mean <- map(1:10, function(x) {
  list("lower" = map(rew_names,
                     ~ list_mean(pluck(extc_sims_higher, x),
                                 y = .x)),
       "higher" = map(rew_names,
                      ~ list_mean(pluck(extc_sims_higher, x),
                                  y = .x,
                                  lower = F)))})

# calculate percentages of remaining sp
# # initial extinction on lower level original web
# sp_remain_lower_org <- map(1:10, function(x) {
#   list("lower" = map(1:3, ~ per_surv(pluck(extc_sims_lower_org_mean, x),
#                                      y = .x,
#                                      original = T)),
#        "higher" = map(1:3, ~ per_surv(pluck(extc_sims_lower_org_mean, x),
#                                       y = .x,
#                                       lower = F,
#                                       original = T)))})
# 
# # initial extinction on higher level original web
# sp_remain_higher_org <- map(1:10, function(x) {
#   list("lower" = map(1:3, ~ per_surv(pluck(extc_sims_higher_org_mean, x),
#                                      y = .x,
#                                      original = T)),
#        "higher" = map(1:3, ~ per_surv(pluck(extc_sims_higher_org_mean, x),
#                                       y = .x,
#                                       lower = F,
#                                       original = T)))})

# initial extinction on lower level
sp_remain_lower <- map(1:10, function(x) {
  list("lower" = map(1:3, ~ per_surv(pluck(extc_sims_lower_mean, x),
                                     y = .x)),
       "higher" = map(1:3, ~ per_surv(pluck(extc_sims_lower_mean, x),
                                      y = .x,
                                      lower = F)))})

# initial extinction on higher level
sp_remain_higher <- map(1:10, function(x) {
  list("lower" = map(1:3, ~ per_surv(pluck(extc_sims_higher_mean, x),
                                     y = .x)),
       "higher" = map(1:3, ~ per_surv(pluck(extc_sims_higher_mean, x),
                                      y = .x,
                                      lower = F)))})

# calculating mean over every web
# initial extincion on lower level
sp_remain_lower <- list("lower" =
                          map(.x = 1:3, function(x) map(.x = 1:3, function(y) {
                            map(1:10, ~ pluck(sp_remain_lower, .x, "lower", y, x))  %>%
                              as.data.table() %>%
                              replace_duplicate() %>%
                              rowMeans()})),
                        "higher" = map(.x = 1:3, function(x) map(.x = 1:3, function(y) {
                          map(1:10, ~ pluck(sp_remain_lower, .x, "higher", y, x))  %>%
                            as.data.table() %>%
                            replace_duplicate() %>%
                            rowMeans()})))


# initial extincion on higher level
sp_remain_higher <- list("lower" =
                           map(.x = 1:3, function(x) map(.x = 1:3, function(y) {
                             map(1:10, ~ pluck(sp_remain_higher, .x, "lower", y, x))  %>%
                               as.data.table() %>%
                               replace_duplicate() %>%
                               rowMeans()})),
                        "higher" = map(.x = 1:3, function(x) map(.x = 1:3, function(y) {
                          map(1:10, ~ pluck(sp_remain_higher, .x, "higher", y, x))  %>%
                            as.data.table() %>%
                            replace_duplicate() %>%
                            rowMeans()})))

# visualize extinction models ----
plot_extc(sp_remain_lower)
plot_extc(sp_remain_higher)

plot_extc_alt(sp_remain_lower, lower = T)
plot_extc_alt(sp_remain_higher)
