# first test using tapnet to simulate network
setwd("~/Documents/Uni/M.sc/Master Thesis/Networks/models/")

library(tapnet)
source("tapnet_helper.R")
source("simnetfromtap_aug.R")
source("simulate_tapnet.R")
library(bipartite)
library(phytools)
# library(vegan)
library(ape)
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

# test simulation -----
# simnet <- simulate_tapnet(nlower = 10, nhigher = 10, ntraits_nopem = 2,
#                           ntraits_pem = 0, abuns = "lognormal", Nobs = 420)
# 
# # visualize phylo
# par(mfrow = c(2,1))
# phytools::plotTree(simnet$trees$low)
# phytools::plotTree(simnet$trees$high)
# par(mfrow = c(1,1))
# 
# # computing phylogenetic distances
# cophenetic(simnet$trees$low)
# cophenetic(simnet$trees$high)
# 
# # overview of network matrices
# head(simnet$networks[[1]]$abuns)
# head(simnet$networks[[1]]$traits)
# # head(simnet$networks[[1]]$pems) # phylogenetic eigenmatrices
# head(simnet$networks[[1]]$I_mat)
# 
# # visualize network
# plotweb(simnet$networks[[1]]$I_mat)
# visweb(simnet$networks[[1]]$I_mat)

# initial simulation ----
# create initial network for community variables
init_sim <- simulate_tapnet_aug(nlower = 20, nhigher = 30, ntraits_nopem = 2,
                                ntraits_pem = 2, abuns = "lognormal", Nobs = 1111,
                                names = sp_names, initial_sim = T)

# set no of webs, nets and simulations
n_nets <- 4
n_sims <- 100

# Create simulated networks with setting contributions ----
source("simnetfromtap_ctrb.R")
#source("clean.R")

sims <- list()

# list w/ ctrb_vecs
ctrb_list <- list("Atl" = c("high", "low", "low"),
                  "aTl" = c("low", "high", "low"),
                  "atL" = c("low", "low", "high"),
                  "ATL" = c("high", "high", "high"))

sims <- map(.x = ctrb_list, ~ simnetfromtap_ctrb(ctrb_vec = .x,
                                                        traits = init_sim$traits_all,
                                                        abuns = init_sim$networks[[1]]$abuns,
                                                        paramsList = init_sim$sim_params,
                                                        pems = init_sim$networks[[1]]$pems,
                                                        tmatch_type_pem = "normal",
                                                        tmatch_type_obs = "normal",
                                                        Nobs = 1111))

#sims <- clean(sims, single = F) # DON'T RUN ! Causes false indexing in extc. sim

# visualize networks ----
map(map(sims, 1), ~ plotweb(.x))

# basic metrics ----
# connectance
connectance <- function(x) {
  sum(x != 0) / (nrow(x) * ncol(x))
}

connectance(sim$networks[[1]]$web)
connectance(sim1web)
connectance(sim2web)
connectance(sim3web)
connectance(sim4web)

# nestedness; nested rank
nestedness(sim$networks[[1]]$web)
nestedness(sim1web)
nestedness(sim2web)
nestedness(sim3web)

# rewiring probabilites ----
# # function to normalize values to range 0,1
# feature_scale <- function(x, method = "simple") {
#   if (method == "simple") {
#     out <- (abs(x)/max(abs(x)))
#   }
#   
#   if (method == "minmax") {
#     out <- (abs(x) - min(abs(x))) / (max(abs(x)) - min(abs(x)))
#   }
#   return(out)
# }
# 
# # pairwise relative abundances
# rew_abund <- init_sim$rew_probs[[1]]$A_mat
# 
# # relative abundances of lower
# rew_abund_low <- sweep(rew_abund, 2, colSums(rew_abund), "/")
# 
# # relative abundances of higher
# rew_abund_high <- sweep(rew_abund, 1, rowSums(rew_abund), "/")
# 
# # traits
# rew_trait <- init_sim$rew_probs[[1]]$T_mat
# 
# # relative traits of lower
# rew_trait_low <- sweep(rew_trait, 2, colSums(rew_trait), "/")
# 
# # relative traits of higher
# rew_trait_high <- sweep(rew_trait, 1, rowSums(rew_trait), "/")
# 
# # phylogeny
# rew_phylo <- init_sim$rew_probs[[1]]$L_mat
# 
# # relative latent traits of lower
# rew_phylo_low <- sweep(rew_phylo, 2, colSums(rew_phylo), "/")
# 
# # relative latent traits of higher
# rew_phylo_high <- sweep(rew_phylo, 1, rowSums(rew_phylo), "/")

# rewiring partner choice ----
  # abundances
  abunds <- init_sim$networks[[1]]$abuns
  
  # traits 
  traits <- init_sim$networks[[1]]$traits
  
  # phylogenetic distances
  library(ape)
  
  phylos <- list("low" = cophenetic.phylo(init_sim$trees$low), 
                    "high" = cophenetic.phylo(init_sim$trees$high))
  
# run extinction models ----

source("rewiring_vizentin-bugoni_2019/Functions/extinction.mod.R")
source("one.second.extinct.mod.R")
source("rewiring_vizentin-bugoni_2019/Functions/calc.mean.one.second.extinct.mod.R")
source("rewiring_vizentin-bugoni_2019/Functions/matrix.p1.R")
source("rewiring_vizentin-bugoni_2019/Functions/IC.R") # calc of 95 percent confidence interval

# fy to simulate extinctions for a nested list of networks; n_sims is no of extc
# simulations, n_nets is no of networks in list, n_webs is no of webs per network
run_extc <- function(web,
                     participant,
                     method,
                     rewiring,
                     partner.choice,
                     interactions,
                     method.rewiring,
                     n_sims,
                     y) {
  
map(web, ~replicate(n_sims, simplify = F, 
                  one.second.extinct.mod_aug(web = pluck(.x, 1), 
                                             participant = participant,
                                             method = method,
                                             rewiring = rewiring,
                                             partner.choice = partner.choice,
                                             interactions = pluck(.x, 2),
                                             method.rewiring = method.rewiring)))
}

# only use Atl, aTl, and atL for extc sims
n_nets <- 3 
sims <- sims[-4] # delete ATL

# original web
extc_sims_lower_org <- map2(.x = list("abund" = abunds,
                                      "trait" = traits,
                                      "phylo" = phylos),
                            .y = c("abund", "trait", "phylo"),
                            ~ replicate(n_sims, simplify = F, 
                                 one.second.extinct.mod_aug(web = init_sim$networks[[1]]$web, 
                                                            participant = "lower",
                                                            method = "random",
                                                            rewiring = T,
                                                            partner.choice = .x,
                                                            interactions = init_sim$networks[[1]]$I_mat,
                                                            method.rewiring = .y)))


# initial extinction on lower level
extc_sims_lower <- map2(.x = list("abund" = abunds,
                                  "trait" = traits,
                                  "phylo" = phylos),
                        .y = c("abund", "trait", "phylo"),
                        ~ run_extc(web = sims,
                                   participant = "lower",
                                   method = "random",
                                   rewiring = T,
                                   partner.choice = .x, 
                                   interactions = sims,
                                   method.rewiring = .y,
                                   n_sims = n_sims))

# initial extiction on higher level
extc_sims_higher <- map2(.x = list("abund" = abunds,
                                   "trait" = traits,
                                   "phylo" = phylos),
                        .y = c("abund", "trait", "phylo"),
                        ~ run_extc(web = sims,
                                   participant = "higher",
                                   method = "random",
                                   rewiring = T,
                                   partner.choice = .x,
                                   interactions = sims,
                                   method.rewiring = .y,
                                   n_sims = n_sims))

# compute means of all simulations
source("list_means.R")

rew_names <- c("abund" = 1, "trait" = 2, "phylo" = 3)

# initial extinction on lower level org
extc_sims_lower_org_mean <- list("lower" = list_mean(extc_sims_lower_org)),
                                 "higher" = map(rew_names, ~ list_mean(extc_sims_lower_org, y = .x, lower = F)))

# initial extinction on lower level
extc_sims_lower_mean <- list("lower" = map(rew_names, ~ list_mean(extc_sims_lower, y = .x)),
                             "higher" = map(rew_names, ~ list_mean(extc_sims_lower, y = .x, lower = F)))

# initial extiction on higher level
extc_sims_higher_mean <- list("lower" = map(rew_names, ~ list_mean(extc_sims_higher, y = .x)),
                             "higher" = map(rew_names, ~ list_mean(extc_sims_higher, y = .x, lower = F)))

# function to calculate percentages of remaining species
per_surv <- function(x, y, lower = T) {
  list_divide <- function(x) {
    x / x[1] * 100
  }
  
  extc_in_network <- ifelse(isTRUE(lower), 1, 2)
  
  map(c("Atl" = 1, "aTl" = 2, "atL" = 3),
      ~ list_divide(pluck(x, extc_in_network, y, .x)))
}


# calculate percentages of remaining sp
# initial extinction on lower level
sp_remain_lower <- list("lower" = map(rew_names, ~ per_surv(extc_sims_lower_mean, y = .x)),
                        "higher" = map(rew_names, ~ per_surv(extc_sims_lower_mean, y = .x, lower = F)))

# initial extinction on higher level
sp_remain_higher <- list("lower" = map(rew_names, ~ per_surv(extc_sims_higher_mean, y = .x)),
                         "higher" = map(rew_names, ~ per_surv(extc_sims_higher_mean, y = .x, lower = F)))

# visualize extinction models ----
library(ggplot2)

plot_extc <- function(x){
  com_vars <- c("Atl" = 1, "aTl" = 2, "atL" = 3)
  map(com_vars, ~ ggplot() + 
        geom_line(aes(rev(pluck(x, 1, 1, .x)), pluck(x, 2, 1, .x),
                      color = "Abundance"), linetype = 1) +
        geom_line(aes(rev(pluck(x, 1, 2, .x)), pluck(x, 2, 2, .x),
                      color = "Traits"), linetype = 2) +
        geom_line(aes(rev(pluck(x, 1, 3, .x)), pluck(x, 2, 3, .x),
                      color = "Phylogeny"), linetype = 3) +
        scale_color_manual(name = "Rewiring Method",
                           values = c("Abundance" = "black",
                                      "Traits" = "firebrick",
                                      "Phylogeny" = "dodgerblue")) +
        labs(x = "plants removed", y = "animals persisting",
             title = paste("Extinction cascade", names(com_vars[.x])))
      )
}

plot_extc(sp_remain_lower)
plot_extc(sp_remain_higher)
