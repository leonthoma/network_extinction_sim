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

### initial simulation ----
# set no of nets and simulations
n_webs <- 1000
n_nets <- 4
n_sims <- 10
ctrbs <- c("Atl" = 1, "aTl" = 2, "atL" = 3, "ATL" = 4)
rew_names <- c("abund", "trait", "phylo")
coextc_thr <- NULL

# create initial network for community variables
init_sim <- map(seq(n_webs), ~ simulate_tapnet_aug(nlower = 40, nhigher = 50, ntraits_nopem = 2,
                                                   ntraits_pem = 2, abuns = "lognormal", Nobs = 1111,
                                                   names = sp_names, initial_sim = T))


# Create simulated networks with specified contributions
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


# Delete sp w/ only dead interactions
clean_sims_web <- map(seq(n_webs), function(x) {
  map(1:4, ~ del_dead_int(pluck(sims, x, .x, "web")))
})

# crop webs to size of smallest network w/o dead interactions
sims_web <- equalize_sp(clean_sims_web)

# create list w/ webs and I_mat
sims <- map(seq(n_webs), function(x) {
  map(seq(1:4), ~ list("web" = pluck(sims_web, x, .x), 
         "I_mat" = pluck(sims, x, .x, "I_mat")))
         })

# rename list
sims <- modify_depth(sims, 1, 
                    ~ set_names(.x, nm = c("Atl",
                                           "aTl",
                                           "atL",
                                           "ATL")))

### Need to also crop interaction mat to size of clean sims !!

# # Create tapnet from tinoco data
# data(Tinoco)
# 
# # create tapnet object of all three networks from tinoco data
# tin_init <- make_tapnet(tree_low = plant_tree,
#                    tree_high = humm_tree,
#                    networks = networks$forest,
#                    abun_low = plant_abun$forest,
#                    abun_high = humm_abun$forest,
#                    traits_low = plant_traits,
#                    traits_high = humm_traits)
# 
# # simulate interactions matrices
# tin_imat <- simnetfromtap(traits = tin_init$networks[[1]]$traits,
#               abuns = tin_init$networks[[1]]$abuns,
#               paramsList = list("lat_low" = rep(1, 19),
#                                 "lat_high" = rep(1, 7),
#                                 "pem_shift" = 0,
#                                 "tmatch_width_pem" = 1,
#                                 "tmatch_width_obs" = rep(1, 4),
#                                 "tmatch_type_pem" = "normal",
#                                 "tmatch_type_obs" = "normal"),
#               pems = tin_init$networks[[1]]$pems,
#               tmatch_type_pem = "normal",
#               tmatch_type_obs = "normal")
# 
# # list of webs and imats for all tinoco networks
# tin <- list("web" = tin_init$networks[[1]]$web,
#             "I_mat" = tin_imat)

### rewiring partner choice ----
# abundances
abunds <- map(seq(n_webs), ~ pluck(init_sim, .x)$networks[[1]]$abuns)
# abunds_tin <- pluck(tin_init$networks, 1, "abuns")

# traits 
traits <- map(seq(n_webs), ~ pluck(init_sim, .x)$networks[[1]]$traits)
# traits_tin <- pluck(tin_init$networks, 1, "traits")

# phylogenetic distances
library(ape)

phylos <- map(seq(n_webs), ~ list("low" = cophenetic.phylo(pluck(init_sim, .x)$trees$low), 
                                  "high" = cophenetic.phylo(pluck(init_sim, .x)$trees$high)))

# phylos_tin <- list("low" = cophenetic.phylo(pluck(tin_init$trees$low)),
#                    "high" = cophenetic.phylo(pluck(tin_init$trees$high)))


### run extinction models ----
  source("rewiring_vizentin-bugoni_2019/Functions/extinction.mod.R")
  source("one.second.extinct.mod.R")
  
  # only use Atl, aTl, and atL for extc sims
  n_nets <- 3
  sims_all <- sims # copy to preserve all contribution importances 
  sims <- map(seq(n_webs), ~ pluck(sims, .x)[-4]) # delete ATL
  
  # initial extinction on lower lever; original
  extc_sims_lower_org <- map(seq(n_webs), function(x) {
    map2(.x = list("abund" = abunds,
                   "trait" = traits,
                   "phylo" = phylos),
         .y = c("abund", "trait", "phylo"),
         ~ replicate(n_sims, simplify = F,
                     one.second.extinct.mod.aug(
                       web = pluck(init_sim, x)$networks[[1]]$web,
                       participant = "lower",
                       method = "random",
                       rewiring = T,
                       abund.partner.choice = pluck(.x, x),
                       trait.partner.choice = pluck(.x, x),
                       phylo.partner.choice = pluck(.x, x),
                       interactions = pluck(init_sim, x)$networks[[1]]$I_mat,
                       method.rewiring = .y, 
                       coextc.thr = coextc_thr)))})
  
  # initial extinction on higher lever; original
  extc_sims_higher_org <- map(seq(n_webs), function(x) {
    map2(.x = list("abund" = abunds,
                   "trait" = traits,
                   "phylo" = phylos),
         .y = c("abund", "trait", "phylo"),
         ~ replicate(n_sims, simplify = F,
                     one.second.extinct.mod.aug(
                       web = pluck(init_sim, x)$networks[[1]]$web,
                       participant = "higher",
                       method = "random",
                       rewiring = T,
                       abund.partner.choice = pluck(.x, x),
                       trait.partner.choice = pluck(.x, x),
                       phylo.partner.choice = pluck(.x, x),
                       interactions = pluck(init_sim, x)$networks[[1]]$I_mat,
                       method.rewiring = .y, 
                       coextc.thr = coextc_thr)))})

  # initial extinction on lower lever; original no rewiring
  
  ## CAN ONLY BE ONE (!!) NETWORK  
  
  extc_sims_lower_org_norew <- map(seq(n_webs), function(x) {
    replicate(n_sims, simplify = F,
              one.second.extinct.mod.aug(
                web = pluck(init_sim, x)$networks[[1]]$web,
                participant = "lower",
                method = "random",
                rewiring = F,
                abund.partner.choice = NULL,
                trait.partner.choice = NULL,
                phylo.partner.choice = NULL,
                interactions = pluck(init_sim, x)$networks[[1]]$I_mat,
                method.rewiring = "NULL", 
                coextc.thr = coextc_thr))})
  
  # initial extinction on lower lever; original no rewiring
  extc_sims_higher_org_norew <- map(seq(n_webs), function(x) {
    replicate(n_sims, simplify = F,
              one.second.extinct.mod.aug(
                web = pluck(init_sim, x)$networks[[1]]$web,
                participant = "higher",
                method = "random",
                rewiring = F,
                abund.partner.choice = NULL,
                trait.partner.choice = NULL,
                phylo.partner.choice = NULL,
                interactions = pluck(init_sim, x)$networks[[1]]$I_mat,
                method.rewiring = "NULL", 
                coextc.thr = coextc_thr))})
  
  # # initial extinction on lower level, abundance driven extinction 
  # extc_sims_lower_abund <- map(seq(n_webs), function(x) {
  #   map2(.x = list("abund" = abunds,
  #                  "trait" = traits,
  #                  "phylo" = phylos),
  #        .y = c("abund", "trait", "phylo"),
  #        ~ run_extc(web = pluck(sims, x),
  #                   participant = "lower",
  #                   method = "abund",
  #                   rewiring = T,
  #                   partner.choice = pluck(.x, x), 
  #                   interactions = pluck(sims, x),
  #                   method.rewiring = .y,
  #                   n_sims = n_sims,
  #                   multiple.webs = T))})
  # 
  # # initial extinction on lower level, abundance driven extinction 
  # extc_sims_higher_abund <- map(seq(n_webs), function(x) {
  #   map2(.x = list("abund" = abunds,
  #                  "trait" = traits,
  #                  "phylo" = phylos),
  #        .y = c("abund", "trait", "phylo"),
  #        ~ run_extc(web = pluck(sims, x),
  #                   participant = "higher",
  #                   method = "abund",
  #                   rewiring = T,
  #                   partner.choice = pluck(.x, x), 
  #                   interactions = pluck(sims, x),
  #                   method.rewiring = .y,
  #                   n_sims = n_sims,
  #                   multiple.webs = T))})
  # 
  # # initial extinction on lower level, abundance driven extinction, no rewiring
  # extc_sims_lower_abund_norew <- map(seq(n_webs), function(x) {
  #   run_extc(web = pluck(sims, x),
  #                   participant = "lower",
  #                   method = "abund",
  #                   rewiring = F,
  #                   partner.choice = NULL, 
  #                   interactions = pluck(sims, x),
  #                   method.rewiring = "NULL",
  #                   n_sims = n_sims,
  #                   multiple.webs = T)})
  # 
  # # initial extinction on lower level, abundance driven extinction, no rewiring
  # extc_sims_higher_abund_norew <- map(seq(n_webs), function(x) {
  #   run_extc(web = pluck(sims, x),
  #                   participant = "higher",
  #                   method = "abund",
  #                   rewiring = F,
  #                   partner.choice = NULL, 
  #                   interactions = pluck(sims, x),
  #                   method.rewiring = "NULL",
  #                   n_sims = n_sims,
  #                   multiple.webs = T)})
  # 
  # # initial extinction on both levels
  # extc_sims_both <- map(seq(n_webs), function(x) {
  #   map2(.x = list("abund" = abunds,
  #                  "trait" = traits,
  #                  "phylo" = phylos),
  #        .y = c("abund", "trait", "phylo"),
  #        ~ run_extc(web = pluck(sims, x),
  #                   participant = "both",
  #                   method = "random",
  #                   rewiring = T,
  #                   partner.choice = pluck(.x, x), 
  #                   interactions = pluck(sims, x),
  #                   method.rewiring = .y,
  #                   n_sims = n_sims,
  #                   multiple.webs = T, 
  #                   coextc_thr = coextc_thr))})
  # 
  # # initial extinction on both levels, no rewiring
  # extc_sims_both_norew <- map(seq(n_webs), function(x) {
  #   run_extc(web = pluck(sims, x),
  #                   participant = "both",
  #                   method = "random",
  #                   rewiring = F,
  #                   partner.choice = NULL, 
  #                   interactions = pluck(sims, x),
  #                   method.rewiring = "NULL",
  #                   n_sims = n_sims,
  #                   multiple.webs = T, 
  #            coextc_thr = coextc_thr)})
  # 
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
                    abund.partner.choice = pluck(.x, x),
                    trait.partner.choice = pluck(.x, x),
                    phylo.partner.choice = pluck(.x, x),
                    interactions = pluck(sims, x),
                    method.rewiring = .y,
                    n_sims = n_sims,
                    multiple.webs = T, 
                    coextc.thr = coextc_thr))})
  
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
                    abund.partner.choice = pluck(.x, x),
                    trait.partner.choice = pluck(.x, x),
                    phylo.partner.choice = pluck(.x, x), 
                    interactions = pluck(sims, x),
                    method.rewiring = .y,
                    n_sims = n_sims,
                    multiple.webs = T, 
                    coextc.thr = coextc_thr))})
  
  # initial extinction on lower level, no rewiring
  extc_sims_lower_norew <- map(seq(n_webs), function(x) {
    run_extc(web = pluck(sims, x),
             participant = "lower",
             method = "random",
             rewiring = F,
             abund.partner.choice = NULL,
             trait.partner.choice = NULL,
             phylo.partner.choice = NULL, 
             interactions = pluck(sims, x),
             method.rewiring = "NULL",
             n_sims = n_sims,
             multiple.webs = T, 
             coextc.thr = coextc_thr)})
  
  # initial extinction on higher level, no rewiring
  extc_sims_higher_norew <- map(seq(n_webs), function(x) {
    run_extc(web = pluck(sims, x),
             participant = "higher",
             method = "random",
             rewiring = F,
             abund.partner.choice = NULL,
             trait.partner.choice = NULL,
             phylo.partner.choice = NULL,
             interactions = pluck(sims, x),
             method.rewiring = "NULL",
             n_sims = n_sims,
             multiple.webs = T, 
             coextc.thr = coextc_thr)})
  
  # # initial extinction on lower level; tinoco
  # extc_sims_lower_tin <- map2(.x = list("abund" = abunds_tin,
  #                                       "trait" = traits_tin,
  #                                       "phylo" = phylos_tin),
  #                             .y = c("abund", "trait", "phylo"),
  #   ~ replicate(n_sims, simplify = F,
  #               one.second.extinct.mod_aug(
  #                 web = pluck(tin, "web"),
  #                 participant = "lower",
  #                 method = "random",
  #                 rewiring = T,
  #                 partner.choice = pluck(.x),
  #                 interactions = pluck(tin, "I_mat"),
  #                 method.rewiring = .y,
  #                 coextc_thr = coextc_thr)))
  # 
  # # initial extinction on higher level; tinoco
  # extc_sims_higher_tin <- map2(.x = list("abund" = abunds_tin,
  #                                        "trait" = traits_tin,
  #                                        "phylo" = phylos_tin),
  #                              .y = c("abund", "trait", "phylo"),
  #        ~ replicate(n_sims, simplify = F,
  #                    one.second.extinct.mod_aug(
  #                      web = pluck(tin, "web"),
  #                      participant = "higher",
  #                      method = "random",
  #                      rewiring = T,
  #                      partner.choice = pluck(.x),
  #                      interactions = pluck(tin, "I_mat"),
  #                      method.rewiring = .y,
  #                      coextc_thr = coextc_thr)))
  # 
  # # initial extinction on lower lever; tinoco no rewiring
  # extc_sims_lower_tin_norew <- replicate(n_sims, simplify = F,
  #              one.second.extinct.mod_aug(web = pluck(tin, "web"),
  #                                         participant = "lower",
  #                                         method = "random",
  #                                         rewiring = F,
  #                                         partner.choice = NULL, 
  #                                         interactions = pluck(tin, "I_mat"),
  #                                         method.rewiring = "NULL",
  #                                         coextc_thr = coextc_thr))
  # 
  # # initial extinction on lower lever; tinoco no rewiring
  # extc_sims_higher_tin_norew <- replicate(n_sims, simplify = F,
  #             one.second.extinct.mod_aug(web = pluck(tin, "web"),
  #                                        participant = "higher",
  #                                        method = "random",
  #                                        rewiring = F,
  #                                        partner.choice = NULL, 
  #                                        interactions = pluck(tin, "I_mat"),
  #                                        method.rewiring = "NULL",
  #                                        coextc_thr = coextc_thr))

### compute means of all simulations ----
  # initial extinction on lower level
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
  
  ### add org norew
  # initial extinction on lower level; org, no rewiring
  extc_sims_lower_org_mean_norew <- map(seq(n_webs), function(x) {
    list("lower" = pluck(extc_sims_lower_org_norew, x) %>% as.data.table(.) %>% 
           select(., contains(c("no", "n.lower"))) %>%
           replace_duplicate(.) %>% rowMeans(.),
         "higher" = pluck(extc_sims_lower_org_norew, x) %>% as.data.table(.) %>% 
           select(., contains(c("no", "n.higher"))) %>%
           replace_duplicate(.) %>% rowMeans(.))})
  
  # initial extinction on lower level; org, no rewiring
  extc_sims_higher_org_mean_norew <- map(seq(n_webs), function(x) {
    list("lower" = pluck(extc_sims_higher_org_norew, x) %>% as.data.table(.) %>% 
           select(., contains(c("no", "n.lower"))) %>%
           replace_duplicate(.) %>% rowMeans(.),
         "higher" = pluck(extc_sims_higher_org_norew, x) %>% as.data.table(.) %>% 
           select(., contains(c("no", "n.higher"))) %>%
           replace_duplicate(.) %>% rowMeans(.))})
  
  # # initial extinction on lower level, abundance driven extinction 
  # extc_sims_lower_abund_mean <- map(seq(n_webs), function(x) {
  #   list("lower" = map(rew_names,
  #                      ~ list_mean(pluck(extc_sims_lower_abund, x),
  #                                  y = .x)),
  #        "higher" = map(rew_names,
  #                       ~ list_mean(pluck(extc_sims_lower_abund, x),
  #                                   y = .x,
  #                                   lower = F)))})
  # 
  # # initial extinction on lower level, abundance driven extinction 
  # extc_sims_higher_abund_mean <- map(seq(n_webs), function(x) {
  #   list("lower" = map(rew_names,
  #                      ~ list_mean(pluck(extc_sims_higher_abund, x),
  #                                  y = .x)),
  #        "higher" = map(rew_names,
  #                       ~ list_mean(pluck(extc_sims_higher_abund, x),
  #                                   y = .x,
  #                                   lower = F)))})
  # 
  # # initial extinction on lower level, abundance driven extinction, no rewiring
  # extc_sims_lower_abund_mean_norew <- map(seq(n_webs), function(x) {
  #   list("lower" = map(1:3,
  #                      ~ list_mean(pluck(extc_sims_lower_abund_norew, x),
  #                                  y = .x,
  #                                  original = T)),
  #        "higher" = map(1:3,
  #                       ~ list_mean(pluck(extc_sims_lower_abund_norew, x),
  #                                   y = .x,
  #                                   lower = F,
  #                                   original = T)))})
  # 
  # # initial extinction on lower level, abundance driven extinction, no rewiring
  # extc_sims_higher_abund_mean_norew <- map(seq(n_webs), function(x) {
  #   list("lower" = map(1:3,
  #                      ~ list_mean(pluck(extc_sims_higher_abund_norew, x),
  #                                  y = .x,
  #                                  original = T)),
  #        "higher" = map(1:3,
  #                       ~ list_mean(pluck(extc_sims_higher_abund_norew, x),
  #                                   y = .x,
  #                                   lower = F,
  #                                   original = T)))})
  # 
  # # initial extinction on both levels
  # extc_sims_both_mean <- map(seq(n_webs), function(x) {
  #   list("lower" = map(rew_names,
  #                      ~ list_mean(pluck(extc_sims_both, x),
  #                                  y = .x)),
  #        "higher" = map(rew_names,
  #                       ~ list_mean(pluck(extc_sims_both, x),
  #                                   y = .x,
  #                                   lower = F)))})
  # 
  # # initial extinction on both levels, no rewiring
  # extc_sims_both_mean_norew <- map(seq(n_webs), function(x) {
  #   list("lower" = map(1:3,
  #                      ~ list_mean(pluck(extc_sims_both_norew, x),
  #                                  y = .x,
  #                                  original = T)),
  #        "higher" = map(1:3,
  #                       ~ list_mean(pluck(extc_sims_both_norew, x),
  #                                   y = .x,
  #                                   lower = F,
  #                                   original = T)))})
  # 
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
  
  # # initial extinction on lower level; tinoco
  # extc_sims_lower_mean_tin <- list("lower" = map(rew_names,
  #                      ~ list_mean(extc_sims_lower_tin,
  #                                  y = .x,
  #                                  original = T)),
  #                      "higher" = map(rew_names,
  #                       ~ list_mean(extc_sims_lower_tin,
  #                                   y = .x,
  #                                   lower = F,
  #                                   original = T)))
  # 
  # # initial extinction on higher level; tinoco
  # extc_sims_higher_mean_tin <- list("lower" = map(rew_names,
  #                      ~ list_mean(extc_sims_higher_tin,
  #                                  y = .x,
  #                                  original = T)),
  #                      "higher" = map(rew_names,
  #                       ~ list_mean(extc_sims_higher_tin,
  #                                   y = .x,
  #                                   lower = F,
  #                                   original = T)))
  # 
  # # initial extinction on lower level; tinoco no rewiring
  # extc_sims_lower_mean_tin_norew <- list("lower" = map(1:3,
  #                      ~ list_mean(extc_sims_lower_tin_norew,
  #                                  y = .x,
  #                                  original = T)),
  #                      "higher" = map(1:3,
  #                       ~ list_mean(extc_sims_lower_tin_norew,
  #                                   y = .x,
  #                                   lower = F,
  #                                   original = T)))
  # 
  # # initial extinction on higher level; tinoco no rewiring
  # extc_sims_higher_mean_tin_norew <- list("lower" = map(1:3,
  #                      ~ list_mean(extc_sims_higher_tin_norew,
  #                                  y = .x,
  #                                  original = T)),
  #                      "higher" = map(1:3,
  #                       ~ list_mean(extc_sims_higher_tin_norew,
  #                                   y = .x,
  #                                   lower = F,
  #                                   original = T)))
  # 
### calculate percentages of remaining sp ----
  # initial extinction on lower level original web
  sp_remain_lower_org <- map(seq(n_webs), function(x) {
    list("lower" = map(1:3, ~ per_surv(pluck(extc_sims_lower_mean_org, x),
                                       y = .x,
                                       original = T)),
         "higher" = map(1:3, ~ per_surv(pluck(extc_sims_lower_mean_org, x),
                                        y = .x,
                                        lower = F,
                                        original = T)))})
  
  # add dropped names
  sp_remain_lower_org <- modify_depth(sp_remain_lower_org, 2, 
                                  ~ set_names(.x, nm = c("abund",
                                                         "trait",
                                                         "phylo")))
  
  # initial extinction on higher level original web
  sp_remain_higher_org <- map(seq(n_webs), function(x) {
    list("lower" = map(1:3, ~ per_surv(pluck(extc_sims_higher_mean_org, x),
                                       y = .x,
                                       original = T)),
         "higher" = map(1:3, ~ per_surv(pluck(extc_sims_higher_mean_org, x),
                                        y = .x,
                                        lower = F,
                                        original = T)))})
  
  # add dropped names
  sp_remain_higher_org <- modify_depth(sp_remain_higher_org, 2, 
                                      ~ set_names(.x, nm = c("abund",
                                                             "trait",
                                                             "phylo")))
  
  # initial extinction on lower level; org, no rewiring
  sp_remain_lower_org_norew <- map(seq(n_webs), function(x) {
    list("lower" = list_divide(pluck(extc_sims_lower_org_mean_norew, x, 1)),
         "higher" = list_divide(pluck(extc_sims_lower_org_mean_norew, x, 2)))})
  
  # initial extinction on higher level; org, no rewiring
  sp_remain_higher_org_norew <- map(seq(n_webs), function(x) {
    list("lower" = list_divide(pluck(extc_sims_higher_org_mean_norew, x, 1)),
         "higher" = list_divide(pluck(extc_sims_higher_org_mean_norew, x, 2)))})
  
 
  # # initial extinction on lower level, abundance driven extinction
  # sp_remain_lower_abund <- map(seq(n_webs), function(x) {
  #   list("lower" = map(1:3, ~ per_surv(pluck(extc_sims_lower_abund_mean, x),
  #                                      y = .x)),
  #        "higher" = map(1:3, ~ per_surv(pluck(extc_sims_lower_abund_mean, x),
  #                                       y = .x,
  #                                       lower = F)))})
  # 
  # # add dropped names
  # sp_remain_lower_abund <- modify_depth(sp_remain_lower_abund, 2, 
  #                              ~ set_names(.x, nm = c("abund",
  #                                                     "trait",
  #                                                     "phylo")))
  # 
  # # initial extinction on higher level, abundance driven extinction
  # sp_remain_higher_abund <- map(seq(n_webs), function(x) {
  #   list("lower" = map(1:3, ~ per_surv(pluck(extc_sims_higher_abund_mean, x),
  #                                      y = .x)),
  #        "higher" = map(1:3, ~ per_surv(pluck(extc_sims_higher_abund_mean, x),
  #                                       y = .x,
  #                                       lower = F)))})
  # 
  # # add dropped names
  # sp_remain_higher_abund <- modify_depth(sp_remain_higher_abund, 2, 
  #                                       ~ set_names(.x, nm = c("abund",
  #                                                              "trait",
  #                                                              "phylo")))
  # 
  # # initial extinction on lower level, abundance driven extinction, no rewiring
  # sp_remain_lower_abund_norew <- map(seq(n_webs), function(x) {
  #   list("lower" = map(1:3,
  #                      ~ per_surv(pluck(extc_sims_lower_abund_mean_norew, x),
  #                                 y = .x,
  #                                 original = T)),
  #        "higher" = map(1:3,
  #                       ~ per_surv(pluck(extc_sims_lower_abund_mean_norew, x),
  #                                  y = .x,
  #                                  lower = F,
  #                                  original = T)))})
  # 
  # # add dropped names
  # sp_remain_lower_abund_norew <- modify_depth(sp_remain_lower_abund_norew, 2, 
  #                                       ~ set_names(.x, nm = c("Atl",
  #                                                              "aTl",
  #                                                              "atL")))
  # 
  # # initial extinction on higher level, abundance driven extinction, no rewiring
  # sp_remain_higher_abund_norew <- map(seq(n_webs), function(x) {
  #   list("lower" = map(1:3,
  #                      ~ per_surv(pluck(extc_sims_higher_abund_mean_norew, x),
  #                                 y = .x,
  #                                 original = T)),
  #        "higher" = map(1:3,
  #                       ~ per_surv(pluck(extc_sims_higher_abund_mean_norew, x),
  #                                  y = .x,
  #                                  lower = F,
  #                                  original = T)))})
  # 
  # # add dropped names
  # sp_remain_higher_abund_norew <- modify_depth(sp_remain_higher_abund_norew, 2, 
  #                                        ~ set_names(.x, nm = c("Atl",
  #                                                               "aTl",
  #                                                               "atL")))
  # 
  # # initial extinction on both levels
  # sp_remain_both <- map(seq(n_webs), function(x) {
  #   list("lower" = map(1:3, ~ per_surv(pluck(extc_sims_both_mean, x),
  #                                      y = .x)),
  #        "higher" = map(1:3, ~ per_surv(pluck(extc_sims_both_mean, x),
  #                                       y = .x,
  #                                       lower = F)))})
  # 
  # # add dropped names
  # sp_remain_both <- modify_depth(sp_remain_both, 2, 
  #                                      ~ set_names(.x, nm = c("abund",
  #                                                             "trait",
  #                                                             "phylo")))
  # 
  # # initial extinction on both levels, no rewiring
  # sp_remain_both_norew <- map(seq(n_webs), function(x) {
  #   list("lower" = map(1:3, ~ per_surv(pluck(extc_sims_both_mean_norew, x),
  #                                      y = .x,
  #                                      original = T)),
  #        "higher" = map(1:3, ~ per_surv(pluck(extc_sims_both_mean_norew, x),
  #                                       y = .x,
  #                                       lower = F,
  #                                       original = T)))})
  # 
  # # add dropped names
  # sp_remain_both_norew <- modify_depth(sp_remain_both_norew, 2, 
  #                                      ~ set_names(.x, nm = c("Atl",
  #                                                             "aTl",
  #                                                             "atL")))
  # 
  # initial extinction on lower level
  sp_remain_lower <- map(seq(n_webs), function(x) {
    list("lower" = map(1:3, ~ per_surv(pluck(extc_sims_lower_mean, x),
                                       y = .x)),
         "higher" = map(1:3, ~ per_surv(pluck(extc_sims_lower_mean, x),
                                        y = .x,
                                        lower = F)))})
  
  # add dropped names
  sp_remain_lower <- modify_depth(sp_remain_lower, 2, 
                                  ~ set_names(.x, nm = c("abund",
                                                         "trait",
                                                         "phylo")))
  
  # initial extinction on higher level
  sp_remain_higher <- map(seq(n_webs), function(x) {
    list("lower" = map(1:3, ~ per_surv(pluck(extc_sims_higher_mean, x),
                                       y = .x)),
         "higher" = map(1:3, ~ per_surv(pluck(extc_sims_higher_mean, x),
                                        y = .x,
                                        lower = F)))})
  
  # add dropped names
  sp_remain_higher <- modify_depth(sp_remain_higher, 2, 
                                  ~ set_names(.x, nm = c("abund",
                                                         "trait",
                                                         "phylo")))
  
  # initial extinction on lower level, no rewiring
  sp_remain_lower_norew <- map(seq(n_webs), function(x) {
    list("lower" = map(1:3, ~ per_surv(pluck(extc_sims_lower_mean_norew, x),
                                       y = .x,
                                       original = T)),
         "higher" = map(1:3, ~ per_surv(pluck(extc_sims_lower_mean_norew, x),
                                        y = .x,
                                        lower = F,
                                        original = T)))})
  
  # add dropped names
  sp_remain_lower_norew <- modify_depth(sp_remain_lower_norew, 2, 
                                         ~ set_names(.x, nm = c("Atl",
                                                                "aTl",
                                                                "atL")))
  
  # initial extinction on higher level, no rewiring
  sp_remain_higher_norew <- map(seq(n_webs), function(x) {
    list("lower" = map(1:3, ~ per_surv(pluck(extc_sims_higher_mean_norew, x),
                                       y = .x,
                                       original = T)),
         "higher" = map(1:3, ~ per_surv(pluck(extc_sims_higher_mean_norew, x),
                                        y = .x,
                                        lower = F,
                                        original = T)))})
  
  # add dropped names
  sp_remain_higher_norew <- modify_depth(sp_remain_higher_norew, 2, 
                                  ~ set_names(.x, nm = c("Atl",
                                                         "aTl",
                                                         "atL")))
  
  # # initial extinction on lower level; tinoco
  # sp_remain_lower_tin <- list("lower" = map(1:3, ~ per_surv(extc_sims_lower_mean_tin,
  #                                      y = .x,
  #                                      original = T)),
  #                             "higher" = map(1:3, ~ per_surv(extc_sims_lower_mean_tin,
  #                                       y = .x,
  #                                       lower = F,
  #                                       original = T)))
  # 
  # # add dropped names
  # sp_remain_lower_tin <- modify_depth(sp_remain_lower_tin, 1,
  #                                     ~ set_names(.x, nm = c("abund",
  #                                                            "trait",
  #                                                            "phylo")))
  # 
  # # initial extinction on higher level; tinoco
  # sp_remain_higher_tin <- list("lower" = map(1:3, ~ per_surv(extc_sims_higher_mean_tin,
  #                                      y = .x,
  #                                      original = T)),
  #                              "higher" = map(1:3, ~ per_surv(extc_sims_higher_mean_tin,
  #                                       y = .x,
  #                                       lower = F,
  #                                       original = T)))
  # 
  # # add dropped names
  # sp_remain_higher_tin <- modify_depth(sp_remain_higher_tin, 1,
  #                                     ~ set_names(.x, nm = c("abund",
  #                                                            "trait",
  #                                                            "phylo")))
  # 
  # # initial extinction on lower level; tinoco no rewiring
  # sp_remain_lower_tin_norew <- list("lower" = map(1:3, ~ per_surv(extc_sims_lower_mean_tin_norew,
  #                                      y = .x,
  #                                      original = T)),
  #                                   "higher" = map(1:3, ~ per_surv(extc_sims_lower_mean_tin_norew,
  #                                       y = .x,
  #                                       lower = F,
  #                                       original = T)))
  # 
  # # add dropped names
  # sp_remain_lower_tin_norew <- modify_depth(sp_remain_lower_tin_norew, 1,
  #                                     ~ set_names(.x, nm = c("abund",
  #                                                            "trait",
  #                                                            "phylo")))
  # 
  # # initial extinction on higher level; tinoco no rewiring
  # sp_remain_higher_tin_norew <- list("lower" = map(1:3, ~ per_surv(extc_sims_higher_mean_tin_norew,
  #                                      y = .x,
  #                                      original = T)),
  #                                    "higher" = map(1:3, ~ per_surv(extc_sims_higher_mean_tin_norew,
  #                                       y = .x,
  #                                       lower = F,
  #                                       original = T)))
  # 
  # # add dropped names
  # sp_remain_higher_tin_norew <- modify_depth(sp_remain_higher_tin_norew, 1,
  #                                      ~ set_names(.x, nm = c("abund",
  #                                                             "trait",
  #                                                             "phylo")))
  # 
### mean over all webs ----
  # initial extincion on lower level original web
  sp_remain_lower_web_mean_org <- list(
    "lower" =
      map(.x = 1:3, function(x) {
        map(seq(n_webs), ~ pluck(sp_remain_lower_org, .x, "lower", x))  %>%
          as.data.table() %>%
          match_lengths() %>%
          rowMeans() %>% 
          replace(., is.nan(.), 0)}),
    "higher" = map(.x = 1:3, function(x){
      map(seq(n_webs), ~ pluck(sp_remain_lower_org, .x, "higher", x))  %>%
        as.data.table() %>%
        match_lengths() %>%
        rowMeans() %>% 
        replace(., is.nan(.), 0)}))
  
  # add dropped names
  sp_remain_lower_web_mean_org <- modify_depth(
    sp_remain_lower_web_mean_org, 1, 
    ~ set_names(.x, nm = c("abund", "trait", "phylo")))
  
  # initial extincion on higher level original web
  sp_remain_higher_web_mean_org <- list(
    "lower" =
      map(.x = 1:3, function(x) {
        map(seq(n_webs), ~ pluck(sp_remain_higher_org, .x, "lower", x))  %>%
          as.data.table() %>%
          match_lengths() %>%
          rowMeans() %>% 
          replace(., is.nan(.), 0)}),
    "higher" = map(.x = 1:3, function(x){
      map(seq(n_webs), ~ pluck(sp_remain_higher_org, .x, "higher", x))  %>%
        as.data.table() %>%
        match_lengths() %>%
        rowMeans() %>% 
        replace(., is.nan(.), 0)}))
  
  # add dropped names
  sp_remain_higher_web_mean_org <- modify_depth(
    sp_remain_higher_web_mean_org, 1, 
    ~ set_names(.x, nm = c("abund", "trait", "phylo")))
  
  # initial extincion on higher level original web no rewiring
  sp_remain_lower_web_mean_org_norew <- list(
    "lower" = map(seq(n_webs), ~ pluck(sp_remain_lower_org_norew, .x, "lower"))  %>%
          as.data.table() %>%
          match_lengths() %>%
          rowMeans() %>% 
      replace(., is.nan(.), 0),
    "higher" = map(seq(n_webs), ~ pluck(sp_remain_lower_org_norew, .x, "higher"))  %>%
        as.data.table() %>%
        match_lengths() %>%
        rowMeans() %>% 
      replace(., is.nan(.), 0))
  
  # initial extincion on higher level original web no rewiring
  sp_remain_higher_web_mean_org_norew <- list(
    "lower" = map(seq(n_webs), ~ pluck(sp_remain_higher_org_norew, .x, "lower"))  %>%
          as.data.table() %>%
          match_lengths() %>%
          rowMeans() %>% 
      replace(., is.nan(.), 0),
    "higher" = map(seq(n_webs), ~ pluck(sp_remain_higher_org_norew, .x, "higher"))  %>%
        as.data.table() %>%
        match_lengths() %>%
        rowMeans() %>% 
      replace(., is.nan(.), 0))
  
  # # initial extincion on lower level, abundance driven extinction
  # sp_remain_lower_web_mean_abund <- list(
  #   "lower" =
  #     map(.x = 1:3, function(x) map(.x = 1:3, function(y) {
  #       map(seq(n_webs), ~ pluck(sp_remain_lower_abund, .x, "lower", x, y))  %>%
  #         as.data.table() %>%
  #         match_lengths() %>%
  #         rowMeans()})),
  #   "higher" = map(.x = 1:3, function(x) map(.x = 1:3, function(y) {
  #     map(seq(n_webs), ~ pluck(sp_remain_lower_abund, .x, "higher", x, y))  %>%
  #       as.data.table() %>%
  #       match_lengths() %>%
  #       rowMeans()})))
  # 
  # # add dropped names
  # sp_remain_lower_web_mean_abund <- modify_depth(
  #   sp_remain_lower_web_mean_abund, 1, 
  #   ~ set_names(.x, nm = c("abund", "trait", "phylo")))
  # 
  # # initial extincion on higher level, abundance driven extinction
  # sp_remain_higher_web_mean_abund <- list(
  #   "lower" =
  #     map(.x = 1:3, function(x) map(.x = 1:3, function(y) {
  #       map(seq(n_webs), ~ pluck(sp_remain_higher_abund, .x, "lower", x, y))  %>%
  #         as.data.table() %>%
  #         match_lengths() %>%
  #         rowMeans()})),
  #   "higher" = map(.x = 1:3, function(x) map(.x = 1:3, function(y) {
  #     map(seq(n_webs), ~ pluck(sp_remain_higher_abund, .x, "higher", x, y))  %>%
  #       as.data.table() %>%
  #       match_lengths() %>%
  #       rowMeans()})))
  # 
  # # add dropped names
  # sp_remain_higher_web_mean_abund <- modify_depth(
  #   sp_remain_higher_web_mean_abund, 1, 
  #   ~ set_names(.x, nm = c("abund", "trait", "phylo")))
  # 
  # # initial extincion on lower level, abundance driven extinction, no rewiring
  # sp_remain_lower_web_mean_abund_norew <- list(
  #   "lower" =
  #     map(.x = 1:3, function(x) {
  #       map(seq(n_webs), ~ pluck(sp_remain_lower_abund_norew, .x, "lower", x))  %>%
  #         as.data.table() %>%
  #         match_lengths() %>%
  #         rowMeans()}),
  #   "higher" = map(.x = 1:3, function(x) {
  #     map(seq(n_webs), ~ pluck(sp_remain_lower_abund_norew, .x, "higher", x))  %>%
  #       as.data.table() %>%
  #       match_lengths() %>%
  #       rowMeans()}))
  # 
  # # add dropped names
  # sp_remain_lower_web_mean_abund_norew <- modify_depth(
  #   sp_remain_lower_web_mean_abund_norew, 1, 
  #   ~ set_names(.x, nm = c("abund", "trait", "phylo")))
  # 
  # # initial extincion on higher level, abundance driven extinction, no rewiring
  # sp_remain_higher_web_mean_abund_norew <- list(
  #   "lower" =
  #     map(.x = 1:3, function(x) {
  #       map(seq(n_webs), ~ pluck(sp_remain_higher_abund_norew, .x, "lower", x))  %>%
  #         as.data.table() %>%
  #         match_lengths() %>%
  #         rowMeans()}),
  #   "higher" = map(.x = 1:3, function(x) {
  #     map(seq(n_webs), ~ pluck(sp_remain_higher_abund_norew, .x, "higher", x))  %>%
  #       as.data.table() %>%
  #       match_lengths() %>%
  #       rowMeans()}))
  # 
  # # add dropped names
  # sp_remain_higher_web_mean_abund_norew <- modify_depth(
  #   sp_remain_higher_web_mean_abund_norew, 1, 
  #   ~ set_names(.x, nm = c("abund", "trait", "phylo")))
  # 
  # # initial extincion on both levels
  # sp_remain_both_web_mean <- list(
  #   "lower" =
  #     map(.x = 1:3, function(x) map(.x = 1:3, function(y) {
  #       map(seq(n_webs), ~ pluck(sp_remain_both, .x, "lower", x, y))  %>%
  #         as.data.table() %>%
  #         match_lengths() %>%
  #         rowMeans()})),
  #   "higher" = map(.x = 1:3, function(x) map(.x = 1:3, function(y) {
  #     map(seq(n_webs), ~ pluck(sp_remain_both, .x, "higher", x, y))  %>%
  #       as.data.table() %>%
  #       match_lengths() %>%
  #       rowMeans()})))
  # 
  # # add dropped names
  # sp_remain_both_web_mean <- modify_depth(
  #   sp_remain_both_web_mean, 1, 
  #   ~ set_names(.x, nm = c("abund", "trait", "phylo")))
  # 
  # # initial extincion on both levels, no rewiring
  # sp_remain_both_web_mean_norew <- list(
  #   "lower" =
  #     map(.x = 1:3, function(x) {
  #       map(seq(n_webs), ~ pluck(sp_remain_both_norew, .x, "lower", x))  %>%
  #         as.data.table() %>%
  #         match_lengths() %>%
  #         rowMeans()}),
  #   "higher" = map(.x = 1:3, function(x) {
  #     map(seq(n_webs), ~ pluck(sp_remain_both_norew, .x, "higher", x))  %>%
  #       as.data.table() %>%
  #       match_lengths() %>%
  #       rowMeans()}))
  # 
  # # add dropped names
  # sp_remain_both_web_mean_norew <- modify_depth(
  #   sp_remain_both_web_mean_norew, 1, 
  #   ~ set_names(.x, nm = c("abund", "trait", "phylo")))
  # 
  # initial extincion on lower level
  sp_remain_lower_web_mean <- list(
    "lower" =
      map(.x = 1:3, function(x) map(.x = 1:3, function(y) {
        map(seq(n_webs), ~ pluck(sp_remain_lower, .x, "lower", x, y))  %>%
          as.data.table() %>%
          match_lengths() %>%
          rowMeans(., na.rm = T) %>% 
          replace(., is.nan(.), 0)})),
    "higher" = map(.x = 1:3, function(x) map(.x = 1:3, function(y) {
      map(seq(n_webs), ~ pluck(sp_remain_lower, .x, "higher", x, y))  %>%
        as.data.table() %>%
        match_lengths() %>%
        rowMeans(., na.rm = T) %>% 
        replace(., is.nan(.), 0)})))
  
  # add dropped names
  sp_remain_lower_web_mean <- modify_depth(
    sp_remain_lower_web_mean, 1, 
    ~ set_names(.x, nm = c("abund", "trait", "phylo")))
  
  sp_remain_lower_web_mean <- modify_depth(
    sp_remain_lower_web_mean, 2, 
    ~ set_names(.x, nm = c("Atl", "aTl", "atL")))
  
  # initial extincion on higher level
  sp_remain_higher_web_mean <- list(
    "lower" =
      map(.x = 1:3, function(x) map(.x = 1:3, function(y) {
        map(seq(n_webs), ~ pluck(sp_remain_higher, .x, "lower", x, y))  %>%
          as.data.table() %>%
          match_lengths() %>%
          rowMeans() %>% 
          replace(., is.nan(.), 0)})),
    "higher" = map(.x = 1:3, function(x) map(.x = 1:3, function(y) {
      map(seq(n_webs), ~ pluck(sp_remain_higher, .x, "higher", x, y))  %>%
        as.data.table() %>%
        match_lengths() %>%
        rowMeans() %>% 
        replace(., is.nan(.), 0)})))
  
  # add dropped names
  sp_remain_higher_web_mean <- modify_depth(
    sp_remain_higher_web_mean, 1, 
    ~ set_names(.x, nm = c("abund", "trait", "phylo")))
  
  sp_remain_higher_web_mean <- modify_depth(
    sp_remain_higher_web_mean, 2, 
    ~ set_names(.x, nm = c("Atl", "aTl", "atL")))
  
  # initial extincion on lower level, no rewiring
  sp_remain_lower_web_mean_norew <-list(
    "lower" =
      map(.x = 1:3, function(x) {
        map(seq(n_webs), ~ pluck(sp_remain_lower_norew, .x, "lower", x))  %>%
          as.data.table() %>%
          match_lengths() %>%
          rowMeans() %>% 
          replace(., is.nan(.), 0)}),
    "higher" = map(.x = 1:3, function(x){
      map(seq(n_webs), ~ pluck(sp_remain_lower_norew, .x, "higher", x))  %>%
        as.data.table() %>%
        match_lengths() %>%
        rowMeans() %>% 
        replace(., is.nan(.), 0)}))
  
  # add dropped names
  sp_remain_lower_web_mean_norew <- modify_depth(
    sp_remain_lower_web_mean_norew, 1, 
    ~ set_names(.x, nm = c("Atl", "aTl", "atL")))
  
  
  # initial extincion on higher level, no rewiring
  sp_remain_higher_web_mean_norew <- list(
    "lower" =
      map(.x = 1:3, function(x) {
        map(seq(n_webs), ~ pluck(sp_remain_higher_norew, .x, "lower", x))  %>%
          as.data.table() %>%
          match_lengths() %>%
          rowMeans() %>% 
          replace(., is.nan(.), 0)}),
    "higher" = map(.x = 1:3, function(x){
      map(seq(n_webs), ~ pluck(sp_remain_higher_norew, .x, "higher", x))  %>%
        as.data.table() %>%
        match_lengths() %>%
        rowMeans() %>% 
        replace(., is.nan(.), 0)}))
  
  # add dropped names
  sp_remain_higher_web_mean_norew <- modify_depth(
    sp_remain_higher_web_mean_norew, 1, 
    ~ set_names(.x, nm = c("Atl", "aTl", "atL")))

  # # initial extincion on lower level; tinoco
  # sp_remain_lower_web_mean_tin <- list(
  #   "lower" = pluck(sp_remain_lower_tin, "lower"),
  #   "higher" = pluck(sp_remain_lower_tin, "higher"))
  # 
  # # initial extincion on lower level; tinoco
  # sp_remain_higher_web_mean_tin <- list(
  #   "lower" = pluck(sp_remain_higher_tin, "lower"),
  #   "higher" = pluck(sp_remain_higher_tin, "higher"))
  # 
  # # initial extincion on lower level; tinoco no rewiring
  # sp_remain_lower_web_mean_tin_norew <- list(
  #   "lower" = pluck(sp_remain_lower_tin_norew, "lower"),
  #   "higher" = pluck(sp_remain_lower_tin_norew, "higher"))
  # 
  # # initial extincion on lower level; tinoco no rewiring
  # sp_remain_higher_web_mean_tin_norew <- list(
  #   "lower" = pluck(sp_remain_higher_tin_norew, "lower"),
  #   "higher" = pluck(sp_remain_higher_tin_norew, "higher"))
  # 
### data wrangling for visualization ----
### get ci of every web
  # initial extinction on lower level original web
  ci_lower_org <- get_ci(x = sp_remain_lower_org,
                         means = sp_remain_lower_web_mean_org,
                         org = T)
  
  # initial extinction on higher level original web
  ci_higher_org <- get_ci(x = sp_remain_higher_org,
                          means = sp_remain_higher_web_mean_org,
                          org = T)
  
  ### add org norew
  # initial extinction on lower level original web no rewiring
  ci_lower_org_norew <- get_ci(x = sp_remain_lower_org_norew,
                         means = sp_remain_lower_web_mean_org_norew,
                         org = T, norew = T)
  
  # initial extinction on higher level original web no rewiring
  ci_higher_org_norew <- get_ci(x = sp_remain_higher_org_norew,
                          means = sp_remain_higher_web_mean_org_norew,
                          org = T, norew = T)
  
  
  # # initial extinction on lower level, abundance driven extinction
  # ci_lower_abund <- get_ci(x = sp_remain_lower_abund,
  #                          means = sp_remain_lower_web_mean_abund)
  # 
  # # initial extinction on lower level, abundance driven extinction
  # ci_higher_abund <- get_ci(x = sp_remain_higher_abund,
  #                           means = sp_remain_higher_web_mean_abund)
  # 
  # 
  # # initial extinction on higher level, abundance driven extinction, no rewiring
  # ci_higher_abund_norew <- get_ci(x = sp_remain_higher_abund_norew,
  #                                 means = sp_remain_higher_web_mean_abund_norew,
  #                                 norew = T)
  # 
  # # initial extinction on lower level, abundance driven extinction, no rewiring
  # ci_lower_abund_norew <- get_ci(x = sp_remain_lower_abund_norew,
  #                                means = sp_remain_lower_web_mean_abund_norew,
  #                                norew = T)
  # 
  # # initial extinction on both levels
  # ci_both <- get_ci(x = sp_remain_both,
  #                   means = sp_remain_both_web_mean)
  # 
  # # initial extinction on both levels, no rewiring
  # ci_both_norew <- get_ci(x = sp_remain_both_norew,
  #                         means = sp_remain_both_web_mean_norew,
  #                         norew = T)
  # 
  # initial extinction on lower level
  ci_lower <- get_ci(x = sp_remain_lower,
                     means = sp_remain_lower_web_mean)
  
  # initial extinction on higher level
  ci_higher <- get_ci(x = sp_remain_higher,
                      means = sp_remain_higher_web_mean)
  
  # initial extinction on lower level no rewiring
  ci_lower_norew <- get_ci(x = sp_remain_lower_norew,
                           means = sp_remain_lower_web_mean_norew,
                           norew = T)
  
  # initial extinction on higher level no rewiring
  ci_higher_norew <- get_ci(x = sp_remain_higher_norew,
                            means = sp_remain_higher_web_mean_norew,
                            norew = T)
  
  # No CIs for tin since there's only 1 web per rewiring scenario
  # # initial extinction on lower level; tinoco
  # ci_lower_tin <- get_ci(x = sp_remain_lower_tin,
  #                        means = sp_remain_lower_web_mean_tin,
  #                        org = T)
  # 
  # # initial extinction on higher level; tinoco
  # ci_higher_tin <- get_ci(x = sp_remain_higher_tin,
  #                        means = sp_remain_higher_web_mean_tin,
  #                        org = T)
  
### reshape sp_remain lists to dfs
  # initial extinction on lower trophic level, original web
  sp_remain_lower_web_mean_df_org <- list_to_df(sp_remain_lower_web_mean_org,
                                                org = T)
  
  # initial extinction on higher trophic level, original web
  sp_remain_higher_web_mean_df_org <- list_to_df(sp_remain_higher_web_mean_org,
                                                 org = T)
  
  ### add org norew
  # initial extinction on lower trophic level, original web no rewiring
  sp_remain_lower_web_mean_df_org_norew <- list_to_df(sp_remain_lower_web_mean_org_norew,
                                                org = T, norew = T)
  
  # initial extinction on higher trophic level, original web no rewiring
  sp_remain_higher_web_mean_df_org_norew <- list_to_df(sp_remain_higher_web_mean_org_norew,
                                                 org = T, norew = T)
  
  # # initial extinction on lower level, abundance driven extinction
  # sp_remain_lower_web_mean_df_abund <- list_to_df(sp_remain_lower_web_mean_abund)
  # 
  # # initial extinction on higher level, abundance driven extinction
  # sp_remain_higher_web_mean_df_abund <- list_to_df(sp_remain_higher_web_mean_abund)
  # 
  # # initial extinction on lower level, abundance driven extinction, no rewiring
  # sp_remain_lower_web_mean_df_abund_norew <- list_to_df(sp_remain_lower_web_mean_abund_norew,
  #                                                       norew = T)
  # 
  # # initial extinction on higher level, abundance driven extinction, no rewiring
  # sp_remain_higher_web_mean_df_abund_norew <- list_to_df(sp_remain_higher_web_mean_abund_norew,
  #                                                        norew = T)
  # 
  # # initial extinction on both levels
  # sp_remain_both_web_mean_df <- list_to_df(sp_remain_both_web_mean)
  # 
  # # initial extinction on both levels, no rewiring
  # sp_remain_both_web_mean_df_norew <- list_to_df(sp_remain_both_web_mean_norew, norew = T)
  
  # initial extinction on lower trophic level
  sp_remain_lower_web_mean_df <- list_to_df(sp_remain_lower_web_mean)
  
  # initial extinction on higher trophic level
  sp_remain_higher_web_mean_df <- list_to_df(sp_remain_higher_web_mean)
  
  # initial extinction on lower trophic level, no rewiring
  sp_remain_lower_web_mean_df_norew <- list_to_df(sp_remain_lower_web_mean_norew,
                                                  norew = T)
  
  # initial extinction on higher trophic level, no rewiring
  sp_remain_higher_web_mean_df_norew <- list_to_df(sp_remain_higher_web_mean_norew,
                                                   norew = T)
  
  # # initial extinction on lower trophic level, tinoco
  # sp_remain_lower_web_mean_tin_df <- list_to_df(sp_remain_lower_web_mean_tin,
  #                                               org = T)
  # 
  # # initial extinction on higher trophic level, tinoco
  # sp_remain_higher_web_mean_tin_df <- list_to_df(sp_remain_higher_web_mean_tin,
  #                                                org = T)
  # 
  # # initial extinction on lower trophic level, tinoco no rewiring
  # sp_remain_lower_web_mean_tin_df_norew <- list_to_df(sp_remain_lower_web_mean_tin_norew,
  #                                                 norew = T, org = T)
  # 
  # # initial extinction on higher trophic level, tinoco no rewiring
  # sp_remain_higher_web_mean_tin_df_norew <- list_to_df(sp_remain_higher_web_mean_tin_norew,
  #                                                  norew = T, org = T)
  
### reshape ci lists to df
  ci_lower_df_org <- list_to_df(ci_lower_org, ci = T, org = T)
  ci_higher_df_org <- list_to_df(ci_higher_org, ci = T, org = T)
 
  ci_lower_df_org_norew <- list_to_df(ci_lower_org_norew, ci = T, org = T, norew = T)
  ci_higher_df_org_norew <- list_to_df(ci_higher_org_norew, ci = T, org = T, norew = T)
  
  # ci_lower_df_abund <- list_to_df(ci_lower_abund, ci = T)
  # ci_higher_df_abund <- list_to_df(ci_higher_abund, ci = T)
  # 
  # ci_lower_df_abund_norew <- list_to_df(ci_lower_abund_norew, ci = T, norew = T)
  # ci_higher_df_abund_norew <- list_to_df(ci_higher_abund_norew, ci = T, norew = T)
  # 
  # ci_both_df <- list_to_df(ci_both, ci = T)
  # ci_both_df_norew <- list_to_df(ci_both_norew, ci = T, norew = T)
  
  ci_lower_df <- list_to_df(ci_lower, ci = T)
  ci_higher_df <- list_to_df(ci_higher, ci = T)
  
  ci_lower_df_norew <- list_to_df(ci_lower_norew, ci = T, norew = T)
  ci_higher_df_norew <- list_to_df(ci_higher_norew, ci = T, norew = T)
  
  
###################  
  
  extc_sims_lower_AT <- map(seq(n_webs), function(x) {
    run_extc(web = pluck(sims, x),
             participant = "lower",
             method = "random",
             rewiring = T,
             abund.partner.choice = pluck(abunds, x), 
             trait.partner.choice = pluck(traits, x),
             phylo.partner.choice = NULL,
             interactions = pluck(sims, x),
             method.rewiring = c("abund", "trait"),
             n_sims = n_sims,
             multiple.webs = T, 
             coextc.thr = coextc_thr)})  
  
  
  extc_sims_lower_mean_AT <- map(seq(n_webs), function(x) {
    list("lower" = map(1:3,
                       ~ list_mean(pluck(extc_sims_lower_AT, x),
                                   y = .x, original = T)),
         "higher" = map(1:3,
                        ~ list_mean(pluck(extc_sims_lower_AT, x),
                                    y = .x,
                                    lower = F, original = T)))})
  
  sp_remain_lower_AT <- map(seq(n_webs), function(x) {
    list("lower" = map(1:3, ~ per_surv(pluck(extc_sims_lower_mean_AT, x),
                                       y = .x,
                                       original = T)),
         "higher" = map(1:3, ~ per_surv(pluck(extc_sims_lower_mean_AT, x),
                                        y = .x,
                                        lower = F,
                                        original = T)))})
  
  # add dropped names
  sp_remain_lower_AT <- modify_depth(sp_remain_lower_AT, 2, 
                                      ~ set_names(.x, nm = c("Atl",
                                                             "aTl",
                                                             "atL")))
  
  sp_remain_lower_web_mean_AT <- list(
    "lower" =
      map(.x = 1:3, function(x) {
        map(seq(n_webs), ~ pluck(sp_remain_lower_AT, .x, "lower", x))  %>%
          as.data.table() %>%
          match_lengths() %>%
          rowMeans()}),
    "higher" = map(.x = 1:3, function(x){
      map(seq(n_webs), ~ pluck(sp_remain_lower_AT, .x, "higher", x))  %>%
        as.data.table() %>%
        match_lengths() %>%
        rowMeans()}))
  
  # add dropped names
  sp_remain_lower_web_mean_AT <- modify_depth(
    sp_remain_lower_web_mean_AT, 1, 
    ~ set_names(.x, nm = c("Atl", "aTl", "atL")))
  
  
  extc_sims_lower_AP <- map(seq(n_webs), function(x) {
    run_extc(web = pluck(sims, x),
             participant = "lower",
             method = "random",
             rewiring = T,
             abund.partner.choice = pluck(abunds, x), 
             trait.partner.choice = NULL,
             phylo.partner.choice = pluck(phylos, x),
             interactions = pluck(sims, x),
             method.rewiring = c("abund", "phylo"),
             n_sims = n_sims,
             multiple.webs = T, 
             coextc.thr = coextc_thr)})  
  
  
  extc_sims_lower_mean_AP <- map(seq(n_webs), function(x) {
    list("lower" = map(1:3,
                       ~ list_mean(pluck(extc_sims_lower_AP, x),
                                   y = .x, original = T)),
         "higher" = map(1:3,
                        ~ list_mean(pluck(extc_sims_lower_AP, x),
                                    y = .x,
                                    lower = F, original = T)))})
  
  sp_remain_lower_AP <- map(seq(n_webs), function(x) {
    list("lower" = map(1:3, ~ per_surv(pluck(extc_sims_lower_mean_AP, x),
                                       y = .x,
                                       original = T)),
         "higher" = map(1:3, ~ per_surv(pluck(extc_sims_lower_mean_AP, x),
                                        y = .x,
                                        lower = F,
                                        original = T)))})
  
  # add dropped names
  sp_remain_lower_AP <- modify_depth(sp_remain_lower_AP, 2, 
                                     ~ set_names(.x, nm = c("Atl",
                                                            "aTl",
                                                            "atL")))
  
  sp_remain_lower_web_mean_AP <- list(
    "lower" =
      map(.x = 1:3, function(x) {
        map(seq(n_webs), ~ pluck(sp_remain_lower_AP, .x, "lower", x))  %>%
          as.data.table() %>%
          match_lengths() %>%
          rowMeans()}),
    "higher" = map(.x = 1:3, function(x){
      map(seq(n_webs), ~ pluck(sp_remain_lower_AP, .x, "higher", x))  %>%
        as.data.table() %>%
        match_lengths() %>%
        rowMeans()}))
  
  # add dropped names
  sp_remain_lower_web_mean_AP <- modify_depth(
    sp_remain_lower_web_mean_AP, 1, 
    ~ set_names(.x, nm = c("Atl", "aTl", "atL")))
  
  