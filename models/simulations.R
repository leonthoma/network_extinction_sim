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
n_sims <- 10
ctrbs <- c("Atl" = 1, "aTl" = 2, "atL" = 3, "atl" = 4, "ATL" = 5)
rew_names <- c("abund", "trait", "phylo")
coextc_thr <- NULL

# create initial network for community variables
init_sim <- map(seq(n_webs), ~ simulate_tapnet_aug(nlower = 40,
                                                   nhigher = 50,
                                                   ntraits_nopem = 2,
                                                   ntraits_pem = 2,
                                                   abuns = "lognormal",
                                                   Nobs = 1111,
                                                   names = sp_names,
                                                   initial_sim = T))


# Create simulated networks with specified contributions
sims <- list()

# list w/ ctrb_vecs
ctrb_list <- list("Atl" = c("high", "low", "low"),
                  "aTl" = c("low", "high", "low"),
                  "atL" = c("low", "low", "high"),
                  "atl" = c("low", "low", "low"),
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

# copy to preserve all improbable interactions 
sims_all <- sims 


# Delete sp w/ only dead interactions
clean_sims_web <- map(seq(n_webs), function(x) {
  map(ctrbs, ~ del_dead_int(pluck(sims, x, .x, "web")))
})

clean_init_sim_web <- map(seq(n_webs), function(x) {
  del_dead_int(pluck(init_sim, x, "networks" , 1, "web"))
})

# crop webs to size of smallest network w/o dead interactions
sims_web <- equalize_sp(clean_sims_web)

init_sim_web <- equalize_sp(clean_init_sim_web, init = T)

# create list w/ webs and I_mat
sims <- map(seq(n_webs), function(x) {
  map(ctrbs, ~ list("web" = pluck(sims_web, x, .x), 
         "I_mat" = pluck(sims, x, .x, "I_mat")))
         })

# rename list
sims <- modify_depth(sims, 1, 
                    ~ set_names(.x, nm = c("Atl",
                                           "aTl",
                                           "atL",
                                           "atl",
                                           "ATL")))

init_sims <- map(seq(n_webs), function(x) {
  list("web" = pluck(init_sim_web, x), 
                    "I_mat" = pluck(init_sim, x, "networks", 1,  "I_mat"))
})

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
# # abundances
# abunds <- map(seq(n_webs), function(x) {
#   list("low" = pluck(init_sim, x)$networks[[1]]$abuns$low[names(pluck(init_sim, x)$networks[[1]]$abuns$low) %in% rownames(pluck(init_sims, x, "web"))],
#        "high" = pluck(init_sim, x)$networks[[1]]$abuns$high[names(pluck(init_sim, x)$networks[[1]]$abuns$high) %in% colnames(pluck(init_sims, x, "web"))])
#   })
# 
# # abunds_tin <- pluck(tin_init$networks, 1, "abuns")
# 
# # traits 
# traits <- map(seq(n_webs), function(x) {
#   list("low" = pluck(init_sim, x)$networks[[1]]$traits$low[rownames(pluck(init_sim, x)$networks[[1]]$traits$low) %in% rownames(pluck(init_sims, x, "web")),],
#        "high" = pluck(init_sim, x)$networks[[1]]$traits$high[rownames(pluck(init_sim, x)$networks[[1]]$traits$high) %in% colnames(pluck(init_sims, x, "web")),])
# })
# # traits_tin <- pluck(tin_init$networks, 1, "traits")
# 
# # phylogenetic distances
# library(ape)
# 
# phylos <- map(seq(n_webs), ~ list("low" = cophenetic.phylo(pluck(init_sim, .x)$trees$low), 
#                                   "high" = cophenetic.phylo(pluck(init_sim, .x)$trees$high)))
# 
# phylos <- map(seq(n_webs), function(x) {
#   list("low" = pluck(phylos, x)$low[rownames(pluck(phylos, x)$low) %in% rownames(pluck(init_sims, x, "web")),
#                                     colnames(pluck(phylos, x)$low) %in% rownames(pluck(init_sims, x, "web"))],
#        "high" = pluck(phylos, x)$high[rownames(pluck(phylos, x)$high) %in% colnames(pluck(init_sims, x, "web")),
#                                      colnames(pluck(phylos, x)$high) %in% colnames(pluck(init_sims, x, "web"))])
# })

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

## original web ----  
  # initial extinction on lower level; original
  extc_sims_lower_org <- map(seq(n_webs), function(x) {
    map2(.x = list("abund" = abunds,
                   "trait" = traits,
                   "phylo" = phylos),
         .y = c("abund", "trait", "phylo"),
         ~ replicate(n_sims, simplify = F,
                     one.second.extinct.mod.aug(
                       web = pluck(init_sims, x, "web"),
                       participant = "lower",
                       method = "random",
                       rewiring = T,
                       abund.partner.choice = pluck(.x, x),
                       trait.partner.choice = pluck(.x, x),
                       phylo.partner.choice = pluck(.x, x),
                       interactions = pluck(init_sims, x, "I_mat"),
                       method.rewiring = .y, 
                       coextc.thr = coextc_thr)))})
  
  # initial extinction on higher level; original
  extc_sims_higher_org <- map(seq(n_webs), function(x) {
    map2(.x = list("abund" = abunds,
                   "trait" = traits,
                   "phylo" = phylos),
         .y = c("abund", "trait", "phylo"),
         ~ replicate(n_sims, simplify = F,
                     one.second.extinct.mod.aug(
                       web = pluck(init_sims, x, "web"),
                       participant = "higher",
                       method = "random",
                       rewiring = T,
                       abund.partner.choice = pluck(.x, x),
                       trait.partner.choice = pluck(.x, x),
                       phylo.partner.choice = pluck(.x, x),
                       interactions = pluck(init_sims, x, "I_mat"),
                       method.rewiring = .y, 
                       coextc.thr = coextc_thr)))})

  # initial extinction on lower level; original no rewiring
  extc_sims_lower_org_norew <- map(seq(n_webs), function(x) {
    replicate(n_sims, simplify = F,
              one.second.extinct.mod.aug(
                web = pluck(init_sims, x, "web"),
                participant = "lower",
                method = "random",
                rewiring = F,
                abund.partner.choice = NULL,
                trait.partner.choice = NULL,
                phylo.partner.choice = NULL,
                interactions = pluck(init_sims, x, "I_mat"),
                method.rewiring = "NULL", 
                coextc.thr = coextc_thr))})
  
  # initial extinction on lower lever; original no rewiring
  extc_sims_higher_org_norew <- map(seq(n_webs), function(x) {
    replicate(n_sims, simplify = F,
              one.second.extinct.mod.aug(
                web = pluck(init_sims, x, "web"),
                participant = "higher",
                method = "random",
                rewiring = F,
                abund.partner.choice = NULL,
                trait.partner.choice = NULL,
                phylo.partner.choice = NULL,
                interactions = pluck(init_sims, x, "I_mat"),
                method.rewiring = "NULL", 
                coextc.thr = coextc_thr))})
  
  # initial extinction on lower level; original AT rewiring
  extc_sims_lower_org_AT <- map(seq(n_webs), function(x) {
    replicate(n_sims, simplify = F,
              one.second.extinct.mod.aug(web = pluck(init_sims, x, "web"),
                                         participant = "lower",
                                         method = "random",
                                         rewiring = T,
                                         abund.partner.choice = pluck(abunds, x),
                                         trait.partner.choice = pluck(traits, x),
                                         phylo.partner.choice = NULL, 
                                         interactions = pluck(init_sims, x, "I_mat"),
                                         method.rewiring = c("abund", "trait"),
                                         coextc.thr = coextc_thr))
  })
  
  # initial extinction on higher level; original AT rewiring
  extc_sims_higher_org_AT <- map(seq(n_webs), function(x) {
    replicate(n_sims, simplify = F,
              one.second.extinct.mod.aug(web = pluck(init_sims, x, "web"),
                                         participant = "higher",
                                         method = "random",
                                         rewiring = T,
                                         abund.partner.choice = pluck(abunds, x),
                                         trait.partner.choice = pluck(traits, x),
                                         phylo.partner.choice = NULL, 
                                         interactions = pluck(init_sims, x, "I_mat"),
                                         method.rewiring = c("abund", "trait"),
                                         coextc.thr = coextc_thr))
  })
  
  # initial extinction on lower level; original AP rewiring
  extc_sims_lower_org_AP <- map(seq(n_webs), function(x) {
    replicate(n_sims, simplify = F,
              one.second.extinct.mod.aug(web = pluck(init_sims, x, "web"),
                                         participant = "lower",
                                         method = "random",
                                         rewiring = T,
                                         abund.partner.choice = pluck(abunds, x),
                                         trait.partner.choice = NULL,
                                         phylo.partner.choice = pluck(phylos, x), 
                                         interactions = pluck(init_sims, x, "I_mat"),
                                         method.rewiring = c("abund", "phylo"),
                                         coextc.thr = coextc_thr))
  })
  
  # initial extinction on higher level; original AP rewiring
  extc_sims_higher_org_AP <- map(seq(n_webs), function(x) {
    replicate(n_sims, simplify = F,
              one.second.extinct.mod.aug(web = pluck(init_sims, x, "web"),
                                         participant = "higher",
                                         method = "random",
                                         rewiring = T,
                                         abund.partner.choice = pluck(abunds, x),
                                         trait.partner.choice = NULL,
                                         phylo.partner.choice = pluck(phylos, x), 
                                         interactions = pluck(init_sims, x, "I_mat"),
                                         method.rewiring = c("abund", "phylo"),
                                         coextc.thr = coextc_thr))
  })
  
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
  
## sim webs ----
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
  
  # initial extinction on lower level, AT
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
  
  # initial extinction on higher level, AT
  extc_sims_higher_AT <- map(seq(n_webs), function(x) {
    run_extc(web = pluck(sims, x),
             participant = "higher",
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
  
  # initial extinction on lower level, AP
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
  
  # initial extinction on higher level, AP
  extc_sims_higher_AP <- map(seq(n_webs), function(x) {
    run_extc(web = pluck(sims, x),
             participant = "higher",
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
  
### calculate remaining species ----  
## original web ----
  sp_remain_lower_org_sims <- map(seq(n_webs), function(x) {
    list("lower" = map(1:3, function(y) { # loop rewiring
      pluck(extc_sims_lower_org, x, y) %>%
        as.data.table(.) %>% 
        select(., contains(c("no", "n.lower"))) %>%
        replace_duplicate() %>% apply(., 2, list_divide)}),
      "higher" = map(1:3, function(y) { # loop rewiring
        pluck(extc_sims_lower_org, x, y) %>%
          as.data.table(.) %>% 
          select(., contains(c("no", "n.higher"))) %>%
          replace_duplicate()%>% apply(., 2, list_divide)}))
  })
  
  sp_remain_higher_org_sims <- map(seq(n_webs), function(x) {
    list("lower" = map(1:3, function(y) { # loop rewiring
      pluck(extc_sims_higher_org, x, y) %>%
        as.data.table(.) %>% 
        select(., contains(c("no", "n.lower"))) %>%
        replace_duplicate() %>% apply(., 2, list_divide)}),
      "higher" = map(1:3, function(y) { # loop rewiring
        pluck(extc_sims_higher_org, x, y) %>%
          as.data.table(.) %>% 
          select(., contains(c("no", "n.higher"))) %>%
          replace_duplicate()%>% apply(., 2, list_divide)}))
  })
  
  sp_remain_lower_org_norew_sims <- map(seq(n_webs), function(x) {
    list("lower" = pluck(extc_sims_lower_org_norew, x) %>%
           as.data.table(.) %>% 
           select(., contains(c("no", "n.lower"))) %>%
           replace_duplicate() %>% apply(., 2, list_divide),
         "higher" = pluck(extc_sims_lower_org_norew, x) %>%
           as.data.table(.) %>% 
           select(., contains(c("no", "n.higher"))) %>%
           replace_duplicate()%>% apply(., 2, list_divide))
  })
  
  sp_remain_higher_org_norew_sims <- map(seq(n_webs), function(x) {
    list("lower" = pluck(extc_sims_higher_org_norew, x) %>%
           as.data.table(.) %>% 
           select(., contains(c("no", "n.lower"))) %>%
           replace_duplicate() %>% apply(., 2, list_divide),
         "higher" = pluck(extc_sims_higher_org_norew, x) %>%
           as.data.table(.) %>% 
           select(., contains(c("no", "n.higher"))) %>%
           replace_duplicate()%>% apply(., 2, list_divide))
  })
  
  sp_remain_lower_org_AT_sims <- map(seq(n_webs), function(x) {
    list("lower" = pluck(extc_sims_lower_org_AT, x) %>%
           as.data.table(.) %>% 
           select(., contains(c("no", "n.lower"))) %>%
           replace_duplicate() %>% apply(., 2, list_divide),
         "higher" = pluck(extc_sims_lower_org_AT, x) %>%
           as.data.table(.) %>% 
           select(., contains(c("no", "n.higher"))) %>%
           replace_duplicate()%>% apply(., 2, list_divide))
  })
  
  sp_remain_higher_org_AT_sims <- map(seq(n_webs), function(x) {
    list("lower" = pluck(extc_sims_higher_org_AT, x) %>%
           as.data.table(.) %>% 
           select(., contains(c("no", "n.lower"))) %>%
           replace_duplicate() %>% apply(., 2, list_divide),
         "higher" = pluck(extc_sims_higher_org_AT, x) %>%
           as.data.table(.) %>% 
           select(., contains(c("no", "n.higher"))) %>%
           replace_duplicate()%>% apply(., 2, list_divide))
  })
  
  sp_remain_lower_org_AP_sims <- map(seq(n_webs), function(x) {
    list("lower" = pluck(extc_sims_lower_org_AP, x) %>%
           as.data.table(.) %>% 
           select(., contains(c("no", "n.lower"))) %>%
           replace_duplicate() %>% apply(., 2, list_divide),
         "higher" = pluck(extc_sims_lower_org_AP, x) %>%
           as.data.table(.) %>% 
           select(., contains(c("no", "n.higher"))) %>%
           replace_duplicate()%>% apply(., 2, list_divide))
  })
  
  sp_remain_higher_org_AP_sims <- map(seq(n_webs), function(x) {
    list("lower" = pluck(extc_sims_higher_org_AP, x) %>%
           as.data.table(.) %>% 
           select(., contains(c("no", "n.lower"))) %>%
           replace_duplicate() %>% apply(., 2, list_divide),
         "higher" = pluck(extc_sims_higher_org_AP, x) %>%
           as.data.table(.) %>% 
           select(., contains(c("no", "n.higher"))) %>%
           replace_duplicate()%>% apply(., 2, list_divide))
  })
  
## sim webs ----  
  sp_remain_lower_sims <- map(seq(n_webs), function(x) {
    list("lower" = map(1:3, function(y) { # loop rewiring
      map(ctrbs, function(z) { # loop com vars
        pluck(extc_sims_lower, x, y, z) %>%
          as.data.table(.) %>% 
          select(., contains(c("no", "n.lower"))) %>%
          replace_duplicate() %>% apply(., 2, list_divide)})}),
      "higher" = map(1:3, function(y) { # loop rewiring
        map(ctrbs, function(z) { # loop com vars
          pluck(extc_sims_lower, x, y, z) %>%
            as.data.table(.) %>% 
            select(., contains(c("no", "n.higher"))) %>%
            replace_duplicate()%>% apply(., 2, list_divide)})}))
  })
  
  sp_remain_higher_sims <- map(seq(n_webs), function(x) {
    list("lower" = map(1:3, function(y) { # loop rewiring
      map(ctrbs, function(z) { # loop com vars
        pluck(extc_sims_higher, x, y, z) %>%
          as.data.table(.) %>% 
          select(., contains(c("no", "n.lower"))) %>%
          replace_duplicate() %>% apply(., 2, list_divide)})}),
      "higher" = map(1:3, function(y) { # loop rewiring
        map(ctrbs, function(z) { # loop com vars
          pluck(extc_sims_higher, x, y, z) %>%
            as.data.table(.) %>% 
            select(., contains(c("no", "n.higher"))) %>%
            replace_duplicate()%>% apply(., 2, list_divide)})}))
  })
  
  sp_remain_lower_norew_sims <- map(seq(n_webs), function(x) {
    list("lower" = map(ctrbs, function(y) { # loop com vars
      pluck(extc_sims_lower_norew, x, y) %>%
        as.data.table(.) %>% 
        select(., contains(c("no", "n.lower"))) %>%
        replace_duplicate() %>% apply(., 2, list_divide)}),
      "higher" = map(ctrbs, function(y) { # loop com vars
        pluck(extc_sims_lower_norew, x, y) %>%
          as.data.table(.) %>% 
          select(., contains(c("no", "n.higher"))) %>%
          replace_duplicate()%>% apply(., 2, list_divide)}))
  })
  
  sp_remain_higher_norew_sims <- map(seq(n_webs), function(x) {
    list("lower" = map(ctrbs, function(y) { # loop com vars
      pluck(extc_sims_higher_norew, x, y) %>%
        as.data.table(.) %>% 
        select(., contains(c("no", "n.lower"))) %>%
        replace_duplicate() %>% apply(., 2, list_divide)}),
      "higher" = map(ctrbs, function(y) { # loop com vars
        pluck(extc_sims_higher_norew, x, y) %>%
          as.data.table(.) %>% 
          select(., contains(c("no", "n.higher"))) %>%
          replace_duplicate()%>% apply(., 2, list_divide)}))
  })
  
  sp_remain_lower_AT_sims <- map(seq(n_webs), function(x) {
    list("lower" = map(ctrbs, function(y) { # loop com vars
      pluck(extc_sims_lower_AT, x, y) %>%
        as.data.table(.) %>% 
        select(., contains(c("no", "n.lower"))) %>%
        replace_duplicate() %>% apply(., 2, list_divide)}),
      "higher" = map(ctrbs, function(y) { # loop com vars
        pluck(extc_sims_lower_AT, x, y) %>%
          as.data.table(.) %>% 
          select(., contains(c("no", "n.higher"))) %>%
          replace_duplicate()%>% apply(., 2, list_divide)}))
  })
  
  sp_remain_higher_AT_sims <- map(seq(n_webs), function(x) {
    list("lower" = map(ctrbs, function(y) { # loop com vars
      pluck(extc_sims_higher_AT, x, y) %>%
        as.data.table(.) %>% 
        select(., contains(c("no", "n.lower"))) %>%
        replace_duplicate() %>% apply(., 2, list_divide)}),
      "higher" = map(ctrbs, function(y) { # loop com vars
        pluck(extc_sims_higher_AT, x, y) %>%
          as.data.table(.) %>% 
          select(., contains(c("no", "n.higher"))) %>%
          replace_duplicate()%>% apply(., 2, list_divide)}))
  })
  
  sp_remain_lower_AP_sims <- map(seq(n_webs), function(x) {
    list("lower" = map(ctrbs, function(y) { # loop com vars
      pluck(extc_sims_lower_AP, x, y) %>%
        as.data.table(.) %>% 
        select(., contains(c("no", "n.lower"))) %>%
        replace_duplicate() %>% apply(., 2, list_divide)}),
      "higher" = map(ctrbs, function(y) { # loop com vars
        pluck(extc_sims_lower_AP, x, y) %>%
          as.data.table(.) %>% 
          select(., contains(c("no", "n.higher"))) %>%
          replace_duplicate()%>% apply(., 2, list_divide)}))
  })
  
  sp_remain_higher_AP_sims <- map(seq(n_webs), function(x) {
    list("lower" = map(ctrbs, function(y) { # loop com vars
      pluck(extc_sims_higher_AP, x, y) %>%
        as.data.table(.) %>% 
        select(., contains(c("no", "n.lower"))) %>%
        replace_duplicate() %>% apply(., 2, list_divide)}),
      "higher" = map(ctrbs, function(y) { # loop com vars
        pluck(extc_sims_higher_AP, x, y) %>%
          as.data.table(.) %>% 
          select(., contains(c("no", "n.higher"))) %>%
          replace_duplicate()%>% apply(., 2, list_divide)}))
  })
  

### create dataframes w/ remaining species ----
## original web ----
  # lower org
  lower_org_sims_df <- list_to_df2(sp_remain_lower_org_sims[[1]], org = T)
  
  lower_org_sims_df_mean <- map(1:3, ~ data.frame("x" = rowMeans(pluck(sp_remain_lower_org_sims,1 ,"lower", .x)),
                                                  "y" = rowMeans(pluck(sp_remain_lower_org_sims, 1, "higher", .x))))
  
  # lower_org_sims_df <- map(seq(n_webs), ~list_to_df2(sp_remain_lower_org_sims[[.x]], org = T))
  # 
  # lower_org_sims_df_mean <- map(seq(n_webs), function(x) {
  #   map(1:3, ~ data.frame("x" = rowMeans(pluck(sp_remain_lower_org_sims, x ,"lower", .x))%>% rev(),
  #                         "y" = rowMeans(pluck(sp_remain_lower_org_sims, x, "higher", .x))))%>% bind_rows(., .id = "id")
  # }) %>% bind_rows(., .id = "web")
  
  # higher org
  higher_org_sims_df <- list_to_df2(sp_remain_higher_org_sims[[1]], org = T)
  
  higher_org_sims_df_mean <- map(1:3, ~ data.frame("x" = rowMeans(pluck(sp_remain_higher_org_sims,1 ,"lower", .x)),
                                                   "y" = rowMeans(pluck(sp_remain_higher_org_sims, 1, "higher", .x))))
  
  # lower org norew
  lower_org_norew_sims_df <- list_to_df2(sp_remain_lower_org_norew_sims[[1]], org = T, norew = T)
  
  lower_org_norew_sims_df_mean <- data.frame("x" = rowMeans(pluck(sp_remain_lower_org_norew_sims, 1, "lower")),
                                             "y" = rowMeans(pluck(sp_remain_lower_org_norew_sims, 1, "higher")))
  
  # higher org norew
  higher_org_norew_sims_df <- list_to_df2(sp_remain_higher_org_norew_sims[[1]], org = T, norew = T)
  
  higher_org_norew_sims_df_mean <- data.frame("x" = rowMeans(pluck(sp_remain_higher_org_norew_sims, 1, "lower")),
                                              "y" = rowMeans(pluck(sp_remain_higher_org_norew_sims, 1, "higher")))
  
  # lower org AT
  lower_org_AT_sims_df <- list_to_df2(sp_remain_lower_org_AT_sims[[1]], org = T, norew = T)
  lower_org_AT_sims_df$id <- "abundtrait"
  
  lower_org_AT_sims_df_mean <- data.frame("x" = rowMeans(pluck(sp_remain_lower_org_AT_sims, 1, "lower")),
                                          "y" = rowMeans(pluck(sp_remain_lower_org_AT_sims, 1, "higher")))
  
  # higher org AT
  higher_org_AT_sims_df <- list_to_df2(sp_remain_higher_org_AT_sims[[1]], org = T, norew = T)
  higher_org_AT_sims_df$id <- "abundtrait"
  
  higher_org_AT_sims_df_mean <- data.frame("x" = rowMeans(pluck(sp_remain_higher_org_AT_sims, 1, "lower")),
                                           "y" = rowMeans(pluck(sp_remain_higher_org_AT_sims, 1, "higher")))
  
  # lower org AP
  lower_org_AP_sims_df <- list_to_df2(sp_remain_lower_org_AP_sims[[1]], org = T, norew = T)
  lower_org_AP_sims_df$id <- "abundphylo"
  
  lower_org_AP_sims_df_mean <- data.frame("x" = rowMeans(pluck(sp_remain_lower_org_AP_sims, 1, "lower")),
                                          "y" = rowMeans(pluck(sp_remain_lower_org_AP_sims, 1, "higher")))
  
  # higher org AP
  higher_org_AP_sims_df <- list_to_df2(sp_remain_higher_org_AP_sims[[1]], org = T, norew = T)
  higher_org_AP_sims_df$id <- "abundphylo"
  
  higher_org_AP_sims_df_mean <- data.frame("x" = rowMeans(pluck(sp_remain_higher_org_AP_sims, 1, "lower")),
                                           "y" = rowMeans(pluck(sp_remain_higher_org_AP_sims, 1, "higher")))
  
## sim webs ----  
  # lower
  lower_sims_df <- list_to_df2(sp_remain_lower_sims[[1]])
  
  lower_sims_df_mean <- map(1:3, function(x) {
    map(ctrbs, function(y) {
      data.frame("x" = rowMeans(pluck(sp_remain_lower_sims,1 ,"lower", x, y)),
                 "y" = rowMeans(pluck(sp_remain_lower_sims, 1, "higher", x, y)))
    })
  })

  # higher
  higher_sims_df <- list_to_df2(sp_remain_higher_sims[[1]])
  
  higher_sims_df_mean <- map(1:3, function(x) {
    map(ctrbs, function(y) {
      data.frame("x" = rowMeans(pluck(sp_remain_higher_sims,1 ,"lower", x, y)),
                 "y" = rowMeans(pluck(sp_remain_higher_sims, 1, "higher", x, y)))
    })
  })
  
  # lower norew
  lower_norew_sims_df <- list_to_df2(sp_remain_lower_norew_sims[[1]], norew = T)
  
  lower_norew_sims_df_mean <- map(ctrbs, ~ data.frame("x" = rowMeans(pluck(sp_remain_lower_norew_sims,1 ,"lower", .x)),
                                                      "y" = rowMeans(pluck(sp_remain_lower_norew_sims, 1, "higher", .x))))
  
  # higher norew
  higher_norew_sims_df <- list_to_df2(sp_remain_higher_norew_sims[[1]], norew = T)
  
  higher_norew_sims_df_mean <- map(ctrbs, ~ data.frame("x" = rowMeans(pluck(sp_remain_higher_norew_sims,1 ,"lower", .x)),
                                                       "y" = rowMeans(pluck(sp_remain_higher_norew_sims, 1, "higher", .x))))
  
  # lower AT
  lower_AT_sims_df <- list_to_df2(sp_remain_lower_AT_sims[[1]], norew = T)
  lower_AT_sims_df$id <- "abundtrait"
  
  lower_AT_sims_df_mean <- map(ctrbs, ~ data.frame("x" = rowMeans(pluck(sp_remain_lower_AT_sims,1 ,"lower", .x)),
                                                   "y" = rowMeans(pluck(sp_remain_lower_AT_sims, 1, "higher", .x))))
  
  # higher AT
  higher_AT_sims_df <- list_to_df2(sp_remain_higher_AT_sims[[1]], norew = T)
  higher_AT_sims_df$id <- "abundtrait"
  
  higher_AT_sims_df_mean <- map(ctrbs, ~ data.frame("x" = rowMeans(pluck(sp_remain_higher_AT_sims,1 ,"lower", .x)),
                                                    "y" = rowMeans(pluck(sp_remain_higher_AT_sims, 1, "higher", .x))))
  
  # lower AP
  lower_AP_sims_df <- list_to_df2(sp_remain_lower_AP_sims[[1]], norew = T)
  lower_AP_sims_df$id <- "abundphylo"
  
  lower_AP_sims_df_mean <- map(ctrbs, ~ data.frame("x" = rowMeans(pluck(sp_remain_lower_AP_sims,1 ,"lower", .x)),
                                                   "y" = rowMeans(pluck(sp_remain_lower_AP_sims, 1, "higher", .x))))
  
  # higher AP
  higher_AP_sims_df <- list_to_df2(sp_remain_higher_AP_sims[[1]], norew = T)
  higher_AP_sims_df$id <- "abundphylo"
  
  higher_AP_sims_df_mean <- map(ctrbs, ~ data.frame("x" = rowMeans(pluck(sp_remain_higher_AP_sims,1 ,"lower", .x)),
                                                    "y" = rowMeans(pluck(sp_remain_higher_AP_sims, 1, "higher", .x))))

### Creating results df for anova, plots, etc.  
## two dimensional shannon entropy
  H2 <- map(seq(n_webs), function(x) {
    map(ctrbs, ~ H2fun(pluck(sims, x, .x, "web"))["H2"])})
  
  H2_org <- map(seq(n_webs), ~ H2fun(pluck(init_sim_web, .x))[["H2"]])
  
## auc values original webs ----
  # lower org
  auc_lower_org <- map(seq(n_webs), function(x) {
    map(1:3, function(y) {
      map(seq(ncol(pluck(sp_remain_lower_org_sims, x, "lower", y))), function(.x) {
        auc(x = pluck(sp_remain_lower_org_sims, x,  "lower", y)[, .x] %>%
              rev(),
            y = pluck(sp_remain_lower_org_sims, x, "higher", y)[, .x])}) %>% 
        unlist %>% 
        data.table("robustness" = .) %>% 
        cbind("rew" = rew_names[y]) %>% 
        cbind("com_vars" = "org") %>% 
        cbind("p_single" = singleton(pluck(init_sim_web, x))) %>% 
        cbind("web_no" = x) %>% 
        cbind("H2" = H2_org[[x]])
    }) %>% bind_rows()
  })
  
  auc_lower_org <- map(seq(n_webs), function(x) {
    bind_cols(pluck(auc_lower_org, x),
              "nrow" = nrow(pluck(init_sim_web, x)) %>% 
                rep(., nrow(pluck(auc_lower_org[[x]]))),
              "ncol" = ncol(pluck(init_sim_web, x)) %>%
                rep(., nrow(pluck(auc_lower_org[[x]]))))
  }) %>% bind_rows()
  
  # higher org
  auc_higher_org <- map(seq(n_webs), function(x) {
    map(1:3, function(y) {
      map(seq(ncol(pluck(sp_remain_higher_org_sims, x, "lower", y))), function(.x) {
        auc(x = pluck(sp_remain_higher_org_sims, x,  "lower", y)[, .x] %>%
              rev(),
            y = pluck(sp_remain_higher_org_sims, x, "higher", y)[, .x])}) %>% 
        unlist %>% 
        data.table("robustness" = .) %>% 
        cbind("rew" = rew_names[y]) %>% 
        cbind("com_vars" = "org") %>% 
        cbind("p_single" = singleton(pluck(init_sim_web, x), lower = F)) %>% 
        cbind("web_no" = x) %>% 
        cbind("H2" = H2_org[[x]])
    }) %>% bind_rows()
  })
  
  auc_higher_org <- map(seq(n_webs), function(x) {
    bind_cols(pluck(auc_higher_org, x),
              "nrow" = nrow(pluck(init_sim_web, x)) %>% 
                rep(., nrow(pluck(auc_higher_org[[x]]))),
              "ncol" = ncol(pluck(init_sim_web, x)) %>%
                rep(., nrow(pluck(auc_higher_org[[x]]))))
  }) %>% bind_rows()
  
  # lower org norew
  auc_lower_org_norew <- map(seq(n_webs), function (x) {
    map(seq(ncol(pluck(sp_remain_lower_org_norew_sims, x, "lower"))), function(y) {
      auc(x = pluck(sp_remain_lower_org_norew_sims, x, "lower")[, y] %>%
            rev(),
          y = pluck(sp_remain_lower_org_norew_sims, x, "higher")[, y])}) %>% unlist %>%
      data.table("robustness" = .) %>% 
      cbind("rew" = "norew") %>%
      cbind("com_vars" = "org") %>% 
      cbind("p_single" = singleton(pluck(init_sim_web, x))) %>% 
      cbind("web_no" = x) %>% 
      cbind("H2" = H2_org[[x]])})
  
  auc_lower_org_norew <- map(seq(n_webs), function(x) {
    bind_cols(pluck(auc_lower_org_norew, x),
              "nrow" = nrow(pluck(init_sim_web, x)),
              "ncol" = ncol(pluck(init_sim_web, x)))
  }) %>% bind_rows()
  
  # higher org norew
  auc_higher_org_norew <- map(seq(n_webs), function (x) {
    map(seq(ncol(pluck(sp_remain_higher_org_norew_sims, x, "lower"))), function(y) {
      auc(x = pluck(sp_remain_higher_org_norew_sims, x, "lower")[, y] %>%
            rev(),
          y = pluck(sp_remain_higher_org_norew_sims, x, "higher")[, y])}) %>% unlist %>%
      data.table("robustness" = .) %>% 
      cbind("rew" = "norew") %>%
      cbind("com_vars" = "org") %>% 
      cbind("p_single" = singleton(pluck(init_sim_web, x), lower = F)) %>% 
      cbind("web_no" = x) %>% 
      cbind("H2" = H2_org[[x]])})
  
  auc_higher_org_norew <- map(seq(n_webs), function(x) {
    bind_cols(pluck(auc_higher_org_norew, x),
              "nrow" = nrow(pluck(init_sim_web, x)),
              "ncol" = ncol(pluck(init_sim_web, x)))
  }) %>% bind_rows()
  
  # lower org AT
  auc_lower_org_AT <- map(seq(n_webs), function (x) {
    map(seq(ncol(pluck(sp_remain_lower_org_AT_sims, x, "lower"))), function(y) {
      auc(x = pluck(sp_remain_lower_org_AT_sims, x, "lower")[, y] %>%
            rev(),
          y = pluck(sp_remain_lower_org_AT_sims, x, "higher")[, y])}) %>% unlist %>%
      data.table("robustness" = .) %>% 
      cbind("rew" = "abundtrait") %>%
      cbind("com_vars" = "org") %>% 
      cbind("p_single" = singleton(pluck(init_sim_web, x))) %>% 
      cbind("web_no" = x) %>% 
      cbind("H2" = H2_org[[x]])})
  
  auc_lower_org_AT <- map(seq(n_webs), function(x) {
    bind_cols(pluck(auc_lower_org_AT, x),
              "nrow" = nrow(pluck(init_sim_web, x)),
              "ncol" = ncol(pluck(init_sim_web, x)))
  }) %>% bind_rows()
  
  # higher org AT
  auc_higher_org_AT <- map(seq(n_webs), function (x) {
    map(seq(ncol(pluck(sp_remain_higher_org_AT_sims, x, "lower"))), function(y) {
      auc(x = pluck(sp_remain_higher_org_AT_sims, x, "lower")[, y] %>%
            rev(),
          y = pluck(sp_remain_higher_org_AT_sims, x, "higher")[, y])}) %>% unlist %>%
      data.table("robustness" = .) %>% 
      cbind("rew" = "abundtrait") %>%
      cbind("com_vars" = "org") %>% 
      cbind("p_single" = singleton(pluck(init_sim_web, x), lower = F)) %>% 
      cbind("web_no" = x) %>% 
      cbind("H2" = H2_org[[x]])})
  
  auc_higher_org_AT <- map(seq(n_webs), function(x) {
    bind_cols(pluck(auc_higher_org_AT, x),
              "nrow" = nrow(pluck(init_sim_web, x)),
              "ncol" = ncol(pluck(init_sim_web, x)))
  }) %>% bind_rows()
  
  # lower org AP
  auc_lower_org_AP <- map(seq(n_webs), function (x) {
    map(seq(ncol(pluck(sp_remain_lower_org_AP_sims, x, "lower"))), function(y) {
      auc(x = pluck(sp_remain_lower_org_AP_sims, x, "lower")[, y] %>%
            rev(),
          y = pluck(sp_remain_lower_org_AP_sims, x, "higher")[, y])}) %>% unlist %>%
      data.table("robustness" = .) %>% 
      cbind("rew" = "abundphylo") %>%
      cbind("com_vars" = "org") %>% 
      cbind("p_single" = singleton(pluck(init_sim_web, x))) %>% 
      cbind("web_no" = x) %>% 
      cbind("H2" = H2_org[[x]])})
  
  auc_lower_org_AP <- map(seq(n_webs), function(x) {
    bind_cols(pluck(auc_lower_org_AP, x),
              "nrow" = nrow(pluck(init_sim_web, x)),
              "ncol" = ncol(pluck(init_sim_web, x)))
  }) %>% bind_rows()
  
  # higher org AP
  auc_higher_org_AP <- map(seq(n_webs), function (x) {
    map(seq(ncol(pluck(sp_remain_higher_org_AP_sims, x, "lower"))), function(y) {
      auc(x = pluck(sp_remain_higher_org_AP_sims, x, "lower")[, y] %>%
            rev(),
          y = pluck(sp_remain_higher_org_AP_sims, x, "higher")[, y])}) %>% unlist %>%
      data.table("robustness" = .) %>% 
      cbind("rew" = "abundphylo") %>%
      cbind("com_vars" = "org") %>% 
      cbind("p_single" = singleton(pluck(init_sim_web, x), lower = F)) %>% 
      cbind("web_no" = x) %>% 
      cbind("H2" = H2_org[[x]])})
  
  auc_higher_org_AP <- map(seq(n_webs), function(x) {
    bind_cols(pluck(auc_higher_org_AP, x),
              "nrow" = nrow(pluck(init_sim_web, x)),
              "ncol" = ncol(pluck(init_sim_web, x)))
  }) %>% bind_rows()
  
## auc values sim webs ----  
  # lower
  auc_lower <- map(seq(n_webs), function(x) {
    map(.x = 1:3, function(y) {
      map(.x = ctrbs, function(z) {
        map(seq(ncol(pluck(sp_remain_lower_sims, x, "lower", y, z))), function(.x) {
          auc(x = pluck(sp_remain_lower_sims, x,  "lower", y, z)[, .x] %>%
                rev(),
              y = pluck(sp_remain_lower_sims, x,  "higher", y, z)[, .x])}) %>% 
          unlist %>% 
          data.table("robustness" = .) %>% 
          cbind("com_vars" = names(ctrbs)[z]) %>% 
          cbind("p_single" = singleton(pluck(sims, x, z, "web"))) %>% 
          cbind("web_no" = x) %>% 
          cbind("H2" = pluck(H2, x, z))
      }) %>% bind_rows() %>% 
        cbind("rew" = rew_names[y])
    }) %>% bind_rows()
  })
  
  auc_lower <- map(seq(n_webs), function(x) {
    bind_cols(pluck(auc_lower, x),
              "nrow" = nrow(pluck(sims, x, 1, "web")) %>%
                rep(., nrow(pluck(auc_lower[[x]]))),
              "ncol" = ncol(pluck(sims, x, 1, "web")) %>% 
                rep(., nrow(pluck(auc_lower[[x]]))))
  }) %>% bind_rows()
  
  # higher
  auc_higher <- map(seq(n_webs), function(x) {
    map(.x = 1:3, function(y) {
      map(.x = ctrbs, function(z) {
        map(seq(ncol(pluck(sp_remain_higher_sims, x, "lower", y, z))), function(.x) {
          auc(x = pluck(sp_remain_higher_sims, x,  "lower", y, z)[, .x] %>%
                rev(),
              y = pluck(sp_remain_higher_sims, x,  "higher", y, z)[, .x])}) %>% 
          unlist %>% 
          data.table("robustness" = .) %>% 
          cbind("com_vars" = names(ctrbs)[z]) %>% 
          cbind("p_single" = singleton(pluck(sims, x, z, "web"), lower = F)) %>% 
          cbind("web_no" = x) %>% 
          cbind("H2" = pluck(H2, x, z))
      }) %>% bind_rows() %>% 
        cbind("rew" = rew_names[y])
    }) %>% bind_rows()
  })
  
  auc_higher <- map(seq(n_webs), function(x) {
    bind_cols(pluck(auc_higher, x),
              "nrow" = nrow(pluck(sims, x, 1, "web")) %>%
                rep(., nrow(pluck(auc_higher[[x]]))),
              "ncol" = ncol(pluck(sims, x, 1, "web")) %>% 
                rep(., nrow(pluck(auc_higher[[x]]))))
  }) %>% bind_rows
  
  # lower norew
  auc_lower_norew <- map(seq(n_webs), function(x) {
    map(.x = ctrbs, function(y) {
      map(seq(ncol(pluck(sp_remain_lower_norew_sims, x, "lower", y))), function(.x) {
        auc(x = pluck(sp_remain_lower_norew_sims, x,  "lower", y)[, .x] %>%
              rev(),
            y = pluck(sp_remain_lower_norew_sims, x,  "higher", y)[, .x])}) %>% 
        unlist %>% 
        data.table("robustness" = .) %>% 
        cbind("com_vars" = names(ctrbs)[y]) %>% 
        cbind("rew" = "norew") %>% 
        cbind("p_single" = singleton(pluck(sims, x, y, "web"))) %>% 
        cbind("web_no" = x) %>% 
        cbind("H2" = pluck(H2, x, y))
    }) %>% bind_rows()
  })
  
  auc_lower_norew <- map(seq(n_webs), function(x) {
    bind_cols(pluck(auc_lower_norew, x),
              "nrow" = nrow(pluck(sims, x, 1, "web")) %>%
                rep(., nrow(pluck(auc_lower_norew[[x]]))),
              "ncol" = ncol(pluck(sims, x, 1, "web")) %>% 
                rep(., nrow(pluck(auc_lower_norew[[x]]))))
  }) %>% bind_rows()
  
  # higher norew
  auc_higher_norew <- map(seq(n_webs), function(x) {
    map(.x = ctrbs, function(y) {
      map(seq(ncol(pluck(sp_remain_higher_norew_sims, x, "lower", y))), function(.x) {
        auc(x = pluck(sp_remain_higher_norew_sims, x,  "lower", y)[, .x] %>%
              rev(),
            y = pluck(sp_remain_higher_norew_sims, x,  "higher", y)[, .x])}) %>% 
        unlist %>% 
        data.table("robustness" = .) %>% 
        cbind("com_vars" = names(ctrbs)[y]) %>% 
        cbind("rew" = "norew") %>% 
        cbind("p_single" = singleton(pluck(sims, x, y, "web"), lower = F)) %>% 
        cbind("web_no" = x) %>% 
        cbind("H2" = pluck(H2, x, y))
    }) %>% bind_rows()
  })
  
  auc_higher_norew <- map(seq(n_webs), function(x) {
    bind_cols(pluck(auc_higher_norew, x),
              "nrow" = nrow(pluck(sims, x, 1, "web")) %>%
                rep(., nrow(pluck(auc_higher_norew[[x]]))),
              "ncol" = ncol(pluck(sims, x, 1, "web")) %>% 
                rep(., nrow(pluck(auc_higher_norew[[x]]))))
  }) %>% bind_rows()
  
  # lower AT
  auc_lower_AT <- map(seq(n_webs), function(x) {
    map(.x = ctrbs, function(y) {
      map(seq(ncol(pluck(sp_remain_lower_AT_sims, x, "lower", y))), function(.x) {
        auc(x = pluck(sp_remain_lower_AT_sims, x,  "lower", y)[, .x] %>%
              rev(),
            y = pluck(sp_remain_lower_AT_sims, x,  "higher", y)[, .x])}) %>% 
        unlist %>% 
        data.table("robustness" = .) %>% 
        cbind("com_vars" = names(ctrbs)[y]) %>% 
        cbind("rew" = "abundtrait") %>% 
        cbind("p_single" = singleton(pluck(sims, x, y, "web"))) %>% 
        cbind("web_no" = x) %>% 
        cbind("H2" = pluck(H2, x, y))
    }) %>% bind_rows()
  })
  
  auc_lower_AT <- map(seq(n_webs), function(x) {
    bind_cols(pluck(auc_lower_AT, x),
              "nrow" = nrow(pluck(sims, x, 1, "web")) %>%
                rep(., nrow(pluck(auc_lower_AT[[x]]))),
              "ncol" = ncol(pluck(sims, x, 1, "web")) %>% 
                rep(., nrow(pluck(auc_lower_AT[[x]]))))
  }) %>% bind_rows()
  
  # higher AT
  auc_higher_AT <- map(seq(n_webs), function(x) {
    map(.x = ctrbs, function(y) {
      map(seq(ncol(pluck(sp_remain_higher_AT_sims, x, "lower", y))), function(.x) {
        auc(x = pluck(sp_remain_higher_AT_sims, x,  "lower", y)[, .x] %>%
              rev(),
            y = pluck(sp_remain_higher_AT_sims, x,  "higher", y)[, .x])}) %>% 
        unlist %>% 
        data.table("robustness" = .) %>% 
        cbind("com_vars" = names(ctrbs)[y]) %>% 
        cbind("rew" = "abundtrait") %>% 
        cbind("p_single" = singleton(pluck(sims, x, y, "web"), lower = F)) %>% 
        cbind("web_no" = x) %>% 
        cbind("H2" = pluck(H2, x, y))
    }) %>% bind_rows()
  })
  
  auc_higher_AT <- map(seq(n_webs), function(x) {
    bind_cols(pluck(auc_higher_AT, x),
              "nrow" = nrow(pluck(sims, x, 1, "web")) %>%
                rep(., nrow(pluck(auc_higher_AT[[x]]))),
              "ncol" = ncol(pluck(sims, x, 1, "web")) %>% 
                rep(., nrow(pluck(auc_higher_AT[[x]]))))
  }) %>% bind_rows()
  
  # lower AP
  auc_lower_AP <- map(seq(n_webs), function(x) {
    map(.x = ctrbs, function(y) {
      map(seq(ncol(pluck(sp_remain_lower_AP_sims, x, "lower", y))), function(.x) {
        auc(x = pluck(sp_remain_lower_AP_sims, x,  "lower", y)[, .x] %>%
              rev(),
            y = pluck(sp_remain_lower_AP_sims, x,  "higher", y)[, .x])}) %>% 
        unlist %>% 
        data.table("robustness" = .) %>% 
        cbind("com_vars" = names(ctrbs)[y]) %>% 
        cbind("rew" = "abundphylo") %>% 
        cbind("p_single" = singleton(pluck(sims, x, y, "web"))) %>% 
        cbind("web_no" = x) %>% 
        cbind("H2" = pluck(H2, x, y))
    }) %>% bind_rows()
  })
  
  auc_lower_AP <- map(seq(n_webs), function(x) {
    bind_cols(pluck(auc_lower_AP, x),
              "nrow" = nrow(pluck(sims, x, 1, "web")) %>%
                rep(., nrow(pluck(auc_lower_AP[[x]]))),
              "ncol" = ncol(pluck(sims, x, 1, "web")) %>% 
                rep(., nrow(pluck(auc_lower_AP[[x]]))))
  }) %>% bind_rows()
  
  # higher AP
  auc_higher_AP <- map(seq(n_webs), function(x) {
    map(.x = ctrbs, function(y) {
      map(seq(ncol(pluck(sp_remain_higher_AP_sims, x, "lower", y))), function(.x) {
        auc(x = pluck(sp_remain_higher_AP_sims, x,  "lower", y)[, .x] %>%
              rev(),
            y = pluck(sp_remain_higher_AP_sims, x,  "higher", y)[, .x])}) %>% 
        unlist %>% 
        data.table("robustness" = .) %>% 
        cbind("com_vars" = names(ctrbs)[y]) %>% 
        cbind("rew" = "abundphylo") %>% 
        cbind("p_single" = singleton(pluck(sims, x, y, "web"), lower = F)) %>% 
        cbind("web_no" = x) %>% 
        cbind("H2" = pluck(H2, x, y))
    }) %>% bind_rows()
  })
  
  auc_higher_AP <- map(seq(n_webs), function(x) {
    bind_cols(pluck(auc_higher_AP, x),
              "nrow" = nrow(pluck(sims, x, 1, "web")) %>%
                rep(., nrow(pluck(auc_higher_AP[[x]]))),
              "ncol" = ncol(pluck(sims, x, 1, "web")) %>% 
                rep(., nrow(pluck(auc_higher_AP[[x]]))))
  }) %>% bind_rows()
  
  
  # auc; initial extinction on lower level all scenarios
  auc_all_lower <- auc_lower %>% bind_rows(auc_lower_norew) %>% 
    bind_rows(auc_lower_org) %>% 
    bind_rows(auc_lower_org_norew) %>%
    bind_rows(auc_lower_AT) %>%
    bind_rows(auc_lower_AP) %>% 
    bind_rows(auc_lower_org_AT) %>%
    bind_rows(auc_lower_org_AP)
  
  # auc; initial extinction on higher level all scenarios
  auc_all_higher <- auc_higher %>% bind_rows(auc_higher_norew) %>% 
    bind_rows(auc_higher_org) %>% 
    bind_rows(auc_higher_org_norew) %>% 
    bind_rows(auc_higher_AT) %>%
    bind_rows(auc_higher_AP) %>% 
    bind_rows(auc_higher_org_AT) %>%
    bind_rows(auc_higher_org_AP)
  
  #  all auc
  auc_all <- bind_rows(cbind(auc_all_lower, "lvl" = "lower"),
  cbind(auc_all_higher, "lvl" = "higher"))
  