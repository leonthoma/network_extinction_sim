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
n_webs <- 101

# create initial network for community variables
init_sim <- map(seq(n_webs), ~ simulate_tapnet_aug(nlower = 20, nhigher = 30, ntraits_nopem = 2,
                                ntraits_pem = 2, abuns = "lognormal", Nobs = 1111,
                                names = sp_names, initial_sim = T))

# set no of nets and simulations
n_nets <- 4
n_sims <- 11
ctrbs <- c("Atl" = 1, "aTl" = 2, "atL" = 3, "ATL" = 4)
rew_names <- c("abund", "trait", "phylo")

# Create simulated networks with setting contributions ----
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


# visualize networks ----
# map(map(sims, 1), ~ plotweb(.x))
# map(seq(n_webs), ~ plotweb(pluck(sims, .x, "atL", "web")))

# basic metrics ----
# connectance
connectance <- function(x) {
  sum(x != 0) / (nrow(x) * ncol(x))
}

# histograms
# empty rows
par(mfrow = c(3,1))
map(1:3, function(x) {
  hist(map(1:100, ~ length(which(rowSums(sims_metrics[[.x]][[x]]$web) == 0))) %>%
       unlist, xlab = paste("empty rows", names(ctrbs)[x]), main = NULL)})

# empty cols
map(1:3, function(x) {
  hist(map(1:100, ~ length(which(colSums(sims_metrics[[.x]][[x]]$web) == 0))) %>%
         unlist, xlab = paste("empty cols", names(ctrbs)[x]), main = NULL)})

# two dimensional shannon entropy
H2 <- map(ctrbs, function(x) {
  map(sims, ~ H2fun(pluck(.x, x, "web"))[[1]])})

H2_org <- map(seq(n_webs), ~ H2fun(init_sim[[.x]]$networks[[1]]$web)[[1]]) %>% 
  unlist() %>%
  data.table("vals" = .) %>% 
  cbind("group" = "org") %>% bind_rows()

H2_grouped <- map(ctrbs, function(x) pluck(H2, x) %>%
                    unlist() %>%
                    data.table("vals" = .) %>% 
                    cbind("group" = names(ctrbs)[x])) %>% bind_rows()

H2_grouped <- bind_rows(H2_org, H2_grouped)

library(ggpubr)
ggboxplot(data = H2_grouped, x = "group", y = "vals", 
          color = "group", palette = "npg",
          order = c("Atl", "aTl", "atL", "ATL", "org"),
          ylab = "H2", xlab = "Contribution importances", legend = "none") 

# anova of all community variables 
cntc <- map(sims, function(x) map(1:4, ~ connectance(pluck(x, .x, "web")))) 

cntc_grouped <- map(ctrbs, function(x) map(cntc, ~ pluck(.x, x)) %>%
      unlist() %>%
      data.table("vals" = .) %>% 
      cbind("group" = names(ctrbs)[x])) %>% bind_rows()

anova(lm(cntc_grouped$vals ~ cntc_grouped$group))
kruskal.test(cntc_grouped$vals ~ cntc_grouped$group)

library(ggpubr)
ggboxplot(data = cntc_grouped, x = "group", y = "vals", 
          color = "group", palette = "npg",
          order = c("Atl", "aTl", "atL", "ATL"),
          ylab = "Connectance", xlab = "Contribution importances", legend = "right")

# # nestedness; nested rank
# map(sims, function(x) map(1:3, ~ nestedness(pluck(x, .x, "web")))) 

# rewiring partner choice ----
  # abundances
  abunds <- map(seq(n_webs), ~ pluck(init_sim, .x)$networks[[1]]$abuns)
  
  # traits 
  traits <- map(seq(n_webs), ~ pluck(init_sim, .x)$networks[[1]]$traits)
  
  # phylogenetic distances
  library(ape)
  
  phylos <- map(seq(n_webs), ~ list("low" = cophenetic.phylo(pluck(init_sim, .x)$trees$low), 
                    "high" = cophenetic.phylo(pluck(init_sim, .x)$trees$high)))
  
# run extinction models ----
source("rewiring_vizentin-bugoni_2019/Functions/extinction.mod.R")
source("one.second.extinct.mod.R")
source("rewiring_vizentin-bugoni_2019/Functions/calc.mean.one.second.extinct.mod.R")
source("rewiring_vizentin-bugoni_2019/Functions/matrix.p1.R")
source("rewiring_vizentin-bugoni_2019/Functions/IC.R") # calc of 95 percent confidence interval

# only use Atl, aTl, and atL for extc sims
n_nets <- 3
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


# compute means of all simulations
rew_names <- c("abund", "trait", "phylo")

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

# calculate percentages of remaining sp
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

# mean over all webs
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
                            map(seq(n_webs), ~ pluck(sp_remain_lower, .x, "lower", x, y))  %>%
                              as.data.table() %>%
                              match_lengths() %>%
                              rowMeans()})),
                        "higher" = map(.x = 1:3, function(x) map(.x = 1:3, function(y) {
                          map(seq(n_webs), ~ pluck(sp_remain_lower, .x, "higher", x, y))  %>%
                            as.data.table() %>%
                            match_lengths() %>%
                            rowMeans()})))


# initial extincion on higher level
sp_remain_higher_web_mean <- list("lower" =
                           map(.x = 1:3, function(x) map(.x = 1:3, function(y) {
                             map(seq(n_webs), ~ pluck(sp_remain_higher, .x, "lower", x, y))  %>%
                               as.data.table() %>%
                               match_lengths() %>%
                               rowMeans()})),
                        "higher" = map(.x = 1:3, function(x) map(.x = 1:3, function(y) {
                          map(seq(n_webs), ~ pluck(sp_remain_higher, .x, "higher", x, y))  %>%
                            as.data.table() %>%
                            match_lengths() %>%
                            rowMeans()})))

# visualize extinction models ----
plot_extc(sp_remain_lower_web_mean, save = F, view = T, lower = T)
plot_extc(sp_remain_higher_web_mean, save = F, view = T)

plot_extc_alt(sp_remain_lower_web_mean,
              sp_remain_lower_web_mean_org,
              sp_remain_lower_web_mean_norew,
              lower = T,
              save = F)
plot_extc_alt(sp_remain_higher_web_mean,
              sp_remain_higher_web_mean_org,
              sp_remain_higher_web_mean_norew,
              save = F)

# visualize deviation of sims
# # get auc of every web original web 
# aucs_org <- map(.x = 1:3, function(x) {
#   map(seq(n_webs), ~ auc(x = pluck(sp_remain_lower_org, .x, "lower", x) %>% 
#                            rev(), 
#                          y = pluck(sp_remain_lower_org, .x, "higher", x)))})
# 
# # get min/max auc original web
# auc_mins_org <- map(.x = 1:3, ~ pluck(aucs, .x) %>% unlist %>% min) # min
# 
# auc_maxs_org <- map(.x = 1:3, ~ pluck(aucs, .x) %>% unlist %>% max) # max
# 
# auc_meds_org <- map(.x = 1:3, ~ pluck(aucs, .x) %>% unlist %>% sort %>% median) # mean
# 
# # get min and max positions of sims based on auc original web
# mins_org <- map(1:3,~ which(pluck(aucs_org, .x) == pluck(auc_mins_org, .x)))
# 
# maxs_org <- map(1:3, ~ which(pluck(aucs_org, .x) == pluck(auc_maxs_org, .x)))
# 
# meds_org <- map(1:3, ~ which(pluck(aucs_org, .x) == pluck(auc_meds_org, .x)))

# get ci of every web 
# initial extinction on lower level
ci_lower <- list("lower" = list("2.5" = conf_int(x = sp_remain_lower,
                                                 mean = sp_remain_lower_web_mean),
                                "97.5" = conf_int(x = sp_remain_lower,
                                                  mean = sp_remain_lower_web_mean,
                                                  lower = F)),
                 "higher" = list("2.5" = conf_int(x = sp_remain_lower,
                                     mean = sp_remain_lower_web_mean,
                                     trph_lvl = "higher"),
                                 "97.5" = conf_int(x = sp_remain_lower,
                                                   mean = sp_remain_lower_web_mean,
                                                   trph_lvl = "higher",
                                                   lower = F)))

# initial extinction on higher level
ci_higher <- list("lower" = list("2.5" = conf_int(x = sp_remain_higher,
                                                 mean = sp_remain_higher_web_mean),
                                "97.5" = conf_int(x = sp_remain_higher,
                                                  mean = sp_remain_higher_web_mean,
                                                  lower = F)),
                 "higher" = list("2.5" = conf_int(x = sp_remain_higher,
                                                  mean = sp_remain_higher_web_mean,
                                                  trph_lvl = "higher"),
                                 "97.5" = conf_int(x = sp_remain_higher,
                                                   mean = sp_remain_higher_web_mean,
                                                   trph_lvl = "higher",
                                                   lower = F)))

# plot; intital extinction on lower trophic level
dev_plots <- map(c("Atl" = 1, "aTl" = 2, "atL" = 3), function(x) ggplot() + 
    # geom_line(aes(rev(pluck(ci_lower, 1, 1, 1, x)),
    #               pluck(ci_lower, 2, 1, 1, x),
    #               color = "Abundance"), linetype = 3, size = .15) + # lower
    # geom_line(aes(rev(pluck(ci_lower, 1, 2, 1, x)),
    #               pluck(ci_lower, 2, 2, 1, x),
    #               color = "Abundance"), linetype = 3, size = .15) + # upper
    #   geom_line(aes(rev(pluck(sp_remain_lower_web_mean, 1, x, 1)),
    #                 pluck(sp_remain_lower_web_mean, 2, x, 1),
    #                color = "Abundance"), linetype = 1, size = .5) + # mean
    # geom_line(aes(rev(pluck(ci_lower, 1, 1, 2, x)),
    #               pluck(ci_lower, 2, 1, 2, x),
    #               color = "Traits"), linetype = 2, size = .15) + # lower
    # geom_line(aes(rev(pluck(ci_lower, 1, 2, 2, x)),
    #               pluck(ci_lower, 2, 2, 2, x),
    #               color = "Traits"), linetype = 2, size = .15) + # upper
    # geom_line(aes(rev(pluck(sp_remain_lower_web_mean, 1, x, 2)),
    #               pluck(sp_remain_lower_web_mean, 2, x, 2),
    #               color = "Traits"), linetype = 1, size = .5) + # mean
    geom_line(aes(rev(pluck(ci_lower, 1, 1, 3, x)),
                  pluck(ci_lower, 2, 1, 3, x),
                  color = "Phylogeny"), linetype = 3, size = .15) + # lower
    geom_line(aes(rev(pluck(ci_lower, 1, 2, 3, x)),
                  pluck(ci_lower, 2, 2, 3, x),
                  color = "Phylogeny"), linetype = 3, size = .15) + # upper
    geom_line(aes(rev(pluck(sp_remain_lower_web_mean, 1, x, 3)),
                  pluck(sp_remain_lower_web_mean, 2, x, 3),
                  color = "Phylogeny"), linetype = 1, size = .5) + # mean
      scale_color_manual(name = "Rewiring Method",
                         values = c("Abundance" = "black",
                                    "Traits" = "firebrick",
                                    "Phylogeny" = "dodgerblue")) +
      labs(x = "plants removed", y = "animals persisting",
           title = paste("Extinction cascade", names(ctrbs)[x]))
      )

dev_plots

map(1:3, ~ ggsave(
  paste("deviation extinction cascade lower trophic level",
             names(c("Atl" = 1, "aTl" = 2, "atL" = 3)[.x])),
  path = paste0(getwd(), "/plot_sink"),
  plot = dev_plots[[.x]],
  device = "pdf",
  width = 1900,
  height = 1205,
  units = "px"))

# Deviance of all networks for one scenario
map(c("Abund" = 1, "Traits" = 2, "Phylo" = 3), ~ ggplot() +
  geom_line(aes(rev(pluck(sp_remain_lower, 1, 1, .x, "atL")),
                pluck(sp_remain_lower, 1, 2, .x, "atL"),
                color = "a"), linetype = 1) +
  geom_line(aes(rev(pluck(sp_remain_lower, 2, 1, .x, "atL")),
                pluck(sp_remain_lower, 2, 2, .x, "atL"),
                color = "b"), linetype = 1) +
  geom_line(aes(rev(pluck(sp_remain_lower, 3, 1, .x, "atL")),
                pluck(sp_remain_lower, 3, 2, .x, "atL"),
                color = "c"), linetype = 1) +
  geom_line(aes(rev(pluck(sp_remain_lower, 4, 1, .x, "atL")),
                pluck(sp_remain_lower, 4, 2, .x, "atL"),
                color = "d"), linetype = 1) +
  geom_line(aes(rev(pluck(sp_remain_lower, 5, 1, .x, "atL")),
                pluck(sp_remain_lower, 5, 2, .x, "atL"),
                color = "e"), linetype = 1) +
  geom_line(aes(rev(pluck(sp_remain_lower, 6, 1, .x, "atL")),
                pluck(sp_remain_lower, 6, 2, .x, "atL"),
                color = "f"), linetype = 1) +
  geom_line(aes(rev(pluck(sp_remain_lower, 7, 1, .x, "atL")),
                pluck(sp_remain_lower, 7, 2, .x, "atL"),
                color = "g"), linetype = 1) +
  geom_line(aes(rev(pluck(sp_remain_lower, 8, 1, .x, "atL")),
                pluck(sp_remain_lower, 8, 2, .x, "atL"),
                color = "h"), linetype = 1) +
  geom_line(aes(rev(pluck(sp_remain_lower, 9, 1, .x, "atL")),
                pluck(sp_remain_lower, 9, 2, .x, "atL"),
                color = "i"), linetype = 1) +
  geom_line(aes(rev(pluck(sp_remain_lower, 10, 1, .x, "atL")),
                pluck(sp_remain_lower, 10, 2, .x, "atL"),
                color = "j"), linetype = 1) +
  geom_line(aes(rev(pluck(sp_remain_lower, 11, 1, .x, "atL")),
                pluck(sp_remain_lower, 11, 2, .x, "atL"),
                color = "k"), linetype = 1) +
  labs(x = "plants removed", y = "animals persisting",
       title = paste("Extinction cascade lower level", rew_names[.x]))  
    )

## Robustness
# lower level extinction org; R_lower
r_lower <- map(1:3, function(x){
    map_dbl(seq(n_webs), function(y) {
      map_dbl(seq(n_sims), function(.x) {
        robustness_aug(pluck(extc_sims_lower_norew, y, x, .x))}) %>% 
        mean}) %>% mean})

# higher level extinction org; R_lower
r_higher <- map(1:3, function(x){
  map_dbl(seq(n_webs), function(y) {
    map_dbl(seq(n_sims), function(.x) {
      robustness_aug(pluck(extc_sims_higher_norew, y, x, .x), lower = F)}) %>% 
      mean}) %>% mean})

# lower level extinction; Rw_lower
rw_lower <- map(1:3, function(x){
  map(1:3, function(y) {
    map_dbl(seq(n_webs), function(z) {
      map_dbl(seq(n_sims), function(.x) {
        robustness_aug(pluck(extc_sims_lower, z, x, y, .x))}) %>% 
      mean}) %>% mean})})

# higher level extinction; Rw_higher
rw_higher <- map(1:3, function(x){
  map(1:3, function(y) {
    map_dbl(seq(n_webs), function(z) {
      map_dbl(seq(n_sims), function(.x) {
        robustness_aug(pluck(extc_sims_higher, z, x, y, .x), lower = F)}) %>% 
        mean}) %>% mean})})

# Effect sizes
d_lower <- map(1:3, function(x) map(1:3, ~ pluck(rw_lower, x, .x) - pluck(r_lower, .x)))
d_higher <- map(1:3, function(x) map(1:3, ~ pluck(rw_higher, x, .x) - pluck(r_higher, .x)))

map(1:3 , ~ ggplot() +
  geom_col(aes(x = "Atl", y = pluck(d_lower, .x,  1), fill = "Atl")) +
  geom_col(aes(x = "aTl", y = pluck(d_lower, .x, 2), fill = "aTl")) +
  geom_col(aes(x = "atL", y = pluck(d_lower, .x, 3), fill = "atL")) +
  scale_fill_manual(name = "Community importances",
                     values = c("Atl" = "black",
                                "aTl" = "firebrick",
                                "atL" = "dodgerblue")) +
  labs(x = "Community importances", y = "Effect size",
       title = paste("Effect size lower level", rew_names[.x])) +
  theme(legend.position = "none"))

map(1:3 , ~ ggplot() +
      geom_col(aes(x = "Atl", y = pluck(d_higher, .x,  1), fill = "Atl")) +
      geom_col(aes(x = "aTl", y = pluck(d_higher, .x, 2), fill = "aTl")) +
      geom_col(aes(x = "atL", y = pluck(d_higher, .x, 3), fill = "atL")) +
      scale_fill_manual(name = "Community importances",
                        values = c("Atl" = "black",
                                   "aTl" = "firebrick",
                                   "atL" = "dodgerblue")) +
      labs(x = "Community importances", y = "Effect size",
           title = paste("Effect size higher level", rew_names[.x]))+
      theme(legend.position = "none"))
