# Calculation & visualization of results
setwd("~/Documents/Uni/M.sc/Master Thesis/Networks/models/")
load("sims.RData") # load simulations

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


# visualize networks ----
# map(map(sims, 1), ~ plotweb(.x))
# map(seq(n_webs), ~ plotweb(pluck(sims, .x, "atL", "web")))

# basic metrics ----
## connectance
# anova of all community variables 
cntc <- map(sims_all, function(x) map(1:4, ~ connectance(pluck(x, .x, "web")))) 

cntc_org <- map(seq(n_webs), ~ connectance(init_sim[[.x]]$networks[[1]]$web)) %>% 
  unlist() %>% data.table("vals" = .) %>% cbind("group" = "atl") %>% bind_rows()

cntc_grouped <- map(ctrbs, function(x) map(cntc, ~ pluck(.x, x)) %>%
                      unlist() %>%
                      data.table("vals" = .) %>% 
                      cbind("group" = names(ctrbs)[x])) %>% bind_rows()
cntc_grouped <- bind_rows(cntc_org, cntc_grouped)

anova(lm(cntc_grouped$vals ~ cntc_grouped$group))
kruskal.test(cntc_grouped$vals ~ cntc_grouped$group)

library(ggpubr)
ggboxplot(data = cntc_grouped, x = "group", y = "vals", 
          #color = "group", palette = "npg",
          order = c("Atl", "aTl", "atL", "ATL", "atl"),
          ylab = "Connectance", xlab = "Contribution importances", legend = "right")+
  theme(text = element_text(size = 20))

## histograms of dead interactions
# empty rows
hist_dead(sims)

# empty cols
hist_dead(sims, lower = F)

## two dimensional shannon entropy
H2 <- map(ctrbs, function(x) {
  map(sims_all, ~ H2fun(pluck(.x, x, "web"))["H2"])})

H2_org <- map(seq(n_webs), ~ H2fun(init_sim[[.x]]$networks[[1]]$web)[[1]]) %>% 
  unlist() %>%
  data.table("vals" = .) %>% 
  cbind("group" = "atl") %>% bind_rows()

H2_grouped <- map(ctrbs, function(x) pluck(H2, x) %>%
                    unlist() %>%
                    data.table("vals" = .) %>% 
                    cbind("group" = names(ctrbs)[x])) %>% bind_rows()

H2_grouped <- bind_rows(H2_org, H2_grouped)

library(ggpubr)
ggboxplot(data = H2_grouped, x = "group", y = "vals", 
          #color = "group", palette = "npg",
          order = c("Atl", "aTl", "atL", "ATL", "atl"),
          ylab = "H2", xlab = "Contribution importances", legend = "none") +
  theme(text = element_text(size = 20))

# visualize extinction models ----

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

# plot extinction sequences w/ CIs
# plot_extc(sp_remain_lower_web_mean, ci_lower, save = F, view = T, lower = T)
# plot_extc(sp_remain_higher_web_mean, ci_higher, save = F, view = T)

# base model; lower 
plot_extc_facet(extc = sp_remain_lower_web_mean_df,
                extc_norew = sp_remain_lower_web_mean_df_norew,
                ci = ci_lower_df,
                ci_norew = ci_lower_df_norew, save = T)

# base model; higher
plot_extc_facet(extc = sp_remain_higher_web_mean_df,
                extc_norew = sp_remain_higher_web_mean_df_norew,
                ci = ci_higher_df,
                ci_norew = ci_higher_df_norew, save = T)

## Extc alg aborts; gets caught in loop when sp w/ lowest abund only has dead interactions
# # abund model; lower
# plot_extc_facet(extc = sp_remain_lower_web_mean_df_abund,
#                 extc_norew = sp_remain_lower_web_mean_df_abund_norew,
#                 ci = ci_lower_df_abund,
#                 ci_norew = ci_lower_df_abund_norew, save = T, abund = T)
# 
# # abund model; higher
# plot_extc_facet(extc = sp_remain_higher_web_mean_df_abund,
#                 extc_norew = sp_remain_higher_web_mean_df_abund_norew,
#                 ci = ci_higher_df_abund,
#                 ci_norew = ci_higher_df_abund_norew, save = T, abund = T)

# both model
plot_extc_facet(extc = sp_remain_both_web_mean_df,
                extc_norew = sp_remain_both_web_mean_df_norew,
                ci = ci_both_df,
                ci_norew = ci_both_df_norew, save = T, both = T)

# org model; lower
plot_extc_facet(extc = sp_remain_lower_web_mean_df_org,
                extc_norew = sp_remain_lower_web_mean_df_org_norew,
                ci = ci_lower_df_org,
                ci_norew = ci_lower_df_org_norew, save = T, org = T)

# org model; higher
plot_extc_facet(extc = sp_remain_higher_web_mean_df_org,
                extc_norew = sp_remain_higher_web_mean_df_org_norew,
                ci = ci_higher_df_org,
                ci_norew = ci_higher_df_org_norew, save = T, org = T)

plot_extc_alt(x = sp_remain_lower_web_mean,
              org = sp_remain_lower_web_mean_org,
              ci = ci_lower, 
              ci_org = ci_lower_org,
              lower = T,
              save = T,
              spaghetti = T)

plot_extc_alt(sp_remain_higher_web_mean,
              sp_remain_higher_web_mean_org,
              ci_higher,
              ci_higher_org,
              save = T, spaghetti = T)

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

# Number of shifts per species
# run extc_sim with shift = T; creates shifts object
map(unique(names(unlist(shifts$lower))), ~unlist(shifts$lower)[names(unlist(shifts$lower)) == .x])
map(unique(names(unlist(shifts$higher))), ~unlist(shifts$higher)[names(unlist(shifts$higher)) == .x])

