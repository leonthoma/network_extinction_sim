# Calculation & visualization of results
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

load("sims.RData") # load simulations

# visualize networks ----
# map(map(sims, 1), ~ plotweb(.x))
# map(seq(n_webs), ~ plotweb(pluck(sims, .x, "atL", "web")))

# basic metrics ----
## connectance
# anova of all community variables 
cntc <- map(sims_all, function(x) map(1:4, ~ connectance(pluck(x, .x, "web")))) 

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

## histograms of dead interactions
# empty rows
par(mfrow = c(3,1))
map(1:3, function(x) {
  hist(map(1:100, ~ length(which(rowSums(sims_metrics[[.x]][[x]]$web) == 0))) %>%
         unlist, xlab = paste("empty rows", names(ctrbs)[x]), main = NULL)})

# empty cols
map(1:3, function(x) {
  hist(map(1:100, ~ length(which(colSums(sims_metrics[[.x]][[x]]$web) == 0))) %>%
         unlist, xlab = paste("empty cols", names(ctrbs)[x]), main = NULL)})

## two dimensional shannon entropy
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
