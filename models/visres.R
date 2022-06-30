# Calculation & visualization of results
setwd("~/Documents/Uni/M.sc/Master Thesis/Networks/models/")
load("sims_NULL.RData") # load simulations

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
library(ggpubr)

# visualize networks ----
# map(map(sims, 1), ~ plotweb(.x))
# map(seq(n_webs), ~ plotweb(pluck(sims, .x, "atL", "web")))

# basic metrics ----
## connectance
# anova of all community variables 
cntc <- map(sims, function(x) map(seq(ctrbs), ~ connectance(pluck(x, .x, "web")))) 

cntc_org <- map(seq(n_webs), ~ connectance(init_sim[[.x]]$networks[[1]]$web)) %>% 
  unlist() %>% data.table("vals" = .) %>% cbind("group" = "atl") %>% bind_rows()

cntc_grouped <- map(ctrbs, function(x) map(cntc, ~ pluck(.x, x)) %>%
                      unlist() %>%
                      data.table("vals" = .) %>% 
                      cbind("group" = names(ctrbs)[x])) %>% bind_rows()
cntc_grouped <- bind_rows(cntc_org, cntc_grouped)

anova(lm(cntc_grouped$vals ~ cntc_grouped$group))
kruskal.test(cntc_grouped$vals ~ cntc_grouped$group)

ggboxplot(data = cntc_grouped, x = "group", y = "vals", 
          #color = "group", palette = "npg",
          order = c("Atl", "aTl", "atL", "ATL", "atl"),
          ylab = "Connectance", xlab = "Contribution importances", legend = "right")+
  theme(text = element_text(size = 20))

## histograms of dead interactions
# empty rows
## needs change in hist_dead if networks prior to deleting dead interactions should
## be used !!!
hist_dead(sims_all)

# empty cols
hist_dead(sims_all, lower = F)

## two dimensional shannon entropy
H2 <- map(seq(n_webs), function(x) {
  map(ctrbs, ~ H2fun(pluck(sims, x, .x, "web"))["H2"])})

H2_org <- map(seq(n_webs), ~ H2fun(pluck(init_sim_web, .x))[[1]]) %>% 
  unlist() %>%
  data.table("vals" = .) %>% 
  cbind("group" = "org") %>% bind_rows()

H2_grouped <- map(seq(n_webs), function(x) pluck(H2, x) %>%
                    unlist() %>%
                    data.table("vals" = .) %>% 
                    cbind("group" = names(ctrbs))) %>% bind_rows()

H2_grouped <- bind_rows(H2_org, H2_grouped)

ggboxplot(data = H2_grouped, x = "group", y = "vals", 
          #color = "group", palette = "npg",
          order = c("Atl", "aTl", "atL", "atl", "org"),
          ylab = "H2", xlab = "Contribution importances", legend = "none") +
  theme(text = element_text(size = 20))

## visualize extinction models ----

# count no of NULL in sims

## adapt for org/norew !!!
# # org
# count_abort(extc_sims_lower_org)
# count_abort(extc_sims_higher_org)
# 
# # org no rewiring
# count_abort(extc_sims_lower_org_norew)
# count_abort(extc_sims_higher_org_norew)

# lower
count_abort(extc_sims_lower)

# higher
count_abort(extc_sims_higher)

# lower no rewiring
count_abort(extc_sims_lower_norew)

# higher no rewiring
count_abort(extc_sims_higher_norew)

# # find sim with lowest aborted runs
# sims_null_lower_org_norew <- map(seq(n_webs), function (z) { # loop over webs
#       map(seq(n_sims), ~ pluck(extc_sims_lower_org_norew, z, .x) %>% is.null
#       ) %>% unlist %>% sum}) 
# 
# min_abort_lower_org_norew <- unlist(sims_null_lower_org_norew) %>% which.min
# 
# 
# sims_null_lower <- map(seq(n_webs), function (z) { # loop over webs
#   map(1:3, function(x) { # loop over rew
#     map(ctrbs, function(y) { # loop over com_vars
#       map(seq(n_sims), ~ pluck(extc_sims_lower, z, x, y, .x) %>% is.null # loop over sims
#       )}) %>% unlist}) %>% unlist %>% sum})
# 
# min_abort_lower <- unlist(sims_null_lower) %>% which.min


## Plot of all sims for one web
# lower org norew
ggplot(aes(rev(x), y), data = lower_org_norew_sims_df) + 
  geom_step(aes(group = sims1)) +
  geom_line(color = "firebrick", data = lower_org_norew_sims_df_mean, size = .8) +
  labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
  theme(aspect.ratio = 1)

# lower org
map(rew_names, function(z) {
  ggplot(aes(rev(x), y),
         data = filter(lower_org_sims_df,
                       id == z)) +
    geom_step(aes(group = sims1)) +
    geom_line(color = "firebrick", data = pluck(lower_org_sims_df_mean,
                                                which(rew_names == z)),
              size = .8) +
    labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
    theme(aspect.ratio = 1)
})

# lower norew
map(names(ctrbs), function(z) {
  ggplot(aes(rev(x), y), data = filter(lower_norew_sims_df,
                                       com_vars == z)) + 
    geom_step(aes(group = sims1)) +
    geom_line(color = "firebrick", data = pluck(lower_norew_sims_df_mean,
                                                z),
              size = .8) +
    labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
    theme(aspect.ratio = 1)
})

# lower AT
map(names(ctrbs), function(z) {
  ggplot(aes(rev(x), y), data = filter(lower_AT_sims_df,
                                       com_vars == z)) + 
    geom_step(aes(group = sims1)) +
    geom_line(color = "firebrick", data = pluck(lower_AT_sims_df_mean,
                                                z),
              size = .8) +
    labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
    theme(aspect.ratio = 1)
})

# lower AP
map(names(ctrbs), function(z) {
  ggplot(aes(rev(x), y), data = filter(lower_AP_sims_df,
                                       com_vars == z)) + 
    geom_step(aes(group = sims1)) +
    geom_line(color = "firebrick", data = pluck(lower_AP_sims_df_mean,
                                                z),
              size = .8) +
    labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
    theme(aspect.ratio = 1)
})

# lower
map(rew_names, function(z) {
  map(names(ctrbs), ~ ggplot(aes(rev(x), y), data = filter(lower_sims_df,
                                                           id == z,
                                                           com_vars == .x)) + 
        geom_step(aes(group = sims1)) +
        geom_line(color = "firebrick", data = pluck(lower_sims_df_mean,
                                                    which(rew_names == z),
                                                    .x),
                  size = .8) +
        labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
        theme(aspect.ratio = 1)
  )})

# higher org norew
ggplot(aes(rev(x), y), data = higher_org_norew_sims_df) + 
  geom_step(aes(group = sims1)) +
  geom_line(color = "firebrick", data = higher_org_norew_sims_df_mean, size = .8) +
  labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
  theme(aspect.ratio = 1)

# higher org
map(rew_names, function(z) {
  ggplot(aes(rev(x), y),
         data = filter(higher_org_sims_df,
                       id == z)) +
    geom_step(aes(group = sims1)) +
    geom_line(color = "firebrick", data = pluck(higher_org_sims_df_mean,
                                                which(rew_names == z)),
              size = .8) +
    labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
    theme(aspect.ratio = 1)
})

# higher norew
map(names(ctrbs), function(z) {
  ggplot(aes(rev(x), y), data = filter(higher_norew_sims_df,
                                       com_vars == z)) + 
    geom_step(aes(group = sims1)) +
    geom_line(color = "firebrick", data = pluck(higher_norew_sims_df_mean,
                                                z),
              size = .8) +
    labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
    theme(aspect.ratio = 1)
})

# higher AT
map(names(ctrbs), function(z) {
  ggplot(aes(rev(x), y), data = filter(higher_AT_sims_df,
                                       com_vars == z)) + 
    geom_step(aes(group = sims1)) +
    geom_line(color = "firebrick", data = pluck(higher_AT_sims_df_mean,
                                                z),
              size = .8) +
    labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
    theme(aspect.ratio = 1)
})

# higher AP
map(names(ctrbs), function(z) {
  ggplot(aes(rev(x), y), data = filter(higher_AP_sims_df,
                                       com_vars == z)) + 
    geom_step(aes(group = sims1)) +
    geom_line(color = "firebrick", data = pluck(higher_AP_sims_df_mean,
                                                z),
              size = .8) +
    labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
    theme(aspect.ratio = 1)
})

# higher
map(rew_names, function(z) {
  map(names(ctrbs), ~ ggplot(aes(rev(x), y), data = filter(higher_sims_df,
                                                           id == z,
                                                           com_vars == .x)) + 
        geom_step(aes(group = sims1)) +
        geom_line(color = "firebrick", data = pluck(higher_sims_df_mean,
                                                    which(rew_names == z),
                                                    .x),
                  size = .8) +
        labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
        theme(aspect.ratio = 1)
  )})

# plot of all mean extc sims

prew <- c("abund" = 3) # set rewiring method to plot
pcm <- "Atl" # set com_var to plot

# by com_vars
ggplot(mapping = aes(rev(x), y)) +
  geom_line(data = pluck(lower_sims_df_mean, which(rew_names == names(prew)), pcm), aes(color = "red")) + # lower
  geom_line(data = pluck(lower_norew_sims_df_mean, pcm), aes(color = "blue")) + # lower norew
  geom_line(data = pluck(lower_org_sims_df_mean, prew), aes(color = "green")) + # lower org
  # geom_line(data = lower_org_norew_sims_df_mean, aes(color = "yellow")) + # lower org norew
  # geom_line(data = pluck(lower_AT_sims_df_mean, pcm), aes(color = "purple")) + # lower AT
  # geom_line(data = pluck(lower_AP_sims_df_mean, pcm), aes(color = "brown")) + # lower AP
  # geom_line(data = lower_org_AT_sims_df_mean, aes(color = "black")) + # lower org AT
  # geom_line(data = lower_org_AP_sims_df_mean, aes(color = "grey")) + # lower org AP
  scale_color_manual(values = c("red",
                                "blue",
                                "green", 
                                "yellow",
                                "purple",
                                "brown",
                                "black",
                                "grey"), labels = c(paste(names(prew), pcm),
                                                 "norew",
                                                 "org",
                                                 "org norew",
                                                 paste("AT", pcm),
                                                 paste("AP", pcm),
                                                 "org AT",
                                                 "org AP")) +
  labs(title = paste(names(prew), pcm), xlab = "Primary extinctions",
       ylab = "Secondary extinctions")
  
### AUC / Robustness
## ANOVA of auc values
# overall
aov_lo <- aov(robustness ~ nrow + ncol + rew + com_vars, data = auc_all_lower) # lower
aov_hi <- aov(robustness ~ nrow + ncol + rew + com_vars, data = auc_all_higher) # higher

# aov_both <- aov(robustness ~ nrow + ncol + rew + com_vars + lvl, data = auc_all) # combined

# use type III test from car package
library(car)
Anova(aov_lo, type = "III")
Anova(aov_hi, type = "III")
# Anova(aov_both, type = "III")

# by com_vars
map(c(names(ctrbs), "org"), ~
      aov(robustness ~ nrow + ncol + rew,
          data = filter(auc_all_lower, com_vars == .x))) # lower
map(c(names(ctrbs), "org"), ~
      aov(robustness ~ nrow + ncol + rew,
          data = filter(auc_all_higher, com_vars == .x))) # lower

# by rew
map(unique(auc_all_lower$rew), ~
      aov(robustness ~ nrow + ncol + com_vars,
          data = filter(auc_all_lower, rew == .x))) # lower
map(unique(auc_all_lower$rew), ~
      aov(robustness ~ nrow + ncol + com_vars,
          data = filter(auc_all_higher, rew == .x))) # lower

# AUC boxplots
# by rewiring
auc_boxplot(auc_all, x.axis = "rew")

auc_rew_plots <- map(c(names(ctrbs), "org"), function(x) {
  auc_boxplot(filter(auc_all, com_vars == x), x.axis = "rew") +
    ggtitle(x)
  })

auc_cv_by_rew <- ggarrange(auc_rew_plots[[1]], auc_rew_plots[[2]], auc_rew_plots[[3]],
          auc_rew_plots[[4]], auc_rew_plots[[5]], auc_rew_plots[[6]], common.legend = T)
ggsave("auc_com_vars_by_rew.pdf",
              path = paste0(getwd(), "/plot_sink"),
              plot = auc_cv_by_rew,
              device = "pdf",
              width = 28.7,
              height = 20,
              units = "cm")

# by com vars
auc_boxplot(auc_all, x.axis = "com_vars")

auc_com_var_plots <- map(c(rew_names,"norew", "abundtrait", "abundphylo"), function(x) {
  auc_boxplot(filter(auc_all, rew == x), x.axis = "com_vars") +
    ggtitle(x)
})

auc_rew_by_cv <- ggarrange(auc_com_var_plots[[1]], auc_com_var_plots[[2]], auc_com_var_plots[[3]],
          auc_com_var_plots[[4]], auc_com_var_plots[[5]],auc_com_var_plots[[6]], common.legend = T)

ggsave("auc_rew_by_com_vars.pdf",
       path = paste0(getwd(), "/plot_sink"),
       plot = auc_cv_by_rew,
       device = "pdf",
       width = 28.7,
       height = 20,
       units = "cm")

# by specific combination
auc_all_plots <- map(.x = c(names(ctrbs), "org"), function(x) { # loop com_vars
  map(.x = c(rew_names, "norew", "abundtrait", "abundphylo"), # loop rewiring methods
 ~ auc_boxplot(auc_all, x.axis = c(x, .x), save = T))
  })

auc_abund_pair <- ggarrange(auc_all_plots[[1]][[1]], auc_all_plots[[1]][[4]], auc_all_plots[[1]][[5]],
          auc_all_plots[[1]][[6]], nrow = 2, ncol = 2, common.legend = T) 

auc_trait_pair <- ggarrange(auc_all_plots[[2]][[2]], auc_all_plots[[2]][[4]], auc_all_plots[[2]][[5]],
                       auc_all_plots[[2]][[6]], nrow = 2, ncol = 2, common.legend = T) 

auc_phylo_pair <- ggarrange(auc_all_plots[[3]][[3]], auc_all_plots[[3]][[4]], auc_all_plots[[3]][[5]],
                            auc_all_plots[[3]][[6]], nrow = 2, ncol = 2, common.legend = T) 






# base model; lower 
plot_extc_facet(extc = sp_remain_lower_web_mean_df,
                extc_norew = sp_remain_lower_web_mean_df_norew,
                ci = ci_lower_df,
                ci_norew = ci_lower_df_norew,
                org = sp_remain_lower_web_mean_df_org,
                org_norew = sp_remain_lower_web_mean_df_org_norew,
                ci_org = ci_lower_df_org,
                ci_org_norew = ci_lower_df_org_norew,
                save = F)

# base model; higher
plot_extc_facet(extc = sp_remain_higher_web_mean_df,
                extc_norew = sp_remain_higher_web_mean_df_norew,
                ci = ci_higher_df,
                ci_norew = ci_higher_df_norew,
                org = sp_remain_higher_web_mean_df_org,
                org_norew = sp_remain_higher_web_mean_df_org_norew,
                ci_org = ci_higher_df_org,
                ci_org_norew = ci_higher_df_org_norew, save = F)

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
                ci_norew = ci_both_df_norew, save = F, both = T)

plot_extc_alt(x = sp_remain_lower_web_mean,
              org = sp_remain_lower_web_mean_org,
              ci = ci_lower, 
              ci_org = ci_lower_org,
              lower = T,
              save = F,
              spaghetti = T, by_com_var = T)

plot_extc_alt(sp_remain_higher_web_mean,
              sp_remain_higher_web_mean_org,
              ci_higher,
              ci_higher_org,
              save = F, spaghetti = T, by_com_var = T)

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
# lower level extinction norew; R_lower
r_lower <- map(1:3, function(x){
  map_dbl(seq(n_webs), function(y) {
    map_dbl(seq(n_sims), function(.x) {
      robustness_aug(pluck(extc_sims_lower_norew, y, x, .x))}) %>% 
      mean}) %>% mean})

# higher level extinction norew; R_lower
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

