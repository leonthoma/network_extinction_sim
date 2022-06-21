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
H2 <- map(seq(n_webs), function(x) {
  map(ctrbs, ~ H2fun(pluck(sims_all, x, .x, "web"))["H2"])})

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



# AUC / Robustness
source("helper_functions.R")
# # lower org; web_mean
# map(.x = 1:3, function(x) {
#    auc(x = pluck(sp_remain_lower_web_mean_org, "lower", x) %>%
#                            rev(),
#                          y = pluck(sp_remain_lower_web_mean_org, "higher", x))})
# 
# # lower org norew; web_mean
# auc(x = pluck(sp_remain_lower_web_mean_org_norew, "lower") %>%
#                            rev(),
#                          y = pluck(sp_remain_lower_web_mean_org_norew, "higher"))

# lower org; all
auc_lower_org <- map(.x = 1:3, function(x) {
  map(seq(n_webs), ~ auc(x = pluck(sp_remain_lower_org, .x, "lower", x) %>%
                           rev(),
                         y = pluck(sp_remain_lower_org, .x, "higher", x))) %>% unlist %>% 
    data.table("robustness" = .) %>% 
    cbind("rew" = rew_names[x]) %>% 
    cbind("com_vars" = "org") 
    }) %>% bind_rows()

auc_lower_org <- bind_cols(auc_lower_org,
                           "nrow" = map(seq(n_webs), function(x) {
                             nrow(pluck(init_sim_web, x))
                             }) %>% unlist() %>% rep(., 3),
                           "ncol" = map(seq(n_webs), function(x) {
                             ncol(pluck(init_sim_web, x))
                           }) %>% unlist() %>% rep(., 3))
# # lower org; all sd
# map(.x = 1:3, function(x) {
#   map(seq(n_webs), ~ auc(x = pluck(sp_remain_lower_org, .x, "lower", x) %>%
#                            rev(),
#                          y = pluck(sp_remain_lower_org, .x, "higher", x))) %>% unlist %>% sd})

# lower org norew; all
auc_lower_org_norew <- map(seq(n_webs), ~ auc(x = pluck(sp_remain_lower_org_norew, .x, "lower") %>%
                         rev(),
                       y = pluck(sp_remain_lower_org_norew, .x, "higher"))) %>% unlist %>% data.table("robustness" = .) %>% 
  cbind("rew" = "norew") %>%
  cbind("com_vars" = "org") %>%
  bind_rows()

auc_lower_org <- bind_cols(auc_lower_org,
                           "nrow" = map(seq(n_webs), function(x) {
                             nrow(pluck(init_sim_web, x))
                           }) %>% unlist(),
                           "ncol" = map(seq(n_webs), function(x) {
                             ncol(pluck(init_sim_web, x))
                           }) %>% unlist())

# # lower org norew; all sd
# map(seq(n_webs), ~ auc(x = pluck(sp_remain_lower_org_norew, .x, "lower") %>%
#       rev(),
#     y = pluck(sp_remain_lower_org_norew, .x, "higher"))) %>% unlist %>% sd

# # higher org; web_mean
# map(.x = 1:3, function(x) {
#   auc(x = pluck(sp_remain_higher_web_mean_org, "lower", x) %>%
#         rev(),
#       y = pluck(sp_remain_higher_web_mean_org, "higher", x))})
# 
# 
# # higher org norew; web_mean
# auc(x = pluck(sp_remain_higher_web_mean_org_norew, "lower") %>%
#       rev(),
#     y = pluck(sp_remain_higher_web_mean_org_norew, "higher"))

# higher org; all
auc_higher_org <- map(.x = 1:3, function(x) {
  map(seq(n_webs), ~ auc(x = pluck(sp_remain_higher_org, .x, "lower", x) %>%
                           rev(),
                         y = pluck(sp_remain_higher_org, .x, "higher", x))) %>% unlist %>%
    data.table("robustness" = .) %>% 
    cbind("rew" = rew_names[x]) %>% 
    cbind("com_vars" = "org")
}) %>% bind_rows()

auc_higher_org <- bind_cols(auc_higher_org,
                           "nrow" = map(seq(n_webs), function(x) {
                             nrow(pluck(init_sim_web, x))
                           }) %>% unlist() %>% rep(., 3),
                           "ncol" = map(seq(n_webs), function(x) {
                             ncol(pluck(init_sim_web, x))
                           }) %>% unlist() %>% rep(., 3))


# # higher org; all sd
# map(.x = 1:3, function(x) {
#   map(seq(n_webs), ~ auc(x = pluck(sp_remain_higher_org, .x, "lower", x) %>%
#                            rev(),
#                          y = pluck(sp_remain_higher_org, .x, "higher", x))) %>% unlist %>% sd})

# higher org norew; all
auc_higher_org_norew <- map(seq(n_webs), ~ auc(x = pluck(sp_remain_higher_org_norew, .x, "lower") %>%
                         rev(),
                       y = pluck(sp_remain_higher_org_norew, .x, "higher"))) %>% unlist %>% 
  data.table("robustness" = .) %>%
  cbind("rew" = "norew") %>%
  cbind("com_vars" = "org") %>%
  bind_rows()

auc_higher_org <- bind_cols(auc_higher_org,
                           "nrow" = map(seq(n_webs), function(x) {
                             nrow(pluck(init_sim_web, x))
                           }) %>% unlist(),
                           "ncol" = map(seq(n_webs), function(x) {
                             ncol(pluck(init_sim_web, x))
                           }) %>% unlist())

# # higher org norew; all sd
# map(seq(n_webs), ~ auc(x = pluck(sp_remain_higher_org_norew, .x, "lower") %>%
#                          rev(),
#                        y = pluck(sp_remain_higher_org_norew, .x, "higher"))) %>% unlist %>% sd

# # lower; web_mean
# map(.x = 1:3, function(x) {
#   map(.x = 1:3, function(z) {
#     auc(x = pluck(sp_remain_lower_web_mean, "lower", z, x) %>%
#                            rev(),
#                          y = pluck(sp_remain_lower_web_mean, "higher", z, x))})
#   })

# lower; all
auc_lower <- map(.x = 1:3, function(x) {
  map(.x = 1:3, function(z) {
    map(seq(n_webs), ~ auc(x = pluck(sp_remain_lower,.x,  "lower", z, x) %>%
                             rev(),
                           y = pluck(sp_remain_lower,.x,  "higher", z, x))) %>% unlist %>% 
      data.table("robustness" = .) %>% 
      cbind("rew" = rew_names[z])
  }) %>% bind_rows() %>% 
    cbind("com_vars" = names(ctrbs[-4])[x])
}) %>% bind_rows()

auc_lower <- bind_cols(auc_lower,
                           "nrow" = map(1:3, function(y) { 
                             map(seq(n_webs), function(x) {
                               nrow(pluck(sims, x, y, "web"))
                           })}) %>% unlist() %>% rep(., 3),
                           "ncol" = map(1:3, function(y) {
                             map(seq(n_webs), function(x) {
                               ncol(pluck(sims, x, y, "web"))
                           })}) %>% unlist() %>% rep(., 3))

# # lower; all sd
# map(.x = 1:3, function(x) {
#   map(.x = 1:3, function(z) {
#     map(seq(n_webs), ~ auc(x = pluck(sp_remain_lower,.x,  "lower", z, x) %>%
#                              rev(),
#                            y = pluck(sp_remain_lower,.x,  "higher", z, x))) %>% unlist %>% sd})})
# 
# # higher; web_mean
# map(.x = 1:3, function(x) {
#   map(.x = 1:3, function(z) {
#     auc(x = pluck(sp_remain_higher_web_mean, "lower", z, x) %>%
#                              rev(),
#                            y = pluck(sp_remain_higher_web_mean, "higher", z, x))})
# })

# higher; all
auc_higher <- map(.x = 1:3, function(x) {
  map(.x = 1:3, function(z) {
    map(seq(n_webs), ~ auc(x = pluck(sp_remain_higher,.x,  "lower", z, x) %>%
                             rev(),
                           y = pluck(sp_remain_higher,.x,  "higher", z, x))) %>% unlist %>%
      data.table("robustness" = .) %>% 
      cbind("rew" = rew_names[z])
  }) %>% bind_rows() %>% 
    cbind("com_vars" = names(ctrbs[-4])[x])
}) %>% bind_rows()

auc_higher <- bind_cols(auc_higher,
                       "nrow" = map(1:3, function(y) { 
                         map(seq(n_webs), function(x) {
                           nrow(pluck(sims, x, y, "web"))
                         })}) %>% unlist() %>% rep(., 3),
                       "ncol" = map(1:3, function(y) {
                         map(seq(n_webs), function(x) {
                           ncol(pluck(sims, x, y, "web"))
                         })}) %>% unlist() %>% rep(., 3))

# # higher; all sd
# map(.x = 1:3, function(x) {
#   map(.x = 1:3, function(z) {
#     map(seq(n_webs), ~ auc(x = pluck(sp_remain_higher,.x,  "lower", z, x) %>%
#                              rev(),
#                            y = pluck(sp_remain_higher,.x,  "higher", z, x))) %>% unlist %>% sd})})

# # lower norew; web_mean
# map(.x = 1:3, function(x) {
#   auc(x = pluck(sp_remain_lower_web_mean_norew, "lower", x) %>%
#                            rev(),
#                          y = pluck(sp_remain_lower_web_mean_norew, "higher", x))})
# 
# # higher norew; web_mean
# map(.x = 1:3, function(x) {
#   auc(x = pluck(sp_remain_higher_web_mean_norew, "lower", x) %>%
#                            rev(),
#                          y = pluck(sp_remain_higher_web_mean_norew, "higher", x))})

# lower norew; all
auc_lower_norew <- map(.x = 1:3, function(x) {
  map(seq(n_webs), ~ auc(x = pluck(sp_remain_lower_norew, .x, "lower", x) %>%
        rev(),
      y = pluck(sp_remain_lower_norew, .x, "higher", x))) %>% unlist %>%
    data.table("robustness" = .) %>%
    cbind("rew" = "norew") %>%
    cbind("com_vars" = rew_names[x])}) %>% 
  bind_rows()

auc_lower <- bind_cols(auc_lower,
                       "nrow" = map(1:3, function(y) { 
                         map(seq(n_webs), function(x) {
                           nrow(pluck(sims, x, y, "web"))
                         })}) %>% unlist(),
                       "ncol" = map(1:3, function(y) {
                         map(seq(n_webs), function(x) {
                           ncol(pluck(sims, x, y, "web"))
                         })}) %>% unlist())
# # lower norew; all sd
# map(.x = 1:3, function(x) {
#   map(seq(n_webs), ~ auc(x = pluck(sp_remain_lower_norew, .x, "lower", x) %>%
#                            rev(),
#                          y = pluck(sp_remain_lower_norew, .x, "higher", x))) %>% unlist %>% sd})

# higher norew; all
auc_higher_norew <- map(.x = 1:3, function(x) {
  map(seq(n_webs), ~ auc(x = pluck(sp_remain_higher_norew, .x, "lower", x) %>%
        rev(),
      y = pluck(sp_remain_higher_norew, .x, "higher", x))) %>% unlist %>%
    data.table("robustness" = .) %>%
    cbind("rew" = "norew") %>%
    cbind("com_vars" = rew_names[x])}) %>% 
  bind_rows()

auc_higher <- bind_cols(auc_higher,
                       "nrow" = map(1:3, function(y) { 
                         map(seq(n_webs), function(x) {
                           nrow(pluck(sims, x, y, "web"))
                         })}) %>% unlist(),
                       "ncol" = map(1:3, function(y) {
                         map(seq(n_webs), function(x) {
                           ncol(pluck(sims, x, y, "web"))
                         })}) %>% unlist())
# # higher norew; all sd
# map(.x = 1:3, function(x) {
#   map(seq(n_webs), ~ auc(x = pluck(sp_remain_higher_norew, .x, "lower", x) %>%
#                            rev(),
#                          y = pluck(sp_remain_higher_norew, .x, "higher", x))) %>% unlist %>% sd})

## differences between all and mean due to jensen inequality

## Anova of Auc values 
# initial extinction on lower level
auc_all_lower <- auc_lower %>% bind_rows(auc_lower_norew) %>% 
  bind_rows(auc_lower_org) %>% 
  bind_rows(auc_lower_org_norew)

# initial extinction on higher level
auc_all_higher <- auc_higher %>% bind_rows(auc_higher_norew) %>% 
  bind_rows(auc_higher_org) %>% 
  bind_rows(auc_higher_org_norew)

aov(robustness ~ nrow + ncol + rew + com_vars, data = auc_all_lower)


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

