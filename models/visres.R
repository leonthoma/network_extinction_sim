# Calculation & visualization of results
setwd("~/Documents/Uni/M.sc/Master Thesis/Networks/models/")
load("sims_NULL_full.RData") # load simulations

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
  unlist() %>% data.table("vals" = .) %>% cbind("group" = "org") %>% bind_rows()

cntc_grouped <- map(ctrbs, function(x) map(cntc, ~ pluck(.x, x)) %>%
                      unlist() %>%
                      data.table("vals" = .) %>% 
                      cbind("group" = names(ctrbs)[x])) %>% bind_rows()
cntc_grouped <- bind_rows(cntc_org, cntc_grouped)

anova(lm(cntc_grouped$vals ~ cntc_grouped$group))
kruskal.test(cntc_grouped$vals ~ cntc_grouped$group)

ggboxplot(data = cntc_grouped, x = "group", y = "vals", 
          #color = "group", palette = "npg",
          order = c("Atl", "aTl", "atL", "atl", "ATL", "org"),
          ylab = "Connectance", xlab = "Contribution importances", legend = "right")+
  theme(text = element_text(size = 20))

## histograms of dead interactions
# empty rows
## needs change in hist_dead if networks prior to deleting dead interactions should
## be used !!!
dead_low <- hist_dead(sims_all, save = F)

# empty cols
dead_high <- hist_dead(sims_all, lower = F, save = F)

ggarrange(dead_low, dead_high, labels = "AUTO")

## two dimensional shannon entropy
H2_org <- map(seq(n_webs), ~ H2_org[[.x]]) %>% 
  unlist() %>%
  data.table("vals" = .) %>% 
  cbind("group" = "org") %>% bind_rows()

H2_grouped <- map(seq(n_webs), function(x) pluck(H2, x) %>%
                    unlist() %>%
                    data.table("vals" = .) %>% 
                    cbind("group" = names(ctrbs))) %>% bind_rows()

H2_grouped <- bind_rows(H2_org, H2_grouped)

h2_plot <- ggboxplot(data = data.frame("vals" = H2_grouped$vals, "group" = H2_grouped$group), x = "group", y = "vals", 
          #color = "group", palette = "npg",
          order = c( "org", "Atl", "aTl", "atL", "atl", "ATL"),
          ylab = "H2'", xlab = "Community variable importances", legend = "none") +
  theme(text = element_text(size = 20), aspect.ratio = 1) +
  scale_x_discrete(labels = c( "Original", "Atl", "aTl", "atL", "atl", "ATL"))


ggsave("h2'_all.pdf",
       path = paste0(getwd(), "/plot_sink"),
       plot = h2_plot,
       device = "pdf",
       width = 28.7,
       height = 20,
       units = "cm")

# Percentage singletons by com_vars
ps_by_cv <- ggboxplot(data = auc_all, x = "lvl", y = "p_single", 
          color = "com_vars", palette = "npg",
          outlier.shape = NA,
          add = "mean",
          add.params = list(size = .2),
          ylab = "Percentage of singletons",
          xlab = "Species removed",
          ylim = c(0, .4)) +
  theme(text = element_text(size = 20), aspect.ratio = 1,
        legend.position = "bottom")  +
  scale_x_discrete(labels = c("Lower level", "Higher level")) +
  scale_color_discrete(name = "Community variable importance",
                       labels = c("Original", "Atl", "aTl", "atL", "atl", "ATL"))

ggsave("ps_by_com_vars.pdf",
       path = paste0(getwd(), "/plot_sink"),
       plot = ps_by_cv,
       device = "pdf",
       width = 28.7,
       height = 20,
       units = "cm")

# Robustness by percent singletons
titles <- c(names(ctrbs), "org")
names(titles) <- c(names(ctrbs), "Original")

auc_all_lower$com_vars <- factor(auc_all_lower$com_vars, 
                                    levels = c("org", "Atl", "aTl", "atL", "atl", "ATL"))

r_by_ps_lower <- group_by(auc_all_lower, web_no, com_vars) %>%
  summarise(m_r = mean(robustness, na.rm = T),
            m_p_s = mean(p_single)) %>%
  ggplot(aes(m_p_s, m_r), data = .) +
  geom_point(aes(color = com_vars)) +
  facet_wrap(vars(com_vars), scale = "free",
             labeller = as_labeller(c("org" = "Original", "Atl" = "Atl", 
                                      "aTl" = "aTl", "atL" = "atL",
                                      "atl" = "atl", "ATL" = "ATL"))) +
  labs(x = "Mean percentage of singletons",
       y = "Robustness") +
  scale_color_manual(values = c("#1b9e77",
                                "#d95f02",
                                "#7570b3",
                                "#e7298a",
                                "#66a61e",
                                "#e6ab02")) +
  theme(text = element_text(size = 20), aspect.ratio = 1,
        legend.position = "none")

ggsave("r_by_ps_lower.pdf",
       path = paste0(getwd(), "/plot_sink"),
       plot = r_by_ps_lower,
       device = "pdf",
       width = 28.7,
       height = 20,
       units = "cm")


auc_all_higher$com_vars <- factor(auc_all_higher$com_vars, 
                                 levels = c("org", "Atl", "aTl", "atL", "atl", "ATL"))

r_by_ps_higher <- group_by(auc_all_higher, web_no, com_vars) %>%
  summarise(m_r = mean(robustness, na.rm = T),
            m_p_s = mean(p_single)) %>%
  ggplot(aes(m_p_s, m_r), data = .) +
  geom_point(aes(color = com_vars)) +
  facet_wrap(vars(com_vars), scale = "free",
             labeller = as_labeller(c("org" = "Original", "Atl" = "Atl", 
                                      "aTl" = "aTl", "atL" = "atL",
                                      "atl" = "atl", "ATL" = "ATL"))) +
  labs(x = "Mean percentage of singletons",
       y = "Robustness") +
  scale_color_manual(values = c("#1b9e77",
                                "#d95f02",
                                "#7570b3",
                                "#e7298a",
                                "#66a61e",
                                "#e6ab02")) +
  theme(text = element_text(size = 20), aspect.ratio = 1,
        legend.position = "none")

ggsave("r_by_ps_higher.pdf",
       path = paste0(getwd(), "/plot_sink"),
       plot = r_by_ps_higher,
       device = "pdf",
       width = 28.7,
       height = 20,
       units = "cm")

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
  theme(aspect.ratio = 1, text = element_text(size = 20))

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
    theme(aspect.ratio = 1, text = element_text(size = 20))
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
    theme(aspect.ratio = 1, text = element_text(size = 20))
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
    theme(aspect.ratio = 1, text = element_text(size = 20))
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
    theme(aspect.ratio = 1, text = element_text(size = 20))
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
        theme(aspect.ratio = 1, text = element_text(size = 20))
  )})

# higher org norew
ggplot(aes(rev(x), y), data = higher_org_norew_sims_df) + 
  geom_step(aes(group = sims1)) +
  geom_line(color = "firebrick", data = higher_org_norew_sims_df_mean, size = .8) +
  labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
  theme(aspect.ratio = 1, text = element_text(size = 20))

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
    theme(aspect.ratio = 1, text = element_text(size = 20))
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
    theme(aspect.ratio = 1, text = element_text(size = 20))
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
    theme(aspect.ratio = 1, text = element_text(size = 20))
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
    theme(aspect.ratio = 1, text = element_text(size = 20))
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
        theme(aspect.ratio = 1, text = element_text(size = 20))
  )})

# com_vars by norew
# lower
cv_norew <- map(names(ctrbs), function(z) {
  ggplot(aes(rev(x), y), data = filter(lower_norew_sims_df,
                                       com_vars == z)) + 
    geom_step(aes(group = sims1)) +
    geom_line(color = "firebrick", data = pluck(lower_norew_sims_df_mean,
                                                z),
              size = .8) +
    labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
    theme(aspect.ratio = 1, text = element_text(size = 15)) +
    ggtitle(z)
})

cv_norew[[6]] <- ggplot(aes(rev(x), y), data = lower_org_norew_sims_df) + 
  geom_step(aes(group = sims1)) +
  geom_line(color = "firebrick", data = lower_org_norew_sims_df_mean, size = .8) +
  labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
  theme(aspect.ratio = 1, text = element_text(size = 15)) +
  ggtitle("Original")

cv_norew_lower_all <- ggarrange(cv_norew[[6]], cv_norew[[1]], cv_norew[[2]],
          cv_norew[[3]], cv_norew[[4]], cv_norew[[5]], nrow  = 2 , ncol = 3)

ggsave("extc_sims_cv_norew_lower.pdf",
       path = paste0(getwd(), "/plot_sink"),
       plot = cv_norew_lower_all,
       device = "pdf",
       width = 28.7,
       height = 20,
       units = "cm")

# higher
cv_norew <- map(names(ctrbs), function(z) {
  ggplot(aes(rev(x), y), data = filter(higher_norew_sims_df,
                                       com_vars == z)) + 
    geom_step(aes(group = sims1)) +
    geom_line(color = "firebrick", data = pluck(higher_norew_sims_df_mean,
                                                z),
              size = .8) +
    labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
    theme(aspect.ratio = 1, text = element_text(size = 15)) +
    ggtitle(z)
})

cv_norew[[6]] <- ggplot(aes(rev(x), y), data = higher_org_norew_sims_df) + 
  geom_step(aes(group = sims1)) +
  geom_line(color = "firebrick", data = higher_org_norew_sims_df_mean, size = .8) +
  labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
  theme(aspect.ratio = 1, text = element_text(size = 15)) +
  ggtitle("Original")

cv_norew_higher_all <- ggarrange(cv_norew[[6]], cv_norew[[1]], cv_norew[[2]],
                                cv_norew[[3]], cv_norew[[4]], cv_norew[[5]], nrow  = 2 , ncol = 3)

ggsave("extc_sims_cv_norew_higher.pdf",
       path = paste0(getwd(), "/plot_sink"),
       plot = cv_norew_higher_all,
       device = "pdf",
       width = 28.7,
       height = 20,
       units = "cm")


# org by rew
# lower
org_rew_titles <- c("abund" =  "Abundance", "trait" = "Trait", "phylo" = "Phylogeny")
org_rew <- map(rew_names, function(z) {
  ggplot(aes(rev(x), y),
         data = filter(lower_org_sims_df,
                       id == z)) +
    geom_step(aes(group = sims1)) +
    geom_line(color = "firebrick", data = pluck(lower_org_sims_df_mean,
                                                which(rew_names == z)),
              size = .8) +
    labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
    theme(aspect.ratio = 1, text = element_text(size = 15)) +
    ggtitle(org_rew_titles[z])
})

org_rew[[4]] <- ggplot(aes(rev(x), y), data = lower_org_AT_sims_df) + 
  geom_step(aes(group = sims1)) +
  geom_line(color = "firebrick", data = lower_org_AT_sims_df_mean,
            size = .8) +
  labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
  theme(aspect.ratio = 1, text = element_text(size = 15)) +
  ggtitle("Abundance x Trait")

org_rew[[5]] <- ggplot(aes(rev(x), y), data = lower_org_AP_sims_df) + 
  geom_step(aes(group = sims1)) +
  geom_line(color = "firebrick", data = lower_org_AT_sims_df_mean,
            size = .8) +
  labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
  theme(aspect.ratio = 1, text = element_text(size = 15)) +
  ggtitle("Abundance x Phylogeny")

org_rew[[6]] <- ggplot(aes(rev(x), y), data = lower_org_norew_sims_df) + 
  geom_step(aes(group = sims1)) +
  geom_line(color = "firebrick", data = lower_org_norew_sims_df_mean, size = .8) +
  labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
  theme(aspect.ratio = 1, text = element_text(size = 15)) +
  ggtitle("No rewiring")

org_rew_lower_all <- ggarrange(org_rew[[6]], org_rew[[1]], org_rew[[2]],
                               org_rew[[3]], org_rew[[4]], org_rew[[5]], nrow  = 2 , ncol = 3)

ggsave("extc_sims_org_rew_lower.pdf",
       path = paste0(getwd(), "/plot_sink"),
       plot = org_rew_lower_all,
       device = "pdf",
       width = 28.7,
       height = 20,
       units = "cm")

# higher
org_rew_titles <- c("abund" =  "Abundance", "trait" = "Trait", "phylo" = "Phylogeny")
org_rew <- map(rew_names, function(z) {
  ggplot(aes(rev(x), y),
         data = filter(higher_org_sims_df,
                       id == z)) +
    geom_step(aes(group = sims1)) +
    geom_line(color = "firebrick", data = pluck(higher_org_sims_df_mean,
                                                which(rew_names == z)),
              size = .8) +
    labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
    theme(aspect.ratio = 1, text = element_text(size = 15)) +
    ggtitle(org_rew_titles[z])
})

org_rew[[4]] <- ggplot(aes(rev(x), y), data = higher_org_AT_sims_df) + 
  geom_step(aes(group = sims1)) +
  geom_line(color = "firebrick", data = higher_org_AT_sims_df_mean,
            size = .8) +
  labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
  theme(aspect.ratio = 1, text = element_text(size = 15)) +
  ggtitle("Abundance x Trait")

org_rew[[5]] <- ggplot(aes(rev(x), y), data = higher_org_AP_sims_df) + 
  geom_step(aes(group = sims1)) +
  geom_line(color = "firebrick", data = higher_org_AT_sims_df_mean,
            size = .8) +
  labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
  theme(aspect.ratio = 1, text = element_text(size = 15)) +
  ggtitle("Abundance x Phylogeny")

org_rew[[6]] <- ggplot(aes(rev(x), y), data = higher_org_norew_sims_df) + 
  geom_step(aes(group = sims1)) +
  geom_line(color = "firebrick", data = higher_org_norew_sims_df_mean, size = .8) +
  labs(x = "Primary extinction [%]", y = "Secondary extinction [%]") +
  theme(aspect.ratio = 1, text = element_text(size = 15)) +
  ggtitle("No rewiring")

org_rew_higher_all <- ggarrange(org_rew[[6]], org_rew[[1]], org_rew[[2]],
                               org_rew[[3]], org_rew[[4]], org_rew[[5]], nrow  = 2 , ncol = 3)

ggsave("extc_sims_org_rew_higher.pdf",
       path = paste0(getwd(), "/plot_sink"),
       plot = org_rew_higher_all,
       device = "pdf",
       width = 28.7,
       height = 20,
       units = "cm")


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
  
# mean number of secondary extinctions per step without last two steps
# by com_vars
# lower norew
map(ctrbs, function(x) {
  map(seq(n_webs), function(y) {
    map(seq(n_sims), function(z) {
      mean(extc_sims_lower_norew[[y]][[x]][[z]][1:(nrow(extc_sims_lower_norew[[y]][[x]][[z]]) - 2), 3])
    }) %>% unlist %>% mean
  }) %>% unlist %>% mean
})

# higher norew
map(ctrbs, function(x) {
  map(seq(n_webs), function(y) {
    map(seq(n_sims), function(z) {
      mean(extc_sims_higher_norew[[y]][[x]][[z]][1:(nrow(extc_sims_higher_norew[[y]][[x]][[z]]) - 2), 2])
    }) %>% unlist %>% mean
  }) %>% unlist %>% mean
})

# lower org norew
m_extc_lower_org_norew <- map(seq(n_webs), function(y) {
    map(seq(n_sims), function(z) {
      mean(extc_sims_lower_org_norew[[y]][[z]][1:(nrow(extc_sims_lower_org_norew[[y]][[z]]) - 2), 3])
    }) %>% unlist %>% mean
  }) %>% unlist %>% data.table("vals" = ., "id" = "org")

# higher org norew
map(seq(n_webs), function(y) {
  map(seq(n_sims), function(z) {
    mean(extc_sims_higher_org_norew[[y]][[z]][1:(nrow(extc_sims_higher_org_norew[[y]][[z]]) - 2), 2])
  }) %>% unlist %>% mean
}) %>% unlist %>% data.table("vals" = ., "id" = "org")


# by rew
# lower org
m_extc_lower_org <- map(rew_names, function(x) {
  map(seq(n_webs), function(y) {
    map(seq(n_sims), function(z) {
      mean(extc_sims_lower_org[[y]][[x]][[z]][1:(nrow(extc_sims_lower_org[[y]][[x]][[z]]) - 2), 3])
    }) %>% unlist %>% mean
  }) %>% unlist %>% data.table("vals" = ., "id" = x)
}) %>% bind_rows



# higher org
map(rew_names, function(x) {
  map(seq(n_webs), function(y) {
    map(seq(n_sims), function(z) {
      mean(extc_sims_higher_org[[y]][[x]][[z]][1:(nrow(extc_sims_higher_org[[y]][[x]][[z]]) - 2), 2])
    }) %>% unlist %>% mean
  }) %>% unlist %>% data.table("vals" = ., "id" = x)
}) %>% bind_rows

# lower org AT
map(seq(n_webs), function(y) {
  map(seq(n_sims), function(z) {
    mean(extc_sims_lower_org_AT[[y]][[z]][1:(nrow(extc_sims_lower_org_AT[[y]][[z]]) - 2), 3])
  }) %>% unlist %>% mean
}) %>% unlist %>% mean

# higher org AT
map(seq(n_webs), function(y) {
  map(seq(n_sims), function(z) {
    mean(extc_sims_higher_org_AT[[y]][[z]][1:(nrow(extc_sims_higher_org_AT[[y]][[z]]) - 2), 2])
  }) %>% unlist %>% mean
}) %>% unlist %>% mean

# lower org AP
map(seq(n_webs), function(y) {
  map(seq(n_sims), function(z) {
    mean(extc_sims_lower_org_AP[[y]][[z]][1:(nrow(extc_sims_lower_org_AP[[y]][[z]]) - 2), 3])
  }) %>% unlist %>% mean
}) %>% unlist %>% mean

# higher org AP
map(seq(n_webs), function(y) {
  map(seq(n_sims), function(z) {
    mean(extc_sims_higher_org_AP[[y]][[z]][1:(nrow(extc_sims_higher_org_AP[[y]][[z]]) - 2), 2])
  }) %>% unlist %>% mean
}) %>% unlist %>% mean


## mean number of extinctions in last two steps
# by com_vars
# lower norew
map(ctrbs, function(x) {
  map(seq(n_webs), function(y) {
    map(seq(n_sims), function(z) {
      mean(extc_sims_lower_norew[[y]][[x]][[z]][(nrow(extc_sims_lower_norew[[y]][[x]][[z]]) - 1):nrow(extc_sims_lower_norew[[y]][[x]][[z]]), 3])
    }) %>% unlist %>% mean
  }) %>% unlist %>% mean
})

# higher norew
map(ctrbs, function(x) {
  map(seq(n_webs), function(y) {
    map(seq(n_sims), function(z) {
      mean(extc_sims_higher_norew[[y]][[x]][[z]][(nrow(extc_sims_higher_norew[[y]][[x]][[z]]) - 1):nrow(extc_sims_higher_norew[[y]][[x]][[z]]), 2])
    }) %>% unlist %>% mean
  }) %>% unlist %>% mean
})

# lower org norew
m_extc_lower_org_norew_last <- map(seq(n_webs), function(y) {
  map(seq(n_sims), function(z) {
    mean(extc_sims_lower_org_norew[[y]][[z]][(nrow(extc_sims_lower_org_norew[[y]][[z]]) - 1):nrow(extc_sims_lower_org_norew[[y]][[z]]), 3])
  }) %>% unlist %>% mean
}) %>% unlist %>% data.table("vals" = ., "id" = "norew")

# higher org norew
m_extc_higher_org_norew_last <- map(seq(n_webs), function(y) {
  map(seq(n_sims), function(z) {
    mean(extc_sims_higher_org_norew[[y]][[z]][(nrow(extc_sims_higher_org_norew[[y]][[z]]) - 1):nrow(extc_sims_higher_org_norew[[y]][[z]]), 2])
  }) %>% unlist %>% mean
}) %>% unlist %>% data.table("vals" = ., "id" = "norew")


# by rew
# lower org
m_extc_lower_org_last <- map(rew_names, function(x) {
  map(seq(n_webs), function(y) {
    map(seq(n_sims), function(z) {
      mean(extc_sims_lower_org[[y]][[x]][[z]][(nrow(extc_sims_lower_org[[y]][[x]][[z]]) - 1):nrow(extc_sims_lower_org[[y]][[x]][[z]]), 3])
    }) %>% unlist %>% mean
  }) %>% unlist %>% data.table("vals" = ., "id" = x)
}) %>% bind_rows

ggboxplot("id", "vals", data = m_extc_lower_org_last)

# higher org
m_extc_higher_org_last <- map(rew_names, function(x) {
  map(seq(n_webs), function(y) {
    map(seq(n_sims), function(z) {
      mean(extc_sims_higher_org[[y]][[x]][[z]][(nrow(extc_sims_higher_org[[y]][[x]][[z]]) - 1):nrow(extc_sims_higher_org[[y]][[x]][[z]]), 2])
    }) %>% unlist %>% mean
  }) %>% unlist %>% data.table("vals" = ., "id" = x)
}) %>% bind_rows

# lower org AT
map(seq(n_webs), function(y) {
  map(seq(n_sims), function(z) {
    mean(extc_sims_lower_org_AT[[y]][[z]][(nrow(extc_sims_lower_org_AT[[y]][[z]]) - 1):nrow(extc_sims_lower_org_AT[[y]][[z]]), 3])
  }) %>% unlist %>% mean
}) %>% unlist %>% mean

# higher org AT
map(seq(n_webs), function(y) {
  map(seq(n_sims), function(z) {
    mean(extc_sims_higher_org_AT[[y]][[z]][(nrow(extc_sims_higher_org_AT[[y]][[z]]) - 1):nrow(extc_sims_higher_org_AT[[y]][[z]]), 2])
  }) %>% unlist %>% mean
}) %>% unlist %>% mean

# lower org AP
map(seq(n_webs), function(y) {
  map(seq(n_sims), function(z) {
    mean(extc_sims_lower_org_AP[[y]][[z]][(nrow(extc_sims_lower_org_AP[[y]][[z]]) - 1):nrow(extc_sims_lower_org_AP[[y]][[z]]), 3])
  }) %>% unlist %>% mean
}) %>% unlist %>% mean

# higher org AP
map(seq(n_webs), function(y) {
  map(seq(n_sims), function(z) {
    mean(extc_sims_higher_org_AP[[y]][[z]][(nrow(extc_sims_higher_org_AP[[y]][[z]]) - 1):nrow(extc_sims_higher_org_AP[[y]][[z]]), 2])
  }) %>% unlist %>% mean
}) %>% unlist %>% mean


m_extc_lower_org_all <- bind_rows(m_extc_lower_org, m_extc_lower_org_norew)
m_extc_higher_org_all <- bind_rows(m_extc_higher_org, m_extc_higher_org_norew)

m_extc_lower_org_last_all <- bind_rows(m_extc_lower_org_last, m_extc_lower_org_norew_last)
m_extc_higher_org_last_all <- bind_rows(m_extc_higher_org_last, m_extc_higher_org_norew_last)

ggboxplot("id", "vals", data = m_extc_lower_org_last_all) 


### AUC / Robustness
## ANOVA of auc values
# overall
aov_lo <- aov(robustness ~ nrow + ncol + com_vars + rew * H2, data = auc_all_lower) # lower
aov_hi <- aov(robustness ~ nrow + ncol + com_vars + rew * H2, data = auc_all_higher) # higher

aov_both <- aov(robustness ~ nrow + ncol + com_vars + lvl + rew * H2, data = auc_all) # combined
summary(aov_both)

# use type III test from car package
library(car)
Anova(aov_lo, type = "III")
Anova(aov_hi, type = "III")
# Anova(aov_both, type = "III")

# by com_vars
aov_cv_low <- map(c(names(ctrbs), "org"), ~
      summary(aov(robustness ~ nrow + ncol + rew *H2,
          data = filter(auc_all_lower, com_vars == .x)))) # lower
aov_cv_high <- map(c(names(ctrbs), "org"), ~
      summary(aov(robustness ~ nrow + ncol + rew * H2,
          data = filter(auc_all_higher, com_vars == .x)))) # higher

## extract R2 for each predictor of each scenario
map(1:6, function(x) {# map scenarios
  map(1:6, ~ aov_cv_low[[x]][[1]]$"Sum Sq"[.x] / sum(aov_cv_low[[x]][[1]]$"Sum Sq"))
})

# by rew
aov_rew_low <- map(unique(auc_all_lower$rew), ~
        summary(aov(robustness ~ nrow + ncol + com_vars + H2,
          data = filter(auc_all_lower, rew == .x)))) # lower
aov_rew_high <- map(unique(auc_all_lower$rew), ~
        summary(aov(robustness ~ nrow + ncol + com_vars + H2,
          data = filter(auc_all_higher, rew == .x)))) # higher

## extract R2 for each predictor of each scenario
map(1:6, function(x) {# map scenarios
  map(1:5, ~ aov_rew_low[[x]][[1]]$"Sum Sq"[.x] / sum(aov_rew_low[[x]][[1]]$"Sum Sq")*100)
})

# AUC boxplots
auc_boxplot(auc_all, x.axis = "rew")

# com_vars by norew
auc_all$com_vars <- factor(auc_all$com_vars, level = c("org", "Atl", "aTl",
                                                       "atL", "atl", "ATL"))

auc_cv_by_norew <- auc_boxplot(auc_all, by = "com_vars", select = "norew") 

ggsave("auc_com_vars_by_norew.pdf",
              path = paste0(getwd(), "/plot_sink"),
              plot = auc_cv_by_norew,
              device = "pdf",
              width = 28.7,
              height = 20,
              units = "cm")

# org by rew
auc_all$rew <- factor(auc_all$rew, level = c("norew", "abund", "trait",
                                                       "phylo", "abundtrait", "abundphylo"))

# hardcoded change to create correct color labels & ylim !!
auc_org_by_cv <- auc_boxplot(auc_all, by = "rew", select = "org") 

ggsave("auc_org_by_cv.pdf",
       path = paste0(getwd(), "/plot_sink"),
       plot = auc_org_by_cv,
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
       plot = auc_rew_by_cv,
       device = "pdf",
       width = 28.7,
       height = 20,
       units = "cm")

# by specific combination
auc_all_plots <- map(.x = c(names(ctrbs), "org"), function(x) { # loop com_vars
  map(.x = c(rew_names, "norew", "abundtrait", "abundphylo"), # loop rewiring methods
 ~ auc_boxplot(auc_all, x.axis = c(x, .x), save = T))
  })

auc_abund_pair <- ggboxplot(y = "robustness",
                            color = "lvl", x = "com_vars",
                            data = filter(auc_all, rew == "abund", com_vars %in% c("Atl", "atl", "org"))) +
  scale_color_manual(name = "Trophic level",
                     values = c("#000000", "#787878")) +
  theme(legend.position = "bottom") +
  labs(x = element_blank()) +
  theme(aspect.ratio = 1) +
  ggtitle("abund/Atl")

auc_trait_pair <- ggboxplot(y = "robustness",
                            color = "lvl", x = "com_vars",
                            data = filter(auc_all, rew == "trait", com_vars %in% c("aTl", "atl", "org"))) +
  scale_color_manual(name = "Trophic level",
                     values = c("#000000", "#787878")) +
  theme(legend.position = "bottom") +
  labs(x = element_blank()) +
  theme(aspect.ratio = 1) +
  ggtitle("trait/aTl")

auc_phylo_pair <- ggboxplot(y = "robustness",
                            color = "lvl", x = "com_vars",
                            data = filter(auc_all, rew == "phylo", com_vars %in% c("atL", "atl", "org"))) +
  scale_color_manual(name = "Trophic level",
                     values = c("#000000", "#787878")) +
  theme(legend.position = "bottom") +
  labs(x = element_blank()) +
  theme(aspect.ratio = 1) +
  ggtitle("phylo/atL")

auc_pairs <- ggarrange(auc_abund_pair, auc_trait_pair, auc_phylo_pair, nrow = 1, common.legend = T)


# calculate mean of every combination
auc_means <- map(.x = c(names(ctrbs), "org"), function(x) { # loop com_vars
  map(.x = c(rew_names, "norew", "abundtrait", "abundphylo"), function(y) { # loop rewiring methods
      lo <- filter(auc_all, rew == y, com_vars == x, lvl == "lower") %>% select(robustness)
      hi <- filter(auc_all, rew == y, com_vars == x, lvl == "higher") %>% select(robustness)
      return(c(mean(hi$robustness, na.rm = T), mean(lo$robustness, na.rm = T)))
        })
})

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

