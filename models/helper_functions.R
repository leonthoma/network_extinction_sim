## Helper functions needed to run full extinction simulation w/ simulated networks

# fy to simulate extinctions for a nested list of networks; n_sims is no of extc
# simulations, n_nets is no of networks in list, n_webs is no of webs per network
run_extc <- function(web,
                     participant = "lower",
                     method = "random",
                     rewiring = F,
                     partner.choice = NULL,
                     interactions = NULL,
                     method.rewiring = "NULL",
                     n_sims = 0,
                     multiple.webs = F,
                     make.bipartite = F
                     ) {
  if (multiple.webs == T) {
   map(.x = c("Atl", "aTl", "atL"),
       ~ replicate(n_sims, simplify = F, 
                   one.second.extinct.mod_aug(web = pluck(web, .x, "web"), 
                                                   participant = participant,
                                                   method = method,
                                                   rewiring = rewiring,
                                                   partner.choice = partner.choice,
                                                   interactions = pluck(interactions, .x, "I_mat"),
                                                   method.rewiring = method.rewiring,
                                                   make.bipartite = make.bipartite)))
  } else {
  map(web, ~replicate(n_sims, simplify = F, 
                      one.second.extinct.mod_aug(web = pluck(.x, 1), 
                                                 participant = participant,
                                                 method = method,
                                                 rewiring = rewiring,
                                                 partner.choice = partner.choice,
                                                 interactions = pluck(.x, 2),
                                                 method.rewiring = method.rewiring,
                                                 make.bipartite = make.bipartite)))
  }
}

# get all simulations from a network, match their lengths by adding zeros to
# shorter dfs and compute the mean 
library(stringr)
library(dplyr)
library(data.table)

list_mean <- function(x, y, lower = T, original = F) {
  
  if (original == T) {
    if (lower == T) {
      # mean of n.lower
      out <- pluck(x, y) %>% as.data.table(.) %>% 
        select(., contains(c("no", "n.lower"))) %>%
        replace_duplicate(.) %>% rowMeans(.)
    } else {
      # mean of n.higher
      out <- pluck(x, y) %>% as.data.table(.) %>% 
        select(., contains(c("no", "n.higher"))) %>%
        replace_duplicate(.) %>% rowMeans(.)
    }  
  } else {
    if (lower == T) {
      # mean of n.lower
      out <- map(c("Atl" = 1, "aTl" = 2, "atL" = 3),
                 ~ pluck(x, y, .x) %>% as.data.table(.) %>% 
                   select(., contains(c("no", "n.lower"))) %>%
                   replace_duplicate(.) %>% rowMeans(.))
    } else {
      # mean of n.higher
      out <- map(c("Atl" = 1, "aTl" = 2, "atL" = 3),
                 ~ pluck(x, y, .x) %>% as.data.table(.) %>% 
                   select(., contains(c("no", "n.higher"))) %>%
                   replace_duplicate(.) %>% rowMeans(.))
    }
  }
  
  return (out)
}

replace_duplicate <- function(x) {
  # get length of each df; add 1 to match no with actual length
  lens <- select(x, contains("no")) %>%
    map_dbl(., ~ max(.x) + 1)  
  max_len <- max(lens) # max length of all dfs
  
  # for shorter dfs replace recycled values with 0s
  out <- map2(.x = select(x, contains("n.")), .y = lens, ~ replace(.x, .y:max_len, values = 0))
  
  return(as.data.table(out) %>% select(., contains("n.")))
}

# calculating mean over every web
web_mean <- function(x) {
  end <- map(x, ~ which(.x == 0)) # get indices of last entries
  max <- unlist(end) %>% max(.) # max length
  out <- map2(x, end,  ~ replace(.x, .y:max, values = 0)) # replace recycled vals
  return(as.data.table(out))
}

# function to calculate percentages of remaining species
per_surv <- function(x, y, lower = T, original = F) {
  list_divide <- function(x) {
    x / x[1] * 100
  }
  
  extc_in_network <- ifelse(isTRUE(lower), 1, 2)
  
  if (original == T) {
    list_divide(pluck(x, extc_in_network, y))
  } else {
    map(c("Atl" = 1, "aTl" = 2, "atL" = 3),
        ~ list_divide(pluck(x, extc_in_network, y, .x)))
  }
}

# functions to plot mean persistance of sp from one trophic level by removal of
# sp from the other trophic level in percent
library(ggplot2)

plot_extc <- function(x, save = F, view = T, lower = F){
  com_vars <- c("Atl" = 1, "aTl" = 2, "atL" = 3)
  ifelse(lower == T, lvl <- "lower", lvl <- "higher")
  
  plots <- map(com_vars, ~ ggplot() + 
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
             title = paste("Extinction cascade", lvl, "trophic level", names(com_vars[.x])))
  )
  
  if (view == T)
    map(com_vars, ~ plot(plots[[.x]]))
  
  if (save == T) {
    map(c("Atl" = 1, "aTl" = 2, "atL" = 3), ~ ggsave(
      paste("extinction cascade", lvl, "trophic level",
            names(com_vars[.x])),
      path = paste0(getwd(), "/plot_sink"),
      plot = plots[[.x]],
      device = "pdf",
      width = 1900,
      height = 1205,
      units = "px"))
  }
}

plot_extc_alt <- function(x, y, z, save = F, view = T, lower = F){
  com_vars <- c("abund" = 1, "trait" = 2, "phylo" = 3)
  ifelse(lower == T, lvl <- "lower", lvl <- "higher")
  
   plots <- map(com_vars, ~ ggplot() + 
          geom_line(aes(rev(pluck(x, 1, .x, 1)), pluck(x, 2, .x, 1),
                        color = "Atl"), linetype = 1) +
          geom_line(aes(rev(pluck(x, 1, .x, 2)), pluck(x, 2, .x, 2),
                        color = "aTl"), linetype = 2) +
          geom_line(aes(rev(pluck(x, 1, .x, 3)), pluck(x, 2, .x, 3),
                        color = "atL"), linetype = 3) +
          geom_line(aes(rev(pluck(y, 1, .x)), pluck(y, 2, .x),
                        color = "original"), linetype = 4) +
          geom_line(aes(rev(pluck(z, 1, .x)), pluck(z, 2, .x),
                        color = "w/o rewiring"), linetype = 5) +
          scale_color_manual(name = "Rewiring Method",
                             values = c("Atl" = "black",
                                        "aTl" = "firebrick",
                                        "atL" = "dodgerblue",
                                        "original" = "burlywood4",
                                        "w/o rewiring" = "seagreen")) +
          labs(x = "plants removed", y = "animals persisting",
               title = paste("Extinction cascade", lvl, "trophic level", names(com_vars[.x])))
    )
   
   if (view == T)
     map(com_vars, ~ plot(plots[[.x]]))
   
   if (save == T) {
     map(c("abund" = 1, "trait" = 2, "phylo" = 3), ~ ggsave(
       paste("alt extinction cascade", lvl, "trophic level",
             names(com_vars[.x])),
       path = paste0(getwd(), "/plot_sink"),
       plot = plots[[.x]],
       device = "pdf",
       width = 1900,
       height = 1205,
       units = "px"))
  }
}

# fy for simulating webs with user specified community var contributions
simnetfromtap_ctrb <- function(traits,
                               abuns,
                               paramsList,
                               pems,
                               tmatch_type_pem,
                               tmatch_type_obs,
                               ctrb_vec = c(NULL),
                               initial_sim = FALSE,
                               Nwebs,
                               Nobs) {
  # create networks
  sim <- simnetfromtap_aug(traits = traits,
                           abuns = abuns,
                           paramsList = paramsList,
                           pems = pems,
                           tmatch_type_pem = tmatch_type_pem,
                           tmatch_type_obs = tmatch_type_obs,
                           ctrb_vec = ctrb_vec,
                           initial_sim = initial_sim)
  # get web from Imat
  simweb <- matrix(rmultinom(1, Nobs, sim$I_mat), nrow = nrow(sim$I_mat),
                   ncol = ncol(sim$I_mat))
  # set sp names
  dimnames(simweb) <- dimnames(sim$I_mat)
  
  return(list("web" = simweb,
              "I_mat" = sim$I_mat))
}

# fy to calculate area under the curve
auc <- function(x, y) {
  sum(diff(x) * head(y, -1) + tail(y, -1))/2
  }

# calculate network robustness; adapted to work with augmented extinction fy
robustness_aug <- function (object, lower = T) 
{
  N <- colSums(object)
  if (isTRUE(lower)) {
    y <- -object[, "ext.higher"]
  }
  else {
    y <- -object[, "ext.lower"]
  }
  y <- (sum(y) - cumsum(y))/sum(y)
  x <- (object[, "no"]/max(object[, "no"]))
  ext.curve <- splinefun(x, y)
  ext.area <- integrate(ext.curve, 0, 1)
  return(as.numeric(ext.area[[1]]))
}
