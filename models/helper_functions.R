## Helper functions needed to run full extinction simulation w/ simulated networks

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
                     multiple.webs = F) {
  if (multiple.webs == T) {
   map(.x = c("Atl", "aTl", "atL"),
       ~ replicate(n_sims, simplify = F, 
                   one.second.extinct.mod_aug(web = pluck(web, .x, "web"), 
                                                   participant = participant,
                                                   method = method,
                                                   rewiring = rewiring,
                                                   partner.choice = partner.choice,
                                                   interactions = pluck(interactions, .x, "I_mat"),
                                                   method.rewiring = method.rewiring)))
  } else {
  map(web, ~replicate(n_sims, simplify = F, 
                      one.second.extinct.mod_aug(web = pluck(.x, 1), 
                                                 participant = participant,
                                                 method = method,
                                                 rewiring = rewiring,
                                                 partner.choice = partner.choice,
                                                 interactions = pluck(.x, 2),
                                                 method.rewiring = method.rewiring)))
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
        select(str_which(colnames(.), c("n.lower"))) %>%
        replace_duplicate(.) %>% rowMeans(.)
    } else {
      # mean of n.higher
      out <- pluck(x, y) %>% as.data.table(.) %>% 
        select(str_which(colnames(.), c("n.higher"))) %>%
        replace_duplicate(.) %>% rowMeans(.)
    }  
  } else {
    if (lower == T) {
      # mean of n.lower
      out <- map(c("Atl" = 1, "aTl" = 2, "atL" = 3),
                 ~ pluck(x, y, .x) %>% as.data.table(.) %>% 
                   select(str_which(colnames(.), c("n.lower"))) %>%
                   replace_duplicate(.) %>% rowMeans(.))
    } else {
      # mean of n.higher
      out <- map(c("Atl" = 1, "aTl" = 2, "atL" = 3),
                 ~ pluck(x, y, .x) %>% as.data.table(.) %>% 
                   select(str_which(colnames(.), c("n.higher"))) %>%
                   replace_duplicate(.) %>% rowMeans(.))
    }
  }
  
  return (out)
}

replace_duplicate <- function(x) {
  
  lens <- map(x, ~ which(.x == 0)) %>% unlist(.) # get length of each df
  max_len <- max(lens) # max length of all dfs
  
  # for shorter dfs replace recycled values with 0s
  out <- map2(.x = x, .y = lens, ~ replace(.x, .y:max_len, values = 0))
  
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

plot_extc_multi <- function(x){
  com_vars <- c("Atl" = 1, "aTl" = 2, "atL" = 3)
  map(1:3, ~ ggplot() + 
        geom_line(aes(rev(pluck(x, 1, .x)), pluck(x, 2, .x),
                      color = "Abundance"), linetype = 1) +
        geom_line(aes(rev(pluck(x, 1, .x)), pluck(x, 2, .x),
                      color = "Traits"), linetype = 2) +
        geom_line(aes(rev(pluck(x, 1, .x)), pluck(x, 2, .x),
                      color = "Phylogeny"), linetype = 3) +
        scale_color_manual(name = "Rewiring Method",
                           values = c("Abundance" = "black",
                                      "Traits" = "firebrick",
                                      "Phylogeny" = "dodgerblue")) +
        labs(x = "plants removed", y = "animals persisting",
             title = paste("Extinction cascade", names(com_vars[.x])))
  )
}

plot_extc_alt <- function(x, lower = F){
  com_vars <- c("abund" = 1, "trait" = 2, "phylo" = 3)
  if (lower == T) {
    map(com_vars, ~ ggplot() + 
          geom_line(aes(rev(pluck(x, 1, .x, 1)), pluck(x, 2, .x, 1),
                        color = "Atl"), linetype = 1) +
          geom_line(aes(rev(pluck(x, 1, .x, 2)), pluck(x, 2, .x, 2),
                        color = "aTl"), linetype = 2) +
          geom_line(aes(rev(pluck(x, 1, .x, 3)), pluck(x, 2, .x, 3),
                        color = "atL"), linetype = 3) +
          geom_line(aes(rev(pluck(sp_remain_lower_org, 1, .x)), pluck(sp_remain_lower_org, 2, .x),
                        color = "original"), linetype = 4) +
          scale_color_manual(name = "Rewiring Method",
                             values = c("Atl" = "black",
                                        "aTl" = "firebrick",
                                        "atL" = "dodgerblue",
                                        "original" = "burlywood4")) +
          labs(x = "plants removed", y = "animals persisting",
               title = paste("Extinction cascade", names(com_vars[.x])))
    )
  } else {
    map(com_vars, ~ ggplot() + 
          geom_line(aes(rev(pluck(x, 1, .x, 1)), pluck(x, 2, .x, 1),
                        color = "Atl"), linetype = 1) +
          geom_line(aes(rev(pluck(x, 1, .x, 2)), pluck(x, 2, .x, 2),
                        color = "aTl"), linetype = 2) +
          geom_line(aes(rev(pluck(x, 1, .x, 3)), pluck(x, 2, .x, 3),
                        color = "atL"), linetype = 3) +
          geom_line(aes(rev(pluck(sp_remain_higher_org, 1, .x)), pluck(sp_remain_higher_org, 2, .x),
                        color = "original"), linetype = 4) +
          scale_color_manual(name = "Rewiring Method",
                             values = c("Atl" = "black",
                                        "aTl" = "firebrick",
                                        "atL" = "dodgerblue",
                                        "original" = "burlywood4")) +
          labs(x = "plants removed", y = "animals persisting",
               title = paste("Extinction cascade", names(com_vars[.x])))
    ) 
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
