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
  
  #out <- match_lengths(out, df = F)
  
  return(out)
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
match_lengths <- function(x, df = T) {
  end <- map(x, ~ which(.x == 0) %>% first()) # get indices of first zero (i.e. last non replicated entry)
  max <- unlist(end) %>% max(.) # max length
  out <- map2(x, end,  ~ replace(.x, .y:max, values = 0)) # replace recycled vals
  if(df)
    out <- as.data.table(out)
  return(out)
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
library(ggpubr)

plot_extc_facet <- function(x, ci, save = F, view = T, norew = F, org = F) {
  # extinction on lower or higher level ?
  lvl_match <- grepl("lower", substitute(x))
  ifelse(lvl_match, lvl <- "lower", lvl <- "higher")
  
  # automate axis labels
  ifelse(lvl == "lower", xlab <- "plants", xlab <- "animals")
  ifelse(lvl == "lower", ylab <- "animals", ylab <- "plants")
  
  # colors for facet boundary boxes
  colors <- c("abund" = "black", "trait" = "firebrick", "phylo" = "dodgerblue")
  
  # create df for ggplot
  df <- cbind(ci,
             select(x, -c("id", "com_vars"))) %>% 
    tidyr::unite(., col = "group", c("id", "com_vars"), remove = F)
  
  
  # create list of labels for facet
  # norew or org simulations used ?
  if(norew | org) {
    if(norew) {
  # set grouping variable to factor
  df$group <- factor(df$group, levels = c("Atl...1",
                                          "aTl...1",
                                          "atL...1"))
  
      id_labs <- c("Atl", "aTl", "atL") 
      names(id_labs) <- c(unique(df$id))

    } else {
      # set grouping variable to factor
      df$group <- factor(df$group, levels = c("abund_...1",
                                              "trait_...1",
                                              "phylo_...1"))
      id_labs <- c("Abundance", "Traits", "Phylogeny")
      names(id_labs) <- c(unique(df$id))
    }
    # set facet layout
    nrow <- 1
    ncol <- 3
    
  } else {
    # set grouping variable to factor
    df$group <- factor(df$group, levels = c("abund_Atl", "abund_aTl", "abund_atL",
                                            "trait_Atl", "trait_aTl", "trait_atL",
                                            "phylo_Atl", "pyhlo_aTl", "phylo_atL"))
    
    id_labs <- map(unique(df$id), function(x) {
      id_labs <- c("Atl", "aTl", "atL") %>%
        set_names(c(paste(x, "Atl", sep = "_"),
                    paste(x, "aTl", sep = "_"),
                    paste(x, "atL", sep = "_")))
    })
    names(id_labs) <- c(unique(df$id))
    
    # set facet layout
    nrow <- 3
    ncol <- 1
  }
  
  # create plots
  plots <- map(unique(df$id), function(z) {
    if(is.list(id_labs)) {
      labs <- pluck(id_labs, z)
    } else {
      labs <- id_labs[z]
      df$title <- labs
    }
    
    if (norew | org) {
      ggplot(filter(df, id == z)) +
        geom_line(aes(rev(x), y, group = group, color = "firebrick")) + # means
        geom_line(aes(rev(x_lower), x_higher, group = group),
                  linetype = 2) + # lower ci
        geom_line(aes(rev(y_lower), y_higher, group = group),
                  linetype = 2) + # higher ci
        labs(x = paste(xlab, "removed"), y = paste(ylab, "persisting")) +
        guides(color = "none") +
        facet_wrap(. ~ title)  +
        theme(strip.text = element_text(),
              strip.background = element_rect(color = colors[z])) # set bb color  
    } else {
      ggplot(filter(df, id == z)) +
        geom_line(aes(rev(x), y, group = group, color = "firebrick")) + # means
        geom_line(aes(rev(x_lower), x_higher, group = group),
                  linetype = 2) + # lower ci
        geom_line(aes(rev(y_lower), y_higher, group = group),
                  linetype = 2) + # higher ci
        labs(x = paste(xlab, "removed"), y = paste(ylab, "persisting")) +
        guides(color = "none") +
        facet_wrap(~ group, labeller = labeller(group = labs))  +
        theme(strip.text = element_text(),
              strip.background = element_rect(color = colors[z])) # set bb color
    }
    })
  
  
  if (view == T) {
    out <- ggarrange(plots[[1]], plots[[2]], plots[[3]],
                     nrow = nrow, ncol = ncol)
    print(out)
  }
  
  if (save == T) {
    ggsave(
      paste("extinction cascade", lvl, "trophic level"),
      path = paste0(getwd(), "/plot_sink"),
      plot = out,
      device = "pdf",
      width = 1900,
      height = 1205,
      units = "px")
  }
}

# DEPRECATED use plot_extc_facet instead
plot_extc <- function(x, ci, save = F, view = T, lower = F){
  com_vars <- c("Atl" = 1, "aTl" = 2, "atL" = 3)
  ifelse(lower == T, lvl <- "lower", lvl <- "higher")
  ifelse(lower == T, xlab <- "plants", xlab <- "animals")
  ifelse(lower == T, ylab <- "animals", ylab <- "plants")

  plots <- map(com_vars, ~ ggplot() +
                 geom_line(aes(rev(pluck(ci, "lower", 1, 1, .x)),
                               pluck(ci, "higher", 1, 1, .x),
                               color = "Abundance"), linetype = 3, size = .3) + # lower
                 geom_line(aes(rev(pluck(ci, "lower", 2, 1, .x)),
                               pluck(ci, "higher", 2, 1, .x),
                               color = "Abundance"), linetype = 3, size = .3) + # upper
                 geom_line(aes(rev(pluck(x, "lower", "abund", .x)),
                               pluck(x, "higher", "abund", .x),
                               color = "Abundance"), linetype = 1, size = .5) + # mean
                 geom_line(aes(rev(pluck(ci, "lower", 1, 2, .x)),
                               pluck(ci, "higher", 1, 2, .x),
                               color = "Traits"), linetype = 2, size = .3) + # lower
                 geom_line(aes(rev(pluck(ci, "lower", 2, 2, .x)),
                               pluck(ci, "higher", 2, 2, .x),
                               color = "Traits"), linetype = 2, size = .3) + # upper
                 geom_line(aes(rev(pluck(x, "lower", "trait", .x)),
                               pluck(x, "higher", "trait", .x),
                               color = "Traits"), linetype = 1, size = .5) + # mean
                 geom_line(aes(rev(pluck(ci, "lower", 1, 3, .x)),
                               pluck(ci, "higher", 1, 3, .x),
                               color = "Phylogeny"), linetype = 3, size = .3) + # lower
                 geom_line(aes(rev(pluck(ci, "lower", 2, 3, .x)),
                               pluck(ci, "higher", 2, 3, .x),
                               color = "Phylogeny"), linetype = 3, size = .3) + # upper
                 geom_line(aes(rev(pluck(x,"lower", "phylo", .x)),
                               pluck(x, "higher", "phylo", .x),
                               color = "Phylogeny"), linetype = 1, size = .5) + # mean
                 scale_color_manual(name = "Rewiring Method",
                                    values = c("Abundance" = "black",
                                               "Traits" = "firebrick",
                                               "Phylogeny" = "dodgerblue")) +
                 labs(x = paste(xlab, "removed"), y = paste(ylab, "persisting"),
                      title = paste("Extinction cascade", names(ctrbs)[.x]))
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
  
  return(plots)
}

plot_extc_alt <- function(x, org, norew, ci, ci_org, ci_norew, save = F, view = T, lower = F){
  com_vars <- c("abund" = 1, "trait" = 2, "phylo" = 3)
  ifelse(lower == T, lvl <- "lower", lvl <- "higher")
  ifelse(lower == T, xlab <- "plants", xlab <- "animals")
  ifelse(lower == T, ylab <- "animals", ylab <- "plants")
  
   plots <- map(com_vars, ~ ggplot() + 
          geom_line(aes(rev(pluck(ci, "lower", 1, .x, 1)),
                        pluck(ci, "higher", 1, .x, 1),
                        color = "Atl"), linetype = 2, size = .3) + # lower
          geom_line(aes(rev(pluck(ci, "lower", 2, .x, 1)),
                        pluck(ci, "higher", 2, .x, 1),
                        color = "Atl"), linetype = 2, size = .3) + # upper
          geom_line(aes(rev(pluck(x, "lower", .x, "Atl")), pluck(x, "higher", .x, "Atl"),
                        color = "Atl"), linetype = 1, size = .5) + # mean
          geom_line(aes(rev(pluck(ci, "lower", 1, .x, 2)),
                        pluck(ci, "higher", 1, .x, 2),
                        color = "aTl"), linetype = 2, size = .3) + # lower
          geom_line(aes(rev(pluck(ci, "lower", 2, .x, 2)),
                        pluck(ci, "higher", 2, .x, 2),
                        color = "aTl"), linetype = 2, size = .3) + # upper
          geom_line(aes(rev(pluck(x, "lower", .x, "aTl")), pluck(x, "higher", .x, "aTl"),
                        color = "aTl"), linetype = 1, size = .5) + # mean
          geom_line(aes(rev(pluck(ci, "lower", 1, .x, 2)),
                        pluck(ci, "higher", 1, .x, 2),
                        color = "atL"), linetype = 2, size = .3) + # lower
          geom_line(aes(rev(pluck(ci, "lower", 2, .x, 2)),
                        pluck(ci, "higher", 2, .x, 2),
                        color = "atL"), linetype = 2, size = .3) + # upper
          geom_line(aes(rev(pluck(x, "lower", .x, "atL")), pluck(x, "higher", .x, "atL"),
                        color = "atL"), linetype = 1, size = .5) + # mean
          geom_line(aes(rev(pluck(ci_org, "lower", 1, .x)),
                        pluck(ci_org, "higher", 1, .x),
                        color = "original"), linetype = 2, size = .3) + # lower
          geom_line(aes(rev(pluck(ci_org, "lower", 2, .x)),
                        pluck(ci_org, "higher", 2, .x),
                        color = "original"), linetype = 2, size = .3) + # upper
          geom_line(aes(rev(pluck(org, "lower", .x)), pluck(org, "higher", .x),
                      color = "original"), linetype = 1, size = .5) + # mean
          geom_line(aes(rev(pluck(ci_norew, "lower", 1, .x)),
                        pluck(ci_norew, "higher", 1, .x),
                        color = "w/o rewiring"), linetype = 2, size = .3) + # lower
          geom_line(aes(rev(pluck(ci_norew, "lower", 2, .x)),
                        pluck(ci_norew, "higher", 2, .x),
                        color = "w/o rewiring"), linetype = 2, size = .3) + # upper
          geom_line(aes(rev(pluck(norew, "lower", .x)), pluck(norew, "higher", .x),
                      color = "w/o rewiring"), linetype = 1, size = .5) + # mean
          scale_color_manual(name = "Contribution",
                             values = c("Atl" = "black",
                                        "aTl" = "firebrick",
                                        "atL" = "dodgerblue",
                                        "original" = "burlywood4",
                                        "w/o rewiring" = "seagreen")) +
          labs(x = paste(xlab, "removed"), y = paste(ylab, "persisting"),
               title = paste("Extinction cascade", lvl, "trophic level", names(com_vars[.x]), "rewiring"))
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

conf_int <- function(x, means, trph_lvl = "lower", lower = T, norew_org = F) {
  if(!(trph_lvl == "lower" | trph_lvl == "higher"))
    stop("Invalid value for trph_lvl ! Specify valid trophic level (lower or higher)")

  # calculate std devs
  if(norew_org) {
    sds <- map(.x = 1:3, function(a) {
      map(seq(n_webs), ~ pluck(x, .x, trph_lvl, a)) %>% as.data.table %>% 
        match_lengths() %>% apply(., 1, sd)})
  } else {
    sds <- map(.x = 1:3, function(a) map(.x = 1:3, function(b) {
      map(seq(n_webs), ~ pluck(x, .x, trph_lvl, a, b)) %>% as.data.table %>% 
        match_lengths() %>% apply(., 1, sd)}))
  }
  
  # calulate margins 
  if(norew_org) {
    margins <- map(1:3, ~ qt(0.975, df = length(pluck(sds, .x)) - 1) * pluck(sds, .x)/sqrt(length(pluck(sds, .x))))
  } else {
    margins <- map(1:3, function(a) {
      map(1:3, ~ qt(0.975, df = length(pluck(sds, a, .x)) - 1) * pluck(sds, a, .x)/sqrt(length(pluck(sds, a, .x))))})
  }
  
  # calculate lower/upper ci
  if (norew_org) {
    if (lower) {
      # out <- map(1:3, function(a) map(1:3, ~ (pluck(means, trph_lvl, a, .x) - pluck(margins, a, .x)) ))
      out <- map(1:3, ~ (pluck(means, trph_lvl, .x) - pluck(margins, .x)) %>% ifelse(. < 0, 0, .))
    } else {
      # out <- map(1:3, function(a) map(1:3, ~ (pluck(means, trph_lvl, a, .x) + pluck(margins, a, .x))))
      out <- map(1:3, ~ (pluck(means, trph_lvl, .x) + pluck(margins, .x)) %>% ifelse(. > 100, 100, .))
    }
  } else {
    if (lower) {
    # out <- map(1:3, function(a) map(1:3, ~ (pluck(means, trph_lvl, a, .x) - pluck(margins, a, .x)) ))
    out <- map(1:3, function(a) map(1:3, ~ (pluck(means, trph_lvl, a, .x) - pluck(margins, a, .x)) %>% ifelse(. < 0, 0, .)))
    } else {
    # out <- map(1:3, function(a) map(1:3, ~ (pluck(means, trph_lvl, a, .x) + pluck(margins, a, .x))))
    out <- map(1:3, function(a) map(1:3, ~ (pluck(means, trph_lvl, a, .x) + pluck(margins, a, .x)) %>% ifelse(. > 100, 100, .)))
    }
  }
  

  return(out)
}

connectance <- function(x) {
  sum(x != 0) / (nrow(x) * ncol(x))
}

# changing lists to df for easier plotting
list_to_df <- function(x, org = F, norew = F, ci = F) {
  ifelse(org|norew, len <- 2, len <- 3) # match list nestedness
  ifelse(org, id <- "org", id <- "norew")
  
  if (ci) {
    if(org|norew) {
      # calculate max lengths of each extc seq, 2.5% ci
      maxlens_loci <- map(1:3, ~ length(pluck(x, "lower", 1, .x))) %>%
        unlist() %>% max()
      
      # 95% ci
      maxlens_hici <- map(1:3, ~ length(pluck(x, "lower", 2, .x))) %>%
        unlist() %>% max()
      
      # calculate lengths of each extc seq, 2.5 % ci
      lens_loci <- map(1:3, ~ length(pluck(x, "lower", 1, .x))) %>% unlist()
      
      # 95% ci
      lens_hici <- map(1:3, ~ length(pluck(x, "lower", 2, .x))) %>% unlist()
      
      # match lengths of lower extc seqs, 2.5 % ci
      len_matched_x_lower_loci <- map(1:3, ~ append(pluck(x, "lower", 1, .x),
                          rep(0, maxlens_loci - lens_loci[.x])))
      
      # 95 % ci
      len_matched_x_lower_hici <- map(1:3, ~ append(pluck(x, "lower", 2, .x),
                          rep(0, maxlens_hici - lens_hici[.x])))
      
      # match length of higher extc seqs, 2.5% ci
      len_matched_x_higher_loci <- map(1:3, ~ append(pluck(x, "higher", 1, .x),
                          rep(0, maxlens_loci - lens_loci[.x])))
      
      # 95% ci
      len_matched_x_higher_hici <- map(1:3, ~ append(pluck(x, "higher", 2, .x),
                          rep(0, maxlens_hici - lens_hici[.x])))
    } else {
      # calculate max lengths of each extc seq, 2.5% ci
      maxlens_loci <- map(1:3, function(y) {
        map(1:len, ~ length(pluck(x, "lower", 1, y, .x))) %>% unlist() %>% max()
      })
      
      # 95% ci
      maxlens_hici <- map(1:3, function(y) {
        map(1:len, ~ length(pluck(x, "lower", 2, y, .x))) %>% unlist() %>% max()
      })
      
      # calculate lengths of each extc seq, 2.5 % ci
      lens_loci <- map(1:3, function(y) {
        map(1:len, ~ length(pluck(x, "lower", 1, y, .x))) %>% unlist()
      })
      
      # 95% ci
      lens_hici <- map(1:3, function(y) {
        map(1:len, ~ length(pluck(x, "lower", 2, y, .x))) %>% unlist()
      })
      
      # match lengths of lower extc seqs, 2.5 % ci
      len_matched_x_lower_loci <- map(1:3, function(y) {
        map(1:len, ~ append(pluck(x, "lower", 1, y, .x),
                          rep(0, maxlens_loci[[y]] - lens_loci[[y]][.x])))
      })
      # 95 % ci
      len_matched_x_lower_hici <- map(1:3, function(y) {
        map(1:len, ~ append(pluck(x, "lower", 2, y, .x),
                          rep(0, maxlens_hici[[y]] - lens_hici[[y]][.x])))
      })
      
      # match length of higher extc seqs, 2.5% ci
      len_matched_x_higher_loci <- map(1:3, function(y) {
        map(1:len, ~ append(pluck(x, "higher", 1, y, .x),
                          rep(0, maxlens_loci[[y]] - lens_loci[[y]][.x])))
      })
      
      # 95% ci
      len_matched_x_higher_hici <- map(1:3, function(y) {
        map(1:len, ~ append(pluck(x, "higher", 2, y, .x),
                          rep(0, maxlens_hici[[y]] - lens_hici[[y]][.x])))
      })
    }
    
    lower <- map(1:3, function(y) {
      cbind(bind_cols(pluck(len_matched_x_lower_loci, y),
                      .id = c(rew_names[y])) %>%
              melt(., id.vars = ".id",
                   value.name = "x_lower",
                   variable.name = "com_vars"),
            bind_cols(pluck(len_matched_x_lower_hici, y),
                      .id = c(rew_names[y])) %>%
              melt(., id.vars = ".id",
                   value.name = "y_lower",
                   variable.name = "com_vars")) %>% 
        select(., -c(4,5)) %>% 
        setNames(., c("id", "com_vars", "x_lower", "y_lower"))}) %>% 
      bind_rows()
    
    higher <- map(1:3, function(y) {
      cbind(bind_cols(pluck(len_matched_x_higher_loci, y),
                          .id = c(rew_names[y])) %>%
                            melt(., id.vars = ".id",
                                 value.name = "x_higher",
                                 variable.name = "com_vars"),
                bind_cols(pluck(len_matched_x_higher_hici, y),
                          .id = c(rew_names[y])) %>%
                  melt(., id.vars = ".id",
                       value.name = "y_higher",
                       variable.name = "com_vars")) %>% 
        select(., -c(4,5)) %>% 
        setNames(., c("id", "com_vars", "x_higher", "y_higher"))}) %>% 
      bind_rows()
    
    out <- cbind(lower, select(higher, -c(1, 2)))
    
  } else {
    if (org|norew) {
        out <- cbind(bind_cols(pluck(x, "lower"), .id = id) %>%
                       melt(., id.vars = ".id",
                            value.name = "x",
                            variable.name = "com_vars"),
                  bind_cols(pluck(x, "higher"), .id = id) %>%
                    melt(., id.vars = ".id",
                         value.name = "y",
                         variable.name = "com_vars")) %>%
          select(., -c(4,5)) %>% 
          setNames(., c("id", "com_vars", "x", "y"))
    } else {
      # calculate max lengths of each extc seq
      maxlens <- map(1:3, function(y) {
        map(1:len, ~ length(pluck(x, "lower", y, .x))) %>% unlist() %>% max()
      })
      
      # calculate lengths of each extc seq
      lens <- map(1:3, function(y) {
        map(1:len, ~ length(pluck(x, "lower", y, .x))) %>% unlist()
      })
      
      # match lengths of lower extc seqs
      len_matched_x_lower <- map(1:3, function(y) {
          map(1:len, ~ append(pluck(x, "lower", y, .x),
                            rep(0, maxlens[[y]] - lens[[y]][.x])))
      })
      
      # match length of higher extc seqs
      len_matched_x_higher <- map(1:3, function(y) {
        map(1:len, ~ append(pluck(x, "higher", y, .x),
                          rep(0, maxlens[[y]] - lens[[y]][.x])))
      })
      
      out <- map(1:3, function(y) {
        cbind(bind_cols(pluck(len_matched_x_lower, y), .id = rew_names[y]) %>%
                    melt(., id.vars = ".id",
                         value.name = "x",
                         variable.name = "com_vars"),
                  bind_cols(pluck(len_matched_x_higher, y), .id = rew_names[y]) %>%
          melt(., id.vars = ".id",
               value.name = "y",
               variable.name = "com_vars"))}) %>% bind_rows(.) %>%
        select(., -c(4,5)) %>% 
        setNames(., c("id", "com_vars", "x", "y"))
    }
  }
  
  # set correct levels for id & com_vars
  levels(out$com_vars) <- c("Atl", "aTl", "atL")
  levels(out$id) <- c("abund", "trait", "phylo")
  
  return(out)
}

get_ci <- function(x, means, org = F, norew = F) {
  # set names
  if (org) {
    out <- list("lower" = list("2.5" = conf_int(x = x,
                                         mean = means,
                                         norew_org = T),
                        "97.5" = conf_int(x = x,
                                          mean = means,
                                          lower = F, norew_org = T)),
         "higher" = list("2.5" = conf_int(x = x,
                                          mean = means,
                                          trph_lvl = "higher", norew_org = T),
                         "97.5" = conf_int(x = x,
                                           mean = means,
                                           trph_lvl = "higher",
                                           lower = F, norew_org = T)))
  
  out <- modify_depth(out, 2, 
                      ~ set_names(.x, nm = c("abund", "trait", "phylo")))
  }
  
  if (norew) {
    out <- list("lower" = list("2.5" = conf_int(x = x,
                                         mean = means,
                                         norew_org = T),
                        "97.5" = conf_int(x = x,
                                          mean = means,
                                          lower = F, norew_org = T)),
         "higher" = list("2.5" = conf_int(x = x,
                                          mean = means,
                                          trph_lvl = "higher", norew_org = T),
                         "97.5" = conf_int(x = x,
                                           mean = means,
                                           trph_lvl = "higher",
                                           lower = F, norew_org = T)))
    
    out <- modify_depth(out, 2, 
                        ~ set_names(.x, nm = c("Atl", "aTl", "atL")))
  }
  
  if (!org & !norew) {
    out <- list("lower" = list("2.5" = conf_int(x = x,
                                         mean = means),
                        "97.5" = conf_int(x = x,
                                          mean = means,
                                          lower = F)),
         "higher" = list("2.5" = conf_int(x = x,
                                          mean = means,
                                          trph_lvl = "higher"),
                         "97.5" = conf_int(x = x,
                                           mean = means,
                                           trph_lvl = "higher",
                                           lower = F)))
    out <- modify_depth(out, 2, 
                              ~ set_names(.x, nm = c("abund", "trait", "phylo")))
    
    out <- modify_depth(out, 3, 
                              ~ set_names(.x, nm = c("Atl", "aTl", "atL")))
  }
  
  return(out)
}
