## Helper functions needed to run full extinction simulation w/ simulated networks

# fy to simulate extinctions for a nested list of networks; n_sims is no of extc
# simulations, n_nets is no of networks in list, n_webs is no of webs per network
run_extc <- function(web,
                     participant = "lower",
                     method = "random",
                     rewiring = F,
                     abund.partner.choice = NULL,
                     trait.partner.choice = NULL,
                     phylo.partner.choice = NULL,
                     interactions = NULL,
                     method.rewiring = "NULL",
                     n_sims = 0,
                     multiple.webs = F,
                     make.bipartite = F,
                     coextc.thr = NULL
                     ) {
  if (multiple.webs == T) {
   map(.x = seq(ctrbs),
       ~ replicate(n_sims, simplify = F, 
                   one.second.extinct.mod.aug(web = pluck(web, .x, "web"), 
                                                   participant = participant,
                                                   method = method,
                                                   rewiring = rewiring,
                                                   abund.partner.choice = abund.partner.choice,
                                                   trait.partner.choice = trait.partner.choice,
                                                   phylo.partner.choice = phylo.partner.choice,
                                                   interactions = pluck(interactions, .x, "I_mat"),
                                                   method.rewiring = method.rewiring,
                                                   make.bipartite = make.bipartite,
                                                   coextc.thr = coextc.thr)))
  } else {
  map(web, ~replicate(n_sims, simplify = F, 
                      one.second.extinct.mod.aug(web = pluck(.x, 1), 
                                                 participant = participant,
                                                 method = method,
                                                 rewiring = rewiring,
                                                 abund.partner.choice = abund.partner.choice,
                                                 trait.partner.choice = trait.partner.choice,
                                                 phylo.partner.choice = phylo.partner.choice,
                                                 interactions = pluck(.x, 2),
                                                 method.rewiring = method.rewiring,
                                                 make.bipartite = make.bipartite,
                                                 coextc.thr = coextc.thr)))
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
        replace_duplicate(.) %>% rowMeans(., na.rm = T)
    } else {
      # mean of n.higher
      out <- pluck(x, y) %>% as.data.table(.) %>% 
        select(., contains(c("no", "n.higher"))) %>%
        replace_duplicate(.) %>% rowMeans(., na.rm = T)
    }  
  } else {
    if (lower == T) {
      # mean of n.lower
      out <- map(c("Atl" = 1, "aTl" = 2, "atL" = 3, "atl" = 4),
                 ~ pluck(x, y, .x) %>% as.data.table(.) %>% 
                   select(., contains(c("no", "n.lower"))) %>%
                   replace_duplicate(.) %>% rowMeans(., na.rm = T))
    } else {
      # mean of n.higher
      out <- map(c("Atl" = 1, "aTl" = 2, "atL" = 3, "atl" = 4),
                 ~ pluck(x, y, .x) %>% as.data.table(.) %>% 
                   select(., contains(c("no", "n.higher"))) %>%
                   replace_duplicate(.) %>% rowMeans(., na.rm = T))
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
  
  # Do nothing if all sims have equal length
  if (all(max_len == lens)) {
    return(as.data.table(x) %>% select(., contains("n.")))
  } else {
  # for shorter dfs replace recycled values with 0s
  out <- map2(.x = select(x, contains("n.")), .y = lens, ~ replace(.x, .y:max_len, values = NA))
  
  return(as.data.table(out) %>% select(., contains("n.")))
  }
}

# calculating mean over every web
match_lengths <- function(x, df = T) {
  end <- map(x, ~ which(.x == 0 | is.nan(.x)) %>% first()) # get indices of first zero (i.e. last non replicated entry)
  max <- unlist(end) %>% max(.) # max length
  out <- map2(x, end,  ~ replace(.x, .y:max, values = NA)) # replace recycled vals
  if(df)
    out <- as.data.table(out)
  return(out)
}

# helper function used in per_surv
list_divide <- function(x) {
  x / x[1] * 100 %>% 
    replace(., is.nan(.), 0)
}

# function to calculate percentages of remaining species
per_surv <- function(x, y, lower = T, original = F) {
  extc_in_network <- ifelse(isTRUE(lower), 1, 2)
  
  if (original == T) {
    list_divide(pluck(x, extc_in_network, y))
  } else {
    map(c("Atl" = 1, "aTl" = 2, "atL" = 3, "atl" = 4),
        ~ list_divide(pluck(x, extc_in_network, y, .x)))
  }
}

# functions to plot mean persistance of sp from one trophic level by removal of
# sp from the other trophic level in percent
library(ggplot2)
library(ggpubr)
library(grid)

plot_extc_facet <- function(extc, ci, extc_norew, ci_norew, org, org_norew, ci_org, ci_org_norew,
                            save = F, view = T, both = F, abund = F) {
  
  # extinction on lower or higher level ?
  lvl_match_extc <- grepl("lower", substitute(extc))
  lvl_match_norew_extc <- grepl("lower", substitute(extc_norew))
  
  lvl_match_ci <- grepl("lower", substitute(ci))
  lvl_match_norew_ci <- grepl("lower", substitute(ci_norew))
  
  ifelse(lvl_match_extc, lvl <- "lower", lvl <- "higher")

  # check if extc/ci & extc_norew/ci_norew match
  if (lvl_match_extc != lvl_match_norew_extc) {
  stop("extc and extc_norew don't match please provide
         data from matching simulations")
  }
  
  if (lvl_match_ci != lvl_match_norew_ci) {
    stop("ci and ci_norew don't match please provide
         data from matching simulations")
  }
  
  # combine dfs 
  extc <- bind_rows(extc, extc_norew)
  ci <- bind_rows(ci, ci_norew)
  org <- bind_rows(org, org_norew)
  ci_org <- bind_rows(ci_org, ci_org_norew)

  # automate axis labels
  ifelse(lvl == "lower", xlab <- "plants", xlab <- "animals")
  ifelse(lvl == "lower", ylab <- "animals", ylab <- "plants")
  
  # colors for facet boundary boxes
  colors <- c("abund" = "firebrick", "trait" = "dodgerblue", "phylo" = "burlywood4")
  
  # create df for ggplot
  df <- cbind(ci,
             select(extc, -c("id", "com_vars"))) %>% 
    tidyr::unite(., col = "group", c("id", "com_vars"), remove = F)
  
  df_org <- cbind(ci_org, select(org, -c("id", "com_vars")))
  
  # set grouping variable to factor
    df$group <- factor(df$group, levels = c("abund_Atl", "abund_aTl", "abund_atL",
                                            "trait_Atl", "trait_aTl", "trait_atL",
                                            "phylo_Atl", "phylo_aTl", "phylo_atL",
                                            "norew_Atl", "norew_aTl", "norew_atL"))
    
  # create list of labels for facet
    id_labs <- map(unique(df$id), function(a) {
      id_labs <- c("Atl", "aTl", "atL") %>%
        set_names(c(paste(a, "Atl", sep = "_"),
                    paste(a, "aTl", sep = "_"),
                    paste(a, "atL", sep = "_")))
    })
    names(id_labs) <- c(unique(df$id))
    
    # set facet layout
    nrow <- 3
    ncol <- 1
    
    # create vector of id to map over
    map_ids <- unique(df$id)[-length(unique(df$id))] # delete last one (norew)
  
  # create plots
  plots <- map(map_ids, function(z) {
    # get facet labels
      labs <- pluck(id_labs, z)
    
      # plot sims
      subplots <- map(c("Atl", "aTl", "atL"), function(a) {
        
        # get facet label for current plot; needs single value for "facet cheat" 
          df$title <- labs[which(unique(df$com_vars) == a)]
        
        ggplot(filter(df, id == z & com_vars == a)) +
          # geom_area(aes(x = rev(x), y = y), fill = "grey90") +
          geom_line(aes(rev(x), y), color = "grey60",
                    data = filter(df, id == "norew" & com_vars == a)) + # means norew
          geom_line(aes(rev(x_lower), x_higher), color = "grey60",
                    linetype = 2,
                    data = filter(df, id == "norew" & com_vars == a)) + # lower ci norew
          geom_line(aes(rev(y_lower), y_higher), color = "grey60",
                    linetype = 2,
                    data = filter(df, id == "norew" & com_vars == a)) + # higher ci norew
          geom_line(aes(rev(x), y, color = "firebrick")) + # means
          geom_line(aes(rev(x_lower), x_higher),
                    linetype = 2) + # lower ci
          geom_line(aes(rev(y_lower), y_higher),
                    linetype = 2) + # higher ci
          labs(x = paste(xlab, "removed"), y = paste(ylab, "persisting")) +
          # labs(x = NULL, y = NULL) +
          guides(color = "none") +
          facet_wrap(. ~ title, labeller = labeller(group = labs))  +
          theme(strip.text = element_text(),
                  strip.background = element_rect(color = colors[z])) # set box color
      })
        # plot org
        
        # get facet label for current plot; needs single value for "facet cheat" 
        df_org$title <- "atl"

        subplots[[4]] <- ggplot(filter(df_org, id == z)) +
          geom_line(aes(rev(x), y), color = "grey80",
                    data = filter(df_org, id == "norew")) + # means norew
          geom_line(aes(rev(x_lower), x_higher), color = "grey80",
                    linetype = 2,
                    data = filter(df_org, id == "norew")) + # lower ci norew
          geom_line(aes(rev(y_lower), y_higher), color = "grey80",
                    linetype = 2,
                    data = filter(df_org, id == "norew")) + # higher ci norew
          geom_line(aes(rev(x), y, color = "firebrick")) + # means
          geom_line(aes(rev(x_lower), x_higher),
                    linetype = 2) + # lower ci
          geom_line(aes(rev(y_lower), y_higher),
                    linetype = 2) + # higher ci
          labs(x = paste(xlab, "removed"), y = paste(ylab, "persisting")) +
          # labs(x = NULL, y = NULL) +
          guides(color = "none") +
          facet_wrap(. ~ title, labeller = labeller(group = labs))  +
          theme(strip.text = element_text(),
                strip.background = element_rect(color = colors[z])) # set box color

      # arrange  
      ggarrange(subplots[[4]], subplots[[1]], subplots[[2]], subplots[[3]],
                nrow = 1)
      })
  
  # view plots    
  if (view) {
    out <- ggarrange(plots[[1]], plots[[2]], plots[[3]],
                     nrow = nrow, ncol = ncol)
    print(out)
  }
  
  # set filenames 
  if (save) {
    name <- paste("extinction_cascade", lvl, "trophic_level", paste0(coextc_thr*100, ".pdf"), sep = "_")
    
    if (both) {
      name <- paste("extinction_cascade_both_trophic_levels", paste0(coextc_thr*100, ".pdf"),
                    sep = "_")
    }
    if (abund) {
      name <- paste("extinction_cascade", lvl, "trophic_level_abund", paste0(coextc_thr*100, ".pdf"),
                    sep = "_")
    }
    ggsave(
      name,
      path = paste0(getwd(), "/plot_sink"),
      plot = out,
      device = "pdf",
      width = 2100,
      height = 1705,
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

plot_extc_alt <- function(x, org, ci, ci_org, save = F, view = T, lower = F, spaghetti = F,
                          by_com_var = F){
  
  com_vars  <- c("Abundance" = 1, "Trait" = 2, "Phylogeny" = 3)
  ifelse(lower == T, lvl <- "lower", lvl <- "higher")
  ifelse(lower == T, xlab <- "plants", xlab <- "animals")
  ifelse(lower == T, ylab <- "animals", ylab <- "plants")
  
  if (by_com_var) {
    if (spaghetti) {
      plots <- map(com_vars, ~ ggplot() + 
                     geom_line(aes(rev(pluck(ci, "lower", 1, 1, .x)),
                                   pluck(ci, "higher", 1, 1, .x),
                                   color = "Abundance"), linetype = 2, size = .3) + # lower
                     geom_line(aes(rev(pluck(ci, "lower", 2, 1, .x)),
                                   pluck(ci, "higher", 2, 1, .x),
                                   color = "Abundance"), linetype = 2, size = .3) + # upper
                     geom_line(aes(rev(pluck(x, "lower", 1, .x)), pluck(x, "higher", 1, .x),
                                   color = "Abundance"), linetype = 1, size = .5) + # mean
                     geom_line(aes(rev(pluck(ci, "lower", 1, 2, .x)),
                                   pluck(ci, "higher", 1, 2, .x),
                                   color = "Trait"), linetype = 2, size = .3) + # lower
                     geom_line(aes(rev(pluck(ci, "lower", 2, 2, .x)),
                                   pluck(ci, "higher", 2, 2, .x),
                                   color = "Trait"), linetype = 2, size = .3) + # upper
                     geom_line(aes(rev(pluck(x, "lower", 2, .x)), pluck(x, "higher", 2, .x),
                                   color = "Trait"), linetype = 1, size = .5) + # mean
                     geom_line(aes(rev(pluck(ci, "lower", 1, 3, .x)),
                                   pluck(ci, "higher", 1, 3, .x),
                                   color = "Phylogeny"), linetype = 2, size = .3) + # lower
                     geom_line(aes(rev(pluck(ci, "lower", 2, 3, .x)),
                                   pluck(ci, "higher", 2, 3, .x),
                                   color = "Phylogeny"), linetype = 2, size = .3) + # upper
                     geom_line(aes(rev(pluck(x, "lower", 3, .x)), pluck(x, "higher", 3, .x),
                                   color = "Phylogeny"), linetype = 1, size = .5) + # mean
                     scale_color_manual(name = "Rewiring",
                                        values = c("Abundance" = "black",
                                                   "Trait" = "firebrick",
                                                   "Phylogeny" = "dodgerblue",
                                                   "original" = "burlywood4")) +
                     labs(x = paste(xlab, "removed"), y = paste(ylab, "persisting"),
                          title = paste("Extinction cascade", lvl, "trophic level", names(com_vars[.x]), "based networks"))
      )
    }
  } else {
    
  if (spaghetti) {
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
                   geom_line(aes(rev(pluck(ci, "lower", 1, .x, 3)),
                                 pluck(ci, "higher", 1, .x, 3),
                                 color = "atL"), linetype = 2, size = .3) + # lower
                   geom_line(aes(rev(pluck(ci, "lower", 2, .x, 3)),
                                 pluck(ci, "higher", 2, .x, 3),
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
                   scale_color_manual(name = "Contribution",
                                      values = c("Atl" = "black",
                                                 "aTl" = "firebrick",
                                                 "atL" = "dodgerblue",
                                                 "original" = "burlywood4")) +
                   labs(x = paste(xlab, "removed"), y = paste(ylab, "persisting"),
                        title = paste("Extinction cascade", lvl, "trophic level", names(com_vars[.x]), "rewiring"))
    )
  } else {
   plots_Atl <- map(com_vars, ~ ggplot() + 
          geom_line(aes(rev(pluck(ci, "lower", 1, .x, 1)),
                        pluck(ci, "higher", 1, .x, 1),
                        color = "Atl"), linetype = 2, size = .3) + # lower
          geom_line(aes(rev(pluck(ci, "lower", 2, .x, 1)),
                        pluck(ci, "higher", 2, .x, 1),
                        color = "Atl"), linetype = 2, size = .3) + # upper
          geom_line(aes(rev(pluck(x, "lower", .x, "Atl")), pluck(x, "higher", .x, "Atl"),
                        color = "Atl"), linetype = 1, size = .5) + # mean
            scale_color_manual(name = "Contribution",
                               values = c("Atl" = "black",
                                          "aTl" = "firebrick",
                                          "atL" = "dodgerblue",
                                          "original" = "burlywood4")) +
            labs(x = paste(xlab, "removed"), y = paste(ylab, "persisting"),
                 # title = paste("Extinction cascade", lvl, "trophic level",
                 #               names(com_vars[.x]), "rewiring")
                 ))
          
          plots_aTl <- map(com_vars, ~ ggplot() + 
            geom_line(aes(rev(pluck(ci, "lower", 1, .x, 2)),
                        pluck(ci, "higher", 1, .x, 2),
                        color = "aTl"), linetype = 2, size = .3) + # lower
          geom_line(aes(rev(pluck(ci, "lower", 2, .x, 2)),
                        pluck(ci, "higher", 2, .x, 2),
                        color = "aTl"), linetype = 2, size = .3) + # upper
          geom_line(aes(rev(pluck(x, "lower", .x, "aTl")), pluck(x, "higher", .x, "aTl"),
                        color = "aTl"), linetype = 1, size = .5) + # mean
              scale_color_manual(name = "Contribution",
                                 values = c("Atl" = "black",
                                            "aTl" = "firebrick",
                                            "atL" = "dodgerblue",
                                            "original" = "burlywood4")) +
              labs(x = paste(xlab, "removed"), y = paste(ylab, "persisting"),
                   # title = paste("Extinction cascade", lvl, "trophic level",
                   #               names(com_vars[.x]), "rewiring")
                   ))
          
          plots_atL <- map(com_vars, ~ ggplot() +     
          geom_line(aes(rev(pluck(ci, "lower", 1, .x, 3)),
                        pluck(ci, "higher", 1, .x, 3),
                        color = "atL"), linetype = 2, size = .3) + # lower
          geom_line(aes(rev(pluck(ci, "lower", 2, .x, 3)),
                        pluck(ci, "higher", 2, .x, 3),
                        color = "atL"), linetype = 2, size = .3) + # upper
          geom_line(aes(rev(pluck(x, "lower", .x, "atL")), pluck(x, "higher", .x, "atL"),
                        color = "atL"), linetype = 1, size = .5) + # mean
                scale_color_manual(name = "Contribution",
                                   values = c("Atl" = "black",
                                              "aTl" = "firebrick",
                                              "atL" = "dodgerblue",
                                              "original" = "burlywood4")) +
                labs(x = paste(xlab, "removed"), y = paste(ylab, "persisting"),
                     # title = paste("Extinction cascade", lvl, "trophic level",
                     #               names(com_vars[.x]), "rewiring")
                     ))
          
          plots_org <- map(com_vars, ~ ggplot() +       
          geom_line(aes(rev(pluck(ci_org, "lower", 1, .x)),
                        pluck(ci_org, "higher", 1, .x),
                        color = "original"), linetype = 2, size = .3) + # lower
          geom_line(aes(rev(pluck(ci_org, "lower", 2, .x)),
                        pluck(ci_org, "higher", 2, .x),
                        color = "original"), linetype = 2, size = .3) + # upper
          geom_line(aes(rev(pluck(org, "lower", .x)), pluck(org, "higher", .x),
                      color = "original"), linetype = 1, size = .5) + # mean
                  scale_color_manual(name = "Contribution",
                                     values = c("Atl" = "black",
                                                "aTl" = "firebrick",
                                                "atL" = "dodgerblue",
                                                "original" = "burlywood4")) +
                  labs(x = paste(xlab, "removed"), y = paste(ylab, "persisting"),
                       # title = paste("Extinction cascade", lvl, "trophic level",
                       #               names(com_vars[.x]), "rewiring")
                       ))
          
          # plots_tin <- map(com_vars, ~ ggplot() +       
          #                    geom_line(aes(rev(pluck(ci_tinoco, "lower", 1, .x)),
          #                                  pluck(ci_tinoco, "higher", 1, .x),
          #                                  color = "original"), linetype = 2, size = .3) + # lower
          #                    geom_line(aes(rev(pluck(ci_tinoco, "lower", 2, .x)),
          #                                  pluck(ci_tinoco, "higher", 2, .x),
          #                                  color = "original"), linetype = 2, size = .3) + # upper
          #                    geom_line(aes(rev(pluck(tinoco, "lower", .x)), pluck(tinoco, "higher", .x),
          #                                  color = "original"), linetype = 1, size = .5) + # mean
          #                    scale_color_manual(name = "Contribution",
          #                                       values = c("Atl" = "black",
          #                                                  "aTl" = "firebrick",
          #                                                  "atL" = "dodgerblue",
          #                                                  "original" = "burlywood4",
          #                                                  "w/o rewiring" = "seagreen")) +
          #                    labs(x = paste(xlab, "removed"), y = paste(ylab, "persisting"),
          #                         # title = paste("Extinction cascade", lvl, "trophic level",
          #                         #               names(com_vars[.x]), "rewiring")
          #                    ))
          # 
          
          # plots <- ggarrange(plots_Atl[[1]], plots_aTl[[1]], plots_atL[[1]],
          #           plots_org[[1]], plots_norew[[1]], plots_tin[[1]], plots_Atl[[2]],
          #           plots_aTl[[2]], plots_atL[[2]], plots_org[[2]],
          #           plots_norew[[2]], plots_tin[[2]], plots_Atl[[3]], plots_aTl[[3]],
          #           plots_atL[[3]], plots_org[[3]], plots_norew[[3]], plots_tin[[3]],
          #           ncol = 6, nrow = 3, common.legend = T, legend = "bottom")
          # 
   
          plots <- ggarrange(plots_Atl[[1]], plots_aTl[[1]], plots_atL[[1]],
                             plots_org[[1]], plots_Atl[[2]],
                             plots_aTl[[2]], plots_atL[[2]], plots_org[[2]],
                             plots_Atl[[3]], plots_aTl[[3]],
                             plots_atL[[3]], plots_org[[3]],
                             ncol = 4, nrow = 3, common.legend = T, legend = "bottom")
          
  }
  }
  
   if (save) {
     if (spaghetti) {
       if (by_com_var) {
         map(c("abund" = 1, "trait" = 2, "phylo" = 3), ~ ggsave(
           paste("alt_extinction_cascade", lvl, "trophic_level",
                 names(com_vars[.x]), coextc_thr*100, "by_com_var.pdf", sep = "_"),
           path = paste0(getwd(), "/plot_sink"),
           plot = plots[[.x]],
           device = "pdf",
           width = 1900,
           height = 1205,
           units = "px"))
       } else {
         map(c("abund" = 1, "trait" = 2, "phylo" = 3), ~ ggsave(
         paste("alt_extinction_cascade", lvl, "trophic_level",
               names(com_vars[.x]), paste0(coextc_thr*100, ".pdf"), sep = "_"),
         path = paste0(getwd(), "/plot_sink"),
         plot = plots[[.x]],
         device = "pdf",
         width = 1900,
         height = 1205,
         units = "px"))
       }
     } else {
       ggsave(
         paste("alt_extinction_cascade", lvl, "trophic_level", paste0(coextc_thr*100, ".pdf"), sep = "_"),
         path = paste0(getwd(), "/plot_sink"),
         plot = plots,
         device = "pdf",
         width = 1900,
         height = 1205,
         units = "px")
     }
  }
   
  if (view){
     return(plots)
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
  (sum(diff(x) * (head(y,-1)+tail(y,-1)))/2)/100
  }

# calculate network robustness; adapted to work with augmented extinction fy
robustness_aug <- function (object, lower = T) {
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

conf_int <- function(x, means, trph_lvl = "lower", lower = T, norew_org = F,
                     tin = F, noreworg = F) {
  if(!(trph_lvl == "lower" | trph_lvl == "higher"))
    stop("Invalid value for trph_lvl ! Specify valid trophic level (lower or higher)")

  # calculate std devs
  if (tin) {
    sds <- map(1:3, ~ pluck(x, trph_lvl, .x) %>% sd)
  } else {
    if (norew_org) {
      sds <- map(.x = 1:3, function(a) {
        map(seq(n_webs), ~ pluck(x, .x, trph_lvl, a)) %>% as.data.table %>% 
          match_lengths() %>% replace(., is.na(.), 0) %>% apply(., 1, sd)})
    } else {
      if (noreworg) {
       sds <- map(seq(n_webs), ~ pluck(x, .x, trph_lvl)) %>% as.data.table %>%
          match_lengths() %>% replace(., is.na(.), 0) %>% apply(., 1, sd)
      } else {
        sds <- map(.x = 1:3, function(a) map(.x = 1:3, function(b) {
          map(seq(n_webs), ~ pluck(x, .x, trph_lvl, a, b)) %>% as.data.table %>% 
            match_lengths() %>% replace(., is.na(.), 0) %>% apply(., 1, sd)}))
        }
    }
  }
  
  
  # calulate margins
  if (tin | noreworg) {
    margins <- qt(0.975, df = length(sds) - 1) * sds/sqrt(length(sds))
  } else {
      if (norew_org) {
        margins <- map(1:3, ~ qt(0.975, df = length(pluck(sds, .x)) - 1) * pluck(sds, .x)/sqrt(length(pluck(sds, .x))))
      } else {
        margins <- map(1:3, function(a) {
          map(1:3, ~ qt(0.975, df = length(pluck(sds, a, .x)) - 1) * pluck(sds, a, .x)/sqrt(length(pluck(sds, a, .x))))})
      }
    }
  
  # calculate lower/upper ci
  if (tin | noreworg) {
    out <- (pluck(means, trph_lvl) - margins) %>% ifelse(. < 0, 0, .)
    } else {
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
  
  ## using with extc sims org, produces df with id set to "org" and com_vars
  ## "abund", "trait", "phylo"
  
  # check if all extc seqs have same length
  if (org|norew) {
    if (ci) {
      if (org & norew) {
        # 2.5 % ci
        lens_loci <- length(pluck(x, "lower", 1)) %>% unlist()
        
        # 95 % ci
        lens_hici <- length(pluck(x, "lower", 2)) %>% unlist()
      } else {
        # calculate lengths of each extc seq, 2.5 % ci
        lens_loci <- map(1:3, ~ length(pluck(x, "lower", 1, .x))) %>% unlist()
        
        # 95% ci
        lens_hici <- map(1:3, ~ length(pluck(x, "lower", 2, .x))) %>% unlist()
      }
    } else {
      if (org & norew) {
        # lens of each extc seq
        lens <- length(pluck(x, "lower")) %>% unlist()
      } else {
        # calculate lengths of each extc seq
        lens <- map(1:3,  ~ length(pluck(x, "lower", .x))) %>% unlist()
      }
    }
  } else {
    if (ci) {
      # calculate lengths of each extc seq, 2.5 % ci
      lens_loci <- map(1:3, function(y) {
        map(1:len, ~ length(pluck(x, "lower", 1, y, .x))) %>% unlist()
      })
      
      # 95% ci
      lens_hici <- map(1:3, function(y) {
        map(1:len, ~ length(pluck(x, "lower", 2, y, .x))) %>% unlist()
      })
    } else {
      # calculate lengths of each extc seq
      lens <- map(1:3, function(y) {
        map(1:3, ~ length(pluck(x, "lower", y, .x))) %>% unlist()
      })
    }
  }
  
  # match lengths of extc seqs if not same
  if (org|norew) {
    if (ci) {
      #  ci for org & norew simulation with different seq lengths
      if (!(all(max(unlist(lens_loci)) == unlist(lens_loci)) |
            all(max(unlist(lens_hici)) == unlist(lens_hici)))) {
        
        if (org & norew) {
          # calculate max lengths of each extc seq, 2.5% ci
          maxlens_loci <- length(pluck(x, "lower", 1)) %>%
            unlist() %>% max()
          
          # 95% ci
          maxlens_hici <- length(pluck(x, "lower", 2)) %>%
            unlist() %>% max()
          
          # match lengths of lower extc seqs, 2.5 % ci
          data_lower_loci <- append(pluck(x, "lower", 1),
                                               rep(0, maxlens_loci - lens_loci))
          
          # 95 % ci
          data_lower_hici <- append(pluck(x, "lower", 2),
                                               rep(0, maxlens_hici - lens_hici))
          
          # match length of higher extc seqs, 2.5% ci
          data_higher_loci <- append(pluck(x, "higher", 1),
                                                rep(0, maxlens_loci - lens_loci))
          
          # 95% ci
          data_higher_hici <- append(pluck(x, "higher", 2),
                                                rep(0, maxlens_hici - lens_hici))
        } else {
          # calculate max lengths of each extc seq, 2.5% ci
          maxlens_loci <- map(1:3, ~ length(pluck(x, "lower", 1, .x))) %>%
            unlist() %>% max()
          
          # 95% ci
          maxlens_hici <- map(1:3, ~ length(pluck(x, "lower", 2, .x))) %>%
            unlist() %>% max()
          
          # match lengths of lower extc seqs, 2.5 % ci
          data_lower_loci <- map(1:3, ~ append(pluck(x, "lower", 1, .x),
                                                        rep(0, maxlens_loci - lens_loci[.x])))
          
          # 95 % ci
          data_lower_hici <- map(1:3, ~ append(pluck(x, "lower", 2, .x),
                                                        rep(0, maxlens_hici - lens_hici[.x])))
          
          # match length of higher extc seqs, 2.5% ci
          data_higher_loci <- map(1:3, ~ append(pluck(x, "higher", 1, .x),
                                                         rep(0, maxlens_loci - lens_loci[.x])))
          
          # 95% ci
          data_higher_hici <- map(1:3, ~ append(pluck(x, "higher", 2, .x),
                                                         rep(0, maxlens_hici - lens_hici[.x])))
        }
      } else {
        if (org & norew) {
          data_lower_loci <- pluck(x, "lower", 1)
        
          data_lower_hici <- pluck(x, "lower", 2)
          
          data_higher_loci <- pluck(x, "higher", 1)
          
          data_higher_hici <- pluck(x, "higher", 2)
        } else {
          data_lower_loci <- map(1:3, ~ pluck(x, "lower", 1, .x))
          
          data_lower_hici <- map(1:3, ~ pluck(x, "lower", 2, .x))
          
          data_higher_loci <- map(1:3, ~ pluck(x, "higher", 1, .x))
          
          data_higher_hici <- map(1:3, ~ pluck(x, "higher", 2, .x))
        }
      }
    } else {
      # extc data for org & norew simulation with different lenghts
      if(!(all(max(unlist(lens)) == unlist(lens)) |
           all(max(unlist(lens)) == unlist(lens)))) {
        if (org &  norew) {
          # calculate max lengths of each extc seq
          maxlens <- length(pluck(x, "lower")) %>% unlist() %>% max()
          
          # match lengths of lower extc seqs
          data_lower <- append(pluck(x, "lower"),
                                          rep(0, maxlens - lens))
          
          # match length of higher extc seqs
          data_higher <- append(pluck(x, "higher"),
                                           rep(0, maxlens - lens))
        } else {
        # calculate max lengths of each extc seq
        maxlens <- map(1:3, ~ length(pluck(x, "lower", .x))) %>% unlist() %>% max()
        
        # match lengths of lower extc seqs
        data_lower <- map(1:3, ~ append(pluck(x, "lower", .x),
                              rep(0, maxlens - lens[.x])))
        
        # match length of higher extc seqs
        data_higher <- map(1:3, ~ append(pluck(x, "higher", .x),
                              rep(0, maxlens - lens[.x])))
        }
      } else {
        if (org & norew) {
          data_lower <- pluck(x, "lower")
          
          data_higher <- pluck(x, "higher")
        } else {
        data_lower <- map(1:3, ~ pluck(x, "lower", .x))
        
        data_higher <- map(1:3, ~ pluck(x, "higher", .x))
        }
      }
    }
  } else {
      if (ci) {
        # ci for non org & norew simulation with different lengths
        if (!(all(max(unlist(lens_loci)) == unlist(lens_loci)) |
              all(max(unlist(lens_hici)) == unlist(lens_hici)))) {
          # calculate max lengths of each extc seq, 2.5% ci
          maxlens_loci <- map(1:3, function(y) {
            map(1:len, ~ length(pluck(x, "lower", 1, y, .x))) %>% unlist() %>% max()
          })
          
          # 95% ci
          maxlens_hici <- map(1:3, function(y) {
            map(1:len, ~ length(pluck(x, "lower", 2, y, .x))) %>% unlist() %>% max()
          })
          
          # match lengths of lower extc seqs, 2.5 % ci
          data_lower_loci <- map(1:3, function(y) {
            map(1:len, ~ append(pluck(x, "lower", 1, y, .x),
                                rep(0, maxlens_loci[[y]] - lens_loci[[y]][.x])))
          })
          # 95 % ci
          data_lower_hici <- map(1:3, function(y) {
            map(1:len, ~ append(pluck(x, "lower", 2, y, .x),
                                rep(0, maxlens_hici[[y]] - lens_hici[[y]][.x])))
          })
          
          # match length of higher extc seqs, 2.5% ci
          data_higher_loci <- map(1:3, function(y) {
            map(1:len, ~ append(pluck(x, "higher", 1, y, .x),
                                rep(0, maxlens_loci[[y]] - lens_loci[[y]][.x])))
          })
          
          # 95% ci
          data_higher_hici <- map(1:3, function(y) {
            map(1:len, ~ append(pluck(x, "higher", 2, y, .x),
                                rep(0, maxlens_hici[[y]] - lens_hici[[y]][.x])))
          })

        } else {
          data_lower_loci <- map(1:3, function(y) {
            map(1:len, ~ pluck(x, "lower", 1, y, .x))})
          
          data_lower_hici <- map(1:3, function(y) {
            map(1:len, ~ pluck(x, "lower", 2, y, .x))})
          
          data_higher_loci <- map(1:3, function(y) {
            map(1:len, ~ pluck(x, "higher", 1, y, .x))})
          
          data_higher_hici <- map(1:3, function(y) {
            map(1:len, ~ pluck(x, "higher", 2, y, .x))})
        }
      } else {
        # extc data for non org & norew simulation with different lengths
        if (!(all(max(unlist(lens)) == unlist(lens)) |
             all(max(unlist(lens)) == unlist(lens)))) {
          # calculate max lengths of each extc seq
          maxlens <- map(1:3, function(y) {
            map(1:len, ~ length(pluck(x, "lower", y, .x))) %>% unlist() %>% max()
          })
        
        # match lengths of lower extc seqs
        data_lower <- map(1:3, function(y) {
          map(1:len, ~ append(pluck(x, "lower", y, .x),
                              rep(0, maxlens[[y]] - lens[[y]][.x])))
        })
        
        # match length of higher extc seqs
        data_higher <- map(1:3, function(y) {
          map(1:len, ~ append(pluck(x, "higher", y, .x),
                              rep(0, maxlens[[y]] - lens[[y]][.x])))
        })

        } else {
          data_lower <- map(1:3, ~ pluck(x, "lower", .x))
          data_higher <- map(1:3, ~ pluck(x, "higher", .x))
        }
      }
    }
  
  if (ci) {
    if (org & norew) {
      # make ci df 
      lower <- cbind(bind_cols(pluck(data_lower_loci),
                        .id = "norew") %>%
                melt(., id.vars = ".id",
                     value.name = "x_lower",
                     variable.name = "com_vars"),
              bind_cols(pluck(data_lower_hici),
                        .id = "norew") %>%
                melt(., id.vars = ".id",
                     value.name = "y_lower",
                     variable.name = "com_vars")) %>% 
          select(., -c(4,5)) %>% 
          setNames(., c("id", "com_vars", "x_lower", "y_lower")) %>% 
        bind_rows()
      
      higher <- cbind(bind_cols(pluck(data_higher_loci),
                        .id = "norew") %>%
                melt(., id.vars = ".id",
                     value.name = "x_higher",
                     variable.name = "com_vars"),
              bind_cols(pluck(data_higher_hici),
                        .id = "norew") %>%
                melt(., id.vars = ".id",
                     value.name = "y_higher",
                     variable.name = "com_vars")) %>% 
          select(., -c(4,5)) %>% 
          setNames(., c("id", "com_vars", "x_higher", "y_higher")) %>% 
        bind_rows()
      
      out <- cbind(lower, select(higher, -c(1, 2)))
    } else {
      # make ci df 
      lower <- map(1:3, function(y) {
        cbind(bind_cols(pluck(data_lower_loci, y),
                        .id = c(rew_names[y])) %>%
                melt(., id.vars = ".id",
                     value.name = "x_lower",
                     variable.name = "com_vars"),
              bind_cols(pluck(data_lower_hici, y),
                        .id = c(rew_names[y])) %>%
                melt(., id.vars = ".id",
                     value.name = "y_lower",
                     variable.name = "com_vars")) %>% 
          select(., -c(4,5)) %>% 
          setNames(., c("id", "com_vars", "x_lower", "y_lower"))}) %>% 
        bind_rows()
      
      higher <- map(1:3, function(y) {
        cbind(bind_cols(pluck(data_higher_loci, y),
                            .id = c(rew_names[y])) %>%
                              melt(., id.vars = ".id",
                                   value.name = "x_higher",
                                   variable.name = "com_vars"),
                  bind_cols(pluck(data_higher_hici, y),
                            .id = c(rew_names[y])) %>%
                    melt(., id.vars = ".id",
                         value.name = "y_higher",
                         variable.name = "com_vars")) %>% 
          select(., -c(4,5)) %>% 
          setNames(., c("id", "com_vars", "x_higher", "y_higher"))}) %>% 
        bind_rows()
      
      out <- cbind(lower, select(higher, -c(1, 2)))
    }
  } else {
    if (org & norew) {
      # make extc df 
      out <- cbind(bind_cols(pluck(data_lower), .id = "norew") %>%
                melt(., id.vars = ".id",
                     value.name = "x",
                     variable.name = "com_vars"),
              bind_cols(pluck(data_higher), .id = "norew") %>%
                melt(., id.vars = ".id",
                     value.name = "y",
                     variable.name = "com_vars")) %>% bind_rows(.) %>%
        select(., -c(4,5)) %>% 
        setNames(., c("id", "com_vars", "x", "y"))
    } else {
      # make extc df 
      out <- map(1:3, function(y) {
        cbind(bind_cols(pluck(data_lower, y), .id = rew_names[y]) %>%
                melt(., id.vars = ".id",
                     value.name = "x",
                     variable.name = "com_vars"),
              bind_cols(pluck(data_higher, y), .id = rew_names[y]) %>%
                melt(., id.vars = ".id",
                     value.name = "y",
                     variable.name = "com_vars"))}) %>% bind_rows(.) %>%
        select(., -c(4,5)) %>% 
        setNames(., c("id", "com_vars", "x", "y"))
    }
  }
  
  if (org & norew) {
    out$com_vars <- "org"
    return (out)
  }
  
  # set correct levels for id & com_vars
  if (org|norew) {
    if (ci) {
      length_each <- max(lens_loci)
    } else {
      length_each <- max(lens)
      }
    } else {
      if (ci) {
        length_each <- max(lens_loci[[1]])
      } else {
        length_each <- max(lens[[1]])
      }
    }
  
  if (org|norew) {
    if (ci) {
      if (org) {
        out$com_vars <- "org"
      } else {
        out$com_vars <- rep(c("Atl", "aTl", "atL"), each = length_each)
        }
    } else {
      if (org) {
        out$com_vars <- "org"
      } else {
        out$com_vars <- rep(c("Atl", "aTl", "atL"), each = length_each) 
        }
      }
    } else {
      levels(out$com_vars) <- c("Atl", "aTl", "atL")
  }
  
  if (org|norew) {
    if (ci) {
      if (org) {
        if (norew) {
          out$id <- rep(c("abund_norew", "trait_norew", "phylo_norew"), each = length_each)
        } else {
        out$id <- rep(c("abund", "trait", "phylo"), each = length_each)
        }
      } else {
        out$id <- "norew"
        }
    } else {
      if (org) {
        if (norew) {
          out$id <- rep(c("abund_norew", "trait_norew", "phylo_norew"), each = length_each)
        } else {
        out$id <- rep(c("abund", "trait", "phylo"), each = length_each)
        }
      } else {
        out$id <- "norew"
        }
      }
    } else {
      levels(out$id) <- c("abund", "trait", "phylo")
    }
  
  return(out)
}

# changing lists to df for easier plotting
list_to_df2 <- function(x, org = F, norew = F) {
  ifelse(org, id <- "org", id <- "norew")
  
  # create df
  if (org | norew) {
    if (norew & !org) {
      # norew
      x_df <- map(ctrbs, ~ reshape2::melt(pluck(x, "lower", .x),
                                  value.name = "x")[, c(2, 3)]  %>% 
              bind_cols("x" =., "com_vars" = names(ctrbs[.x])) %>% 
              bind_cols("x" = .,  "id" = "norew")) %>% 
          bind_rows
      y_df <- map(ctrbs, ~ reshape2::melt(pluck(x, "higher", .x),
                                  value.name = "y")[, c(2, 3)]  %>% 
              bind_cols("y" =., "com_vars" = names(ctrbs[.x])) %>%
              bind_cols("y" = .,  "id" = "norew")) %>% 
          bind_rows
      
      out <- cbind(x_df[, c(1,2)], y_df)
      colnames(out) <- c("sims1", "x", "sims2", "y", "com_vars", "id")
    } else {
      if (org & norew) {
        # org + norew
        x_df <- reshape2::melt(pluck(x, "lower"),
                                          value.name = "x")[, c(2, 3)]  %>% 
                      bind_cols("x" =., "com_vars" = "org") %>% 
                      bind_cols("x" = .,  "id" = "norew")
        y_df <- reshape2::melt(pluck(x, "higher"),
                                          value.name = "y")[, c(2, 3)]  %>% 
                      bind_cols("y" =., "com_vars" = "org") %>%
                      bind_cols("y" = .,  "id" = "norew")
        
        out <- cbind(x_df[, c(1,2)], y_df)
        colnames(out) <- c("sims1", "x", "sims2", "y", "com_vars", "id")
      } else {
      # org
        x_df <- map(1:3, ~ reshape2::melt(pluck(x, "lower", .x),
                                          value.name = "x")[, c(2, 3)]  %>% 
                      bind_cols("x" =., "com_vars" = "org") %>% 
                      bind_cols("x" = .,  "id" = rew_names[.x])) %>% 
          bind_rows
        y_df <- map(1:3, ~ reshape2::melt(pluck(x, "higher", .x),
                                          value.name = "y")[, c(2, 3)]  %>% 
                      bind_cols("y" =., "com_vars" = "org") %>%
                      bind_cols("y" = .,  "id" = rew_names[.x])) %>% 
          bind_rows
        
        out <- cbind(x_df[, c(1,2)], y_df)
        colnames(out) <- c("sims1", "x", "sims2", "y", "com_vars", "id")
      }
    }
  } else {
    # lower/higher
    x_df <- map(1:3, function(y) {
      map(ctrbs, ~ reshape2::melt(pluck(x, "lower", y, .x),
                                value.name = "x")[, c(2, 3)]  %>% 
            bind_cols("x" =., "com_vars" = names(ctrbs[.x])) %>% 
            bind_cols("x" = .,  "id" = rew_names[y])) %>% 
        bind_rows}) %>% 
      bind_rows
    y_df <- map(1:3, function(y) {
      map(ctrbs, ~ reshape2::melt(pluck(x, "higher", y, .x),
                                value.name = "y")[, c(2, 3)]  %>% 
            bind_cols("y" =., "com_vars" = names(ctrbs[.x])) %>%
            bind_cols("y" = .,  "id" = rew_names[y])) %>% 
        bind_rows}) %>% 
      bind_rows
    
    out <- cbind(x_df[, c(1,2)], y_df)
    colnames(out) <- c("sims1", "x", "sims2", "y", "com_vars", "id")
  }
  return(out)
}

get_ci <- function(x, means, org = F, norew = F, tin = F) {
  if (org & norew) {
    out <- list("lower" = list("2.5" = conf_int(x = x,
                                                mean = means,
                                                tin = tin,
                                                noreworg = T),
                               "97.5" = conf_int(x = x,
                                                 mean = means,
                                                 lower = F,
                                                 tin = tin,
                                                 noreworg = T)),
                "higher" = list("2.5" = conf_int(x = x,
                                                 mean = means,
                                                 trph_lvl = "higher",
                                                 tin = tin,
                                                 noreworg = T),
                                "97.5" = conf_int(x = x,
                                                  mean = means,
                                                  trph_lvl = "higher",
                                                  tin = tin,
                                                  noreworg = T)))
    
    return(out) 
  }
  
  # set names
  if (org) {
    out <- list("lower" = list("2.5" = conf_int(x = x,
                                         mean = means,
                                         norew_org = T,
                                         tin = tin),
                        "97.5" = conf_int(x = x,
                                          mean = means,
                                          lower = F,
                                          norew_org = T,
                                          tin = tin)),
         "higher" = list("2.5" = conf_int(x = x,
                                          mean = means,
                                          trph_lvl = "higher",
                                          norew_org = T,
                                          tin = tin),
                         "97.5" = conf_int(x = x,
                                           mean = means,
                                           trph_lvl = "higher",
                                           lower = F,
                                           norew_org = T,
                                           tin = tin)))
  
  out <- modify_depth(out, 2, 
                      ~ set_names(.x, nm = c("abund", "trait", "phylo")))
  }
  
  if (norew) {
    out <- list("lower" = list("2.5" = conf_int(x = x,
                                         mean = means,
                                         norew_org = T,
                                         tin = tin),
                        "97.5" = conf_int(x = x,
                                          mean = means,
                                          lower = F,
                                          norew_org = T,
                                          tin = tin)),
         "higher" = list("2.5" = conf_int(x = x,
                                          mean = means,
                                          trph_lvl = "higher",
                                          norew_org = T,
                                          tin = tin),
                         "97.5" = conf_int(x = x,
                                           mean = means,
                                           trph_lvl = "higher",
                                           lower = F,
                                           norew_org = T,
                                           tin = tin)))
    
    out <- modify_depth(out, 2, 
                        ~ set_names(.x, nm = c("Atl", "aTl", "atL")))
  }
  
  if (!org & !norew) {
    out <- list("lower" = list("2.5" = conf_int(x = x,
                                         mean = means,
                                         tin = tin),
                        "97.5" = conf_int(x = x,
                                          mean = means,
                                          lower = F,
                                          tin = tin)),
         "higher" = list("2.5" = conf_int(x = x,
                                          mean = means,
                                          trph_lvl = "higher",
                                          tin = tin),
                         "97.5" = conf_int(x = x,
                                           mean = means,
                                           trph_lvl = "higher",
                                           lower = F,
                                           tin = tin)))
    out <- modify_depth(out, 2, 
                              ~ set_names(.x, nm = c("abund", "trait", "phylo")))
    
    out <- modify_depth(out, 3, 
                              ~ set_names(.x, nm = c("Atl", "aTl", "atL")))
  }
  
  return(out)
}

hist_dead <- function(x, lower = T) {
  
  com_vars <- c("Atl" = 1, "aTl" = 2, "atL" = 3, "atl" = 4)
  
  ifelse(lower, title <- "Lower trophic level",
         title <- "Higher trophic level")
  
  
  # get data
  if (lower) {
    # sims
    df <- map_dfc(com_vars, function(y) {
      map(seq(n_webs), ~ length(which(rowSums(pluck(x, .x, y, "web")) == 0))) %>%
        unlist})
    
    # org
    df$org <- map(seq(n_webs), function(y) {
      length(which(rowSums(pluck(init_sim_web, y)) == 0))}) %>%
      unlist
    
  } else {
    # sims
    df <- map_dfc(com_vars, function(y) {
      map(seq(n_webs), ~ length(which(colSums(pluck(x, .x, y, "web")) == 0))) %>%
        unlist})
    
    # org
    df$org <- map(seq(n_webs), function(y) {
      length(which(colSums(pluck(init_sim_web, y)) == 0))}) %>%
      unlist
  }
  
  ggplot(df) +
    geom_boxplot(aes(x = names(com_vars)[1], y = Atl)) +
    geom_boxplot(aes(x = names(com_vars)[2], y = aTl)) +
    geom_boxplot(aes(x = names(com_vars)[3], y = atL)) +
    geom_boxplot(aes(x = names(com_vars)[4], y = atl)) +
    geom_boxplot(aes(x = "org", y = org)) +
    scale_x_discrete(limits = c("Atl", "aTl", "atL", "atl", "org")) +
    labs(y = "No. of dead interactions", x = "Community variable importances") +
    theme_classic() +
    theme(text = element_text(size = 20)) 
}

# calculate number of aborted simulations
count_abort <- function(sim) {
  map(1:3, function(x) { # loop over rewiring methods
    map(ctrbs, function(y) { # loop over community variables
      map(seq(n_webs), function (z) { # loop over webs
        map(seq(n_sims), function(.x) { # loop over sims
          pluck(sim, z, x, y, .x) %>% is.null
        }) %>% unlist
      }) %>% unlist %>% which(. == T) %>% length
    })
  })
}

equalize_sp <- function(sims, init = F) {
  if (init) {
    min_low <- map(seq(n_webs), function(x) {
      nrow(sims[[x]]) %>% unlist %>% min
    })
    
    min_high <- map(seq(n_webs), function(x) {
      ncol(sims[[x]]) %>% unlist %>% min
    })
    
    # draw random sp according to lowest no of species from each network
    sp_low <- map(seq(n_webs), function(x) {
      rownames(sims[[x]]) %>%
            unlist %>%
            sample(., size = min_low[[x]])
      })
    
    sp_high <- map(seq(n_webs), function(x) {
      colnames(sims[[x]]) %>%
            unlist %>%
            sample(., size = min_high[[x]])
      })
    
    # crop networks to equal size
    out <- map(seq(n_webs), function(x) {
      sims[[x]][sp_low[[x]], sp_high[[x]]]
    })
  } else {
  # minimum no of species per trophic level
  min_low <- map(seq(n_webs), function(x) {
    map(seq(ctrbs), ~ nrow(sims[[x]][[.x]])) %>% unlist %>% min
  })
  
  min_high <- map(seq(n_webs), function(x) {
    map(seq(ctrbs), ~ ncol(sims[[x]][[.x]])) %>% unlist %>% min
  })
  
  # draw random sp according to lowest no of species from each network
  sp_low <- map(seq(n_webs), function(x) {
    map(seq(ctrbs), ~ rownames(sims[[x]][[.x]]) %>%
          unlist %>%
          sample(., size = min_low[[x]])
    )})
  
  sp_high <- map(seq(n_webs), function(x) {
    map(seq(ctrbs), ~ colnames(sims[[x]][[.x]]) %>%
          unlist %>%
          sample(., size = min_high[[x]])
    )})
  
  # crop networks to equal size
  out <- map(seq(n_webs), function(x) {
    map(seq(ctrbs), ~ sims[[x]][[.x]][sp_low[[x]][[.x]], sp_high[[x]][[.x]]])
  })
  }
  
  return(out)
}

# fx to delete species that have no interactions (i.e. all values in row/col are zero)
del_dead_int <- function(x) {
  # filter sp w/o any interactions
  sp_low <- which(rowSums(x) != 0L)
  sp_high <- which(colSums(x) != 0L)
  
  # drop sp w/o any interactions from web
  out <- x[sp_low, sp_high, drop = F]
}

auc_boxplot <- function(x, x.axis, log = F, save = F) {
  # extinction on lower or higher level ?
  lvl_match_extc <- grepl("lower", substitute(x))
  
  ifelse(lvl_match_extc, lvl <- "lower", lvl <- "higher")
  
  if (is.null(x.axis) | !any(c(rew_names,
                               names(ctrbs),
                               "rew",
                               "com_vars",
                               "org") %in% x.axis)) {
    stop("invalid x axis ! Please provide a character vector of a variable that
    should be used as x axis. Either use rew or com_vars or a specific 
    value/combination from among them")
  }
  
  # should y axis be logarithmic ?
  if (log) {
    scale <- "log10"
  } else {
    scale <- "none"
  }
  
  # filter data according to input
  if (length(x.axis) <= 1) {
    if (x.axis == "rew") {
      aucs <- x
      x_axis <- x.axis
      x_lab <- unique(x$rew)
    }
    
    if (x.axis == "com_vars") {
      aucs <- x
      x_axis <- x.axis
      x_lab <- unique(x$com_vars)
    }
    
    if (any(x.axis == c(rew_names, "norew", "abundtrait", "abundphylo"))){
      aucs <- filter(x, rew == x.axis)
      x_axis <- "com_vars"
      x_lab <- unique(x$com_vars)
    }
    
    if (any(x.axis == c(names(ctrbs), "org"))){
      aucs <- filter(x, com_vars == x.axis)
      x_axis <- "rew"
      x_lab <- unique(x$rew)
    }
  } else {
    if (length(x.axis) > 2)
      stop("Too many values in x.axis. Please provide a character vector of 
           lenght one or two")
    
    rew_choice <- ifelse(any(x.axis[1] == rew_names), x.axis[1], x.axis[2])
    com_vars_choice <- ifelse(any(x.axis[1] == c(names(ctrbs), "org")),
                              x.axis[1], x.axis[2])
    
    if (any(rew_choice == rew_names) & any(com_vars_choice == rew_names) |
        any(rew_choice == c(names(ctrbs), "org")) & any(com_vars_choice == c(names(ctrbs), "org"))) {
      stop("Invalid combination of x.axis values. Did you provide two rew/com_vars values ?")
    }
    
    aucs <- filter(x, com_vars == com_vars_choice, rew == rew_choice)
    x_axis <- NULL
    x_lab <- paste(rew_choice, com_vars_choice)
    
  }
  
  if (length(x.axis) <= 1) {
    out <- ggboxplot(y = "robustness",
              x = x_axis,
              color = "lvl",
              data = aucs,
              order = x_lab) +
      scale_color_manual(name = "Trophic level",
                         values = c("#000000", "#787878")) +
      yscale(scale) +
      theme(legend.position = "bottom") +
      labs(x = element_blank()) +
      theme(aspect.ratio = 1)
  } else {
    out <- ggboxplot(y = "robustness",
              x = x_axis,
              color = "lvl",
              data = aucs,
              order = x_lab) +
      scale_color_manual(name = "Trophic level",
                         values = c("#000000", "#787878")) +
      scale_x_discrete(labels = paste(rew_choice, com_vars_choice)) +
      yscale(scale) +
      theme(legend.position = "bottom") +
      labs(x = element_blank()) +
      theme(aspect.ratio = 1)
  }
  
  if (save) {
    # generate names
    if (length(x.axis) <= 1) {
        size <- c(14.3, 20.5)
      if (x.axis == "rew") {
        name <- paste("auc", x.axis, "all_com_vars.pdf", sep = "_")
      }
      
      if (x.axis == "com_vars") {
        name <- paste("auc", x.axis, "all_rew.pdf", sep = "_")
      }
      
      if (any(x.axis == c(rew_names, "norew", "abundtrait", "abundphylo"))) {
        name <- c(paste("auc", x.axis, x_axis, sep = "_"), ".pdf")
      }
      
      if (any(x.axis == c(names(ctrbs), "org"))) {
        name <- c(paste("auc", x.axis, x_axis, sep = "_"), ".pdf")
      }
      
    } else {
      size <- c(14.3, 10)
      name <- c(paste("auc", x.axis[1], x.axis[2], sep = "_"), ".pdf")
    }
  
  ggsave(name,
    path = paste0(getwd(), "/plot_sink"),
    plot = out,
    device = "pdf",
    width = size[1],
    height = size[2],
    units = "cm")
  }
  return(out)
}

# same as auc_boxplot but x.axis shows higher/lower
auc_boxplot2 <- function(x, by, log = F) {
  # extinction on lower or higher level ?
  lvl_match_extc <- grepl("lower", substitute(x))
  
  ifelse(lvl_match_extc, lvl <- "lower", lvl <- "higher")
  # 
  # if (is.null(x.axis) | !any(c(rew_names,
  #                              names(ctrbs),
  #                              "rew",
  #                              "com_vars",
  #                              "org") %in% x.axis)) {
  #   stop("invalid x axis ! Please provide a character vector of a variable that
  #   should be used as x axis. Either use rew or com_vars or a specific 
  #   value/combination from among them")
  # }
  
  # should y axis be logarithmic ?
  if (log) {
    scale <- "log10"
  } else {
    scale <- "none"
  }
  
  # filter data according to input
  if (length(by) <= 1) {
    if (by == "rew") {
      aucs <- x
      x_axis <- x.axis
      x_lab <- unique(x$rew)
    }
    
    if (x.axis == "com_vars") {
      aucs <- x
      x_axis <- x.axis
      x_lab <- unique(x$com_vars)
    }
    
    if (any(x.axis == rew_names)){
      aucs <- filter(x, rew == x.axis)
      x_axis <- "com_vars"
      x_lab <- unique(x$com_vars)
    }
    
    if (any(x.axis == c(names(ctrbs), "org"))){
      aucs <- filter(x, com_vars == x.axis)
      x_axis <- "rew"
      x_lab <- unique(x$rew)
    }
  } else {
    if (length(x.axis) > 2)
      stop("Too many values in x.axis. Please provide a character vector of 
           lenght one or two")
    
    rew_choice <- ifelse(any(x.axis[1] == rew_names), x.axis[1], x.axis[2])
    com_vars_choice <- ifelse(any(x.axis[1] == c(names(ctrbs), "org")),
                              x.axis[1], x.axis[2])
    
    if (any(rew_choice == rew_names) & any(com_vars_choice == rew_names) |
        any(rew_choice == c(names(ctrbs), "org")) & any(com_vars_choice == c(names(ctrbs), "org"))) {
      stop("Invalid combination of x.axis values. Did you provide two rew/com_vars values ?")
    }
    
    aucs <- filter(x, com_vars == com_vars_choice, rew == rew_choice)
    x_axis <- NULL
    x_lab <- paste(rew_choice, com_vars_choice)
    
  }
  
  if (length(x.axis) <= 1) {
    ggboxplot(y = "robustness",
              x = lvl,
              color = "by",
              data = aucs,
              order = x_lab,) +
      scale_color_manual(name = by,
                         values = c("#000000", "#787878")) +
      yscale(scale) +
      theme(legend.position = "bottom") +
      labs(x = element_blank()) +
      theme(aspect.ratio = 1)
  } else {
    ggboxplot(y = "robustness",
              x = "lvl", color = "rew",
              data = filter(auc_all, com_vars == "org")) +
      theme(legend.position = "bottom") +
      labs(x = element_blank()) +
      theme(aspect.ratio = 1)
  }
}