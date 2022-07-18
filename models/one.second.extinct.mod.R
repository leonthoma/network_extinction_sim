# Small change in the internal function "one.second.extinct" of the bipartite package. The change allows rewiring of interations after each step of extinction
# The arguments are the same of "one.second.extinct" function and additional arguments are described below.
# New arguments:
# rewiring - Logical argument to specify if allow rewiring (default rewiring = FALSE).
# probabilities.rewiring - A matrix with probabilities of rewiring, must be the same dimensions of the web (i.e. network). See section Methods in Vizentin-Bugoni et al. [in review] for details). This matrix is required in step ii of framework (default probabilities.rewiring = NULL).
# probabilities.rewiring2 - A matrix with probabilities of rewiring, must be the same dimensions of web. See section Methods in Vizentin-Bugoni et al. [in review] for details. This matrix is required in step iii of framework (default probabilities.rewiring2 = NULL).
# method.rewiring = Type of method used to trial rewiring, partial match to "one.try.single.partner", "multiple.trials.single.partner", "multiple.trials.multiples.partners", "one.try.each.partner" and "multiple.trials.each.partner". See section Methods in Vizentin-Bugoni et al. [in review] for details; (default method.rewiring = "one.try.single.partner").
# 
one.second.extinct.mod.aug <- function(web,
                                       participant = "higher",
                                       method = "abun",
                                       ext.row = NULL,
                                       ext.col = NULL, 
                                       rewiring = FALSE,
                                       abund.partner.choice = NULL,
                                       trait.partner.choice = NULL,
                                       phylo.partner.choice = NULL,
                                       interactions = NULL,
                                       method.rewiring = "NULL",
                                       make.bipartite = F,
                                       adapt = F,
                                       vis.steps = F,
                                       shift = F,
                                       coextc.thr = NULL,
                                       terminate_early = T,
                                       seed = F) {
  i <- 1L
  j <- 1L
  retry <- 1L # set counter for max retries for empty m2
  abort <- F # abort rewiring if extc sp only has dead interactions
  skip <- F # skip updating dead if extc sp only has dead interactions
  fail <- F # prevent Nobs from being 0
  
  dead <- matrix(nrow = 0, ncol = 5L)
  
  # filter sp w/o any interactions
  sp_low <- which(rowSums(web) != 0L)
  sp_high <- which(colSums(web) != 0L)

  # drop sp w/o any interactions from web
  m2 <- web[sp_low, sp_high, drop = F]
  
  # # set networks
  # m2 <- web
  
  # # match web and interactions
  # abund.partner.choice$low <- abund.partner.choice$low[names(abund.partner.choice$low) %in% rownames(web)]
  # abund.partner.choice$high <- abund.partner.choice$high[names(abund.partner.choice$high) %in% colnames(web)]
  # 
  # trait.partner.choice$low <- trait.partner.choice$low[rownames(trait.partner.choice$low) %in% rownames(web), ]
  # trait.partner.choice$high <- trait.partner.choice$high[rownames(trait.partner.choice$high) %in% colnames(web), ]
  # 
  # phylo.partner.choice$low <- phylo.partner.choice$low[rownames(phylo.partner.choice$low) %in% rownames(web),
  #                                                      rownames(phylo.partner.choice$low) %in% rownames(web)]
  # phylo.partner.choice$high <- phylo.partner.choice$high[colnames(phylo.partner.choice$high) %in% colnames(web),
  #                                                        colnames(phylo.partner.choice$high) %in% colnames(web)]

  # track no of shifts
  if (shift) {
    # create dfs to track number of shifts
    shifts_lower <- vector(mode = "list")
    
    shifts_higher <- vector(mode = "list")
  }
  
  dead <- cbind(0, 0, 0, nrow(m2), ncol(m2)) # correct for sp w/o any interactions
  
  colnames(dead) <- c("no", "ext.lower", "ext.higher", "n.lower", "n.higher")
  
  # create dummy if participant option both is used
  if (participant == "both") {
    partis <- participant
  } else {
    partis <- "single"
  }
  
  # set initial choice
  abund_choice_low <- NULL
  trait_choice_low <- NULL
  phylo_choice_low <- NULL
  
  abund_choice_high <- NULL
  trait_choice_high <- NULL
  phylo_choice_high <- NULL

  
  # determine coextinction of species if coextinction threshold is provided
  if (!is.null(coextc.thr)) {
    coextc <- T
  } else {
    coextc <- F
  }
  
  # adapatability
  if (is.list(adapt)) {
  adapt_list <- list("low" = adapt$low,
                     "high" = adapt$high)
  } else {
    adapt_list <- list("low" = rep(.5, nrow(m2)),
                  "high" = rep(.5, ncol(m2)))
  }

  if(!any(c("NULL", "abund", "phylo", "trait") %in% method.rewiring)) {
    stop("\n Invalid rewiring method. Choose either NULL, abund, phylo or trait \n")
  }
  if(rewiring){
    if(any(web%%1!=0)){
      stop("\n If rewiring is TRUE the web must must contain only integers \n")
    }
    if(is.null(rownames(web)) | is.null(colnames(web))){
      stop("\n If rewiring is TRUE the web must contain rownames and colnames\n")
    }
  }
  repeat {
    ext_temp <- extinction.mod(m2, participant = participant, method = method, ext.row = ext.row, ext.col = ext.col)
    # extinction.mod returns list w/ plants (rows; "rexcl") that lost partners and
    # animals (cols; "cexcl") that lost partners as well as network ("web")
    
    if (sum(ext_temp$web) == 0) {
      fail <- T
      Nobs_fail <- sum(m2)
    }
    
    if (fail) {
      Nobs <- Nobs_fail
    } else {
        Nobs <- sum(ext_temp$web) # no of obs interaction
      }

    # handle edge case where extc. sp has only dead interactions
    if (sum(ext_temp$rexcl, ext_temp$cexcl) == 0) {
      j <- 1
      warning("Extinct sp only has dead interactions, trying to find alternative")
      # should extinction simulations be terminated when trying to select valid
      # extinct species too often
      while (j <= 3) {
        # add max number of iterations before break; if still not valid stop sim and restart with initial imat
        ext_temp <- extinction.mod(m2, participant = participant, method = method, ext.row = ext.row, ext.col = ext.col)
          
        j <- j + 1
  
        if (sum(ext_temp$rexcl, ext_temp$cexcl) != 0){
          break
        }
          
        if (j == 3) {
          warning("Failed to find alternative, skipping current simulation step(", i, ")")
          abort <- T
          # skip <- T
          break
          }
        }
      }
    
    if (!abort
        ) {
      if (rewiring) {
        if (!is.null(ext_temp$rexcl)){ # Plant is extinct, looking for new interaction partners of birds that interacted with lost plant
          sp_ext <- rownames(ext_temp$rexcl) # name of extc. plant
          sp_try_rewiring <- colnames(ext_temp$rexcl)[which(ext_temp$rexcl>0)] # number & position of interaction partners (higher trophic level)
          
          # overwrite participant if using option both
          if (partis == "both")
            participant <- "lower"
          
          # Choice of rewiring partner
          if (any(method.rewiring == "abund")) {
            if (i == 1 | is.null(abund_choice_low)) {
              abund_choice_low <- abund.partner.choice$low # get values
              }
  
            abund_choice_low <- abund_choice_low[which(names(abund_choice_low) %in% rownames(interactions)), drop = F] # drop sp w/o any interactions
            
            #abund_sp_ext_idx <- which(names(abund_choice_low) %in% sp_ext) # get idx of sp.ext (NOT THE SAME IN ext.temp$web !!)
            
            # get sp w/ highest abund
            abund_rew_partner <- names(which.max(abund_choice_low[! names(abund_choice_low) %in% sp_ext, drop = F]))
            
            #abund_rew_partner <- which(rownames(m2) %in% sp_high_abund)
            abund_choice_low_tmp <- abund_choice_low[! names(abund_choice_low) %in% sp_ext, drop = F] # delete extc sp
          }
          
          # choose rewiring partner based on most similar trait
          if (any(method.rewiring == "trait")) {
            if (i == 1 | is.null(trait_choice_low)) {
              trait_choice_low <- trait.partner.choice$low # get values
            }
            
            trait_choice_low <- trait_choice_low[which(rownames(trait_choice_low) %in% rownames(interactions)), , drop = F] # drop sp w/o any interactions
    
            #trait_sp_ext_idx <- which(rownames(trait_choice_low) %in% sp_ext) # get idx of sp.ext (NOT THE SAME IN ext.temp$web !!)
            
            # calculate euclidean distances
            trait_dist <- as.matrix(dist(trait_choice_low, diag = T, upper = T)) 
            
            # find sp w/ smallest trait dist
            trait_rew_partner <- names(which(sort(trait_dist[, sp_ext, drop = F])[2] == trait_dist[, sp_ext])) # use 2nd lowest since sp.ext is lowest
              
            # match interacting sp of closest relative(check.int) to remaining sp (sp.rem)
            #trait_rew_partner <- which(rownames(m2) %in% sp_closest)
            trait_choice_low_tmp <- trait_choice_low[! rownames(trait_choice_low) %in% sp_ext, , drop = F] # delete extc sp
          }
          
          # choose rewiring partner based on closest phylogenetic distance
          if (any(method.rewiring == "phylo")) {
            if (i == 1 | is.null(phylo_choice_low)) {
              phylo_choice_low <- phylo.partner.choice$low # get values
            }
            
            phylo_choice_low <- phylo_choice_low[which(rownames(phylo_choice_low) %in% rownames(interactions)), which(colnames(phylo_choice_low) %in% rownames(interactions)), drop = F] # drop sp w/o any interactions
            
            #phylo_sp_ext_idx <- which(rownames(phylo_choice_low) %in% sp_ext) # get idx of sp.ext (NOT THE SAME IN ext.temp$web !!)
            
            # get closest relative
            sp_close_rel_tmp <- names(which(sort(phylo_choice_low[, sp_ext])[2] == phylo_choice_low[, sp_ext])) # use 2nd lowest since sp.ext is lowest
            
            # if multiple species have same distance, randomly choose one
            if (length(sp_close_rel_tmp) != 1){
              phylo_rew_partner <- sample(sp_close_rel_tmp, 1)
            } else {
              phylo_rew_partner <- sp_close_rel_tmp
            }
            
            #phylo_rew_partner <- which(rownames(m2) %in% sp_close_rel)
            phylo_choice_low_tmp <- phylo_choice_low[! rownames(phylo_choice_low) %in% sp_ext, ! rownames(phylo_choice_low) %in% sp_ext, drop = F] # delete extc sp
          }
          
          # Shift interaction probability from old interaction (sp.ext + sp.try.rewiring) to new interaction (rew.partner + sp.try.rewiring)
          #sp_ext_idx <- which(rownames(interactions) %in% sp_ext) # get idx for i_mat
          
          # Update I_mat
          if (any(method.rewiring == "abund")) {
          abund_prob_int_old <- interactions[sp_ext, sp_try_rewiring] # get old interaction vals
          interactions_abund <- interactions[abund_rew_partner, sp_try_rewiring] + (.5 * abund_prob_int_old) # set new interactions vals  
          }  
          
          if (any(method.rewiring == "trait")) {
          trait_prob_int_old <- interactions[sp_ext, sp_try_rewiring] # get old interaction vals
          interactions_trait <- interactions[trait_rew_partner, sp_try_rewiring] + (.5 * trait_prob_int_old) # set new interactions vals
          }
          
          if (any(method.rewiring == "phylo")) {
          phylo_prob_int_old <- interactions[sp_ext, sp_try_rewiring] # get old interaction vals
          interactions_phylo <- interactions[phylo_rew_partner, sp_try_rewiring] + (.5 * phylo_prob_int_old) # set new interactions vals
          }
          
          if (length(method.rewiring) == 1) {
            if (any(method.rewiring == "abund")) {
              interactions[abund_rew_partner, sp_try_rewiring] <- interactions_abund
            }
            if (any(method.rewiring == "trait")) {
              interactions[trait_rew_partner, sp_try_rewiring] <- interactions_trait
            }
            if (any(method.rewiring == "phylo")) {
              interactions[phylo_rew_partner, sp_try_rewiring] <- interactions_phylo
            }
          } else {
            if (all(method.rewiring == c("abund", "trait")) | all(method.rewiring == c("trait", "abund"))) {
              if (abund_rew_partner == trait_rew_partner) {
                interactions[abund_rew_partner, sp_try_rewiring] <- interactions_abund + interactions_trait
              } else {
                interactions[abund_rew_partner, sp_try_rewiring] <- interactions_abund
                interactions[trait_rew_partner, sp_try_rewiring] <- interactions_trait
              }
            }
            if (all(method.rewiring == c("abund", "phylo")) | all(method.rewiring == c("phylo", "abund"))) {
              if (abund_rew_partner == phylo_rew_partner) {
                interactions[abund_rew_partner, sp_try_rewiring] <- interactions_abund + interactions_phylo
              } else {
                interactions[abund_rew_partner, sp_try_rewiring] <- interactions_abund
                interactions[phylo_rew_partner, sp_try_rewiring] <- interactions_phylo
              }
            }
            
            # if (trait_rew_partner == phylo_rew_partner) {
            #   interactions[trait_rew_partner, sp_try_rewiring] <- interactions_trait + interactions_phylo
            # }
            
          }
        }
        
        if(!is.null(ext_temp$cexcl)){ # Bird is extinct, looking for new interaction partners of plants that interacted with lost bird
          sp_ext <- colnames(ext_temp$cexcl)  # name of extc. bird
          sp_try_rewiring <- rownames(ext_temp$cexcl)[which(ext_temp$cexcl>0)] # number & position of interaction partners (lower trophic level)
          
          # overwrite participant if using option both
          if (partis == "both")
            participant <- "higher"
          
          # Choice of rewiring partner
          if (any(method.rewiring == "abund")) {
            if (i == 1 | is.null(abund_choice_high)) {
              abund_choice_high <- abund.partner.choice$high # get values
            }
            
            abund_choice_high <- abund_choice_high[which(names(abund_choice_high) %in% colnames(interactions)), drop = F] # drop sp w/o any interactions
            
            #abund_sp_ext_idx <- which(names(abund_choice_high) %in% sp_ext) # get idx of sp.ext (NOT THE SAME IN ext.temp$web !!)
            
            # get sp w/ highest abund
            abund_rew_partner <- names(which.max(abund_choice_high[! names(abund_choice_high) %in% sp_ext, drop = F]))
            
            #abund_rew_partner <- which(colnames(m2) %in% sp_high_abund)
            abund_choice_high_tmp <- abund_choice_high[! names(abund_choice_high) %in% sp_ext, drop = F] # delete extc sp
          }
          
          # choose rewiring partner based on most similar trait
          if (any(method.rewiring == "trait")) {
            if (i == 1 | is.null(trait_choice_high)) {
              trait_choice_high <- trait.partner.choice$high # get values
            }
            
            trait_choice_high <- trait_choice_high[which(rownames(trait_choice_high) %in% colnames(interactions)), , drop = F] # drop sp w/o any interactions
            
            #trait_sp_ext_idx <- which(rownames(trait_choice_high) %in% sp_ext) # get idx of sp.ext (NOT THE SAME IN ext.temp$web !!)
            
            # calculate euclidean distances
            trait_dist <- as.matrix(dist(trait_choice_high, diag = T, upper = T)) 
            
            # find sp w/ smallest trait dist
            trait_rew_partner <- names(which(sort(trait_dist[, sp_ext, drop = F])[2] == trait_dist[, sp_ext])) # use 2nd lowest since sp.ext is lowest
            
            # match interacting sp of closest relative(check.int) to remaining sp (sp.rem)
            #trait_rew_partner <- which(colnames(m2) %in% sp_closest)
            trait_choice_high_tmp <- trait_choice_high[! rownames(trait_choice_high) %in% sp_ext, , drop = F] # delete extc sp
          }
          
          # choose rewiring partner based on closest phylogenetic distance
          if (any(method.rewiring == "phylo")) {
            if (i == 1 | is.null(phylo_choice_high)) {
              phylo_choice_high <- phylo.partner.choice$high # get values
            }
            
            phylo_choice_high <- phylo_choice_high[which(rownames(phylo_choice_high) %in% colnames(interactions)), which(colnames(phylo_choice_high) %in% colnames(interactions)), drop = F] # drop sp w/o any interactions
            
            
            #phylo_sp_ext_idx <- which(rownames(phylo_choice_high) %in% sp_ext) # get idx of sp.ext (NOT THE SAME IN ext.temp$web !!)
            
            # get closest relative
            sp_close_rel_tmp <- names(which(sort(phylo_choice_high[, sp_ext])[2] == phylo_choice_high[, sp_ext])) # use 2nd lowest since sp.ext is lowest 
            
            # if multiple species have same distance, randomly choose one
            if (length(sp_close_rel_tmp) != 1){
              phylo_rew_partner <- sample(sp_close_rel_tmp, 1)
            } else {
              phylo_rew_partner <- sp_close_rel_tmp
            }
            
            #phylo_rew_partner <- which(colnames(m2) %in% sp_close_rel)
            phylo_choice_high_tmp <- phylo_choice_high[! rownames(phylo_choice_high) %in% sp_ext, ! rownames(phylo_choice_high) %in% sp_ext, drop = F] # delete extc sp
          }
        
        # Shift interaction probability from old interaction (sp.ext + sp.try.rewiring) to new interaction (rew.partner + sp.try.rewiring)
        #sp_ext_idx <- which(colnames(interactions) %in% sp_ext) # get idx for i_mat
        
        # Update I_mat
        if (any(method.rewiring == "abund")) {
        abund_prob_int_old <- interactions[sp_try_rewiring, sp_ext] # get old interaction vals
        interactions_abund <- interactions[sp_try_rewiring, abund_rew_partner] + (.5 * abund_prob_int_old) # set new interactions vals
        }
        
        if (any(method.rewiring == "trait")) {
          trait_prob_int_old <- interactions[sp_try_rewiring, sp_ext] # get old interaction vals
          interactions_trait <- interactions[sp_try_rewiring, trait_rew_partner] + (.5 * trait_prob_int_old) # set new interactions vals
        }
        
        if (any(method.rewiring == "phylo")) {
          phylo_prob_int_old <- interactions[sp_try_rewiring, sp_ext] # get old interaction vals
          interactions_phylo <- interactions[sp_try_rewiring, phylo_rew_partner] + (.5 * phylo_prob_int_old) # set new interactions vals
        }
        
        if (length(method.rewiring) == 1) {
          if (any(method.rewiring == "abund")) {
            interactions[sp_try_rewiring, abund_rew_partner] <- interactions_abund
          }
          if (any(method.rewiring == "trait")) {
            interactions[sp_try_rewiring, trait_rew_partner] <- interactions_trait
          }
          if (any(method.rewiring == "phylo")) {
            interactions[sp_try_rewiring, phylo_rew_partner] <- interactions_phylo
          }
        } else {
          if (all(method.rewiring == c("abund", "trait")) | all(method.rewiring == c("trait", "abund"))) {
            if (abund_rew_partner == trait_rew_partner) {
              interactions[sp_try_rewiring, abund_rew_partner] <- interactions_abund + interactions_trait
            } else {
              interactions[sp_try_rewiring, abund_rew_partner] <- interactions_abund
              interactions[sp_try_rewiring, trait_rew_partner] <- interactions_trait
            }
          }
          
          if (all(method.rewiring == c("abund", "phylo")) | all(method.rewiring == c("phylo", "abund"))) {
            if (abund_rew_partner == phylo_rew_partner) {
              interactions[sp_try_rewiring, abund_rew_partner] <- interactions_abund + interactions_phylo
            } else {
              interactions[sp_try_rewiring, abund_rew_partner] <- interactions_abund
              interactions[sp_try_rewiring, phylo_rew_partner] <- interactions_phylo
            }
          }
          
          # if (trait_rew_partner == phylo_rew_partner) {
          #   interactions[sp_try_rewiring, trait_rew_partner] <- interactions_trait + interactions_phylo
          # }
          
        }
        
        }
      }
    }
    # removed rows and cols
    if (coextc) {
      if (nrow(ext_temp$web) >= 2L & ncol(ext_temp$web) >= 2L) {
        # sort ext.temp web and m2
        ext_temp$web <- ext_temp$web[sort(rownames(ext_temp$web)),
                     sort(colnames(ext_temp$web))]
        
        m2 <- m2[sort(rownames(m2)), sort(colnames(m2))]
        
        # percentage of interactions remaining
        
        ## !!! FIX If both are 0 -> NAN If m2 is 0 -> Inf FIX !!!!
        coext_r <- which(rowSums(ext_temp$web)/ rowSums(m2) <= coextc.thr)
        coext_c <- which(colSums(ext_temp$web)/ colSums(m2) <= coextc.thr)
      
        # set interaction to 0 if coextinction threshold was exceeded
        if (!length(coext_r) == 0) {
          ext_temp$web[coext_r, ] <- 0
        }
        
        if (!length(coext_c) == 0) {
        ext_temp$web[, coext_c] <- 0
        }
      }
    } 
    
    rem_r_c <- empty(ext_temp$web, count = T) 
    
    
    # get extc. sp if no rewiring
    if (!rewiring) {
      if (participant == "lower") {
        sp_ext <- rownames(ext_temp$rexcl)
      } else {
        sp_ext <- colnames(ext_temp$cexcl)
      }
    }
    
    # check if choice is updated correctly
    if(rewiring) {
      # if (abort) {
      #   if(!is.null(abund_choice_low) & !is.null(trait_choice_low))
      #     abund_choice_low <- abund_choice_low # don't delete sp.ext if it had only dead interactions
      #     trait_choice_low <- trait_choice_low # don't delete sp.ext if it had only dead interactions
      #     phylo_choice_low <- phylo_choice_low # don't delete sp.ext if it had only dead interactions
      #   if(!is.null(choice_high))
      #     choice_high <- choice_high 
      #   } else {
          if (!is.null(abund_choice_low))
            abund_choice_low <- abund_choice_low_tmp # delete sp.ext
          
          if (!is.null(trait_choice_low))  
          trait_choice_low <- trait_choice_low_tmp # delete sp.ext
          
          if (!is.null(phylo_choice_low))  
            phylo_choice_low <- phylo_choice_low_tmp # delete sp.ext
          
          if(!is.null(abund_choice_high))
            abund_choice_high <- abund_choice_high_tmp
          
          if(!is.null(trait_choice_high))
            trait_choice_high <- trait_choice_high_tmp
          
          if(!is.null(phylo_choice_high))
            phylo_choice_high <- phylo_choice_high_tmp
        # }
      }
    
    # remaining sp
    rem_low <- as.vector(attributes(rem_r_c)$dimnames[[1]])
    rem_high <- as.vector(attributes(rem_r_c)$dimnames[[2]])
  
    # compare dead interaction sp with extc sp
    if (abort){ # if sp.ext has only dead interactions delete no species
      ext_low <- 0L
      ext_high <- 0L
      skip <- T
    } else {
      if (i > 1) {
        if (length(dead_low) != 0) {
          ext_low <- length(rownames(m2[-which(rownames(m2) %in% names(dead_low)), , drop = F])) - length(rem_low)
        } else {
          ext_low <- attributes(rem_r_c)$empty[1]
        }
        if (length(dead_high) != 0) {
          ext_high <- length(colnames(m2[, -which(colnames(m2) %in% names(dead_high)), drop = F])) - length(rem_high)
        } else {
          ext_high <- attributes(rem_r_c)$empty[2]
        }
      } else {
        ext_low <- attributes(rem_r_c)$empty[1]
        ext_high <- attributes(rem_r_c)$empty[2]
      }
    }
    
    # track no of extc and no of remaining sp, skip if only dead interactions
    if (!skip) {
      if (i == 1) {
        dead <- rbind(dead, c(i, ext_low, ext_high,
                              nrow(m2) - ext_low,
                              ncol(m2) - ext_high))  
      } else {
        dead <- rbind(dead, c(i, ext_low, ext_high,
                              nrow(m2) - ext_low,
                              ncol(m2) - ext_high))
      }
    }
    
    if (vis.steps) {
      par(mfrow = c(2, 1))
      plotweb(m2)
      visweb(m2, text = "interaction", type = "nested")
    }
  
    # update I_mat
    if (!abort) {
      if (participant == "lower") {
        interactions <- interactions[-which(rownames(interactions) == sp_ext), , drop = F]
      }
      
      if (participant == "higher") {
        interactions <- interactions[, -which(colnames(interactions) == sp_ext), drop = F]
      }
    } else {
      retry <- retry + 1L
      if (terminate_early) {
        if (retry > 3L) {
          warning("Max no. of retries for dead interactions reached, aborting")
          break
        }
      }
    }
    
    if (!abort) {
    drop_rows <- which(!(rownames(interactions) %in% rownames(rem_r_c)))
    drop_cols <- which(!(colnames(interactions) %in% colnames(rem_r_c)))
    } else {
      drop_rows <- numeric(0L)
      drop_cols <- numeric(0L)
    }

    if (i == 1L) {
      dead_low <- numeric(0L)
      dead_high <- numeric(0L)
    } # set initial values
    
    retain_low <- which(rownames(interactions) %in% names(dead_low))
    retain_high <- which(colnames(interactions) %in% names(dead_high))
    
    if (!length(retain_low) == 0L) {
      drop_rows <- drop_rows[-which(drop_rows %in% retain_low), drop = F]
    } # retain dead interactions
    
    if (!length(drop_rows) == 0L) {
    interactions <- interactions[-drop_rows, , drop = F]
    }

    if (!length(retain_high) == 0L) {
      drop_cols <- drop_cols[-which(drop_cols %in% retain_high), drop = F]
    } # retain dead interactions
    
    if (!length(drop_cols) == 0L) {
      interactions <- interactions[, -drop_cols, drop = F]
    }
    
    # break if n.lower/higher is 0
    if ((tail(dead, 1)[, "n.lower"] == 0 | tail(dead, 1)[, "n.higher"] == 0)) {
      dead <- dead[-nrow(dead), , drop = F]
      break
      }
    
    # get current partners to track shifts
    if (shift) {
      # get partners of lower level sp
      lower_partners_old <- map(seq(nrow(m2)), ~ names(which(m2[.x, ] >= 1) == T))
      names(lower_partners_old) <- rownames(m2)
      
      # get partners of higher level sp
      higher_partners_old <- map(seq(ncol(m2)), ~ names(which(m2[, .x] >= 1) == T))
      names(higher_partners_old) <- colnames(m2)
    }
    
    # renormalize interactions
    interactions <- interactions/sum(interactions)

    # calculate new web w/ updated o_mat
    if (nrow(interactions) >= 2L | ncol(interactions) >= 2L) {
      if (seed) {
        set.seed(420)
      }
      m2 <- matrix(rmultinom(1, Nobs, interactions), nrow = nrow(interactions),
                  ncol = ncol(interactions))
      dimnames(m2) <- dimnames(interactions)
      
      # get dead interactions
      dead_low <- which(rowSums(m2) == 0L)
      dead_high <- which(colSums(m2) == 0L)
      
    } else {
      break
    }
    
    # update partners for tracking shifts
    if (shift) {
      # get partners of lower level sp
      lower_partners_new <- map(seq(nrow(m2)), ~ names(which(m2[.x, ] >= 1) == T))
      names(lower_partners_new) <- rownames(m2)
      
      # get partners of higher level sp
      higher_partners_new <- map(seq(ncol(m2)), ~ names(which(m2[, .x] >= 1) == T))
      names(higher_partners_new) <- colnames(m2)
      
      # new interaction partners 
      added_lower <- map(names(lower_partners_old), function (x) {
        length(which(lower_partners_new[[x]] %in% lower_partners_old[[x]] == F))
      })
      names(added_lower) <- names(lower_partners_old)
      
      added_higher <- map(names(higher_partners_old), function(x) {
        length(which(higher_partners_new[[x]] %in% higher_partners_old[[x]] == F))
        })
      names(added_higher) <- names(higher_partners_old)
      
      shifts_lower[[i]] <- unlist(added_lower)
      shifts_higher[[i]] <- unlist(added_higher)
      
    }
    
    if (participant == "lower" & NROW(m2) < 2L) 
      break
    if (participant == "higher" & NCOL(m2) < 2L) 
      break
    if (participant == "both" & min(dim(m2)) < 2L) 
      break
    if (any(dim(ext_temp$web) == 1L)) 
      break
    
    # "Skip" iteration if extc sp only had dead interactions 
    if (!skip) {
      i <- i + 1L
    }
    
    # reset
    skip <- F
    abort <- F
    
    # reset participant overwrite; only used when option is both
    if (partis == "both")
      participant <- "both" 
  }

  # Add last line of dead
  if (abort & terminate_early) {
    dead2 <- NULL
  } else {
    dead2 <- rbind(dead, c(NROW(dead), tail(dead, 1)[, "n.lower"], tail(dead, 1)[, "n.higher"], 0L, 0L)) # add last line of dead (i.e. all sp extc)
  }

  if (participant == "lower" & method == "degree") {
    if (length(table(dead[, 2])) > 1) 
      dead2[, 2] <- 1
  }
  
  # if (nrow(dead) + 1 != nrow(dead2)) 
  #   stop("PANIC! Something went wrong with the extinct sequence! Please contact the author to fix this!!")
  # 
  if (make.bipartite & !abort) {
    out <- dead2
    class(out) <- "bipartite"
    attr(out, "exterminated") <- c("both", "lower", "higher")[pmatch(participant, c("both", "lower", "higher"))]
    attr(out, "exterminated")
  } else {
    out <- dead2
  }
  if (shift) {
    shifts <<- list("lower" = shifts_lower, "higher" = shifts_higher)
  }
  return(out)
}

