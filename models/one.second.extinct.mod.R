# Small change in the internal function "one.second.extinct" of the bipartite package. The change allows rewiring of interations after each step of extinction
# The arguments are the same of "one.second.extinct" function and additional arguments are described below.
# New arguments:
# rewiring - Logical argument to specify if allow rewiring (default rewiring = FALSE).
# probabilities.rewiring - A matrix with probabilities of rewiring, must be the same dimensions of the web (i.e. network). See section Methods in Vizentin-Bugoni et al. [in review] for details). This matrix is required in step ii of framework (default probabilities.rewiring = NULL).
# probabilities.rewiring2 - A matrix with probabilities of rewiring, must be the same dimensions of web. See section Methods in Vizentin-Bugoni et al. [in review] for details. This matrix is required in step iii of framework (default probabilities.rewiring2 = NULL).
# method.rewiring = Type of method used to trial rewiring, partial match to "one.try.single.partner", "multiple.trials.single.partner", "multiple.trials.multiples.partners", "one.try.each.partner" and "multiple.trials.each.partner". See section Methods in Vizentin-Bugoni et al. [in review] for details; (default method.rewiring = "one.try.single.partner").
# 
# method <-  "random"
# participant <- "higher"
# partner.choice <- phylos
# method.rewiring <- "phylo"
# interactions <- sims$aTl$I_mat
# web <- sims$aTl$web

# one.second.extinct.mod_aug(web = web,
#                            participant = participant,
#                            method = "random",
#                            rewiring = T,
#                            partner.choice = partner.choice,
#                            interactions = interactions,
#                            method.rewiring = method.rewiring)



one.second.extinct.mod_aug <- function(web,
                                       participant = "higher",
                                       method = "abun",
                                       ext.row = NULL,
                                       ext.col = NULL, 
                                       rewiring = FALSE,
                                       partner.choice,
                                       interactions,
                                       method.rewiring = NULL,
                                       make.bipartite = F,
                                       adapt = T,
                                       vis.steps = F) {
  i <- 1L
  j <- 1L
  retry <- 1L # set counter for max retries for empty m2
  abort <- F # abort rewiring if extc sp only has dead interactions
  skip <- F # skip updating dead if extc sp only has dead interactions
  
  dead <- matrix(nrow = 0, ncol = 5L)
  
  # filter sp w/o any interactions
  sp_low <- which(rowSums(web) != 0L)
  sp_high <- which(colSums(web) != 0L)

  # drop sp w/o any interactions from web
  m2 <- web[sp_low, sp_high]
  
  ## Debugging
  # m2 <- web
  
  dead <- cbind(0, 0, 0, nrow(m2), ncol(m2)) # correct for sp w/o any interactions
  
  colnames(dead) <- c("no", "ext.lower", "ext.higher", "n.lower", "n.higher")
  
  
  # get interaction vals & drop sp w/o any interactions
  # interactions <- interactions[which(rownames(interactions) %in% names(sp_low)),
  #                              which(colnames(interactions) %in% names(sp_high))]

  # adapatability
  if (adapt == T) {
  adapt_list <- list("low" = runif(nrow(m2), 0, 1),
                "high" = runif(ncol(m2), 0, 1))
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
    ext.temp <- extinction.mod(m2, participant = participant, method = method, ext.row = ext.row, ext.col = ext.col)
    # extinction.mod returns list w/ plants (rows; "rexcl") that lost partners and
    # animals (cols; "cexcl") that lost partners as well as network ("web")

    # handle edge case where extc. sp has only dead interactions
    if (sum(ext.temp$rexcl, ext.temp$cexcl) == 0) {
      j <- 1
      warning("Extinct sp only has dead interactions, trying to find alternative")
      while (j <= 3) {

        # add max number of iterations before break; if still not valid stop sim and restart with initial imat
        ext.temp <- extinction.mod(m2, participant = participant, method = method, ext.row = ext.row, ext.col = ext.col)
        
        j <- j + 1

        if (sum(ext.temp$rexcl, ext.temp$cexcl) != 0)
          break
      }
      if (j > 3) {
        warning("Failed to find alternative, skipping current simulation step(", i, ")")
         abort <- T
        # skip <- T
      }
    }
    if (!isTRUE(abort)) {
      if (rewiring) {
        if(!is.null(ext.temp$rexcl)){ # Plant is extinct, looking for new interaction partners of birds that interacted with lost plant
          sp.ext <- rownames(ext.temp$rexcl) # name of extc. plant
          sp.try.rewiring <- which(ext.temp$rexcl>0) # number & position of interaction partners (higher trophic level)
          Nobs <- sum(m2)# - sum(m2[sp.ext, ]) # no of obs interactions
          
          # Choice of rewiring partner
          if (method.rewiring == "abund") {
            if (i == 1) {
              choice <- partner.choice$low # get values
              }
  
            choice <- choice[which(names(choice) %in% rownames(interactions))] # drop sp w/o any interactions
            
            sp.ext.idx <- which(names(choice) %in% sp.ext) # get idx of sp.ext (NOT THE SAME IN ext.temp$web !!)
            
            # get sp w/ highest abund
            sp.high.abund <- names(which.max(choice[-sp.ext.idx]))
            
            rew.partner <- which(rownames(m2) %in% sp.high.abund)
            choice_tmp <- choice[-sp.ext.idx] # delete extc sp
          }
          
          # choose rewiring partner based on most similar trait
          if (method.rewiring == "trait") {
            if (i == 1) {
              choice <- partner.choice$low # get values
            }
            
            choice <- choice[which(rownames(choice) %in% rownames(interactions)), , drop = F] # drop sp w/o any interactions
    
            sp.ext.idx <- which(rownames(choice) %in% sp.ext) # get idx of sp.ext (NOT THE SAME IN ext.temp$web !!)
            
            # calculate euclidean distances
            trait.dist <- as.matrix(dist(choice, diag = T, upper = T)) 
            
            # find sp w/ smallest trait dist
            sp.closest <- names(which.min(trait.dist[-sp.ext.idx, sp.ext.idx]))
              
            # match interacting sp of closest relative(check.int) to remaining sp (sp.rem)
            rew.partner <- which(rownames(m2) %in% sp.closest)
            choice_tmp <- choice[-sp.ext.idx, ] # delete extc sp
          }
          
          # choose rewiring partner based on closest phylogenetic distance
          if (method.rewiring == "phylo") {
            if (i == 1) {
              choice <- partner.choice$low # get values
            }
            
            choice <- choice[which(rownames(choice) %in% rownames(interactions)), which(colnames(choice) %in% rownames(interactions)), drop = F] # drop sp w/o any interactions
            
            sp.ext.idx <- which(rownames(choice) %in% sp.ext) # get idx of sp.ext (NOT THE SAME IN ext.temp$web !!)
            
            # get closest relative
            close.rel <-  min(choice[-sp.ext.idx, sp.ext.idx]) # get distance
            sp.close.rel.tmp <- names(which(choice[, sp.ext.idx] == close.rel)) # closest sp
            
            # if multiple species have same distance, randomly choose one
            if (length(sp.close.rel.tmp) != 1){
              sp.close.rel <- sample(sp.close.rel.tmp, 1)
            } else {
              sp.close.rel <- sp.close.rel.tmp
            }
            
            rew.partner <- which(rownames(m2) %in% sp.close.rel)
            choice_tmp <- choice[-sp.ext.idx, -sp.ext.idx, drop = F] # delete extc sp
          }
          
          # Shift interaction probability from old interaction (sp.ext + sp.try.rewiring) to new interaction (rew.partner + sp.try.rewiring)
          sp.ext.idx <- which(rownames(interactions) %in% sp.ext) # get idx for i_mat
          
          # Update I_mat
            prob.int.old <- interactions[sp.ext.idx, sp.try.rewiring] # get old interaction vals
            interactions[rew.partner, sp.try.rewiring] <- interactions[rew.partner, sp.try.rewiring] + (adapt_list$high[sp.try.rewiring] * prob.int.old) # set new interactions vals
          }
        
        if(!is.null(ext.temp$cexcl)){ # Bird is extinct, looking for new interaction partners of plants that interacted with lost bird
          sp.ext <- colnames(ext.temp$cexcl)  # name of extc. bird
          sp.try.rewiring <- which(ext.temp$cexcl>0) # number & position of interaction partners (lower trophic level)
          Nobs <- sum(m2)# - sum(m2[, sp.ext]) # no of obs interactions
          
          # Choice of rewiring partner
          if (method.rewiring == "abund") {
            if (i == 1) {
              choice <- partner.choice$high # get values
            }
            
            choice <- choice[which(names(choice) %in% colnames(interactions))] # drop sp w/o any interactions
            
            sp.ext.idx <- which(names(choice) %in% sp.ext) # get idx of sp.ext (NOT THE SAME IN ext.temp$web !!)
            
            # get sp w/ highest abund
            sp.high.abund <- names(which.max(choice[-sp.ext.idx]))
            
            rew.partner <- which(colnames(m2) %in% sp.high.abund)
            choice_tmp <- choice[-sp.ext.idx] # delete extc sp
          }
          
          # choose rewiring partner based on most similar trait
          if (method.rewiring == "trait") {
            if (i == 1) {
              choice <- partner.choice$high # get values
            }
            
            choice <- choice[which(rownames(choice) %in% colnames(interactions)), , drop = F] # drop sp w/o any interactions
            
            sp.ext.idx <- which(rownames(choice) %in% sp.ext) # get idx of sp.ext (NOT THE SAME IN ext.temp$web !!)
            
            # calculate euclidean distances
            trait.dist <- as.matrix(dist(choice, diag = T, upper = T)) 
            
            # find sp w/ smallest trait dist
            sp.closest <- names(which.min(trait.dist[-sp.ext.idx, sp.ext.idx]))
            
            # match interacting sp of closest relative(check.int) to remaining sp (sp.rem)
            rew.partner <- which(colnames(m2) %in% sp.closest)
            choice_tmp <- choice[-sp.ext.idx, ] # delete extc sp
          }
          
          # choose rewiring partner based on closest phylogenetic distance
          if (method.rewiring == "phylo") {
            if (i == 1) {
              choice <- partner.choice$high # get values
            }
            
            choice <- choice[which(rownames(choice) %in% colnames(interactions)), which(colnames(choice) %in% colnames(interactions)), drop = F] # drop sp w/o any interactions
            
            sp.ext.idx <- which(rownames(choice) %in% sp.ext) # get idx of sp.ext (NOT THE SAME IN ext.temp$web !!)
            
            # get closest relative
            close.rel <-  min(choice[-sp.ext.idx, sp.ext.idx]) # get distance
            sp.close.rel.tmp <- names(which(choice[, sp.ext.idx] == close.rel)) # closest sp
            
            # if multiple species have same distance, randomly choose one
            if (length(sp.close.rel.tmp) != 1){
              sp.close.rel <- sample(sp.close.rel.tmp, 1)
            } else {
              sp.close.rel <- sp.close.rel.tmp
            }
            
            rew.partner <- which(colnames(m2) %in% sp.close.rel)
            choice_tmp <- choice[-sp.ext.idx, -sp.ext.idx, drop = F] # delete extc sp
          }
        
        # Shift interaction probability from old interaction (sp.ext + sp.try.rewiring) to new interaction (rew.partner + sp.try.rewiring)
        sp.ext.idx <- which(colnames(interactions) %in% sp.ext) # get idx for i_mat
        
        # Update I_mat
        prob.int.old <- interactions[sp.try.rewiring, sp.ext.idx] # get old interaction vals
        interactions[sp.try.rewiring, rew.partner] <- interactions[sp.try.rewiring, rew.partner] + (adapt_list$low[sp.try.rewiring] * prob.int.old) # set new interactions vals
        }
      }
    }
    # removed rows and cols
    rem_r_c <- empty(ext.temp$web, count = T)
    
    # check if choice is updated correctly
    if (isTRUE(abort)) {
      choice <- choice # don't delete sp.ext if it had only dead interactions
      } else {
        choice <- choice_tmp # delete sp.ext
        }
    
    # remaining sp
    rem_low <- as.vector(attributes(rem_r_c)$dimnames[[1]])
    rem_high <- as.vector(attributes(rem_r_c)$dimnames[[2]])
  
    # compare dead interaction sp with extc sp
    if (isTRUE(abort)){ # if sp.ext has only dead interactions delete no species
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
    if (!isTRUE(skip)) {
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
    
    if (vis.steps == T) {
      plotweb(m2)
    }
  
    # update I_mat
    if (!isTRUE(abort)) {
      if (participant == "lower") {
        interactions <- interactions[-which(rownames(interactions) == sp.ext), , drop = F]
      }
      
      if (participant == "higher") {
        interactions <- interactions[, -which(colnames(interactions) == sp.ext), drop = F]
      }
    } else {
      retry <- retry + 1L
      if (retry >= 3L) {
        warning("Max no. of retries for dead interactions reached, aborting")
        break
      }
    }
    
    if (!isTRUE(abort)) {
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
      drop_rows <- drop_rows[-which(drop_rows %in% retain_low)]
    } # retain dead interactions
    
    if (!length(drop_rows) == 0L) {
    interactions <- interactions[-drop_rows, , drop = F]
    }

    if (!length(retain_high) == 0L) {
      drop_cols <- drop_cols[-which(drop_cols %in% retain_high)]
    } # retain dead interactions
    
    if (!length(drop_cols) == 0L) {
      interactions <- interactions[, -drop_cols, drop = F]
    }
    
    # break if n.lower/higher is 0
    if ((tail(dead, 1)[, "n.lower"] == 0 | tail(dead, 1)[, "n.higher"] == 0)) {
      dead <- dead[-nrow(dead),]
      break
      }
    
    # calculate new web w/ updated I_mat
    if ((nrow(interactions) >= 2L | ncol(interactions) >= 2L)) {
      m2 <- matrix(rmultinom(1, Nobs, interactions), nrow = nrow(interactions),
                  ncol = ncol(interactions))
      dimnames(m2) <- dimnames(interactions)
      
      # get dead interactions
      dead_low <- which(rowSums(m2) == 0L)
      dead_high <- which(colSums(m2) == 0L)
      
    } else {
      break
    }
    
    if (participant == "lower" & NROW(m2) < 2L) 
      break
    if (participant == "higher" & NCOL(m2) < 2L) 
      break
    if (participant == "both" & min(dim(m2)) < 2L) 
      break
    if (any(dim(ext.temp$web) == 1L)) 
      break
    if (method == "external") {
      ext.col[ext.col > ext.col[1L]] <- ext.col[ext.col > ext.col[1]] - 1L
      ext.row[ext.row > ext.row[1L]] <- ext.row[ext.row > ext.row[1]] - 1L
      ext.row <- ext.row[-1L]
      ext.col <- ext.col[-1L]
    }
    
    # "Skip" iteration if extc sp only had dead interactions 
    if (!isTRUE(skip)) {
      i <- i + 1L
    }
    
    # reset
    skip <- F
    abort <- F
    
  }

  # Add last line of dead
  dead2 <- rbind(dead, c(NROW(dead), nrow(interactions), ncol(interactions), 0L, 0L)) # add last line of dead (i.e. all sp extc)

  if (participant == "lower" & method == "degree") {
    if (length(table(dead[, 2])) > 1) 
      dead2[, 2] <- 1
  }
  
  # if (nrow(dead) + 1 != nrow(dead2)) 
  #   stop("PANIC! Something went wrong with the extinct sequence! Please contact the author to fix this!!")
  # 
  if (make.bipartite == T) {
    out <- list(dead2, ext.temp$web)
    class(out) <- "bipartite"
    attr(out, "exterminated") <- c("both", "lower", "higher")[pmatch(participant, c("both", "lower", "higher"))]
    attr(out, "exterminated")
  } else {
    out <- dead2
  }
  return(out)
}
