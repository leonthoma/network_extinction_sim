# Small change in the internal function "one.second.extinct" of the bipartite package. The change allows rewiring of interations after each step of extinction
# The arguments are the same of "one.second.extinct" function and additional arguments are described below.
# New arguments:
# rewiring - Logical argument to specify if allow rewiring (default rewiring = FALSE).
# probabilities.rewiring - A matrix with probabilities of rewiring, must be the same dimensions of the web (i.e. network). See section Methods in Vizentin-Bugoni et al. [in review] for details). This matrix is required in step ii of framework (default probabilities.rewiring = NULL).
# probabilities.rewiring2 - A matrix with probabilities of rewiring, must be the same dimensions of web. See section Methods in Vizentin-Bugoni et al. [in review] for details. This matrix is required in step iii of framework (default probabilities.rewiring2 = NULL).
# method.rewiring = Type of method used to trial rewiring, partial match to "one.try.single.partner", "multiple.trials.single.partner", "multiple.trials.multiples.partners", "one.try.each.partner" and "multiple.trials.each.partner". See section Methods in Vizentin-Bugoni et al. [in review] for details; (default method.rewiring = "one.try.single.partner").
# 
# participant <- "lower"
# partner.choice <- choice_abund
# method.rewiring <- "abund"
# interactions <- sims$aTl$I_mat
# web <- sims$aTl$web

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
  dead <- matrix(nrow = 0, ncol = 5)
  
  # filter sp w/o any interactions
  sp_low <- which(rowSums(web) != 0)
  sp_high <- which(colSums(web) != 0)
  
  # drop sp w/o any interactions from web
  m2 <- web[sp_low, sp_high]
  
  dead <- cbind(0, 0, 0, nrow(m2), ncol(m2)) # correct for sp w/o any interactions
  
  colnames(dead) <- c("no", "ext.lower", "ext.higher", "n.lower", "n.higher")
  
  i <- 1
  
  # get interaction vals & drop sp w/o any interactions
  interactions <- interactions[which(rownames(interactions) %in% names(sp_low)),
                               which(colnames(interactions) %in% names(sp_high))]
  
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
    n_lo <- which(rowSums(m2) != 0) # initial number of sp in low trophic level
    n_hi <- which(colSums(m2) != 0) # initial number of sp in high trophic level
    
    if(rewiring){
      if(!is.null(ext.temp$rexcl)){ # Plant is extinct, looking for new interaction partners of birds that interacted with lost plant
        sp.ext <- rownames(ext.temp$rexcl) # name of extc. plant
        sp.try.rewiring <- which(ext.temp$rexcl>0) # number & position of interaction partners (higher trophic level)
        Nobs <- sum(m2)# - sum(m2[sp.ext, ]) # no of obs interactions
        
        # Choice of rewiring partner
        if (method.rewiring == "abund") {
          if (i == 1) {
            choice_abund <- partner.choice$low # get values
            }

          choice_abund <- choice_abund[which(names(choice_abund) %in% names(n_lo))] # drop sp w/o any interactions
          
          sp.ext.idx <- which(names(choice_abund) %in% sp.ext) # get idx of sp.ext (NOT THE SAME IN ext.temp$web !!)
          
          # get sp w/ highest abund
          sp.high.abund <- names(which.max(choice_abund[-sp.ext.idx]))
          
          rew.partner <- which(rownames(m2) %in% sp.high.abund)
          choice_abund <- choice_abund[-sp.ext.idx] # delete extc sp
        }
        
        # choose rewiring partner based on most similar trait
        if (method.rewiring == "trait") {
          if (i == 1) {
            choice_trait <- partner.choice$low # get values
          }
          
          choice_trait <- choice_trait[which(rownames(choice_trait) %in% names(n_lo)), ] # drop sp w/o any interactions
  
          sp.ext.idx <- which(rownames(choice_trait) %in% sp.ext) # get idx of sp.ext (NOT THE SAME IN ext.temp$web !!)
          
          # calculate euclidean distances
          trait.dist <- as.matrix(dist(choice_trait, diag = T, upper = T)) 
          
          # find sp w/ smallest trait dist
          sp.closest <- names(which.min(trait.dist[-sp.ext.idx, sp.ext.idx]))
            
          # match interacting sp of closest relative(check.int) to remaining sp (sp.rem)
          rew.partner <- which(rownames(m2) %in% sp.closest)
          choice_trait <- choice_trait[-sp.ext.idx, ] # delete extc sp
        }
        
        # choose rewiring partner based on closest phylogenetic distance
        if (method.rewiring == "phylo") {
          if (i == 1) {
            choice_phylo <- partner.choice$low # get values
          }
          
          choice_phylo <- choice_phylo[which(rownames(choice_phylo) %in% names(n_lo)), which(colnames(choice_phylo) %in% names(n_lo))] # drop sp w/o any interactions
          
          sp.ext.idx <- which(rownames(choice_phylo) %in% sp.ext) # get idx of sp.ext (NOT THE SAME IN ext.temp$web !!)
          if (length(sp.ext.idx) == 0) {
            warning("something went wrong")
            print(sp.ext)
            print(rownames(choice_phylo))
            print(i)
            print(names(n_lo))}
          
          # get closest relative
          close.rel <-  min(choice_phylo[-sp.ext.idx, sp.ext.idx]) # get distance
          sp.close.rel.tmp <- names(which(choice_phylo[, sp.ext.idx] == close.rel)) # closest sp
          
          # if multiple species have same distance, randomly choose one
          if (length(sp.close.rel.tmp) != 1){
            sp.close.rel <- sample(sp.close.rel.tmp, 1)
          } else {
            sp.close.rel <- sp.close.rel.tmp
          }
          
          rew.partner <- which(rownames(m2) %in% sp.close.rel)
          choice_phylo <- choice_phylo[-sp.ext.idx, -sp.ext.idx] # delete extc sp
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
            choice_abund <- partner.choice$high # get values
          }
          
          choice_abund <- choice_abund[which(names(choice_abund) %in% names(n_hi))] # drop sp w/o any interactions
          
          sp.ext.idx <- which(names(choice_abund) %in% sp.ext) # get idx of sp.ext (NOT THE SAME IN ext.temp$web !!)
          
          # get sp w/ highest abund
          sp.high.abund <- names(which.max(choice_abund[-sp.ext.idx]))
          
          rew.partner <- which(colnames(m2) %in% sp.high.abund)
          choice_abund <- choice_abund[-sp.ext.idx] # delete extc sp
        }
        
        # choose rewiring partner based on most similar trait
        if (method.rewiring == "trait") {
          if (i == 1) {
            choice_trait <- partner.choice$high # get values
          }
          
          choice_trait <- choice_trait[which(rownames(choice_trait) %in% names(n_hi)), ] # drop sp w/o any interactions
          
          sp.ext.idx <- which(rownames(choice_trait) %in% sp.ext) # get idx of sp.ext (NOT THE SAME IN ext.temp$web !!)
          
          # calculate euclidean distances
          trait.dist <- as.matrix(dist(choice_trait, diag = T, upper = T)) 
          
          # find sp w/ smallest trait dist
          sp.closest <- names(which.min(trait.dist[-sp.ext.idx, sp.ext.idx]))
          
          # match interacting sp of closest relative(check.int) to remaining sp (sp.rem)
          rew.partner <- which(colnames(m2) %in% sp.closest)
          choice_trait <- choice_trait[-sp.ext.idx, ] # delete extc sp
        }
        
        # choose rewiring partner based on closest phylogenetic distance
        if (method.rewiring == "phylo") {
          if (i == 1) {
            choice_phylo <- partner.choice$high # get values
          }
          
          choice_phylo <- choice_phylo[which(rownames(choice_phylo) %in% names(n_hi)), which(colnames(choice_phylo) %in% names(n_hi))] # drop sp w/o any interactions
          
          sp.ext.idx <- which(rownames(choice_phylo) %in% sp.ext) # get idx of sp.ext (NOT THE SAME IN ext.temp$web !!)
          
          # get closest relative
          close.rel <-  min(choice_phylo[-sp.ext.idx, sp.ext.idx]) # get distance
          sp.close.rel.tmp <- names(which(choice_phylo[, sp.ext.idx] == close.rel)) # closest sp
          
          # if multiple species have same distance, randomly choose one
          if (length(sp.close.rel.tmp) != 1){
            sp.close.rel <- sample(sp.close.rel.tmp, 1)
          } else {
            sp.close.rel <- sp.close.rel.tmp
          }
          
          rew.partner <- which(colnames(m2) %in% sp.close.rel)
          choice_phylo <- choice_phylo[-sp.ext.idx, -sp.ext.idx] # delete extc sp
        }
      
      # Shift interaction probability from old interaction (sp.ext + sp.try.rewiring) to new interaction (rew.partner + sp.try.rewiring)
      sp.ext.idx <- which(colnames(interactions) %in% sp.ext) # get idx for i_mat
      
      # Update I_mat
      prob.int.old <- interactions[sp.try.rewiring, sp.ext.idx] # get old interaction vals
      interactions[sp.try.rewiring, rew.partner] <- interactions[sp.try.rewiring, rew.partner] + (adapt_list$low[sp.try.rewiring] * prob.int.old) # set new interactions vals
    }
  }
    # removed rows and cols
    rem_r_c <- empty(ext.temp$web, count = T)
    
    # track no of extc and no of remaining sp
    dead <- rbind(dead, c(i, attributes(rem_r_c)$empty,
                          nrow(m2) - as.vector(attributes(rem_r_c)$empty[1]),
                          ncol(m2) - as.vector(attributes(rem_r_c)$empty[2])))
    
    if (participant == "lower" & NROW(m2) < 2) 
      break
    if (participant == "higher" & NCOL(m2) < 2) 
      break
    if (participant == "both" & min(dim(m2)) < 2) 
      break
    if (any(dim(ext.temp$web) == 1)) 
      break
    if (method == "external") {
      ext.col[ext.col > ext.col[1]] <- ext.col[ext.col > ext.col[1]] - 1
      ext.row[ext.row > ext.row[1]] <- ext.row[ext.row > ext.row[1]] - 1
      ext.row <- ext.row[-1]
      ext.col <- ext.col[-1]
    }
    
    if (vis.steps == T) {
      plotweb(m2)
    }
    
    i <- i + 1
    
    # update I_mat
    drop_rows <- which(!(rownames(m2) %in% rownames(rem_r_c)))
    drop_cols <- which(!(colnames(m2) %in% colnames(rem_r_c)))
    
    if (!length(drop_rows) == 0) {
      if (class(interactions) == "numeric") {
        interactions <- interactions[-drop_rows] # edge case for last col
      } else {
    interactions <- interactions[-drop_rows, ]
      }
    }
    
    if (!length(drop_cols) == 0) {
      if (class(interactions) == "numeric") {
        interactions <- interactions[-drop_cols] # edge case for last row
      } else {
      interactions <- interactions[, -drop_cols]
      }
    }
    
    # calculate new web w/ updated I_mat
    if (class(interactions) != "numeric") {
      m2 <- matrix(rmultinom(1, Nobs, interactions), nrow = nrow(interactions),
                  ncol = ncol(interactions))
      dimnames(m2) <- dimnames(interactions)
    } else {
      break
    }
# 
#     # filter sp w/o any interactions
#     sp_low_m2 <- which(rowSums(m2) != 0)
#     sp_hi_m2 <- which(colSums(m2) != 0)
# 
#     # drop sp w/o any interactions
#     m2 <- m2[sp_low_m2, sp_hi_m2]
    
  }
  dead2 <- rbind(dead, c(NROW(dead), NROW(rem_r_c), NCOL(rem_r_c), 0, 0)) # add last line of dead (i.e. all sp extc)
  if (participant == "lower" & method == "degree") {
    if (length(table(dead[, 2])) > 1) 
      dead2[, 2] <- 1
  }
  
  if (nrow(dead) + 1 != nrow(dead2)) 
    stop("PANIC! Something went wrong with the extinct sequence! Please contact the author to fix this!!")
  
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
