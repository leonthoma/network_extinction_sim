# Small change in the internal function "one.second.extinct" of the bipartite package. The change allows rewiring of interations after each step of extinction
# The arguments are the same of "one.second.extinct" function and additional arguments are described below.
# New arguments:
# rewiring - Logical argument to specify if allow rewiring (default rewiring = FALSE).
# probabilities.rewiring1 - A matrix with probabilities of rewiring, must be the same dimensions of the web (i.e. network). See section Methods in Vizentin-Bugoni et al. [in review] for details). This matrix is required in step ii of framework (default probabilities.rewiring1 = NULL).
# probabilities.rewiring2 - A matrix with probabilities of rewiring, must be the same dimensions of web. See section Methods in Vizentin-Bugoni et al. [in review] for details. This matrix is required in step iii of framework (default probabilities.rewiring2 = NULL).
# method.rewiring = Type of method used to trial rewiring, partial match to "one.try.single.partner", "multiple.trials.single.partner", "multiple.trials.multiples.partners", "one.try.each.partner" and "multiple.trials.each.partner". See section Methods in Vizentin-Bugoni et al. [in review] for details; (default method.rewiring = "one.try.single.partner").
one.second.extinct.mod_aug <- function(web, participant = "higher", method = "abun", ext.row = NULL, ext.col = NULL, 
                                   rewiring = FALSE, probabilities.rewiring1 = NULL, probabilities.rewiring2 = NULL,
                                   method.rewiring = "one.try.single.partner") {
  dead <- matrix(nrow = 0, ncol = 3)
  colnames(dead) <- c("no", "ext.lower", "ext.higher")
  m2 <- web
  i <- 1
  METHOD.REWIRING = c("one.try.single.partner", "multiple.trials.single.partner", "multiple.trials.multiples.partners", "one.try.each.partner", "multiple.trials.each.partner")
  method.rewiring <- pmatch(method.rewiring, METHOD.REWIRING)
  if (length(method.rewiring) > 1) {
    stop("\n Only one argument is accepted in method.rewiring \n")
  }
  if (is.na(method.rewiring)) {
    stop("\n Invalid method.rewiring \n")
  }
  if(method.rewiring == 4 | method.rewiring == 5){
    keep.trying <- TRUE
  } else {
    keep.trying <- FALSE
  }
  method.rewiring <- ifelse(method.rewiring == 4, 1, ifelse(method.rewiring == 5, 2, method.rewiring))
  if(rewiring){
    if(any(web%%1!=0)){
      stop("\n If rewiring is TRUE the web must must contain only integers \n")
    }
    if(is.null(rownames(web)) | is.null(colnames(web))){
      stop("\n If rewiring is TRUE the web must must rownames and colnames\n")
    }
    if(is.null(probabilities.rewiring1) | is.null(probabilities.rewiring2)){
      stop("\n If rewiring is TRUE probabilities.rewiring1 and probabilities.rewiring1 must not be NULL\n")
    }
  }
  repeat {
    ext.temp <- extinction.mod(m2, participant = participant, method = method, ext.row = ext.row, ext.col = ext.col)
    # extinction.mod returns list w/ plants (rows; "rexcl") that lost partners and
    # animals (cols; "cexcl") that lost partners as well as network ("web")
    if(rewiring){
      if(!is.null(ext.temp$rexcl)){ # Plant is extinct, looking for new interaction partners of birds that interacted with lost plant
        sp.ext <- rownames(ext.temp$rexcl) # name of extc. plant
        sp.try.rewiring <- which(ext.temp$rexcl>0) # number & position of possible interaction partners (birds) for rewiring
        sp.surv <- seq_len(nrow(ext.temp$web)) # seq w/ all plant species
        sp.surv <- sp.surv[-1*which(rownames(ext.temp$web) %in% sp.ext)] # survived plant species
        for(jj in sp.try.rewiring){
          sp.surv.temp <- sp.surv
          go <- TRUE
          m <- 0
          if(method.rewiring == 1 | method.rewiring == 3){
            trials <- 1
          } else {
            trials <- ext.temp$rexcl[1, jj]
          }
          while (go) {
            m <- m+1
            # choose rewiring partner based on highest abundance, irrespective
            # if a interaction was formerly observed
            sp.add <- which.max(probabilities.rewiring1[sp.surv.temp, jj]) # get name of sp
            sp.surv.prob2 <- max(probabilities.rewiring1[sp.surv.temp, jj]) # get rew prob
            
            # original code to determine rewiring partner
            # sp.surv.prob1 <- probabilities.rewiring1[sp.surv.temp, jj] # probs of rewiring to a potential partner
            # sp.add <- sample(as.character(sp.surv.temp), 1, prob = sp.surv.prob1) # randomly choosing new partner for birds that interacted w/ lost plant
            # sp.surv.prob2 <- probabilities.rewiring2[as.numeric(sp.add), jj] # prob of rewiring based on rewiring factors
            n.add <- rbinom(1, trials, sp.surv.prob2) # binomial trail to determine rewiring success
            if(n.add>0){
              ext.temp$web[as.numeric(sp.add), jj] <- ext.temp$web[as.numeric(sp.add), jj]+n.add # update interaction freq in interaction matrix
            }
            if(method.rewiring == 1 | method.rewiring == 2){
              if(!keep.trying){
                go <- FALSE  
              } else{
                if((method.rewiring == 1 & n.add>0) | (method.rewiring == 2 & n.add == trials)){ # if rewiring failed w/ single try or multipe tries and trails left
                  go <- FALSE
                } else{
                  sp.surv.temp <- sp.surv.temp[-1*which(sp.surv.temp%in% sp.add)] # remove chosen partner with failed new interaction from possible partners
                  trials <- trials-n.add # update no of trails
                  if(length(sp.surv.temp)<1){
                    go <- FALSE
                  }
                }
              }
            } else {
              if(ext.temp$rexcl[1, jj] == m){
                go <- FALSE
              }
            }
          }
        }
      }
      if(!is.null(ext.temp$cexcl)){ # Bird is extinct, looking for new interaction partners of plants that interacted with lost bird
        sp.ext <- colnames(ext.temp$cexcl)  # name of extc. bird
        sp.try.rewiring <- which(ext.temp$cexcl>0) # number & position of interaction partners (plants) for possible rewiring
        sp.surv <- seq_len(ncol(ext.temp$web))
        sp.surv <- sp.surv[-1*which(colnames(ext.temp$web) %in% sp.ext)] # survived bird species
        for(ii in sp.try.rewiring){
          sp.surv.temp <- sp.surv 
          go <- TRUE
          m <- 0
          if(method.rewiring == 1 | method.rewiring == 3){ # define no of trails
            trials <- 1
          } else {
            trials <- ext.temp$cexcl[ii, 1]
          }
          while (go) {
            m <- m+1
            # choose rewiring partner based on highest abundance, irrespective
            # if a interaction was formerly observed
            sp.add <- which.max(probabilities.rewiring1[sp.surv.temp, ii]) # get name of sp
            sp.surv.prob2 <- max(probabilities.rewiring1[sp.surv.temp, ii]) # get rew prob
            
            # original code to determine rewiring partner
            # sp.surv.prob1 <- probabilities.rewiring1[ii, sp.surv.temp] # probs of rewiring to a potential partner (either random or anbundance)
            # sp.add <- sample(as.character(sp.surv.temp), 1, prob = sp.surv.prob1) # randomly choosing new partner for plants that interacted w/ lost bird
            # sp.surv.prob2 <- probabilities.rewiring2[ii, as.numeric(sp.add)] # prob of rewiring based on chosen rewiring factor (e.g. morphology)
            n.add <- rbinom(1, trials, sp.surv.prob2) # binomial trail to determine rewiring success
            if(n.add>0){
              ext.temp$web[ii, as.numeric(sp.add)] <- ext.temp$web[ii, as.numeric(sp.add)]+n.add # update interaction freq in interaction matrix
            }
            if(method.rewiring == 1 | method.rewiring == 2){ 
              if(!keep.trying){
                go <- FALSE  
              } else{
                if((method.rewiring == 1 & n.add>0) | (method.rewiring == 2 & n.add == trials)){ # if rewiring failed w/ single try or multipe tries and trails left
                  go <- FALSE 
                } else{
                  sp.surv.temp <- sp.surv.temp[-1*which(sp.surv.temp%in% sp.add)] # remove chosen partner with failed new interaction from possible partners
                  trials <- trials-n.add # update no of trails
                  if(length(sp.surv.temp)<1){
                    go <- FALSE
                  }
                }
              }
            } else { # only run when using multiple.trials.multiples.partners
              if(ext.temp$cexcl[ii, 1] == m){
                go <- FALSE
              }
            }
          }
        }
      }
    }
    n <- ext.temp$web
    dead <- rbind(dead, c(i, attributes(m2 <- empty(n, count = TRUE))$empty))
    if (participant == "lower" & NROW(m2) < 2) 
      break
    if (participant == "higher" & NCOL(m2) < 2) 
      break
    if (participant == "both" & min(dim(m2)) < 2) 
      break
    if (any(dim(n) == 1)) 
      break
    if (method == "external") {
      ext.col[ext.col > ext.col[1]] <- ext.col[ext.col > ext.col[1]] - 1
      ext.row[ext.row > ext.row[1]] <- ext.row[ext.row > ext.row[1]] - 1
      ext.row <- ext.row[-1]
      ext.col <- ext.col[-1]
    }
    i <- i + 1
  }
  dead2 <- rbind(dead, c(NROW(dead) + 1, NROW(m2), NCOL(m2)))
  if (participant == "lower" & method == "degree") {
    if (length(table(dead[, 2])) > 1) 
      dead2[, 2] <- 1
  }
  if (nrow(dead) + 1 != nrow(dead2)) 
    stop("PANIC! Something went wrong with the extinct sequence! Please contact the author to fix this!!")
  if (participant == "lower") 
    supposed.length <- NROW(dead2)
  if (participant == "higher") 
    supposed.length <- NCOL(dead2)
  if (participant == "both") 
    supposed.length <- NROW(dead2)
  if (NROW(dead2) != supposed.length) {
    missing <- supposed.length - NROW(dead2)
    addit1 <- (NROW(dead2) + 1):(NROW(dead2) + missing)
    addit2n3 <- rep(0, times = missing)
    dead2 <- rbind(dead2, as.matrix(data.frame(addit1, addit2n3, addit2n3)))
  }
  out <- list(dead2, ext.temp$web)
  class(out) <- "bipartite"
  attr(out, "exterminated") <- c("both", "lower", "higher")[pmatch(participant, c("both", "lower", "higher"))]
  attr(out, "exterminated")
  return(out)
}
