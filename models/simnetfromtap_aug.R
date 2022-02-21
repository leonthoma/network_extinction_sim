# Augmented simnetfromtap fy from tapnet; allows to set contribution (i.e.
# importance of community structure variables (abundance, latent traits, traits)
# individually);
# For each variable (abundance, latent traits, and traits, respectively) in 
# ctrb_vec specify either "high" or "low".

simnetfromtap_aug <- function(traits,
                              abuns,
                              paramsList,
                              pems,
                              tmatch_type_pem,
                              tmatch_type_obs,
                              ctrb_vec = c(NULL),
                              initial_sim = FALSE)
{ if (!is.null(ctrb_vec)){
  if (length(ctrb_vec) != 3) warning("ctrb_vec has to be a character vector of length 3, otherwise default imporatance of community variables is used")
  }
  if (!is.null(traits$low)) 
    traits$low <- traits$low[order(rownames(traits$low)), 
                             , drop = FALSE]
  if (!is.null(traits$high)) 
    traits$high <- traits$high[order(rownames(traits$high)), 
                               , drop = FALSE]
  abuns$low <- abuns$low[order(names(abuns$low))]
  abuns$high <- abuns$high[order(names(abuns$high))]
  pems$low <- pems$low[order(rownames(pems$low)), , drop = FALSE]
  pems$high <- pems$high[order(rownames(pems$high)), , drop = FALSE]
  if (is.null(traits$high) | is.null(traits$low)) {
    T_mat <- matrix(1, nrow = length(abuns$low), ncol = length(abuns$high))
    T_mat <- T_mat/sum(T_mat)
  }
  else {
    T_mat <- tmatch(t(outer(traits$high[, 1], traits$low[, 1], "-")),
                    type = tmatch_type_obs[1], width = paramsList[[5]][1])
    rownames(T_mat) <- rownames(traits$low)
    colnames(T_mat) <- rownames(traits$high)
    T_mat <- T_mat/sum(T_mat)
    if (ncol(traits$low) > 1) {
      if (length(tmatch_type_obs) == 1) 
        tmatch_type_obs <- rep(tmatch_type_obs, ncol(traits$low))
      if (length(paramsList[[5]]) == 1) 
        paramsList[[5]] <- rep(paramsList[[5]], ncol(traits$low))
      for (i in 2:ncol(traits$low)) {
        T_mat_next <- tmatch(t(outer(traits$high[, i], traits$low[, i], "-")),
                             type = tmatch_type_obs[i], 
                             width = paramsList[[5]][i])
        rownames(T_mat) <- rownames(traits$low)
        colnames(T_mat) <- rownames(traits$high)
        T_mat_next <- T_mat_next/sum(T_mat_next)
        T_mat <- T_mat * T_mat_next
      }
    }
  }
  nspec_low <- length(abuns$low)
  nspec_high <- length(abuns$high)
  a_mat_low <- matrix(rep(paramsList[[1]], nspec_low), nrow = nspec_low, 
                      byrow = TRUE)
  a_mat_high <- matrix(rep(paramsList[[2]], nspec_high), nrow = nspec_high, 
                       byrow = TRUE)
  lat_low <- as.vector(scale(rowSums(a_mat_low * pems[[1]])))
  lat_high <- as.vector(scale(rowSums(a_mat_high * pems[[2]]))) + 
    paramsList[[3]]
  L_mat <- tmatch(t(outer(lat_high, lat_low, "-")), type = tmatch_type_pem, 
                  width = paramsList[[4]])
  rownames(L_mat) <- rownames(pems[[1]])
  colnames(L_mat) <- rownames(pems[[2]])
  L_mat <- L_mat/sum(L_mat)
  A_mat <- as.matrix(abuns$low) %*% t(abuns$high)
  
  # Get contributions
  # if no contributions are specified use 1 as exponent for each community var
  if (is.null(ctrb_vec)) {
    A_ctrb <- 1
    L_ctrb <- 1
    T_ctrb <- 1
  }
  
  # if contributions are specified choose val according to setting
  if (!is.null(ctrb_vec)) {
    A_ctrb <- switch(ctrb_vec[[1]], "low" = 1.9, "high" = .1)
    L_ctrb <- switch(ctrb_vec[[2]], "low" = 1.9, "high" = .1)
    T_ctrb <- switch(ctrb_vec[[3]], "low" = 1.9, "high" = .1)
  }
  
  A_mat <- A_mat^A_ctrb # set contribution
  A_mat <- A_mat/sum(A_mat)
  if (is.null(paramsList[["delta"]])) {
    L_mat <- L_mat^L_ctrb # set contribution
    T_mat <- T_mat^T_ctrb # set contribution
    LT <- (L_mat * T_mat)/sum(L_mat * T_mat)
  }
  else {
    delta <- paramsList[["delta"]]
    L_mat <- L_mat^L_ctrb # set contribution
    T_mat <- T_mat^T_ctrb # set contribution
    LT <- (L_mat * T_mat)^(plogis(delta))/sum((L_mat * T_mat)^(plogis(delta)))
  }
  I_mat <- A_mat * LT
  I_mat <- I_mat/sum(I_mat)
  
  if (initial_sim == TRUE) {
    return(list("A_mat" = A_mat, "T_mat" = T_mat,
                "L_mat" = L_mat, "I_mat" = I_mat))
  } else {
    return(list("I_mat" = I_mat))
  }
}

tmatch <- function(delta_t, # Vector of pairwise trait differences (higher - lower)
                   type = "normal", # Trait matching function
                   width = 1, # Width parameter of trait matching function,
                   shift = 0, # shift parameter (optimum trait distance)
                   err = 1E-5 # "baseline" probability of match, even if traits do not match at all
){# Calculate interaction probabilities based on trait matching
  # shift
  delta_t <- delta_t + shift
  
  # lognormal distribution with mode shifted to zero:
  if(type == "shiftlnorm") out <- dlnorm(delta_t + width, meanlog = log(width) + 1) + err
  
  # normal distribution:
  if(type == "normal") out <- dnorm(delta_t, mean = 0, sd = width) + err
  
  return(out)
}
