# first test using tapnet to simulate network

library(tapnet)
library(phytools)
library(bipartite)
library(vegan)

# simulation
simnet <- simulate_tapnet(nlower = 10, nhigher = 10, ntraits_nopem = 2,
                          ntraits_pem = 0, abuns = "lognormal", Nobs = 420)

# visualize phylo
par(mfrow = c(2,1))
phytools::plotTree(simnet$trees$low)
phytools::plotTree(simnet$trees$high)
par(mfrow = c(1,1))

# overview of network matrices
head(simnet$networks[[1]]$abuns)
head(simnet$networks[[1]]$traits)
# head(simnet$networks[[1]]$pems) # phylogenetic eigenmatrices
head(simnet$networks[[1]]$I_mat)

# visualize network
plotweb(simnet$networks[[1]]$I_mat)
visweb(simnet$networks[[1]]$I_mat)

# altering simnetfromtap fy
simnetfromtap1 <- function (traits,
                            abuns,
                            paramsList,
                            pems,
                            tmatch_type_pem,
                            tmatch_type_obs,
                            ctrb_list) 
{
  stopifnot(length(ctrb_list) == 3)
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
  # Set contributions
  A_ctrb <- switch(ctrb_list[[1]], "low" = 1.9, "high" = .1)
  L_ctrb <- switch(ctrb_list[[2]], "low" = 1.9, "high" = .1)
  T_ctrb <- switch(ctrb_list[[3]], "low" = 1.9, "high" = .1)
  
  A_mat <- A_mat^A_ctrb 
  A_mat <- A_mat/sum(A_mat)
  if (is.null(paramsList[["delta"]])) {
    L_mat <- L_mat^1.9
    T_mat <- T_mat^.1
    LT <- (L_mat * T_mat)/sum(L_mat * T_mat)
  }
  else {
    delta <- paramsList[["delta"]]
    L_mat <- L_mat^1.9
    T_mat <- T_mat^.1
    LT <- (L_mat * T_mat)^(plogis(delta))/sum((L_mat * T_mat)^(plogis(delta)))
  }
  I_mat <- A_mat * LT
  I_mat <- I_mat/sum(I_mat)
  #return(I_mat)
  return(list(I_mat, L_mat, T_mat))
}

# initial simulation
sim <- simulate_tapnet(nlower = 10, nhigher = 10, ntraits_nopem = 2,
                       ntraits_pem = 0, abuns = "lognormal", Nobs = 420)
dimnames(sim$networks[[1]]$web) <- dimnames(sim$networks[[1]]$I_mat)
# use parms form initial sim

# abun
sim1 <- simnetfromtap1(traits = sim$traits_all,
                       abuns = sim$networks[[1]]$abuns,
                       paramsList = sim$sim_params,
                       pems = sim$networks[[1]]$pems,
                       tmatch_type_pem = "normal",
                       tmatch_type_obs = "normal",
                       ctrb_list = c("high", "low", "low"))
sim1web <- matrix(rmultinom(1, 420, sim1[[1]]), nrow = nrow(sim1[[1]]),
                  ncol = ncol(sim1[[1]]))
dimnames(sim1web) <- dimnames(sim$networks[[1]]$I_mat)

# phylo
sim2 <- simnetfromtap1(traits = sim$traits_all,
                       abuns = sim$networks[[1]]$abuns,
                       paramsList = sim$sim_params,
                       pems = sim$networks[[1]]$pems,
                       tmatch_type_pem = "normal",
                       tmatch_type_obs = "normal",
                       ctrb_list = c("low", "high", "low"))
sim2web <- matrix(rmultinom(1, 420, sim2[[1]]), nrow = nrow(sim2[[1]]),
                  ncol = ncol(sim2[[1]]))
dimnames(sim2web) <- dimnames(sim$networks[[1]]$I_mat)

# traits
sim3 <- simnetfromtap1(traits = sim$traits_all,
                       abuns = sim$networks[[1]]$abuns,
                       paramsList = sim$sim_params,
                       pems = sim$networks[[1]]$pems,
                       tmatch_type_pem = "normal",
                       tmatch_type_obs = "normal",
                       ctrb_list = c("low", "low", "high"))
sim3web <- matrix(rmultinom(1, 420, sim3[[1]]), nrow = nrow(sim3[[1]]),
                  ncol = ncol(sim3[[1]]))
dimnames(sim3web) <- dimnames(sim$networks[[1]]$I_mat)

# all
sim4 <- simnetfromtap1(traits = sim$traits_all,
                       abuns = sim$networks[[1]]$abuns,
                       paramsList = sim$sim_params,
                       pems = sim$networks[[1]]$pems,
                       tmatch_type_pem = "normal",
                       tmatch_type_obs = "normal",
                       ctrb_list = c("high", "high", "high"))
sim4web <- matrix(rmultinom(1, 420, sim4[[1]]), nrow = nrow(sim4[[1]]),
                  ncol = ncol(sim4[[1]]))
dimnames(sim4web) <- dimnames(sim$networks[[1]]$I_mat)

# visualize
plotweb(sim$networks[[1]]$web)
plotweb(sim1web)
plotweb(sim2web)
plotweb(sim3web)
plotweb(sim4web)

# basic metrics
# connectance
connectance <- function(web) {
  sum(web != 0) / (nrow(web) * ncol(web))
}

connectance(sim$networks[[1]]$web)
connectance(sim1web)
connectance(sim2web)
connectance(sim3web)
connectance(sim4web)

# nestedness; nested rank
nestedness(sim$networks[[1]]$web)
nestedness(sim1web)
nestedness(sim2web)
nestedness(sim3web)

