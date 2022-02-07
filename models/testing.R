# first test using tapnet to simulate network
setwd("~/Documents/Uni/M.sc/Master Thesis/Networks/models/")

library(tapnet)
source("tapnet_helper.R")
source("simnetfromtap_aug.R")
source("./tapnet/tapnet/R/simulate_tapnet.R")
library(phytools)
library(bipartite)
library(vegan)
library(ape)
library(MPSEM)

# read species names and make df 
animals <- read.csv("animal-names.csv", colClasses = c("NULL", "character"),
                    stringsAsFactors = F)
plants <- read.csv("flower-and-plant-names.csv", stringsAsFactors = F)
sp_names <- data.frame("plants" = plants,
                       "animals" = animals[1:nrow(plants),],
                       stringsAsFactors = F)

# simulation -----
simnet <- simulate_tapnet(nlower = 10, nhigher = 10, ntraits_nopem = 2,
                          ntraits_pem = 0, abuns = "lognormal", Nobs = 420)

# visualize phylo
par(mfrow = c(2,1))
phytools::plotTree(simnet$trees$low)
phytools::plotTree(simnet$trees$high)
par(mfrow = c(1,1))

# computing phylogenetic distances
cophenetic(simnet$trees$low)
cophenetic(simnet$trees$high)

# overview of network matrices
head(simnet$networks[[1]]$abuns)
head(simnet$networks[[1]]$traits)
# head(simnet$networks[[1]]$pems) # phylogenetic eigenmatrices
head(simnet$networks[[1]]$I_mat)

# visualize network
plotweb(simnet$networks[[1]]$I_mat)
visweb(simnet$networks[[1]]$I_mat)

# initial simulation ----

# fy for deleting rows/cols without interactions; if multiple networks shall be
# used pass list of lists (NOT vector i.e. c())
clean <- function(x) {
  if (is.null(names(x))) {
    n_nets <- length(x) # no of networks when multiple networks are passed
  } else {
    n_nets <- 1 # no of networks when single network is passed
  }
  
  out <- vector(mode = "list", length = n_nets) # list for webs of all networks
  names(out) <- paste("Net", 1:n_nets)
  
  # iterate over all networks
  for (i in seq(n_nets)) {
    network <- list # list for webs of each network
    
    # calculate no of webs in x
    if (is.null(names(x)))
      n_webs <- length(x[[i]]) # for list of networks
    else
      n_webs <- length(x) # for single network
    
    
    # iterate over all web 
    for (j in seq(n_webs)) {
      if (is.null(names(x))) {
        web <- as.data.frame(x[[i]][j]) # get web for list of networks
        colnames(web) <- colnames(init_sim$networks[[1]]$web) # manually set names
      } else {
      web <- x[[j]] # get web for single network
      }
      
      # sp w/o any interactions
      no_int_low <- which(rowSums(web) == 0) # lower trophic level
      no_int_high <- which(colSums(web) == 0) # higher trophic level
  
      # exclude if rows/cols w/ no interactions occur
      if (!length(no_int_low) == 0)
        network[[j]] <- web[-no_int_low, ]
      if (!length(no_int_high) == 0)
        network[[j]] <- web[, -no_int_high]
    }
    # if every sp has at least one interaction return intial network
    if (length(no_int_low) == 0 & length(no_int_high) == 0)  
      out[[i]] <- x[[i]]
    else out[[i]] <- network
  }
  return(out)
}

# create initial network for community variables
init_sim <- simulate_tapnet_aug(nlower = 26, nhigher = 14, ntraits_nopem = 4,
                       ntraits_pem = 2, abuns = "lognormal", Nobs = 1111,
                       names = sp_names)

source("simnetfromtap_ctrb.R")

# Atl
sim1 <- simnetfromtap_ctrb(traits = sim$traits_all,
                          abuns = sim$networks[[1]]$abuns,
                          paramsList = sim$sim_params,
                          pems = sim$networks[[1]]$pems,
                          tmatch_type_pem = "normal",
                          tmatch_type_obs = "normal",
                          ctrb_vec = c("high", "low", "low"),
                          Nwebs = 1,
                          Nobs = 1111)
clean(sim1)

# sim1 <- simnetfromtap_aug(traits = sim$traits_all,
#                        abuns = sim$networks[[1]]$abuns,
#                        paramsList = sim$sim_params,
#                        pems = sim$networks[[1]]$pems,
#                        tmatch_type_pem = "normal",
#                        tmatch_type_obs = "normal",
#                        ctrb_vec = c("high", "low", "low"))
# sim1web <- matrix(rmultinom(1, 1111, sim1), nrow = nrow(sim1),
#                   ncol = ncol(sim1))
# dimnames(sim1web) <- dimnames(sim$networks[[1]]$I_mat)

# atL
sim2 <- simnetfromtap_aug(traits = sim$traits_all,
                       abuns = sim$networks[[1]]$abuns,
                       paramsList = sim$sim_params,
                       pems = sim$networks[[1]]$pems,
                       tmatch_type_pem = "normal",
                       tmatch_type_obs = "normal",
                       ctrb_vec = c("low", "high", "low"))
sim2web <- matrix(rmultinom(1, 1111, sim2), nrow = nrow(sim2),
                  ncol = ncol(sim2))
dimnames(sim2web) <- dimnames(sim$networks[[1]]$I_mat)

# aTl
sim3 <- simnetfromtap_aug(traits = sim$traits_all,
                       abuns = sim$networks[[1]]$abuns,
                       paramsList = sim$sim_params,
                       pems = sim$networks[[1]]$pems,
                       tmatch_type_pem = "normal",
                       tmatch_type_obs = "normal",
                       ctrb_vec = c("low", "low", "high"))
sim3web <- matrix(rmultinom(1, 1111, sim3), nrow = nrow(sim3),
                  ncol = ncol(sim3))
dimnames(sim3web) <- dimnames(sim$networks[[1]]$I_mat)

# ATL
sim4 <- simnetfromtap_aug(traits = sim$traits_all,
                       abuns = sim$networks[[1]]$abuns,
                       paramsList = sim$sim_params,
                       pems = sim$networks[[1]]$pems,
                       tmatch_type_pem = "normal",
                       tmatch_type_obs = "normal",
                       ctrb_vec = c("high", "high", "high"))
sim4web <- matrix(rmultinom(1, 1111, sim4), nrow = nrow(sim4),
                  ncol = ncol(sim4))
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

# extinction model

source("rewiring_vizentin-bugoni_2019/Functions/extinction.mod.R")
source("rewiring_vizentin-bugoni_2019/Functions/one.second.extinct.mod.R")
source("rewiring_vizentin-bugoni_2019/Functions/calc.mean.one.second.extinct.mod.R")
source("rewiring_vizentin-bugoni_2019/Functions/matrix.p1.R")
source("rewiring_vizentin-bugoni_2019/Functions/IC.R") # calc of 95 percent confidence interval

# rewiring probabilites
# pairwise relative abundances
rew_abund <- sim$networks[[1]]$abuns$low[-no_int_low] %*%
  t(sim$networks[[1]]$abuns$high[-no_int_high])
row.names(rew_abund) <- sim$trees[[1]]$tip.label[-no_int_low]

# relative abundances of lower
rew_abund_low <- sweep(rew_abund, 2, colSums(rew_abund), "/")

# relative abundances of higher
rew_abund_high <- sweep(rew_abund, 1, rowSums(rew_abund), "/")

# run extinction models
# original
extc_sim <- one.second.extinct.mod_aug(web = network,
                                       participant = "lower",
                                       method = "random", 
                                       rewiring = T, 
                                       probabilities.rewiring1 = rew_abund_low, 
                                       probabilities.rewiring2 = rew_abund,
                                       method.rewiring = "one.try.single.partner")

# Atl
extc_sim_Atl <- one.second.extinct.mod_aug(web = sim1web,
                                       participant = "lower",
                                       method = "random", 
                                       rewiring = T, 
                                       probabilities.rewiring1 = rew_abund_low, 
                                       probabilities.rewiring2 = rew_abund,
                                       method.rewiring = "one.try.single.partner")

# atL
extc_sim_atL <- one.second.extinct.mod_aug(web = sim2web,
                                           participant = "lower",
                                           method = "random", 
                                           rewiring = T, 
                                           probabilities.rewiring1 = rew_abund_low, 
                                           probabilities.rewiring2 = rew_abund,
                                           method.rewiring = "one.try.single.partner")

# aTl
extc_sim_aTl <- one.second.extinct.mod_aug(web = sim3web,
                                           participant = "lower",
                                           method = "random", 
                                           rewiring = T, 
                                           probabilities.rewiring1 = rew_abund_low, 
                                           probabilities.rewiring2 = rew_abund,
                                           method.rewiring = "one.try.single.partner")

# ATL
extc_sim_ATL <- one.second.extinct.mod_aug(web = sim4web,
                                           participant = "lower",
                                           method = "random", 
                                           rewiring = T, 
                                           probabilities.rewiring1 = rew_abund_low, 
                                           probabilities.rewiring2 = rew_abund,
                                           method.rewiring = "one.try.single.partner")

# calculate percentages of remaining sp
extc_lo_per <- extc_sim[[1]][, 4]/nrow(network)*100
extc_hi_per <- extc_sim[[1]][, 5]/ncol(network)*100

extc_lo_per_Atl <- extc_sim_Atl[[1]][, 4]/nrow(sim1web)*100
extc_hi_per_Atl <- extc_sim_Atl[[1]][, 5]/ncol(sim1web)*100

extc_lo_per_atL <- extc_sim_atL[[1]][, 4]/nrow(sim2web)*100
extc_hi_per_atL <- extc_sim_atL[[1]][, 5]/ncol(sim2web)*100

extc_lo_per_aTl <- extc_sim_aTl[[1]][, 4]/nrow(sim3web)*100
extc_hi_per_aTl <- extc_sim_aTl[[1]][, 5]/ncol(sim3web)*100

extc_lo_per_ATL <- extc_sim_ATL[[1]][, 4]/nrow(sim4web)*100
extc_hi_per_ATL <- extc_sim_ATL[[1]][, 5]/ncol(sim4web)*100

# library(ggplot2)
#  
# ggplot() +
#   geom_line(aes(1-extc_lo_per, extc_hi_per)) +
#   geom_line(aes(1-extc_lo_per_Atl, extc_hi_per_Atl), col = "firebrick") +
#   geom_line(aes(1-extc_lo_per_atL, extc_hi_per_atL), col = "dodgerblue") +
#   geom_line(aes(1-extc_lo_per_aTl, extc_hi_per_aTl), col = "darkseagreen") +
#   scale_color_manual(labels = c("Original", "Atl", "atL", "aTl")) +
#   # geom_line(aes(1-extc_lo_per_ATL, extc_hi_per_ATL), col = "sienna") +
#   labs(x = "Plants removed [%]",
#        y = "Animals persiting [%]")

# plot
plot(rev(extc_lo_per), extc_hi_per, type = "l", ylab = "animals persisting",
     xlab = "plants removed", main = "Extinction cascade w/ abundance rewiring")
lines(rev(extc_lo_per_Atl), extc_hi_per_Atl, col = "firebrick", lty = 2) # Atl
lines(rev(extc_lo_per_atL), extc_hi_per_atL, col = "dodgerblue", lty = 3) # atL
lines(rev(extc_lo_per_aTl), extc_hi_per_aTl, col = "darkseagreen", lty = 4) # aTl
#lines(1-extc_lo_per_ATL, extc_hi_per_ATL, col = "sienna") # ATL
legend("bottomleft",
       title = "Community var importances",
       legend = c("Original","Atl", "atL", "aTl"),
       col = c("black", "firebrick", "dodgerblue", "darkseagreen"),
       lty = 1:4)

