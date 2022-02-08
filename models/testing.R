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
library(purrr)
library(MPSEM)

# read species names and make df 
animals <- read.csv("animal-names.csv", colClasses = c("NULL", "character"),
                    stringsAsFactors = F)
plants <- read.csv("flower-and-plant-names.csv", stringsAsFactors = F)
sp_names <- data.frame("plants" = plants,
                       "animals" = animals[1:nrow(plants),],
                       stringsAsFactors = F)

# test simulation -----
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
# create initial network for community variables
init_sim <- simulate_tapnet_aug(nlower = 26, nhigher = 14, ntraits_nopem = 4,
                       ntraits_pem = 2, abuns = "lognormal", Nobs = 1111,
                       names = sp_names)


# Create simulated networks with setting contributions ----
source("simnetfromtap_ctrb.R")
source("clean.R")

sims <- list()

# list w/ ctrb_vecs
ctrb_list <- list("Atl" = c("high", "low", "low"),
                  "aTl" = c("low", "high", "low"),
                  "atL" = c("low", "low", "high"),
                  "ATL" = c("high", "high", "high"))

sims <- purrr::map(.x = ctrb_list, ~ simnetfromtap_ctrb(ctrb_vec = .x,
  traits = init_sim$traits_all,
  abuns = init_sim$networks[[1]]$abuns,
  paramsList = init_sim$sim_params,
  pems = init_sim$networks[[1]]$pems,
  tmatch_type_pem = "normal",
  tmatch_type_obs = "normal",
  Nwebs = 2,
  Nobs = 1111))

sims <- clean(sims, single = F)

# visualize ----
for (i in seq(sims)) {
  for (j in seq(sims[[i]])) {
    plotweb(sims[[i]][[j]])
  }
}

# basic metrics ----
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

# extinction model ----

source("rewiring_vizentin-bugoni_2019/Functions/extinction.mod.R")
source("rewiring_vizentin-bugoni_2019/Functions/one.second.extinct.mod.R")
source("rewiring_vizentin-bugoni_2019/Functions/calc.mean.one.second.extinct.mod.R")
source("rewiring_vizentin-bugoni_2019/Functions/matrix.p1.R")
source("rewiring_vizentin-bugoni_2019/Functions/IC.R") # calc of 95 percent confidence interval

  # rewiring probabilites
# pairwise relative abundances
rew_abund <- init_sim$networks[[1]]$abuns$low %*%
  t(init_sim$networks[[1]]$abuns$high)
row.names(rew_abund) <- init_sim$trees[[1]]$tip.label

# relative abundances of lower
rew_abund_low <- sweep(rew_abund, 2, colSums(rew_abund), "/")

# relative abundances of higher
rew_abund_high <- sweep(rew_abund, 1, rowSums(rew_abund), "/")

# fy to simulate extinctions for a nested list of networks; n_sims is no of extc
# simulations, n_nets is no of networks in list, n_webs is no of webs per network
run_extc <- function(web,
                     participant,
                     method,
                     rewiring,
                     probabilities.rewiring1,
                     probabilities.rewiring2,
                     method.rewiring,
                     n_sims,
                     n_nets,
                     n_webs) {
  
  # list for all models
  out <- vector(mode = "list", length = n_nets)
  names(out) <- paste("Net", 1:n_nets)
  
  # iterate over no of nets
  for (i in seq(n_nets)) {
    # iterate over no of webs
    out_temp <- vector(mode = "list", length = n_webs)
    names(out_temp) <- paste("Web", 1:n_webs)
    
    for (j in seq(n_webs)) {
      web_temp <- web[[i]][[j]] # get web
      
      # model extinctions
      res <- replicate(n_sims, simplify = F, 
                       one.second.extinct.mod_aug(web = web_temp,
                                                  participant = participant,
                                                  method = method,
                                                  rewiring = rewiring,
                                                  probabilities.rewiring1 = probabilities.rewiring1,
                                                  probabilities.rewiring2 = probabilities.rewiring2,
                                                  method.rewiring = method.rewiring))
    out_temp[[j]] <- res # add extc sim of webs
      }
  out[[i]] <- out_temp # add extc sim of nets
  }
  return(out)
}

  # run extinction models 
extc_sims <- run_extc(web = sims,
                      participant = "lower",
                      method = "random",
                      rewiring = T,
                      probabilities.rewiring1 = rew_abund_low,
                      probabilities.rewiring2 = rew_abund,
                      method.rewiring = "one.try.single.partner",
                      n_sims = 2,
                      n_nets = 4,
                      n_webs = 2)

# function to calculate percentages of remaining species
per_surv <- function(extc_sims,
                     n_sims,
                     n_nets,
                     n_webs) {
  
  out <- vector(mode = "list", length = n_nets) # set list for lo/hi percentages  
  names(out) <- paste("Net", 1:n_nets)
  
  for(i in seq(n_nets)) {
    
    web_temp <- vector(mode = "list", length = n_webs)
    names(web_temp) <- paste("Web", 1:n_webs)
    
    for(j in seq(n_webs)) {
      
      sim_temp <- vector(mode = "list", length = n_sims)
      names(sim_temp) <- paste("Sim", 1:n_sims)
      
      for(k in seq(n_sims)) {
        lo <- extc_sims[[i]][[j]][[k]][[1]][, 4]/nrow(sims[[i]][[j]])*100
        hi <- extc_sims[[i]][[j]][[k]][[1]][, 5]/ncol(sims[[i]][[j]])*100
        
        sim_temp[[k]] <- data.frame("low" = lo, "high" = hi)
      }
      web_temp[[j]] <- sim_temp
    }
    out[[i]] <- web_temp
  }
  return(out)
}

  # calculate percentages of remaining sp
sp_remain <- per_surv(extc_sims, n_sims = 2, n_nets = 4, n_webs = 2)

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

