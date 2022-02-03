# first test using tapnet to simulate network
setwd("~/Documents/Uni/M.sc/Master Thesis/Networks/models/")

#library(tapnet)
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

# simulation
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

# initial simulation
sim <- simulate_tapnet_aug(nlower = 26, nhigher = 14, ntraits_nopem = 4,
                       ntraits_pem = 2, abuns = "lognormal", Nobs = 1111,
                       names = sp_names)
dimnames(sim$networks[[1]]$web) <- dimnames(sim$networks[[1]]$I_mat)
# use parms form initial sim

# abun
sim1 <- simnetfromtap_aug(traits = sim$traits_all,
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

# extinction model

source("rewiring_vizentin-bugoni_2019/Functions/extinction.mod.R")
source("rewiring_vizentin-bugoni_2019/Functions/one.second.extinct.mod.R")
source("rewiring_vizentin-bugoni_2019/Functions/calc.mean.one.second.extinct.mod.R")
source("rewiring_vizentin-bugoni_2019/Functions/matrix.p1.R")
source("rewiring_vizentin-bugoni_2019/Functions/IC.R") # calc of 95 percent confidence interval

# remove sp w/o interactions from web
no_int_low <- which(rowSums(sim$networks[[1]]$web) == 0)
no_int_high <- which(colSums(sim$networks[[1]]$web) == 0)

network <- sim$networks[[1]]$web[-no_int_low, -no_int_high]

# rewiring probabilites
# pairwise relative abundances
rew_abund <- sim$networks[[1]]$abuns$low[-no_int_low] %*%
  t(sim$networks[[1]]$abuns$high[-no_int_high])
row.names(rew_abund) <- sim$trees[[1]]$tip.label[-no_int_low]

# relative abundances of lower
rew_abund_low <- sweep(rew_abund, 2, colSums(rew_abund), "/")

# relative abundances of higher
rew_abund_high <- sweep(rew_abund, 1, rowSums(rew_abund), "/")

# run extiction model
extc_sim <- one.second.extinct.mod_aug(web = network,
                                       participant = "lower",
                                       method = "random", 
                                       rewiring = T, 
                                       probabilities.rewiring1 = rew_abund_low, 
                                       probabilities.rewiring2 = rew_abund,
                                       method.rewiring = "one.try.single.partner")

# calculate percentages of remaining sp
extc_lo_per <- extc_sim[[1]][, 4]/nrow(network)
extc_hi_per <- extc_sim[[1]][, 5]/ncol(network)

plot(1-extc_lo_per, extc_hi_per, type = "l", xlab = "animals persisting",
     ylab = "plants removed")
