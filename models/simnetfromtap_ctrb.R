# fy for simulating webs with user specified community var contributions
simnetfromtap_ctrb <- function(traits,
                               abuns,
                               paramsList,
                               pems,
                               tmatch_type_pem,
                               tmatch_type_obs,
                               ctrb_vec = c(NULL),
                               Nwebs,
                               Nobs) {
  out <- vector(mode = "list", length = Nwebs)
  names(out) <- paste("Web", 1:Nwebs)
  
  for (i in 1:Nwebs) {
    sim <- simnetfromtap_aug(traits = traits,
                             abuns = abuns,
                             paramsList = paramsList,
                             pems = pems,
                             tmatch_type_pem = tmatch_type_pem,
                             tmatch_type_obs = tmatch_type_obs,
                             ctrb_vec = ctrb_vec)
    simweb <- matrix(rmultinom(1, Nobs, sim), nrow = nrow(sim),
                     ncol = ncol(sim))
    dimnames(simweb) <- dimnames(init_sim$networks[[1]]$I_mat)
    
    out[[i]] <- simweb
  }
  return(out)
}