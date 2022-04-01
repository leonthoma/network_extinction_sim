# fy for simulating webs with user specified community var contributions
simnetfromtap_ctrb <- function(traits,
                               abuns,
                               paramsList,
                               pems,
                               tmatch_type_pem,
                               tmatch_type_obs,
                               ctrb_vec = c(NULL),
                               initial_sim = FALSE,
                               Nwebs,
                               Nobs) {
  # create networks
  sim <- simnetfromtap_aug(traits = traits,
                           abuns = abuns,
                           paramsList = paramsList,
                           pems = pems,
                           tmatch_type_pem = tmatch_type_pem,
                           tmatch_type_obs = tmatch_type_obs,
                           ctrb_vec = ctrb_vec,
                           initial_sim = initial_sim)
  # get web from Imat
  simweb <- matrix(rmultinom(1, Nobs, sim$I_mat), nrow = nrow(sim$I_mat),
                   ncol = ncol(sim$I_mat))
  # set sp names
  dimnames(simweb) <- dimnames(init_sim$networks[[1]]$I_mat)
  
  return(list("web" = simweb,
              "I_mat" = sim$I_mat))
}
