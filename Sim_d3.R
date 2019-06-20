source("Models_paper.R")
source("AuxillaryFunctions.R")
source("SLIDEfunctions.R")
# Necessary functions for JIVE
library(r.jive)
source("JIVEfunctions.R")
# Necessary functions for GFA
library(CCAGFA)
source("GFAfunctions.R")

Nrep=100
library(doParallel)
library(doRNG)
nworkers <- detectCores()
cl <- makePSOCKcluster(nworkers)
registerDoParallel(cl)

# Our approach
set.seed(2940765)
results_d3_struct <- foreach(i=1:Nrep, .packages = c("r.jive", "CCAGFA")) %dorng% {
  # generate the model
  data <- generateModel2(n = 100, pvec = c(100,100,100), snr = 1, orthogonalV = T)
  # Prepare output
  output <- list(slide = NA, jive = NA, gfa = NA, onestep = NA, sca_best = NA)
  # apply SLIDE
  output$SLIDE <- compareBCV_vs_Best_structure_refit_center(data$X, pvec = data$pvec, trueU = data$U, trueV = data$V, nl = 50, center = T)
  # apply one step PCA
  output$onestep <- apply_onestep_PCA_d3(data$X, pvec = data$pvec, trueU = data$U, trueV = data$V, r0 = 2, rvec = c(2,2,2), rshared = 2)
  # apply JIVE
  output$jive = applyJIVE_d3(data$X, pvec = data$pvec, trueU = data$U, trueV = data$V, orthIndiv = F)
  # apply Klami et al. with K being equal to true choice
  output$gfa = applyGFA(data$X, pvec = data$pvec, trueU = data$U, trueV = data$V, K = 14)
  # apply SCA_best
  # Scale the original dataset
  out_s <- standardizeX(data$X, data$pvec, center = T)
  X <- out_s$X
  svec <- out_s$svec
  pcum <- cumsum(data$pvec)
  d <- length(data$pvec)
  
  # Solve for each lambda
  out <- solve_optim1_seq(X = X, pvec = data$pvec, n_lambda = 50, reduced = F)
  lambda_all <- out$lambda
  
  ## Calculate weighted residual error for each lambda
  errors <- getAllErrors(out, data$U, data$V, data$X, data$pvec, out_s$norms)
  best_id <- which.min(errors)
  lambda_best <- lambda_all[best_id]
  
  # Calculate ranks for each lambda
  if (length(data$pvec) == 2){
    ranks <- getAllranks_d2(out, data$pvec)
  }else{
    ranks <- getAllranks_d3(out, data$pvec)
  }
  
  # end output
  output$sca_best = list(lambda_all = lambda_all, error = errors, ranks = ranks, best_error = errors[best_id], best_ranks = ranks[[best_id]])
  
  output
}
save(results_d3_struct, file = "Results_d3_struct.Rda")


