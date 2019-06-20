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

getAllErrors <- function(out, trueU, trueV, X, pvec, norms){
  nl = length(out$lambda)
  errors = rep(NA, nl)
  for (l in 1:nl){
    distUV = 0
    Vadj <- out$param[[l]]$V
    for (i in 1:d){
      if (i == 1){
        index <- c(1:pvec[1])
      }else{
        index <- c((pcum[i-1]+1):pcum[i])
      }
      Vadj[index,] <- out$param[[l]]$V[index,]*sqrt(norms[i])
      # Calculate Frobenius norm squared of the signal part
      signal = sum(tcrossprod(trueU,trueV[index,])^2)
      # If nonzero, calculate error relative to signal otherwise just directly
      if (signal > 0){
        distUV <- distUV + sum((tcrossprod(trueU,trueV[index,])-tcrossprod(out$param[[l]]$U, Vadj[index,]))^2)/sum(tcrossprod(trueU,trueV[index,])^2)
      }else{
        distUV <- distUV + sum((tcrossprod(out$param[[l]]$U, Vadj[index,]))^2)
      }
      
    }
    errors[l] = distUV
  }
  errors
}

set.seed(192375)
results_d2_AC <- foreach(i=1:Nrep, .packages = c("r.jive", "CCAGFA")) %dorng% {
  # generate the model
  data <- generateModel1_AC(n = 100, pvec = c(25,150), c = 0.8, snr = 1, orthogonalV = T)
  # Prepare output
  output <- list(slide = NA, jive = NA, jive_best = NA, onestep = NA, gfa = NA, sca_best = NA)
  # apply SLIDE
  output$slide = compareBCV_vs_Best_structure_refit_center(data$X, pvec = data$pvec, trueU = data$U, trueV = data$V, nl = 50, center = T)
  # apply JIVE
  output$jive = applyJIVE(data$X, pvec = data$pvec, trueU = data$U, trueV = data$V, orthIndiv = F)
  # apply JIVE with true ranks for JIVE model
  output$jive_best = applyJIVE_given(data$X, pvec = data$pvec, trueU = data$U, trueV = data$V, r0 = 0, rvec = c(1,1))
  # apply one step PCA
  output$onestep = apply_onestep_PCA(data$X, pvec = data$pvec, trueU = data$U, trueV = data$V, r0 = 0, rvec = c(1,1))
  # apply Klami et al. with K being equal to true choice
  output$gfa = applyGFA(data$X, pvec = data$pvec, trueU = data$U, trueV = data$V, K = 2)
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
save(results_d2_AC, file = "Results_d2_AC.Rda")



