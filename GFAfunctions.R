# Wrapper for Group Factor Analysis method by Klami et al. 

# Need to specify K - this is the trickiest part as the results depend on the value of K
applyGFA <- function(X, pvec, trueU, trueV, K = NULL){
  pcum <- cumsum(pvec)
  D <- length(pvec)
  
  # Standardize the data in the same way as other methods
  out <- standardizeX(X, pvec, center = T)
  
  # Transfer X into list form needed for GFA
  data <- list()
  for (d in 1:D){
    # Index the measurements corresponding to dth dataset
    if (d == 1){
      index <- c(1:pvec[1])
    }else{
      index <- c((pcum[d-1]+1):pcum[d])
    }
    # Allocate to corresponding list entry, columns are subjects for JIVE
    data[[d]] <- out$X[,index]
  }
  norms <- out$norms
  
  # Apply GFA
  opts = getDefaultOpts()
  if (is.null(K)){ # Set the number of components to be minimum (n, p_i)
    K = min(nrow(X), min(pvec))
  }
  out_gfa = GFAexperiment(Y = data, K = K, opts = opts, Nrep=10)
  out_gfa_trim = GFAtrim(out_gfa)

  # Calculate the ranks
  structure = out_gfa_trim$active
  
  # Calculate the error using trueU and true V
  distUV = 0
  for (d in 1:D){
    estimate = (out_gfa_trim$Z %*% t(out_gfa_trim$W[[d]]))*sqrt(norms[d])
    # Index dth dataset
    if (d == 1){
      index <- c(1:pvec[1])
    }else{
      index <- c((pcum[d-1]+1):pcum[d])
    }
    # Calculate Frobenius norm squared of the signal part
    signal = sum(tcrossprod(trueU,trueV[index,])^2)
    # If nonzero, calculate error relative to signal otherwise just directly
    if (signal > 0){
      distUV <- distUV + sum((tcrossprod(trueU,trueV[index,])-estimate)^2)/sum(tcrossprod(trueU,trueV[index,])^2)
    }else{
      distUV <- distUV + sum(estimate^2)
    }
  }
  return(list(error = distUV, structure = structure))
}