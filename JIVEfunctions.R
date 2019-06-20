# Wrapper for JIVE function
# X - n x p concatenated data matrix
# pvec - values p_1,....,p_d corresponding to the number of measurements within each data type
# trueU, trueV - true parameters from the simulated model
applyJIVE <- function(X, pvec, trueU, trueV, orthIndiv = T){
  pcum <- cumsum(pvec)
  D <- length(pvec)

  # Standardize the data in the same way as in RankLassoR
  out <- standardizeX(data$X, data$pvec, center = T)
  X <- out$X
  norms <- out$norms

  # Transfer X into list form needed for JIVE
  data <- list()
  for (d in 1:D){
    # Index the measurements corresponding to dth dataset
    if (d == 1){
      index <- c(1:pvec[1])
    }else{
      index <- c((pcum[d-1]+1):pcum[d])
    }
    # Allocate to corresponding list entry, columns are subjects for JIVE
    data[[d]] <- t(X[,index])
  }
  # Apply JIVE, permutation method for rank selection, orthogonality between individual structures is  enforced, est uses SVD for dimension reduction, scale is not used to the results agree with rank approach
  jive.out <- jive(data, method = "perm", scale = F, orthIndiv = orthIndiv, est = T)
  #jive.out <- jive_mine(data, method = "perm", scale = F, orthIndiv = T, est = T)
  # Figure out the error using trueU and true V
  distUV = 0
  for (d in 1:D){
    # Adjust for scale
    estimate = t((jive.out$joint[[d]]+jive.out$individual[[d]])*sqrt(norms[d]))
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
  return(list(error = distUV, r0 = jive.out$rankJ, r1 = jive.out$rankA[1], r2 = jive.out$rankA[2], converged = jive.out$converged))
}

# Wrapper for JIVE function
# X - n x p concatenated data matrix
# pvec - values p_1,....,p_d corresponding to the number of measurements within each data type
# trueU, trueV - true parameters from the simulated model
applyJIVE_d3 <- function(X, pvec, trueU, trueV, orthIndiv = T){
  pcum <- cumsum(pvec)
  D <- length(pvec)

  # Standardize the data in the same way ad RankLassoR
  out <- standardizeX(data$X, data$pvec, center = T)
  X <- out$X
  norms <- out$norms

  # Transfer X into list form needed for JIVE
  data <- list()
  for (d in 1:D){
    # Index the measurements corresponding to dth dataset
    if (d == 1){
      index <- c(1:pvec[1])
    }else{
      index <- c((pcum[d-1]+1):pcum[d])
    }
    # Allocate to corresponding list entry, columns are subjects for JIVE
    data[[d]] <- t(X[,index])
  }
  # Apply JIVE, permutation method for rank selection, orthogonality between individual structures is  enforced, est uses SVD for dimension reduction, scale is not used to the results agree with rank approach
  jive.out <- jive(data, method = "perm", scale = F, orthIndiv = orthIndiv, est = T)
  #jive.out <- jive_mine(data, method = "perm", scale = F, orthIndiv = T, est = T)
  # Figure out the error using trueU and true V
  distUV = 0
  for (d in 1:D){
    # Adjust for scale
    estimate = t((jive.out$joint[[d]]+jive.out$individual[[d]])*sqrt(norms[d]))
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
  return(list(error = distUV, r0 = jive.out$rankJ, r1 = jive.out$rankA[1], r2 = jive.out$rankA[2],r3 = jive.out$rankA[3], converged = jive.out$converged))
}


# Wrapper for JIVE function - given ranks
# X - n x p concatenated data matrix
# pvec - values p_1,....,p_d corresponding to the number of measurements within each data type
# trueU, trueV - true parameters from the simulated model
applyJIVE_given <- function(X, pvec, trueU, trueV, r0 = 2, rvec = c(2,2), orthIndiv = T){
  pcum <- cumsum(pvec)
  D <- length(pvec)

  # Standardize the data in the same way ad RankLassoR
  out <- standardizeX(data$X, data$pvec, center = T)
  X <- out$X
  norms <- out$norms

  # Transfer X into list form needed for JIVE
  data <- list()
  for (d in 1:D){
    # Index the measurements corresponding to dth dataset
    if (d == 1){
      index <- c(1:pvec[1])
    }else{
      index <- c((pcum[d-1]+1):pcum[d])
    }
    # Allocate to corresponding list entry, columns are subjects for JIVE
    data[[d]] <- t(X[,index])
  }
  # Apply JIVE, permutation method for rank selection, orthogonality between individual structures is  enforced, est uses SVD for dimension reduction, scale is not used to the results agree with rank approach
  jive.out <- jive(data, method = "given", rankJ = r0, rankA = rvec, scale = F, est = T, orthIndiv = orthIndiv)
  #jive.out <- jive_mine2(data, method = "given", rankJ = r0, rankA = rvec, scale = F, orthIndiv = T, est = T)
  # Figure out the error using trueU and true V
  distUV = 0
  for (d in 1:D){
    # Adjust for scale
    estimate = t((jive.out$joint[[d]] + jive.out$individual[[d]])*sqrt(norms[d]))
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
  return(list(error = distUV, r0 = jive.out$rankJ, r1 = jive.out$rankA[1], r2 = jive.out$rankA[2], converged = jive.out$converged))
}

# Wrapper for ADHOC approach: PCA on concatenated data, then PCA on the residuals
# X - n x p concatenated data matrix
# pvec - values p_1,....,p_d corresponding to the number of measurements within each data type
# trueU, trueV - true parameters from the simulated model
# r0 - joint rank to consider for one step PCA, default value is true rank
# rvec - individual ranks to consider for one step PCA, default values are true ranks
apply_onestep_PCA <- function(X, pvec, trueU, trueV, r0 = 2, rvec = c(2,2)){
  # Standardize the data in the same way ad RankLassoR
  out <- standardizeX(data$X, data$pvec, center = T)
  X <- out$X
  norms <- out$norms

  pcum <- cumsum(pvec)
  D <- length(pvec)

  # Apply one step PCA on concatenated dataset to get the joint structure out
  if (r0 > 0){
    X_svd <- svd(X)
    # Get first r0 components for joint structure
    joint <- X_svd$u[,1:r0]%*%diag(X_svd$d[1:r0])%*%t(X_svd$v[, 1:r0])
  }else{
    joint <- matrix(0, nrow(X), ncol(X))
  }
  # Form residual dataset
  X_res <- X - joint

  # Apply one step PCA to the residuals, and evaluate the error
  distUV = 0
  for (d in 1:D){
    # Index dth dataset
    if (d == 1){
      index <- c(1:pvec[1])
    }else{
      index <- c((pcum[d-1]+1):pcum[d])
    }
    # Apply one-step PCA to the residual matrix
    Xi_res_svd <- svd(X_res[, index])
    # Get first rvec[d] components for individual structure
    if (rvec[d] > 1){
      individual <- Xi_res_svd$u[,1:rvec[d]]%*%diag(Xi_res_svd$d[1:rvec[d]])%*%t(Xi_res_svd$v[, 1:rvec[d]])
    }else if (rvec[d]==1){
      individual <- Xi_res_svd$u[,1]%*%t(Xi_res_svd$v[,1]) * Xi_res_svd$d[1]
    }else{
      individual <- matrix(0, nrow(X), pvec[d])
    }

    # Adjust for scale, and calculate the error
    estimate = (joint[, index] + individual)*sqrt(norms[d])
    # Calculate Frobenius norm squared of the signal part
    signal = sum(tcrossprod(trueU,trueV[index,])^2)
    # If nonzero, calculate error relative to signal otherwise just directly
    if (signal > 0){
      distUV <- distUV + sum((tcrossprod(trueU,trueV[index,])-estimate)^2)/sum(tcrossprod(trueU,trueV[index,])^2)
    }else{
      distUV <- distUV + sum(estimate^2)
    }
  }
  return(list(error = distUV, r0 = r0, r1 = rvec[1], r2 = rvec[2]))
}


# Wrapper for ADHOC approach: PCA on concatenated data, then PCA on the residuals
# X - n x p concatenated data matrix
# pvec - values p_1,....,p_d corresponding to the number of measurements within each data type
# trueU, trueV - true parameters from the simulated model
# r0 - joint rank to consider for one step PCA, default value is true rank
# rvec - individual ranks to consider for one step PCA, default values are true ranks
# rshared - semi-shared ranks between two datasets out of three
apply_onestep_PCA_d3 <- function(X, pvec, trueU, trueV, r0 = 2, rvec = c(2,2,2), rshared = 2){
  # Standardize the data in the same way ad RankLassoR
  out <- standardizeX(data$X, data$pvec, center = T)
  X <- out$X
  norms <- out$norms

  pcum <- cumsum(pvec)
  D <- length(pvec)

  # Apply one step PCA on concatenated dataset to get the joint structure out
  X_svd <- svd(X)
  # Get first r0 components for joint structure
  joint <- X_svd$u[,1:r0]%*%diag(X_svd$d[1:r0])%*%t(X_svd$v[, 1:r0])

  # Form residual dataset
  X_res <- X - joint

  # Apply PCA on the shared parts across two datasets, order 12, 13, 23
  # Shared 12
  Xi_res_svd <- svd(X_res[, 1:pcum[2]])
  # Get first rvec[d] components for individual structure
  shared12 <- Xi_res_svd$u[,1:rshared]%*%diag(Xi_res_svd$d[1:rshared])%*%t(Xi_res_svd$v[, 1:rshared])
  X_res[, 1:pcum[2]] <- X_res[, 1:pcum[2]] - shared12
  # Shared 13
  Xi_res_svd <- svd(X_res[, c(1:pvec[1], (pcum[2]+1):pcum[3])])
  # Get first rvec[d] components for individual structure
  shared13 <- Xi_res_svd$u[,1:rshared]%*%diag(Xi_res_svd$d[1:rshared])%*%t(Xi_res_svd$v[, 1:rshared])
  X_res[, c(1:pvec[1], (pcum[2]+1):pcum[3])] <- X_res[, c(1:pvec[1], (pcum[2]+1):pcum[3])] - shared13
  # Shared 23
  Xi_res_svd <- svd(X_res[, (pvec[1]+1):pcum[3]])
  # Get first rvec[d] components for individual structure
  shared23 <- Xi_res_svd$u[,1:rshared]%*%diag(Xi_res_svd$d[1:rshared])%*%t(Xi_res_svd$v[, 1:rshared])
  X_res[, (pvec[1]+1):pcum[3]] <- X_res[, (pvec[1]+1):pcum[3]] - shared23

  # Apply one step PCA to the residuals, and evaluate the error
  distUV = 0
  for (d in 1:D){
    # Index dth dataset
    if (d == 1){
      index <- c(1:pvec[1])
    }else{
      index <- c((pcum[d-1]+1):pcum[d])
    }
    # Apply one-step PCA to the residual matrix
    Xi_res_svd <- svd(X_res[, index])
    # Get first rvec[d] components for individual structure
    individual <- Xi_res_svd$u[,1:rvec[d]]%*%diag(Xi_res_svd$d[1:rvec[d]])%*%t(Xi_res_svd$v[, 1:rvec[d]])

    # Adjust for scale, and calculate the error
    if (d == 1){
      estimate = (joint[, index] + shared12[, index] + shared13[, index] + individual)*sqrt(norms[d])
    }else if (d == 2){
      estimate = (joint[, index] + shared12[, index] + shared23[, 1:pvec[2]] + individual)*sqrt(norms[d])
    }else{
      # d = 3
      estimate = (joint[, index]  + shared13[, (pvec[1]+1):(pvec[1]+pvec[3])] + shared23[, (pvec[2]+1):(pvec[2]+pvec[3])] + individual)*sqrt(norms[d])
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
  return(list(error = distUV, r0 = r0, r1 = rvec[1], r2 = rvec[2]))
}


