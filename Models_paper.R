# Models setup for SLIDE paper
# Last modified by Irina Gaynanova on 04/26/2019

######################################################################################
# 2 datasets
# n - sample size
# pvec - sizes for each dataset
# cvec - vector of relative scaling
# snr - signal to noise ratio
# orthogonalV - whether the loadings should be generated to be orthogonal to each other, true by default
#
# both U and V are genereated using U[0,1] with subsequent orthogonalization
# noise E is generated using N(0, sigma_d^2) separately for each dataset
# ranks are set as r0=r1=r2=2
#
#########################################################################
generateModel1 <- function(n = 100, pvec = c(25,25), cvec = c(1,1), snr = 1, orthogonalV = T){
  # Fix number of datasets and all the ranks
  D <- 2
  # joint rank
  r0 <- 2
  # individual ranks
  rvec <- c(2,2)
  # Total rank
  r_total <- sum(rvec) + r0
  # Total number of variables
  p_total <- sum(pvec)

  # Space for the whole concatenated X
  X <- matrix(0, n, p_total)

  # Generate all the scores
  U = matrix(runif(n*r_total), n, r_total)
  # Column-center U so that the overall mean is zero
  U = scale(U, scale = F)

  # Orthogonalize the scores
  U <- qr.Q(qr(U))

  # Generate all the loadings
  Vtmp <- list()
  # Generate dataset-specific loadings
  for (d in 1:D){
    Vtmp[[d]] <- matrix(runif(pvec[d]*(r0+rvec[d])), pvec[d],(r0+rvec[d]))
    # orthogonalize the loadings
    Vtmp[[d]] <- qr.Q(qr(Vtmp[[d]]))
  }
  # Put everything together into one V
  V <- matrix(0, nrow = p_total, ncol = r_total)
  for (d in 1:D){
    # Row index for dth dataset
    if (d == 1){
      index <- c(1:pvec[1])
    }else{
      index <- c((pvec[1]+1):sum(pvec))
    }
    # Column index for dth dataset
    if (d==1){
      col_indexd = c(1:(r0+rvec[1]))
    }else{
      col_indexd = c(1:r0,(r0+rvec[1]+1):(r0+sum(rvec)))
    }
    V[index,col_indexd]=Vtmp[[d]]
  }
  # Orthogonalize V
  if (orthogonalV == T){
    V <- qr.Q(qr(V))
  }
  # Add the scale between joint and individual by fixing singular values
  sjoint <- c(1.5,1.3)
  s1ind <- c(1, 0.8)
  s2ind <- c(1, 0.7)
  V <- V%*%diag(c(sjoint,s1ind,s2ind))

  # Add the dataset-specific scale
  for (d in 1:D){
    # Row index
    if (d == 1){
      index <- c(1:pvec[1])
    }else{
      index <- c((pvec[1]+1):sum(pvec))
    }
    V[index,] <- cvec[d]*V[index,]
  }

  # Form X
  X <- tcrossprod(U, V)

  # Add noise
  if (snr != 0){
    if ((cvec[1]==cvec[2])&(pvec[1]==pvec[2])){
      # Datasets on the same scale
      sigma <- sqrt(sum(X^2) / (n * sum(pvec) * snr))
      # Add the noise part0.3*
      X <- X + matrix(rnorm(n * p_total, sd = sigma), n, p_total)
      # Combine joint loadings and individual loadings into one matrix V
      return(list(X = X, pvec = pvec, r0 = r0, rvec = rvec, sigma = sigma, U = U, V = V, cvec = cvec, snr = snr))
    }else{
      # Datasets have different scales or different dimensions
      index1 <- 1:pvec[1]
      sigma1 <- sqrt(sum(X[,index1]^2) / (n * pvec[1] * snr))
      X[,index1] <- X[,index1] + matrix(rnorm(n * pvec[1], sd = sigma1), n, pvec[1])
      index2 <- (pvec[1]+1):sum(pvec)
      sigma2 <- sqrt(sum(X[,index2]^2) / (n * pvec[2] * snr))
      X[,index2] <- X[,index2] + matrix(rnorm(n * pvec[2], sd = sigma2), n, pvec[2])
      # Combine joint loadings and individual loadings into one matrix V
      return(list(X = X, pvec = pvec, r0 = r0, rvec = rvec, sigma1 = sigma1, sigma2 = sigma2, U = U, V = V, cvec = cvec, snr = snr))
    }
  }else{
    # No noise
    # Combine joint loadings and individual loadings into one matrix V
    return(list(X = X, pvec = pvec, r0 = r0, rvec = rvec, sigma = 0, U = U, V = V, cvec = cvec, snr = snr))
  }
}

######################################################################################
# 2 datasets
# n - sample size
# pvec - sizes for each dataset
# cvec - vector of relative scaling
# snr - signal to noise ratio
# orthogonalV - whether the loadings should be generated to be orthogonal to each other, true by default
#
# both U and V are genereated using U[0,1] with subsequent orthogonalization
# noise E is generated using N(0, sigma_d^2) separately for each dataset
# ranks are set as r0=1, r1=r2=1 for SLIDE
#########################################################################
generateModel1_AC <- function(n = 100, pvec = c(25,150), c = 0.8, snr = 1, orthogonalV = T){
  # Fix number of datasets and all the ranks
  D <- 2
  # joint rank
  r0 <- 0
  # individual ranks
  rvec <- c(1,1)
  # Total rank
  r_total <- sum(rvec) + r0
  # Total number of variables
  p_total <- sum(pvec)

  # Space for the whole concatenated X
  X <- matrix(0, n, p_total)

  # Generate all the scores
  U = matrix(runif(n*r_total), n, r_total)
  # Column-center U so that the overall mean is zero
  U = scale(U, scale = F)

  # Orthogonalize the scores
  U <- qr.Q(qr(U))
  # Create the scores with pre-specified inner product
  u1 <- U[,1]
  tmp <- U[,2]
  alpha <- c/(c+ sqrt(1-c^2))
  u2 <- (alpha*u1 + (1-alpha)*tmp)/sqrt(alpha^2 + (1-alpha)^2)
  U <- cbind(u1, u2)

  # Generate all the loadings
  Vtmp <- list()
  # Generate dataset-specific loadings
  for (d in 1:D){
    Vtmp[[d]] <- matrix(runif(pvec[d]), pvec[d],1)
  }
  # Put everything together into one V
  V <- matrix(0, nrow = p_total, ncol = r_total)
  for (d in 1:D){
    # Row index for dth dataset
    if (d == 1){
      index <- c(1:pvec[1])
    }else{
      index <- c((pvec[1]+1):sum(pvec))
    }
    # Column index for dth dataset
    if (d==1){
      col_indexd = c(1:(r0+rvec[1]))
    }else{
      col_indexd = c((r0+rvec[1]+1):(r0+sum(rvec)))
    }
    V[index,col_indexd]=Vtmp[[d]]
  }
  # Orthogonalize V
  if (orthogonalV == T){
    V <- qr.Q(qr(V))
  }
  # Add the scale between joint and individual by fixing singular values

  # Form X
  X <- tcrossprod(U, V)

  # Add noise
  if (snr != 0){
    if (pvec[1]==pvec[2]){
      # Datasets on the same scale
      sigma <- sqrt(sum(X^2) / (n * sum(pvec) * snr))
      # Add the noise part0.3*
      X <- X + matrix(rnorm(n * p_total, sd = sigma), n, p_total)
      # Combine joint loadings and individual loadings into one matrix V
      return(list(X = X, pvec = pvec, r0 = r0, rvec = rvec, sigma = sigma, U = U, V = V, c = c, snr = snr))
    }else{
      # Datasets have different scales or different dimensions
      index1 <- 1:pvec[1]
      sigma1 <- sqrt(sum(X[,index1]^2) / (n * pvec[1] * snr))
      X[,index1] <- X[,index1] + matrix(rnorm(n * pvec[1], sd = sigma1), n, pvec[1])
      index2 <- (pvec[1]+1):sum(pvec)
      sigma2 <- sqrt(sum(X[,index2]^2) / (n * pvec[2] * snr))
      X[,index2] <- X[,index2] + matrix(rnorm(n * pvec[2], sd = sigma2), n, pvec[2])
      # Combine joint loadings and individual loadings into one matrix V
      return(list(X = X, pvec = pvec, r0 = r0, rvec = rvec, sigma1 = sigma1, sigma2 = sigma2, U = U, V = V, c = c, snr = snr))
    }
  }else{
    # No noise
    # Combine joint loadings and individual loadings into one matrix V
    return(list(X = X, pvec = pvec, r0 = r0, rvec = rvec, sigma = 0, U = U, V = V, c = c, snr = snr))
  }
}

######################################################################################
# 3 datasets
# n - sample size
# pvec - sizes for each dataset
# snr - signal to noise ratio
# orthogonalV - whether the loadings should be generated to be orthogonal to each other, true by default
#
# all ranks are set to 2
#########################################################################
generateModel2 <- function(n = 100, pvec = c(100,100,100), snr = 1, orthogonalV = T){
  p_total <- sum(pvec)
  pcum <- cumsum(pvec)
  D <- 3
  r0 <- 2
  rmat <- matrix(2,3,3)
  cvec <- c(1,1,1)
  r_total <- r0 + sum(rmat[upper.tri(rmat,diag=T)])

  # generate noise part of X
  X <- matrix(0, n, p_total)

  # generate all the scores
  U <- matrix(runif(n*r_total), n, r_total)
  # Column-center U so that the overall mean is zero
  U <- scale(U, scale = F)

  # orthogonalize the scores
  U <- qr.Q(qr(U))

  # Generate all the loadings
  Vtmp <- list()
  # generate dataset-specific Vs
  for (d in 1:D){
    Vtmp[[d]] <- matrix(runif(pvec[d]*(r0+sum(rmat[d,]))), pvec[d],(r0+sum(rmat[d,])))
    # orthogonalize the loadings
    Vtmp[[d]] <- qr.Q(qr(Vtmp[[d]]))
  }
  # Put everything together into one V
  V <- matrix(0, nrow = p_total, ncol = r_total)
  for (d in 1:D){
    # Row index
    if (d == 1){
      index <- c(1:pvec[1])
    }else{
      index <- c((pcum[d-1]+1):pcum[d])
    }

    # Column index, order is r0, r12, r13, r23, r1, r2, r3
    if (d==1){
      col_indexd = c(1:(r0+rmat[1,2]+rmat[1,3]),(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]))
    }else if (d==2){
      col_indexd = c(1:(r0+rmat[1,2]),(r0+rmat[1,2]+rmat[1,3]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]), (r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+rmat[2,2]))
    }else{
      col_indexd = c(1:r0,(r0+rmat[1,2]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]),(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+rmat[2,2]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+rmat[2,2]+rmat[3,3]))
    }
    V[index,col_indexd]=Vtmp[[d]]
  }
  # Orthogonallize V
  if (orthogonalV == T){
    V <- qr.Q(qr(V))
  }
  # Add the scale between joint, semi-joint and individual by fixing singular values
  # Column index, order is r0, r12, r13, r23, r1, r2, r3
  sjoint <- c(1.5, 1.3)
  s12 <- c(1, 0.8)
  s13 <- c(1, 0.7)
  s23 <- c(1, 0.5)
  s1ind <- c(1.2, 0.5)
  s2ind <- c(0.9, 0.8)
  s3ind <- c(0.5, 0.4)
  V <- V%*%diag(c(sjoint,s12, s13, s23, s1ind, s2ind, s3ind))


  # Add the dataset-specific scale
#   for (d in 1:D){
#     # Row index
#     if (d == 1){
#       index <- c(1:pvec[1])
#     }else{
#       index <- c((pcum[d-1]+1):pcum[d])
#     }
#     V[index,] <- cvec[d]*V[index,]
#   }

  # Form X
  X <- tcrossprod(U, V)
  if (snr != 0){
    # Datasets have different scales or different dimensions
    sigmavec = rep(0, D)
    for (d in 1:D){
      # Row index for dth dataset
      if (d == 1){
        index <- c(1:pvec[1])
      }else{
        index <- c((pcum[d-1]+1):pcum[d])
      }
      sigmavec[d] <- sqrt(sum(X[,index]^2) / (n * pvec[d] * snr))
      X[,index] <- X[,index] + matrix(rnorm(n * pvec[d], sd = sigmavec[d]), n, pvec[d])
    }
    # Combine joint loadings and individual loadings into one matrix V
    return(list(X = X, pvec = pvec, r0 = r0, rmat = rmat, sigmavec = sigmavec, U = U, V = V, snr = snr))
  }else{
    sigma <- 0
    # Combine joint loadings and individual loadings into one matrix V
    return(list(X = X, pvec = pvec, r0 = r0, rmat=rmat, sigma = sigma, U = U, V = V, snr = snr))
  }
}



######################################################################################
# 3 datasets, one dataset is purely noise
# n - sample size
# pvec - sizes for each dataset
# cvec - vector of relative scaling, currently only support equal scaling (because of how noise is added)
# snr - signal to noise ratio
# orthogonalV - whether the loadings should be generated to be orthogonal to each other, true by default
#
# same setup as generateModel2, but 3rd dataset has only noise
#########################################################################
generateModel2_noise <- function(n = 100, pvec = c(100,100,100), snr = 1, orthogonalV = T){
  p_total <- sum(pvec)
  pcum <- cumsum(pvec)
  D <- 3
  r0 <- 2
  rmat <- matrix(2,3,3)
  cvec <- c(1,1,1)
  r_total <- r0 + sum(rmat[upper.tri(rmat,diag=T)])
  
  # generate noise part of X
  X <- matrix(0, n, p_total)
  
  # generate all the scores
  U <- matrix(runif(n*r_total), n, r_total)
  # Column-center U so that the overall mean is zero
  U <- scale(U, scale = F)
  
  # orthogonalize the scores
  U <- qr.Q(qr(U))
  
  # Generate all the loadings
  Vtmp <- list()
  # generate dataset-specific Vs
  for (d in 1:D){
    Vtmp[[d]] <- matrix(runif(pvec[d]*(r0+sum(rmat[d,]))), pvec[d],(r0+sum(rmat[d,])))
    # orthogonalize the loadings
    Vtmp[[d]] <- qr.Q(qr(Vtmp[[d]]))
  }
  # Put everything together into one V
  V <- matrix(0, nrow = p_total, ncol = r_total)
  for (d in 1:D){
    # Row index
    if (d == 1){
      index <- c(1:pvec[1])
    }else{
      index <- c((pcum[d-1]+1):pcum[d])
    }
    
    # Column index, order is r0, r12, r13, r23, r1, r2, r3
    if (d==1){
      col_indexd = c(1:(r0+rmat[1,2]+rmat[1,3]),(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]))
    }else if (d==2){
      col_indexd = c(1:(r0+rmat[1,2]),(r0+rmat[1,2]+rmat[1,3]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]), (r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+rmat[2,2]))
    }else{
      col_indexd = c(1:r0,(r0+rmat[1,2]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]),(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+rmat[2,2]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]+rmat[1,1]+rmat[2,2]+rmat[3,3]))
    }
    V[index,col_indexd]=Vtmp[[d]]
  }
  # Orthogonallize V
  if (orthogonalV == T){
    V <- qr.Q(qr(V))
  }
  # Add the scale between joint, semi-joint and individual by fixing singular values
  # Column index, order is r0, r12, r13, r23, r1, r2, r3
  sjoint <- c(1.5, 1.3)
  s12 <- c(1, 0.8)
  s13 <- c(1, 0.7)
  s23 <- c(1, 0.5)
  s1ind <- c(1.2, 0.5)
  s2ind <- c(0.9, 0.8)
  s3ind <- c(0.5, 0.4)
  V <- V%*%diag(c(sjoint,s12, s13, s23, s1ind, s2ind, s3ind))
  
  # Get rid of signal part in the 3rd dataset
  index <- c((pcum[2]+1):pcum[3])
  V[index, ] <- 0
  
  # Form X
  X <- tcrossprod(U, V)
  if (snr != 0){
    # Datasets have different scales or different dimensions
    sigmavec = rep(0, D)
    for (d in 1:D){
      # Row index for dth dataset
      if (d == 1){
        index <- c(1:pvec[1])
      }else{
        index <- c((pcum[d-1]+1):pcum[d])
      }
      # Add noise, but for the last one need to do something special
      if (d != 3){
        sigmavec[d] <- sqrt(sum(X[,index]^2) / (n * pvec[d] * snr))
        X[,index] <- X[,index] + matrix(rnorm(n * pvec[d], sd = sigmavec[d]), n, pvec[d])
      }else{
        # 3d dataset, purely noise
        sigmavec[d] <- mean(sigmavec[1:2])
        X[,index] <- matrix(rnorm(n * pvec[d], sd = sigmavec[d]), n, pvec[d])
      }  
    }
    # Combine joint loadings and individual loadings into one matrix V
    return(list(X = X, pvec = pvec, r0 = r0, rmat = rmat, sigmavec = sigmavec, U = U, V = V, snr = snr))
  }else{
    sigma <- 0
    # Combine joint loadings and individual loadings into one matrix V
    return(list(X = X, pvec = pvec, r0 = r0, rmat=rmat, sigma = sigma, U = U, V = V, snr = snr))
  }
}


######################################################################################
# 3 datasets
# n - sample size
# pvec - sizes for each dataset
# snr - signal to noise ratio
# orthogonalV - whether the loadings should be generated to be orthogonal to each other, true by default
#
# ranks are set as r0=10, r12=r13=r23=1
#########################################################################
generateModel2_perturb <- function(n = 100, pvec = c(100,100,100), snr = 1, orthogonalV = T){
  p_total <- sum(pvec)
  pcum <- cumsum(pvec)
  D <- 3
  r0 <- 10 # shared rank 10
  rmat <- matrix(1,3,3) # each partially-shared structure has rank 1, no individual
  diag(rmat) <- 0
  r_total <- r0 + sum(rmat[upper.tri(rmat,diag=T)])
  
  # generate noise part of X
  X <- matrix(0, n, p_total)
  
  # generate all the scores
  U <- matrix(runif(n*r_total), n, r_total)
  # Column-center U so that the overall mean is zero
  U <- scale(U, scale = F)
  
  # orthogonalize the scores
  U <- qr.Q(qr(U))
  
  # Generate all the loadings
  Vtmp <- list()
  # generate dataset-specific Vs
  for (d in 1:D){
    Vtmp[[d]] <- matrix(runif(pvec[d]*(r0+sum(rmat[d,]))), pvec[d],(r0+sum(rmat[d,])))
    # orthogonalize the loadings
    Vtmp[[d]] <- qr.Q(qr(Vtmp[[d]]))
  }
  # Put everything together into one V
  V <- matrix(0, nrow = p_total, ncol = r_total)
  for (d in 1:D){
    # Row index
    if (d == 1){
      index <- c(1:pvec[1])
    }else{
      index <- c((pcum[d-1]+1):pcum[d])
    }
    
    # Column index, order is r0, r12, r13, r23, r1, r2, r3
    if (d==1){
      col_indexd = c(1:(r0+rmat[1,2]+rmat[1,3]))
    }else if (d==2){
      col_indexd = c(1:(r0+rmat[1,2]),(r0+rmat[1,2]+rmat[1,3]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]))
    }else{
      col_indexd = c(1:r0,(r0+rmat[1,2]+1):(r0+rmat[1,2]+rmat[1,3]+rmat[2,3]))
    }
    V[index,col_indexd]=Vtmp[[d]]
  }
  # Orthogonallize V
  if (orthogonalV == T){
    V <- qr.Q(qr(V))
  }
  # Add the scale between joint and semi-joint by drawing singular values
  sjoint <- runif(r0, min = 1, max =  1.5)
  spartial <- runif(sum(rmat[upper.tri(rmat,diag=T)]), min = 0.9, max = 1.2)
  V <- V%*%diag(c(sjoint, spartial))
  
  # Form X
  X <- tcrossprod(U, V)
  if (snr != 0){
    # Datasets have different scales or different dimensions
    sigmavec = rep(0, D)
    for (d in 1:D){
      # Row index for dth dataset
      if (d == 1){
        index <- c(1:pvec[1])
      }else{
        index <- c((pcum[d-1]+1):pcum[d])
      }
      sigmavec[d] <- sqrt(sum(X[,index]^2) / (n * pvec[d] * snr))
      X[,index] <- X[,index] + matrix(rnorm(n * pvec[d], sd = sigmavec[d]), n, pvec[d])
    }
    # Combine joint loadings and individual loadings into one matrix V
    return(list(X = X, pvec = pvec, r0 = r0, rmat = rmat, sigmavec = sigmavec, U = U, V = V, snr = snr))
  }else{
    sigma <- 0
    # Combine joint loadings and individual loadings into one matrix V
    return(list(X = X, pvec = pvec, r0 = r0, rmat=rmat, sigma = sigma, U = U, V = V, snr = snr))
  }
}