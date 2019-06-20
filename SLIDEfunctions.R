# Functions for group lasso penalty for determination of rank paths, and exact algorithm for fitting particular structure
# Last updated: February 21st, 2017

# Fit SLIDE model with structure specified in pattern using starting U and V
est_givenranks_v3<- function(X, pvec, pattern, U, V, eps = 1e-6, k_max = 1000){
  d <- length(pvec)
  pcum <- cumsum(pvec)
  r <- nrow(pattern)
  error <- 1000
  k <- 1

  while((k < k_max)&(error[k] > eps)){
    k <- k + 1
    UVold <- tcrossprod(U,V)
    # Refit U
    # Perform SVD on XV
    XV_svd <- svd(X%*%V)
    U <- tcrossprod(XV_svd$u, XV_svd$v)

    # Refit V
    for (i in 1:d){
      # Index the measurements corresponding to ith dataset
      if (i == 1){
        index <- c(1:pvec[1])
      }else{
        index <- c((pcum[i-1]+1):pcum[i])
      }
      # Identify columns in V that are present for the dataset i
      nonzero <- c(1:r)[pattern[,i] == 1]
      if (length(nonzero) > 0){
        V[index, nonzero] <- crossprod(X[,index], U[,nonzero])
      }
    }
    # Calculate the difference due to refitting
    error[k] <- sum((tcrossprod(U,V)-UVold)^2)
  }
  if (k == k_max){
    warning(paste("Algorithm didn't converge in ", k_max, " iterations!", sep=""))
  }
  return(list(U=U, V=V, error=error))
}

# Fit SLIDE model with structure specified in pattern
est_givenranks_v4<- function(X, pvec, pattern, eps = 1e-6, k_max = 1000, U = NULL){
  d <- length(pvec)
  pcum <- cumsum(pvec)
  r <- nrow(pattern)
  error <- 1000
  k <- 1

  # Initialize U
  if (is.null(U)){
    outs <- svd(X, nu = r)
    U <- outs$u
  }else if ((ncol(U)!=r)|(nrow(U)!=nrow(X))){
    stop("Supplied dimensions of Uinit are wrong")
  }

  # Initialize V
  V <- matrix(0, sum(pvec), r)
  for (i in 1:d){
    # Index the measurements corresponding to ith dataset
    if (i == 1){
      index <- c(1:pvec[1])
    }else{
      index <- c((pcum[i-1]+1):pcum[i])
    }
    # Identify columns in V that are present for the dataset i
    nonzero <- c(1:r)[pattern[,i] == 1]
    if (length(nonzero) >0 ){
      V[index, nonzero] <- crossprod(X[,index], U[,nonzero])
    }
  }

  out <- est_givenranks_v3(X, pvec, pattern, U, V, eps = eps, k_max = k_max)

  return(out)
}

# Refit hatV based on hatU and X
# X - standardized dataset
# pvec - values p_1,....,p_d corresponding to the number of measurements within each data type
# hatU - estimated score matrix
# hatV - estimated loadings matrix
VRefit <- function(X, pvec, hatU, hatV){
  d <- length(pvec)
  pcum <- cumsum(pvec)
  V <- hatV
  for (i in 1:d){
    # Index the measurements corresponding to ith dataset
    if (i == 1){
      index <- c(1:pvec[1])
    }else{
      index <- c((pcum[i-1]+1):pcum[i])
    }
    # Identify columns in V that are present for the dataset i
    if (ncol(hatV) > 1){
      nonzero <- c(1:ncol(hatV))[colSums(hatV[index,]^2)!=0]
      if (length(nonzero) > 0){
        # Refit onlu nonzero columns
        V[index, nonzero] <- crossprod(X[,index], hatU[,nonzero])
      }
    }else if(sum(abs(hatV[index,]))!=0)
      V[index,] <- crossprod(X[,index], hatU)
  }
  return(V)
}

# Block-coordinate update of V for group-lasso rank problem
# X - n x p concatenated data matrix, each dataset has been already standardized
# U - n x r current scores matrix
# lambda - nonnegative tuning parameter
# pvec - values p_1,....,p_d corresponding to the number of measurements within each data type
updateVoptim1 <- function(XU, lambda, pvec){
  d <- length(pvec)
  #r <- ncol(U)
  #p <- sum(pvec)
  #V <- matrix(0, p, r)
  pcum <- cumsum(pvec)
  #XU <- crossprod(X, U)
  V <- matrix(0, nrow(XU), ncol(XU))
  for (i in 1:d){
    # Index the measurements corresponding to ith dataset
    if (i == 1){
      index <- c(1:pvec[1])
    }else{
      index <- c((pcum[i-1]+1):pcum[i])
    }
    # Calculate norm of each column in XU
    norms <- sqrt(colSums(XU[index, ]^2))
    # Perform soft-thresholding for each column - only those that have norms > lambda
    indexr = (norms - 1e-13) > lambda
    if (sum(indexr) > 1){
      V[index, indexr] =  XU[index, indexr] %*% diag(1 - lambda/norms[indexr])
    }else if (sum(indexr) == 1){
      V[index, indexr] =  XU[index, indexr] * (1 - lambda/norms[indexr])
    }
  }
  return(V)
}

# Value of the objective function for the current iterate
# Assume X has already been standardized
evaluatef_optim1 <- function(XU, V, lambda, pvec){
  #f <- sum((X-tcrossprod(U,V))^2)/2
  f <- sum(V^2)/2 - sum(diag(crossprod(V, XU)))
  d <- length(pvec)
  pcum <- cumsum(pvec)
  for (i in 1:d){
    # Index the measurements corresponding to ith dataset
    if (i == 1){
      index <- c(1:pvec[1])
    }else{
      index <- c((pcum[i-1]+1):pcum[i])
    }
    f <- f + lambda * sum(sqrt(colSums(V[index,]^2)))
  }
  return(f)
}

# Iterative algorithm for group-lasso rank problem
# X - n x p concatenated data matrix, standardized
# lambda - nonnegative tuning parameter
# pvec - values p_1,....,p_d corresponding to the number of measurements within each data type
# k_max - maximal number of iterations allowed
# eps - convergence tolerance as measured by the difference in objective function values
# reduced - if true, the rank of U is allowed to decrease with iterations (by taking only nonzero columns of V)
# rank_total - starting rank
solve_optim1 <- function (X, lambda, pvec, k_max = 1000, eps = 1e-6, reduced = F, rank_total = NULL, Uinit = NULL){
  n <- nrow(X)
  p <- ncol(X)
  if (is.null(Uinit)){
    # Perform svd of X to initialize U and V
    X_svd <- svd(X)
    if (is.null(rank_total)){
      rank_total <- sum(X_svd$d > eps)
    }else{
      if (rank_total > length(X_svd$d)){
        stop("Specified rank is larger than the rank of X")
      }
    }
    U <- X_svd$u[ , 1:rank_total]
    V <- X_svd$v[ , 1:rank_total] %*% diag(X_svd$d[1:rank_total])
  }else{
    U = Uinit
    V = crossprod(X, U)
    rank_total = ncol(U)
  }
  # Current function value
  f <- evaluatef_optim1(crossprod(X,U), V, lambda, pvec)
  k <- 1
  error <- 100
  while ((k < k_max)&(error[k] > eps)){
    k <- k + 1
    # Update V
    V <- updateVoptim1(crossprod(X,U), lambda, pvec)
    # Update U
    if (reduced){
      # Eliminate the columns of V that are completely zero
      index_nonzero <- (1:ncol(V))[colSums(abs(V))!=0]
      if (length(index_nonzero)>0){
        V <- V[,index_nonzero, drop = FALSE]
      }else{
        # V is completely zero
        return(list(U=U,V=as.matrix(V),k=k,error=error,f=0))
      }
      # Perform SVD on XV
      XV_svd <- svd(X%*%V)
      # Update U
      U <- tcrossprod(XV_svd$u, XV_svd$v)
    }else{
      # No reduction, full update of U
      # Find out nonzero columns of V
      nonzero = c(1:ncol(V))[colSums(abs(V)) > 0]
      if (length(nonzero) == 0){
        # V is exactly zero
        f[k] = 0
        error[k] <- f[k-1] - f[k]
        return(list(U = U, V = V, Vrefit = V, k = k, error = error, f = f))
      }else{
        # Do full SVD again as QR doesn't work as expected
        XV_svd <- svd(X%*%V)
        U <- tcrossprod(XV_svd$u, XV_svd$v)
      }
    }

    # Current function value
    f[k] <- evaluatef_optim1(crossprod(X,U), V, lambda, pvec)

    # Current difference in function values
    error[k] <- f[k-1] - f[k]
  }
  if (k == k_max){
    warning(paste("Algorithm didn't converge in ", k_max, " iterations!", sep=""))
  }
  Vnew <- VRefit(X, pvec, U, V)
  return(list(U=U, V=as.matrix(V), Vrefit = as.matrix(Vnew), k=k, error=error, f=f))
}

# Iterative algorithm for group-lasso rank problem for the range of lambda values
# X - n x p concatenated data matrix, standardized
# lambdavec - vector of nonnegative tuning parameters
# pvec - values p_1,....,p_d corresponding to the number of measurements within each data type
# k_max - maximal number of iterations allowed
# eps - convergence tolerance as measured by the difference in objective function values
# reduced - if true, the rank of U is allowed to decrease with iterations (by taking only nonzero columns of V)
# rank_total - starting rank
solve_optim1_seq <- function(X, lambda_seq = NULL, pvec, k_max = 1000, eps = 1e-6, reduced = F, rank_total = NULL, n_lambda = 50, lambda_max = 1, lambda_min = 0.01){
  # Check for user supplied lambda_sequence
  if (!is.null(lambda_seq)){
    # Truncate the sequnce by lambda_max
    lambda_seq <- lambda_seq[lambda_seq <= lambda_max]
    # Order the sequence from largest to smallest
    lambda_seq <- sort(lambda_seq, decreasing = T)
    # Update n_lambda value
    n_lambda <- length(lambda_seq)
  }else{
    # generate the sequence of tuning parameters
    lambda_seq <- exp(seq(log(lambda_max), log(lambda_min), length.out=n_lambda))
  }

  # Generate the same Uinit for each, so no need to repeat SVD
  # Perform svd of X to initialize U and V
  X_svd <- svd(X)
  if (is.null(rank_total)){
    rank_total <- sum(X_svd$d > eps)
  }else{
    rank_total = min(rank_total, sum(X_svd$d > eps))
  }
  Uinit <- X_svd$u[ , 1:rank_total]

  # Solve for each lambda
  param = list()
  for (l in 1:n_lambda){
    # Solve group lasso problem
    out <- solve_optim1(X, lambda = lambda_seq[l], pvec = pvec, k_max = k_max, eps = eps, reduced = reduced, rank_total = rank_total, Uinit = Uinit)
    # Update starting points for next time - may not want to do this
    param[[l]] <- list(U = out$U, V = out$V, Vrefit = out$Vrefit, fmin = min(out$f))
  }
  return(list(param = param, lambda = lambda_seq))
}



# Bi-cross-validation approach adapted from Owen and Perry
# X - n x p concatenated data matrix, standardized
# pvec - values p_1,....,p_d corresponding to the number of measurements within each data type
# k_max - maximal number of iterations allowed
# eps - convergence tolerance as measured by the difference in objective function values
# reduced - if true, the rank of U is allowed to decrease with iterations (by taking only nonzero columns of V)
# rank_total - starting rank
# structure_list - list of rank combinations to consider
# n_structure - default number of structures to consider
#
# This one first fits the whole range of structures, and then determines which one is the best
# Centering - need to take into account the error on a dataset in a different way, currently don't consider centering or scaling
bcv_optim1_structure_centering <- function(X, pvec, structure_list, n_fold = 3, p_fold=3, k_max = 1000, eps = 1e-6, center = F){
  d <- length(pvec)
  n <- nrow(X)
  p <- ncol(X)
  pcum <- cumsum(pvec)

  # Get length of structure list
  n_structure = length(structure_list)

  # Get the folds
  # For samples
  fold_id_n <- sample(rep(seq_len(n_fold), length.out = n))
  # For measurements
  fold_id_p <- c()
  for (i in 1:d){
    fold_id_p <- c(fold_id_p, sample(rep(seq_len(p_fold), length.out = pvec[i])))
  }

  # Save prediction error
  error <- matrix(0, n_fold*p_fold, n_structure)
  iter <- 0
  for (k in 1:n_fold){
    for (j in 1:p_fold){
      iter <- iter+1
      # Update the pvec
      pvec_foldj = pvec
      for (i in 1:d){
        # Index the measurements corresponding to ith dataset
        if (i == 1){
          index <- c(1:pvec[1])
        }else{
          index <- c((pcum[i-1]+1):pcum[i])
        }
        pvec_foldj[i] = sum(fold_id_p[index] != j)
      }
      pcum_fold = cumsum(pvec_foldj)
      # Standardize each dataset
      out_s <- standardizeX(X = X[fold_id_n != k, fold_id_p != j], pvec = pvec_foldj, center = center)
      Xfold <- out_s$X

      # Consider all structure from the list, fit the model, evaluate BCV error
      for (l in 1:n_structure){
        # Find U and V based on the given structure
        out <- est_givenranks_v4(X = Xfold, pvec = pvec_foldj, pattern = structure_list[[l]], k_max = k_max, eps = eps)
        Vadj <- out$V
        for (i in 1:d){
          if (i == 1){
            index <- c(1:pvec_foldj[1])
          }else{
            index <- c((pcum_fold[i-1]+1):pcum_fold[i])
          }
          Vadj[index,] <- out$V[index,]*sqrt(out_s$norms[i])
        }
        # Perform SVD on the output
        out_svd <- svd(tcrossprod(out$U, Vadj))
        if (max(out_svd$d) < eps){
          error[iter, l] <- d
        }else{
          # because of est_givenranks, only have nonzero columns already
          if (center == F){
            # No centering was done
            tildeV = crossprod(X[fold_id_n != k, fold_id_p == j], out$U)
            tildeU = X[fold_id_n == k, fold_id_p != j] %*% Vadj %*% solve(crossprod(Vadj))
            # Calculate error
            error[iter, l] = sum((X[fold_id_n == k, fold_id_p == j] - tcrossprod(tildeU, tildeV))^2)
          }else{
            n_newf = sum(fold_id_n != k) # number of samples in the submatrix for model fitting
            # Centering was done, stored in out_s$Xmean
            tildeV = cbind(colMeans(X[fold_id_n != k, fold_id_p == j])*sqrt(n_newf), crossprod(X[fold_id_n != k, fold_id_p == j], out$U))
            Vnew = cbind(sqrt(n_newf)*out_s$Xmean, Vadj)
            # Adjust for potential singularity here
            Vsvd = svd(Vnew)
            Vproj = Vsvd$u[, Vsvd$d > eps] %*% diag(1/Vsvd$d[Vsvd$d > eps]) %*% t(Vsvd$v[, Vsvd$d > eps])
            # This may be singular, need to take that into account by taking ginv
            tildeU = X[fold_id_n == k, fold_id_p != j] %*% Vproj
            pcum_notfold <- cumsum(pvec) - pcum_fold
            for (i in 1:d){
              if (i == 1){
                index <- c(1:(pvec[1]-pvec_foldj[1]))
              }else{
                index <- c((pcum_notfold[i-1]+1):pcum_notfold[i])
              }
              error[iter, l] = error[iter,l] + sum((X[fold_id_n == k, which(fold_id_p == j)[index]] - tcrossprod(tildeU, tildeV[index,]))^2)/sum(scale(X[fold_id_n == k, which(fold_id_p == j)[index]],scale=F)^2)
            }
            #error[iter, l] = sum((X[fold_id_n == k, fold_id_p == j] - tcrossprod(tildeU, tildeV))^2)
          }
        }
      } # end for lamda_seq
    } # end for p folds
  } # end for n folds
  # Calculate average prediction error for each tuning parameter
  error_mean <- colSums(error)
  # Find largest lambda corresponding to minimal average error, here assume that lambda_seq is sorted from largest to smallest
  id_min <- min(c(1:n_structure)[error_mean == min(error_mean)])
  structure_min <- structure_list[[id_min]]
  return(list(structure_list = structure_list, structure_min = structure_min, error_mean = error_mean, error = error))
}



# Iterative algorithm for group-lasso rank problem for the range of lambda values
# X - n x p concatenated data matrix, standardized
# lambdavec - vector of nonnegative tuning parameters
# pvec - values p_1,....,p_d corresponding to the number of measurements within each data type
# k_max - maximal number of iterations allowed
# eps - convergence tolerance as measured by the difference in objective function values
# reduced - if true, the rank of U is allowed to decrease with iterations (by taking only nonzero columns of V)
# rank_total - starting rank
# n_start - additional starting points
solve_optim1_seq_restarts <- function (X, lambda_seq = NULL, pvec, k_max = 1000, eps = 1e-6, reduced = F,  rank_total = NULL, n_lambda = 50, n_start = 10, lambda_max = 1, lambda_min = 0.01){
  require(pracma)
  # Check for user supplied lambda_sequence
  if (!is.null(lambda_seq)){
    # Truncate the sequnce by lambda_max
    lambda_seq <- lambda_seq[lambda_seq <= lambda_max]
    # Order the sequence from largest to smallest
    lambda_seq <- sort(lambda_seq, decreasing = T)
    # Update n_lambda value
    n_lambda <- length(lambda_seq)
  }else{
    # generate the sequence of tuning parameters
    lambda_seq <- exp(seq(log(lambda_max), log(lambda_min), length.out=n_lambda))
  }

  # Generate starting points
  ################
  Ustart <- list()
  # Take svd of X
  Xsvd <- svd(X)
  rankX <- sum(Xsvd$d > eps)
  if (is.null(rank_total)){
    rank_total = rankX
  }else{
    rank_total = min(rank_total, rankX)
  }
  # Initialize U with left singular vectors of concatenated X
  Ustart[[1]] = Xsvd$u[,1:rank_total]
  # Initialize U with left singular vectors of each individual X_i
  d <- length(pvec)
  pcum <- cumsum(pvec)
  n <- nrow(X)
  for (i in 1:d){
    # Index the measurements corresponding to ith dataset
    if (i == 1){
      index <- c(1:pvec[1])
    }else{
      index <- c((pcum[i-1]+1):pcum[i])
    }
    # Perform svd on X_i
    Xi_svd <- svd(X[, index])
    dim_i <- length(Xi_svd$d)
    if (dim_i >= rank_total){
      Ustart[[i+1]] <- Xi_svd$u[,1:rank_total]
    }else{
      # Change this part
      U0 <- Xi_svd$u
      U1 <- svd(diag(n)-tcrossprod(U0))$u[, 1:(rank_total - dim_i)]
      Ustart[[i+1]] <- cbind(U0, U1)
    }
  }
  # Initialize U with left random orthogonal matrices, use n_start points
  for (i in 1:n_start){
    Ustart[[1+d+i]] <- randortho(n)[, 1:rank_total]
  }

  # Solve for each lambda
  param = list()
  startp = rep(0, n_lambda)
  startf = rep(0, n_lambda)
  for (l in 1:n_lambda){
    f = d/2 # This is the value of objective function at V=0
    # Try each starting point
    for (s in 1:(d + 1 + n_start)){
      # Solve group lasso problem
      out <- solve_optim1(X, lambda = lambda_seq[l], pvec = pvec, k_max = k_max, eps = eps, reduced = reduced, rank_total = rank_total, Uinit = Ustart[[s]])
      # Check whether obtained function value is better than the current stored
      if (min(out$f) <= f){
        f <- min(out$f)
        param[[l]] <- out
        # Keep track of which starting point resulted in obtained minimum
        startp[l] <- s
        startf[l] <- f
      }
    }
  }
  return(list(param = param, startp = startp, startf = startf, lambda = lambda_seq))
}



###################
# Calculate structure for each value of tuning parameter in out_solve, any number of datasets
# Return: r times d structure matrix
# Need to take into account both unique and non-unique structures, get rid of all zero structures
# ratio_max - do not consider structures that assign individual rank larger than 2/3 of the total possible rank
get_structure_v2 <- function(out_solve, pvec, eps = 1e-6, ratio_max = 2/3){
  n_lambda <- length(out_solve$lambda)
  n <- nrow(out_solve$param[[1]]$U)
  # Keeps track of which structures were considered, i.e. if s_used = 1, then the structure is used, if s_used = 0, then the structure is not used
  s_used <- rep(0, n_lambda)
  d <- length(pvec)
  pcum <- cumsum(pvec)
  Slist <- list()
  iter <- 0
  # Sold keeps track of the previous structure so that the comparison can be made
  Sold <- NULL
  rall <- ncol(out_solve$param[[1]]$V)
  for (t in 1:n_lambda){
    # Calculate the number of nonzero columns in V
    r <- sum(colSums(abs(out_solve$param[[t]]$V))>0)
    if (r > 0){
      if (r > ratio_max * min(n, sum(pvec))){
        # Total rank is larger than 2/3 of possible rank
        return(list(Slist = Slist, s_used = s_used))
      }
      # Figure out current structure
      S <- matrix(0, rall, d)
      for (j in 1:d){
        if (j == 1){
          index <- c(1:pvec[1])
        }else{
          index <- c((pcum[j-1]+1):pcum[j])
        }
      S[,j] = (colSums(abs(out_solve$param[[t]]$V[index,])) > 0)
       # Check whether the maximal rank is exceeded
        if (sum(S[,j]) > ratio_max * min(n, pvec[j])){
          return(list(Slist = Slist, s_used = s_used))
        }
      }
      # Reduce to only have non-zero columns
      S = S[rowSums(S) > 0, , drop = FALSE]
      # Compare with what is already there to see whether it is new or not
      if (is.null(Sold)){
        # Have not seen this structure before
        iter = iter + 1
        s_used[t] = 1
        Slist[[iter]] <- S
        Sold <- S
      }else if (nrow(Sold)!=nrow(S)){
          # Have not seen this structure before
          iter = iter + 1
          s_used[t] = 1
          Slist[[iter]] <- S
          Sold <- S
      }else if (max(abs(Sold-S))!=0){
        # Have not seen this structure before
        iter = iter + 1
        s_used[t] = 1
        Slist[[iter]] <- S
        Sold <- S
      }
    }
  }#end for lambda_seq
  return(list(Slist = Slist, s_used = s_used))
}


######################################################

# Function for simulations to compare BCV versus best, use alternative version of BCV (full residual)
# X - original dataset
# trueU, trueV- true structure
# pvec - values p_1,....,p_d corresponding to the number of measurements within each data type
# nl - length of grid to consider for BCV
compareBCV_vs_Best_structure_refit_center <- function(X, pvec, trueU, trueV, nl = 30, n_fold = 3, p_fold = 3, k_max = 1000, eps = 1e-6, reduced = F, rank_total = NULL, center = T){
  # Scale the original dataset
  out_s <- standardizeX(X, pvec, center = center)
  X <- out_s$X
  svec <- out_s$svec
  pcum <- cumsum(pvec)
  d <- length(pvec)

  # Solve for each lambda
  out <- solve_optim1_seq(X = X, pvec = pvec, n_lambda = nl, reduced = reduced, lambda_max = max(svec), lambda_min = 0.05)
  lambda_all <- out$lambda

  ## Calculate weighted residual error for each lambda
  #errors <- getAllErrors(out, trueU, trueV, X, pvec, out_s$norms)
  #best_id <- which.min(errors)
  #lambda_best <- lambda_all[best_id]

  # Calculate ranks for each lambda
  if (length(pvec) == 2){
    ranks <- getAllranks_d2(out, pvec)
  }else{
    ranks <- getAllranks_d3(out, pvec)
  }

  # Get all the structures
  out_struct <- get_structure_v2(out, pvec)

  if (sum(out_struct$s_used)==0){
    stop("Something went terribly wrong")
  }

  errors = rep(0, length(out_struct$Slist))
  for (l in 1:length(out_struct$Slist)){
    distUV = 0
    outV <- est_givenranks_v4(X = X, pvec = pvec, pattern = out_struct$Slist[[l]], k_max = k_max, eps = eps)
    Vadj <- outV$V
    for (i in 1:d){
      if (i == 1){
        index <- c(1:pvec[1])
      }else{
        index <- c((pcum[i-1]+1):pcum[i])
      }
      Vadj[index,] <- outV$V[index,]*sqrt(out_s$norms[i])
      # Calculate Frobenius norm squared of the signal part
      signal = sum(tcrossprod(trueU,trueV[index,])^2)
      # If nonzero, calculate error relative to signal otherwise just directly
      if (signal > 0){
        distUV <- distUV + sum((tcrossprod(trueU,trueV[index,])-tcrossprod(outV$U, Vadj[index,]))^2)/sum(tcrossprod(trueU,trueV[index,])^2)
      }else{
        distUV <- distUV + sum((tcrossprod(outV$U, Vadj[index,]))^2)
      }
    }
    errors[l] = distUV
  }

  best_id_used <- which(errors == min(errors))
  # out of all the lambdas, which one does it correspond to
  best_id <- which(out_struct$s_used == 1)[best_id_used]

  # Find the best tuning parameter using bcv from the structures above
  outbcv <- bcv_optim1_structure_centering(X, pvec = data$pvec, structure_list = out_struct$Slist, n_fold = n_fold, p_fold=p_fold, k_max = k_max, eps = eps, center = center)
  # out of the structures considered, which one was picked
  bcv_id_used <- which(outbcv$error_mean == min(outbcv$error_mean))
  # out of all the lambdas, which one does it correspond to
  bcv_id <- which(out_struct$s_used == 1)[bcv_id_used]

  return(list(best = list(lambda = lambda_all[best_id], error = errors[best_id_used], ranks = ranks[[best_id]]), bcv = list(lambda = lambda_all[bcv_id], error = errors[bcv_id_used], ranks = ranks[[bcv_id]])))
}




