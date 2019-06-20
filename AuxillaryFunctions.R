# Function to scale the data
# X - original dataset, unstandardized
# pvec - values p_1,....,p_d corresponding to the number of measurements within each data type
# center  - T/F, whether centering should be performed before scaling
standardizeX <- function(X, pvec, center = F){
  d <- length(pvec)
  norms <- rep(0,d)
  svec <- rep(0,d)
  pcum <- cumsum(pvec)
  # Center each column
  if (center){
    Xmean <- colMeans(X)
    X <- X - matrix(Xmean, nrow(X), ncol(X), byrow = T) 
  }else{
    Xmean <- rep(0, ncol(X))
  }
  for (i in 1:d){
    if (i == 1){
      index <- c(1:pvec[1])
    }else{
      index <- c((pcum[i-1]+1):pcum[i])
    }
    norms[i] <- sum(X[,index]^2)
    # Scale the dataset
    X[, index] <- X[,index]/sqrt(norms[i])
    # Calculate largest singular value
    svec[i] <- max(svd(X[, index])$d)
  }
  return(list(X = X, svec = svec, norms = norms, Xmean = Xmean))
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

###################
# Calculate residul errors for each value of tuning parameter in out_solve (using already refitted V)
getAllErrors <- function(out_solve, U, V, X, pvec, norms){
  n_lambda <- length(out_solve$lambda)
  distUV <- rep(0, n_lambda)
  d <- length(pvec)
  pcum <- cumsum(pvec)
  for (l in 1:n_lambda){
    # Adjust estimated V to correct for scaling
    Vadj <- out_solve$param[[l]]$Vrefit
    for (i in 1:d){
      if (i == 1){
        index <- c(1:pvec[1])
      }else{
        index <- c((pcum[i-1]+1):pcum[i])
      }
      Vadj[index,] <- Vadj[index,]*sqrt(norms[i])
      # Calculate all corresponding distances
      distUV[l] <- distUV[l] + sum((tcrossprod(U,V[index,])-tcrossprod(out_solve$param[[l]]$U,Vadj[index,]))^2)/sum(tcrossprod(U,V[index,])^2)
    } # end for d
  }# end for lamda_seq
  return(distUV)
}

# Calculate all ranks for each value of tuning parameter in out_solve, 2 datasets
getAllranks_d2 <- function(out_solve, pvec, eps = 1e-6){
  n_lambda <- length(out_solve$lambda)
  d <- length(pvec)
  pcum <- cumsum(pvec)
  ranks_list <- list()
  ranks_m <- matrix(0,d,d)
  for (t in 1:n_lambda){
    # Calculate all corresponding ranks
    for (j in 1:d){
      if (j == 1){
        index <- c(1:pvec[1])
      }else{
        index <- c((pcum[j-1]+1):pcum[j])
      }
      if (ncol(out_solve$param[[t]]$V)>1){
        ranks_m[j,j] = sum(colSums(out_solve$param[[t]]$V[index,]^2)>eps)
      }else{
        ranks_m[j,j]= (sum(out_solve$param[[t]]$V[index,]^2)>eps)
      }
      if (j < d){
        for (l in (j+1):d){
          index2 <- c((pcum[l-1]+1):pcum[l])
          if (ncol(out_solve$param[[t]]$V)>1){
            ranks_m[j,l] = sum((colSums(out_solve$param[[t]]$V[index,]^2)>eps)&(colSums(out_solve$param[[t]]$V[index2,]^2)>eps))
          }else{
            ranks_m[j,l] =(sum(out_solve$param[[t]]$V[index,]^2)>eps)&(sum(out_solve$param[[t]]$V[index2,]^2)>eps)
          }
        }
      }
    }
    ranks_list[[t]] = list(r0 = ranks_m[1,2], r1 = ranks_m[1,1]-ranks_m[1,2], r2 = ranks_m[2,2]-ranks_m[1,2])
  } # end for lamda_seq
  return(ranks_list)
}

# Calculate all ranks for each value of tuning parameter in out_solve, 3 datasets
getAllranks_d3 <- function(out_solve, pvec, eps = 1e-6){
  n_lambda <- length(out_solve$lambda)
  d <- length(pvec)
  pcum <- cumsum(pvec)
  ranks_list <- list()
  ranks_m <- matrix(0,d,d)
  for (t in 1:n_lambda){
    # Calculate all corresponding ranks
    for (j in 1:d){
      if (j == 1){
        index <- c(1:pvec[1])
      }else{
        index <- c((pcum[j-1]+1):pcum[j])
      }
      # Total rank for thar particular dataset
      if (ncol(out_solve$param[[t]]$V)>1){
        ranks_m[j,j] = sum(colSums(out_solve$param[[t]]$V[index,]^2)>eps)
      }else{
        ranks_m[j,j]= (sum(out_solve$param[[t]]$V[index,]^2)>eps)
      }
      # How much is shared between any two
      if (j < d){
        for (l in (j+1):d){
          index2 <- c((pcum[l-1]+1):pcum[l])
          if (ncol(out_solve$param[[t]]$V)>1){
            ranks_m[j,l] = sum((colSums(out_solve$param[[t]]$V[index,]^2)>eps)&(colSums(out_solve$param[[t]]$V[index2,]^2)>eps))
          }else{
            # only one column in V
            ranks_m[j,l] =as.integer((sum(out_solve$param[[t]]$V[index,]^2)>eps)&(sum(out_solve$param[[t]]$V[index2,]^2)>eps))
          }
        }
      }
      # Joint rank
      if (ncol(out_solve$param[[t]]$V)>1){
        r0 = sum((colSums(out_solve$param[[t]]$V[1:pvec[1],]^2)>eps)&(colSums(out_solve$param[[t]]$V[(pcum[1]+1):pcum[2],]^2)>eps)&(colSums(out_solve$param[[t]]$V[(pcum[2]+1):pcum[3],]^2)>eps))
      }else{
        r0 = as.integer((sum(out_solve$param[[t]]$V[1:pvec[1],]^2)>eps)&(sum(out_solve$param[[t]]$V[(pcum[1]+1):pcum[2],]^2)>eps)&(sum(out_solve$param[[t]]$V[(pcum[2]+1):pcum[3],]^2)>eps))
      }
    }
    ranks_list[[t]] = list(r0 = r0, r12 = ranks_m[1,2]-r0, r13 = ranks_m[1,3]-r0, r23 = ranks_m[2,3]-r0, r1= ranks_m[1,1]-ranks_m[1,3]-ranks_m[1,2]+r0,r2 = ranks_m[2,2]-ranks_m[1,2]-ranks_m[2,3]+r0, r3=ranks_m[3,3]-ranks_m[1,3]-ranks_m[2,3]+r0)
  } # end for lamda_seq
  return(ranks_list)
}

# Calculate structure for each value of tuning parameter in out_solve, any number of datasets
get_structure <- function(out_solve, pvec, eps = 1e-6){
  n_lambda <- length(out_solve$lambda)
  d <- length(pvec)
  pcum <- cumsum(pvec)
  Slist <- list()
  for (t in 1:n_lambda){
    # Calculate all corresponding ranks
    r <- ncol(out_solve$param[[t]]$V)
    S <- matrix(0, d, r)
    for (j in 1:d){
      if (j == 1){
        index <- c(1:pvec[1])
      }else{
        index <- c((pcum[j-1]+1):pcum[j])
      }
      S[j, ] = (colSums(abs(out_solve$param[[t]]$V[index,])) > 0)
    }
    Slist[[t]] <- S
  }#end for lambda_seq
  return(Slist)
}

# Calculate residul errors for each value of tuning parameter in out_solve
getAllErrors_exact <- function(out_solve, U, V, X, pvec, norms){
  n_lambda <- length(out_solve$lambda)
  distUV <- rep(0, n_lambda)
  d <- length(pvec)
  pcum <- cumsum(pvec)
  Vadj <- V
  # Adjust true V for scale to make the fair comparison
  for (i in 1:d){
    if (i == 1){
      index <- c(1:pvec[1])
    }else{
      index <- c((pcum[i-1]+1):pcum[i])
    }
    Vadj[index,]=V[index,]/sqrt(norms[i])
  }
  for (l in 1:n_lambda){
    # Adjust estimated V to correct for scaling
    Vadj <- out_solve$param[[l]]$V
    for (i in 1:d){
      if (i == 1){
        index <- c(1:pvec[1])
      }else{
        index <- c((pcum[i-1]+1):pcum[i])
      }
      Vadj[index,] <- Vadj[index,]*sqrt(norms[i])
      # Calculate all corresponding distances
      distUV[l] <- distUV[l] + sum((tcrossprod(U,V[index,])-tcrossprod(out_solve$param[[l]]$U,Vadj[index,]))^2)/sum(tcrossprod(U,V[index,])^2)
    } # end for d
  }# end for lamda_seq
  return(distUV)
}

## To get value of r0, can use for example
# r0seq = sapply(Slist, function(x) sum(colSums(x) == nrow(x)))
## To get values of r12
# r12seq = sapply(Slist, funciton(x) sum((colSums(x[1:2,]) == 2)&(colSums(x[-c(1:2),]) == 0)))