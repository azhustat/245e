# implementation of the BNK

ll <- function(Theta) {
  sum.ind <- which(A==1, arr.ind=T)
  logLik <- sum(log(Theta[sum.ind])) - sum(Theta)
  return(logLik)
}

bnk <- function(A, n, K, maxit = 100, thresh = 0.1) {
  # random initialization
  temp <- matrix(runif(n*K), n, K)
  theta.old <- temp/rowSums(temp)
  Theta.old <- theta.old %*% t(theta.old)
  L.old <- ll(Theta.old)
  L.old
  q <- array(0, dim = c(n, n, K))
  
  for (repid in 1:20) {
    # update q
    for (i in 1:n) {
      for (j in 1:n) {
        for (z in 1:K) {
          q[i, j, z] <- theta.old[i, z]*theta.old[j, z]/Theta.old[i, j]
        }
      }
    }
    
    # update thetas
    theta.new <- theta.old
    for (z in 1:K) {
      theta.new[,z] <- rowSums(A*q[,,z])/sqrt(sum(A*q[,,z]))
    }
    
    Theta.new <- theta.new %*% t(theta.new)
    # look at L and see if it converged.
    L.new <- ll(Theta.new)
    L.old
    L.new
    
    converge.check <- abs(L.new - L.old) < thresh
    converge.check
    # take appropriate actions
    if (converge.check){
      cat("converged at L = ",L.new,"\n")
      break
    } else if (repid == maxit) {
      cat("did not converge after ", maxit, " repetitions \n")
    } else {
      L.old <- L.new
      theta.old <- theta.new
      Theta.old <- Theta.new
    }
  }
  
  theta.new <- theta.new/rowSums(theta.new)
  result <- list(theta = theta.new, loglik = L.new, q = q)
  return(result)
}