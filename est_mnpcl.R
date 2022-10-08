est.mnpcl <- function(X, obs_freq, lambda, pivec, n_level, theta, itmax, tol_logL, ...) {
#
### parameter estimation: EM procedure ###
#
#
p <- ncol(X)
n <- nrow(X)
g <- length(pivec)
logL_t <- rep(NA,itmax)
limlogL_t <- rep(NA,itmax)
limlogL <- 1.e-1000
#limlogL_t[1:2] <- rep(1.e-1000,2)
#
#
init_para <- list(obs_freq=obs_freq,lambda=lambda,pivec=pivec,n_level=n_level,theta=theta)
init_para$logL <- do.call("logL.mnpcl", init_para)
#
#
for (niter in 1 : itmax) {
#  
########################## tau estimation ##################################
#cat("EM iteration =",niter,"\n")
#  

  tpivec <- array(as.numeric(unlist(init_para$pivec)), dim = c(g,1))
  ttheta <- array(as.numeric(unlist(init_para$theta)), dim = c(max(n_level),max(n_level),p,p,g))
  tau_array <- array(rep(NA,prod(n*g*p^2)),c(n,g,p,p))
  for (j in 1:(p-1)) {
    for (k in (j+1):p) {
      x <- X[,c(j,k)]
      tau_array[1:n,1:g,j,k] <- tau.mnpcl(x, j, k, g, tpivec, ttheta) 
      if (any(is.nan(tau_array[1:n,1:g,j,k]))) {
        FIT <- c(init_para$pivec, init_para$theta)
        return(FIT)
      }
    }
  }
  

##############################################################################  

  FIT <- do.call('mstep.mnpcl', c(list(X=X, tau_array=tau_array), init_para))
  FIT$logL <- do.call('logL.mnpcl', c(list(obs_freq=obs_freq,lambda=lambda,n_level=n_level), FIT))
  if ((FIT$logL == -Inf) | is.na(FIT$logL)){
    FIT <- paste('Log likelihood computed after the ', niter,
      'th iteration of the M-step is not finite. Try smaller lambda.', sep='')
    class(FIT) <- "error"
    return(FIT)
  }

###################  Aitken acceleration criterion ####################
  logL_t[niter] <- FIT$logL
  if (niter >=3) {
    a_t <- (logL_t[niter]-logL_t[niter-1])/(logL_t[niter-1]-logL_t[niter-2])
    limlogL <- logL_t[niter-1]+(logL_t[niter]-logL_t[niter-1])/(1-a_t)
    accel <- (limlogL - logL_t[niter])
    if (accel > 0 & accel < tol_logL) {
      break
    }
  }
  init_para <- c(list(obs_freq=obs_freq, lambda=lambda, n_level=n_level), FIT)
}

class(FIT) <- "mnpcl"
return(FIT)
}
