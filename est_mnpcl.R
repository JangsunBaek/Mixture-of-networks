est.mnpcl <- function(X, obs_freq, lambda, pivec, n_level, theta, itmax, tol_logL, ...) {
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
init_para$logL <- do.call("loglike.clnmmj", init_para)

#if ((class(init_para$logL) == "try-error") ||
#      (class(init_para$logL) == 'character')) {
#
#  FIT <- paste('in computing the log-likelihood before EM-steps')
#  class(FIT) <- "error"
#  return(FIT)
#}
#cat("init_para$logL", init_para$logL, "\n")
#cat("\n", "theta[,,3,4,1]=","\n", theta[,,3,4,1],"\n","theta[,,3,4,2]=","\n",theta[,,3,4,2],"\n")
for (niter in 1 : itmax) {
  
########################## tau estimation ##################################
#cat("EM iteration =",niter,"\n")
#  

  tpivec <- array(as.numeric(unlist(init_para$pivec)), dim = c(g,1))
  ttheta <- array(as.numeric(unlist(init_para$theta)), dim = c(max(n_level),max(n_level),p,p,g))
  tau_array <- array(rep(NA,prod(n*g*p^2)),c(n,g,p,p))
  for (j in 1:(p-1)) {
    for (k in (j+1):p) {
      x <- X[,c(j,k)]
      tau_array[1:n,1:g,j,k] <- tau.clnmmj(x, j, k, g, tpivec, ttheta) 
      if (any(is.nan(tau_array[1:n,1:g,j,k]))) {
        FIT <- c(init_para$pivec, init_para$theta)
        return(FIT)
      }
    }
  }
  
#  cat("\n","mstep_tau =", t(tau_array[1:10,1:g,1,2]),"\n" )
  
##############################################################################  

  FIT <- do.call('mstep.clnmmj', c(list(X=X, tau_array=tau_array), init_para))


    
#  cat("FIT_pivec =", "\n", FIT$pivec,"\n" )
#  cat("\n", "FIT_theta_marg[,1,1]=","\n", FIT$theta_marg[,1,1],"\n","FIT_theta_marg[,1,2]=","\n", FIT$theta_marg[,1,2],"\n")
#  cat("\n", "FIT_theta[,,3,4,1]=","\n", FIT$theta[,,3,4,1],"\n","FIT_theta[,,3,4,2]=","\n", FIT$theta[,,3,4,2],"\n")

#  if (class(FIT) == 'error') {
#
#    FIT <- paste('in ', niter,
#                   'iteration of the M-step:',
#                   FIT)
#    class(FIT) <- "error"
#    return(FIT)
#  }

  FIT$logL <- do.call('loglike.clnmmj', c(list(obs_freq=obs_freq,lambda=lambda,n_level=n_level), FIT))
#  FIT$logL <- try(do.call('ploglike.clnmm', FIT))

#  if ((class(FIT$logL) == "try-error") ||
#        (class(FIT$logL) == 'character')) {
#
#    FIT <- paste('in computing the log-likelihood after the ', niter,
#                   'th the M-step', FIT$logL, sep = '')
#    class(FIT) <- "error"
#    return(FIT)
#  }
#
  if ((FIT$logL == -Inf) | is.na(FIT$logL)){

    FIT <- paste('Log likelihood computed after the ', niter,
      'th iteration of the M-step is not finite. Try smaller lambda.', sep='')
    class(FIT) <- "error"
    return(FIT)
  }

#cat("niter =", niter, "FIT$logL =", FIT$logL, "\n" )
  
###################  Aitken acceleration criterion ####################
  logL_t[niter] <- FIT$logL
  if (niter >=3) {
    a_t <- (logL_t[niter]-logL_t[niter-1])/(logL_t[niter-1]-logL_t[niter-2])
#cat("niter=",niter,"","logL_t[niter-1]=", logL_t[niter-1],"logL_t[niter-2]=",logL_t[niter-2], "\n")
#    limlogL_t[niter] <- logL_t[niter-1]+(logL_t[niter]-logL_t[niter-1])/(1-a_t)
#    accel <- (limlogL_t[niter] - limlogL_t[niter-1])/abs(limlogL_t[niter-1])
    limlogL <- logL_t[niter-1]+(logL_t[niter]-logL_t[niter-1])/(1-a_t)
    accel <- (limlogL - logL_t[niter])
#cat("accel=",accel,"\n")    
    if (accel > 0 & accel < tol_logL) {
      break
    }
  }
  init_para <- c(list(obs_freq=obs_freq, lambda=lambda, n_level=n_level), FIT)
}

class(FIT) <- "clnmm"
return(FIT)
}
