tau.mnpcl <- function(x, j, k, g, pivec, theta, ...) {
# 
# tau_{i(j,k)}^{g}: nxg matrix
# x: x_jk 
# 
#
  n <- nrow(x)
  tau <- matrix(rep(0,n*g),c(n,g))
  log_theta <- rep(0,n)
  #
  for (h in 1:g) {
    for (i in 1:n) {
#cat("j=",j,"\n","k=",k,"\n")      
      
      if (theta[x[i,1],x[i,2],j,k,h]==0) {
        log_theta[i] <- -Inf
      }
      log_theta[i] <- log(theta[x[i,1],x[i,2],j,k,h])
#cat("i=",i,"\n")
#cat("log_theta_cond[i]=",log_theta_cond[i], "\n") 
    }

    tau[, h] <- log_theta
#cat("\n",tau[,h])
    
  }
  pitau <- sweep (tau, 2, log(pivec), '+')
  pitaumax <- apply (pitau, 1, max)
  pitau <- sweep (pitau, 1, pitaumax, '-')
  pitau <- exp(pitau)
  tau <- sweep(pitau, 1, rowSums(pitau), '/')

#cat("\n",t(tau))  
  
  return(tau)
}