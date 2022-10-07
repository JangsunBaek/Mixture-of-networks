logL.mnpcl <- function(lambda, n_level, obs_freq, pivec, theta, ...) {
#
# composite penalized log-likelihood
#
# obs_freq[max(n_level),max(n_level),p,p]: observed frequency in contingency tables
#  
#
  eps <- 1.e-5
  p <- length(n_level)
  g <- length(pivec)

  cloglike <- 0
  for (j in 1:(p-1)) {
    for (k in (j+1):p) {
      for (l in 1 : n_level[j]) {
        for (m in 1: n_level[k]) {

#cat("j=",j, "k=",k,"l=",l,"m=",m,"\n")                    
          
          
          cellprob <- sum(exp(log(pivec)+log(theta[l,m,j,k,])))
          diff <- log(theta[l,m,j,k,])-log(eps)
#cat("diff=",diff,"\n") 
          sumlogtheta <- sum(diff[diff>0])
          log_cellprob <- log(cellprob)
#cat("sumlogtheta=",sumlogtheta,"\n") 
#cat("log_cellprob=",log(cellprob),"\n")          
          cloglike <- cloglike + obs_freq[l,m,j,k]*log_cellprob - lambda*sumlogtheta
#cat("cloglike=",cloglike,"\n")
#cat("obsfreq=",obs_freq[l,m,j,k],"\n")          
        }
      }
    }
  }
#cat("cloglike=",cloglike,"\n")
  return(cloglike)
}

