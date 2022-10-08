mstep.mnpcl <- function (X, lambda, pivec, n_level, theta, tau_array, ...) {
#
# X: nxp
 tol <- 1.e-10
 eps <- 1.e-5  # tiny probability for the empty cell
 itmax <- 5
 g <- length(pivec)
 p <- ncol(X)
 n <- nrow(X)

 ###### pivec estimation #####
  pivec <- rep(0,g)
  for (i in 1:g) {
    sp <- 0
    for (j in 1:(p-1)) {
      for (k in (j+1):p) {
        sp <- sp + sum(tau_array[,i,j,k])
#cat("\n","tau_i=",tau_array[,i,j,k],"\n")
      }
    }
    pivec[i] <- sp*2/(n*p*(p-1))
  }

##### Estimate of theta #####

  for (i in 1 : g) {
    for (j in 1:(p-1)) {
      for (k in (j+1):p) {
        x <- X[,c(j,k)]
        num <- matrix(rep(0, n_level[j]*n_level[k]), c(n_level[j],n_level[k]))
        level_xj <- 1:n_level[j]
        level_xk <- 1:n_level[k]
        thetaseq <- 1:(n_level[j]*n_level[k])
        midx_eps <- matrix(rep(0, n_level[j]*n_level[k]), c(n_level[j],n_level[k]))
        tp_theta <-  theta[level_xj,level_xk,j,k,i] 
        tb_theta <- tp_theta
        
        for ( t in 1:itmax) {
          for (l in level_xj) {
            for (m in level_xk) {
              st <- sum(tau_array[x[,1]==l&x[,2]==m,i,j,k])
              num[l,m] <- st - lambda*max(0,sign(tp_theta[l,m] - eps))
              if (tp_theta[l,m]<=eps | num[l,m]<=0) {
                midx_eps[l,m] <- 1
              }
            }
          }
          idx_eps <- which(midx_eps == 1)
          loc<-arrayInd(idx_eps,dim(midx_eps))
          n_idx <- nrow(loc)
          tp_theta[loc] <- eps
          tp_theta[thetaseq[-idx_eps]] <- num[-idx_eps]/sum(num[-idx_eps])*(1-n_idx*eps)
          if (sum(abs(tb_theta - tp_theta)) < tol) {
            theta[level_xj,level_xk,j,k,i] <- tp_theta
            break
          }
          tb_theta <- tp_theta
        }
        theta[level_xj,level_xk,j,k,i] <- tp_theta
       
      }
    }
  }
 
model <- list(pivec = pivec, theta=theta)
return(model)
}
