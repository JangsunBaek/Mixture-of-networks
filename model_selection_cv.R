model_selection_cv <- function(x,n_level,k_cv,min_lambda,max_lambda,inc_lambda,min_G,
                            max_G, initmax, itmax, tol_logL) {
#
# model selection for real data
# x: real data matrix (nxp) 
# k_cv: fold cardinal number for cv  
# g : true number of group of simulation data 
# n_level : vector of number of response level of the variables
# theta_marg1: marginal prob. of x_1: n_level[1]xg
# theta_cond : conditional prob. array:dim=(max_n_level,max_n_level,j,k,g)=c(3,3,4,4,2)  
# min_lambda  
# max_lambda 
# inc_lambda : increment of trial lambda  
# max_G : max trial number of groups  
#  
#  install.packages("klaR")
#  library(klaR)
#
#  install.packages("poLCA") for LCA
#  library("poLCA")  
#    eps <- 1.e-5
    n <- nrow(x)
    p <- ncol(x)
    n_sample <- 10 ### number of simulation sample
    G <- min_G:max_G
    Lambda <- seq(min_lambda,max_lambda, by=inc_lambda)
    n_Lambda <- length(Lambda)
    opt_para <- matrix(rep(0,3*2), nrow=2)   ## optimal g and lambda for cCLC and cv
    rmax=max(n_level)
    obs_freq <- array(rep(0,prod(rmax^2*p^2)),c(rmax,rmax,p,p))
    obs_freqcv <- array(rep(0,prod(rmax^2*p^2)),c(rmax,rmax,p,p))
    cCLC <- matrix(rep(0,(max_G-min_G+1)*n_Lambda), nrow=(max_G-min_G+1))  ## cCLC
    pcllcv <- array(rep(0,(max_G-min_G+1)*n_Lambda*k_cv), dim=c(max_G-min_G+1,n_Lambda,k_cv))  ## predictive composite log-likelihood for corresponding CV data set
    pcll <- matrix(rep(0,(max_G-min_G+1)*n_Lambda), nrow=c(max_G-min_G+1,n_Lambda))
   
    
### predictive composite log-likelihood ####
#
#    x <- x[sample(1:n),]  ## suffle data
    n_cv <- floor(n/k_cv)  ## sample size of cv data
    n_tr <- n - n_cv      ## sample size of training data
    for (s in min_G:max_G) {
      for (t in 1:n_Lambda) {
        # separation of k-fold cv data and training data 
        for (i in 1:k_cv) {
          idx_cv <- ((i-1)*n_cv+1):(i*n_cv)
          x_cv <- x[idx_cv,]
          x_tr <- x[-idx_cv,]      
          obs_freqcv <- array(rep(0,prod(rmax^2*p^2)),c(rmax,rmax,p,p))
          
# Transform the training data and CV data to the observed frequency array

         for (j in 1:(p-1)) {
           for (k in (j+1):p) {
             Xjk_cv <- x_cv[,c(j,k)]
             for (m in 1:n_cv) {
               obs_freqcv[Xjk_cv[m,1],Xjk_cv[m,2],j,k] <- obs_freqcv[Xjk_cv[m,1],Xjk_cv[m,2],j,k] +1
             }
           }
         }
         y_pcll<-mnpcl(x_tr, g=s, n_level, lambda=Lambda[t], initmax, itmax, tol_logL, init_para = NULL)
         pcllcv[s-min_G+1,t,i] <- logLpred.mnpcl(n_level, obs_freqcv, y_pcll$pivec, y_pcll$theta)
        }
        pcll[s-min_G+1,t] <- mean(pcllcv[s-min_G+1,t,])
cat("g=",s," ", "lambda=",Lambda[t]," ", "pcll=",pcll[s-min_G+1,t],"\n")  
      }
    }

######### optimal g and lambda  ############
#
      occur2=which.max(pcll) 
      loc2=arrayInd(occur2,dim(pcll)) 
      opt_para[2,] <- c(G[loc2[1,1]],Lambda[loc2[1,2]],pcll[loc2[1,1],loc2[1,2]])
      opt_para
      out <- list(opt_para=opt_para, pcll=pcll)
}