mnpcl <- function(X, g, n_level, lambda, initmax, itmax, tol_logL, init_para = NULL) {
  # Mixture of networks based on penalized composite likelihood
  #
  # install.packages("klaR") for K-modes 
  # library(klaR)
  # install.packages("poLCA") for LCA
  # library("poLCA")
  #
  #
  # input: 
  #    X: categorical data (nxp)  
  #    g: number of components
  #    n_level: vector of the number of response levels for variables: (r_1, r_2, ..., r_p)
  #    lambda: tuning parameter of penalty (must be less than max(obs_freq)),  
  #            where obs_freq[max(n_level),max(n_level),p,p] is the observed frequency in contingency tables
  #    initmax: number of initializations
  #    itmax: maximum number of EM iterations
  #    tol_logL: tolerance of convergence
  #
  #    
  # output:
  #    pivec: vector of pi_i
  #    logL: loglikelihood
  #    cluster: vector of clustered labels in (1,2,..., g)
  #    theta: bivariate joint prob. penalized mle
  #    theta_tilde: bivariate joint prob. mle
  #    thetacond_tilde: conditional prob. mle
  #

  eps <- 1.e-5  # tiny probability for the empty cell
  #
  if (!is.matrix(X)){
    X <- as.matrix(X)
  }

  p <- ncol(X)    # p >1
  n <- nrow(X)
  warn_msg <- NULL
  maxLOGL <- -Inf

  
  #######################################################
  # Transform the data X to the observed frequency array
  #
  rmax=max(n_level)
  obs_freq <- array(rep(0,prod(rmax^2*p^2)),c(rmax,rmax,p,p))
  for (j in 1:(p-1)) {
    for (k in (j+1):p) {
      Xjk <- X[,c(j,k)]
      for (i in 1:n) {
        obs_freq[Xjk[i,1],Xjk[i,2],j,k] <- obs_freq[Xjk[i,1],Xjk[i,2],j,k] +1
      }
    }
  }
  
  
  #######################################################
  # Transform the data X to the Adjacency matrix data
  
  A <- A_array(X,n_level)
  
  ###################################################
  # initializations initmax times
  #
  for (s in 1:initmax) {
    cat("initialization = ", s, "\n")  
    
    ### initial parameter estimation for pivec, w_g^(r,s), theta_g^r and theta_g^(r,s) for components
    ### 1) initial parameter estimation method: K-modes algorithm ###
    # install.packages("klaR")
    # library(klaR)  
    #  test<-FALSE
    #  while(test==FALSE){
    #    y <- kmodes(X, g, iter.max = 10, weighted = FALSE)
    #    init_label <- y$cluster
    #    cl_size<-y$size
    #    pivec <- cl_size/n  
    #    if(length(cl_size)==g) test<-TRUE
    #    else test<-FALSE
    #    }
    #
    ### 2) initial parameter estimation method: LCA algorithm  ###
    # install.packages("poLCA")
    # library("poLCA")

    df_sim<-data.frame(X)
    number<-c(1:p)
    col_label<-paste0('x',number)
    colnames(df_sim) <- col_label 
    x_var<-paste0(col_label,collapse=',')
    code_line<-paste0("f","<-cbind(",x_var,') ~ 1')
    eval(parse(text=code_line))

    test<-FALSE
    model.frame(f, df_sim, na.action = NULL)
    
    while(test==FALSE){
      yy <- poLCA(f, df_sim, nclass = g, verbose=FALSE)
      init_label <- yy$predclass
      cl_size <- as.numeric(table(init_label))
      pivec <- yy$P
      if(length(cl_size)==g) test<-TRUE
      else test<-FALSE
    }
    ###################################################

    
    theta <- array(rep(NA,prod(max(n_level)^2*p^2*g)),dim=c(max(n_level),max(n_level),p,p,g))
    #
    #
    for (i in 1:g){
      index <- (init_label==i)
      Xg <- X[index, byrow=T]
      Ag <- A[ , , , , index]
      
      for (j in 1:(p-1)) {
        for (k in (j+1):p) {
          if (length(index[index==TRUE])==1) {
            jointprob <- apply(Ag[1:n_level[j],1:n_level[k],j,k],c(1,2),sum)/cl_size[i] + eps
          }
          else{jointprob <- apply(Ag[1:n_level[j],1:n_level[k],j,k,],c(1,2),sum)/cl_size[i] + eps}
          theta[1:n_level[j],1:n_level[k],j,k,i]  <-  jointprob/sum(jointprob)
        }
      }
    }
    
    #
    # Fit a model using initial parameter estimates from pairwise contingency tables
    #
    # EM steps
    
    estd_model <- est.mnpcl(X=X, obs_freq=obs_freq, lambda, pivec=pivec, n_level=n_level, theta=theta, 
                             itmax = itmax, tol_logL = tol_logL)
    
    if ((class(estd_model) == "error")) {
      cat(estd_model)
      return(estd_model)
    }
    
    # keep the model with highest log-likelihood
    if ((class(estd_model) == "mnpcl")) {
      if (estd_model$logL > maxLOGL) {
        Hmodel <- estd_model
        maxLOGL <- Hmodel$logL
      }
    }
  }
  
  ###############################################################################
  # clustering
  ###############################################################################
  
  Sijrs <- array(rep(0,n*g*p*p),c(n,g,p,p))
  Sij <- matrix(rep(0,n*g), c(n,g))
  log_theta <- rep(NA, n)
  #
  for (j in 1 : g) {
    for (r in 1:(p-1)) {
      for (s in (r+1):p) {
        #####  find location of the observation in (r,s)th adjacency table
        occur=which(A[,,r,s,]==1) 
        loc=arrayInd(occur,dim(A[,,r,s,]))   
        
        for (k in 1:n) {
          log_theta[k] <- log(Hmodel$theta[loc[k,1],loc[k,2],r,s,j])
        }
        Sijrs[, j, r, s] <- log_theta
      }
    }
    for (i in 1:n) {
      for (r in 1:(p-1)){
        for (s in (r+1):p) {
          Sij[i,j] <- Sij[i,j] + Sijrs[i,j,r,s]
        }
      }
    }
  }
  
  Sij <- sweep(exp(Sij),2,Hmodel$pivec, '*')
  Hmodel$cluster <- rep(NA, n)
  
  if(any(rowSums(Sij)==0)==FALSE) {
    Sij <- sweep (Sij, 1, rowSums(Sij), '/')
    Hmodel$cluster <- apply(Sij, 1, which.max) 
  }
  else{
    for (i in 1:n) {
      if (sum(Sij[i,]) == 0) {
        Hmodel$cluster[i] <- sample(1:g, 1) 
      }
      Hmodel$cluster[i] <- which.max(Sij[i,])
    }
  }

  ################################################################################    
  # estimates of lower triangular elements of joint prob., margianl prob., 
  # and conditional prob.
  #######################
  # marginal pmf of X_j of group i: diag(Hmodel$theta(1:n_level[j],1:n_level[j],j,j,i))  
  thetacond_tilde <- array(rep(NA,prod(max(n_level)^2*p^2*g)),dim=c(max(n_level),max(n_level),p,p,g))
  theta_tilde <- array(rep(NA,prod(max(n_level)^2*p^2*g)),dim=c(max(n_level),max(n_level),p,p,g))
  
  
  cl_size <- as.numeric(table(Hmodel$cluster))
  if(length(cl_size) < g) {
    Hmodel <- paste('The number of clustered groups is less than the ', g, sep='')
    class(Hmodel) <- "error"
    return(Hmodel)
  }
  
  ### lower triangle penalized theta estimates
  for (i in 1:g){
    for (j in 1:(p-1)) {
      for (k in (j+1):p) {
        Hmodel$theta[1:n_level[k],1:n_level[j],k,j,i]  <- t(Hmodel$theta[1:n_level[j],1:n_level[k],j,k,i])
      }
    }
  }
  
  ### estimates of marginal prob. and conditional prob.  
  for (i in 1:g){
    index <- (Hmodel$cluster==i)
    Xg <- X[index, byrow=T]
    n_g <- length(index[index==TRUE])
    cat("n_g =", n_g, "\n")
    
    rmax=max(n_level)
    obs_freq <- array(rep(0,prod(rmax^2*p^2)),c(rmax,rmax,p,p))
    for (j in 1:(p-1)) {
      for (k in (j+1):p) {
        Xjk <- Xg[,c(j,k)]
        for (m in 1: n_g) {
          obs_freq[Xjk[m,1],Xjk[m,2],j,k] <- obs_freq[Xjk[m,1],Xjk[m,2],j,k] +1
        }
        
        jointprob <- obs_freq[1:n_level[j],1:n_level[k],j,k]/n_g
        theta_tilde[1:n_level[j],1:n_level[k],j,k,i] <- jointprob/sum(jointprob)
        theta_tilde[1:n_level[k],1:n_level[j],k,j,i] <- t(theta_tilde[1:n_level[j],1:n_level[k],j,k,i])
      }
    }
    
    for (j in 1:p) {
      if (j < p) {
        for (l in 1:n_level[j]) {
          theta_tilde[l,l,j,j,i] <- sum(obs_freq[l,1:n_level[j+1],j,j+1])/n_g # marginal prob. for X_1 - X_{p-1}
        }
      }
      else {
        for (l in 1:n_level[j]) {
          theta_tilde[l,l,p,p,i] <- sum(obs_freq[1:n_level[j-1],l,j-1,j])/n_g # marginal prob. for X_p
        }
      }
      for (k in (1:p)) {
        thetacond_tilde[1:n_level[j],1:n_level[k],j,k,i] <-  sweep(theta_tilde[1:n_level[j],1:n_level[k],j,k,i],1,apply(theta_tilde[1:n_level[j],1:n_level[k],j,k,i],1,sum),"/")
      }
    }
  }
  
  Hmodel$theta_tilde <-theta_tilde
  Hmodel$thetacond_tilde <- thetacond_tilde
  
  ################################################################################
  ## modefied BIC
  ###############
  # number of effective(non-close zero) theta(j,k): 
  #de <- 0
  #for (j in 1:(p-1)) {
  #  for (k in (j+1):p) {
  #    d_non <- length(Hmodel$theta[1:n_level[j],1:n_level[k],j,k,] <= eps)
  #    de <- de + n_level[j]*n_level[k] - 1 - d_non
  #  }
  #}
  # total number of effective parameters
  #d <- (g-1) + de 
  #Hmodel$BIC <- -2*Hmodel$logL + d * log(n)
  
  ################################################################################  
  Hmodel$call <- match.call()
  return(Hmodel)
}  
