sim_data_generation<-function(n,p,g,n_level,theta_marg1,theta_cond) {
  #
  # function for generating simulation data
  #
  # input: p: dimension of x
  #        n_level: response level vector of p-dimensional x: C(2,3,3,2)
  #        theta_marg1: marginal prob. of x_1: n_level[1]xg
  #        theta_cond: conditional prob. array:dim=(max_n_level,max_n_level,j,k,g)=c(3,3,4,4,2)
  #
  # output: x:array (nxpxg)
  #
  x <- array(rep(NA,prod(n*p*g)),dim=c(n,p,g))
  for (k in 1:g){
    x1 <- sample(1:n_level[1],n,replace=TRUE,prob=theta_marg1[,k])
    x[,1,k] <-x1   
    for (j in 2:p) {
      x[,j,k] <- next_data_generation(x[,(j-1),k],1:n_level[j],theta_cond[1:n_level[j-1],1:n_level[j],(j-1),j,k])
    }
  }
  x
}

