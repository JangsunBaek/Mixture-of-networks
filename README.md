# Mixture of networks
 Mixture of networks for clustering categorical data
 
R codes for Mixture of networks based on penalized composite likelihood

##### Important input data preparation for the first version of R code ##### 

Categorical response of the jth categorical variable, x_j, of an observation in data set must be an integer in 1 to r_j, j = 1, ..., p. That is, categorical responses of each variable should be transformed to sequential integers starting from 1. For example, if the response of X_1, X_2, X_3 belongs to two, three, and two categories, respectively, then x_1 is either 1 or 2, x_2 is either 1, 2, or 3, and x_3 is either 1 or 2. Therefore, the input data of (x_1, x_2, x_3) would be (1, 1, 1), (1, 2, 1), (2, 3, 1), (2, 3, 2), but neither (1, 4, 6) nor (3, 1, 3). Future versions of R code will accept any type of reponses of categorical variables.       

#################################################################

1. Packages needed for Mixture of networks

> install.packages("poLCA")  # LCA for obtaining initial parameters
> library("poLCA")

2. R codes 

2.1 “mnpcl”: clustering program
  #
  # mnpcl(X, g, n_level, lambda, initmax, itmax, tol_logL, init_para = NULL) {
  #
  # input: 
  #    X: categorical data matrix (nxp)  
  #    g: number of components
  #    n_level: vector of the number of response levels for variables: (r_1, r_2, ..., r_p)
  #    lambda: tuning parameter of penalty (must be less than min(obs_freq)),  where obs_freq[max(n_level),max(n_level),p,p] is the observed frequency in contingency tables
  #    initmax: number of initializations
  #    itmax: maximum number of EM iterations
  #    tol_logL: tolerance of convergence
  #
  #    
  # output(list):
  #    $pivec: vector of pi_i
  #    $logL: loglikelihood
  #    $cluster: vector of clustered labels in (1,2,..., g)
  #    $theta: bivariate joint prob. penalized mle
  #    $theta_tilde: bivariate joint prob. mle
  #    $thetacond_tilde: conditional prob. mle

2.2 “est_mnpcl”: parameter estimation using EM procedure
2.3 “mstep_mnpcl”: M-step
2.4 “tau_mnpcl”: posterior prob. estimation
2.5 “logL_mnpcl”: penalized composite log-likelihood
2.6 “sim_data_generation”: synthetic data generation
2.7 “error.rate”: clustering error

3. Zoo dataset, true class labels, and number of response levels 

3.1 xnew_zoo: csv file of shuffled Zoo data
3.2 zoo_level: csv file of the vector of the number of response levels  
3.3 zoonew_label: csv file of the class labels of xnew_zoo

4. Examples

4.1 Synthetic data
4.1.1 Generation of synthetic data 

# theta_cond: conditional probability of X_2 given X_1, X_3 given X_2, and X_4 given X_3 for two groups
> theta_marg1<-matrix(c(.01,.99,.99,.01),c(2,2)) # marginal probability of X_1 for two groups
> xx<-sim_data_generation(100,4,2,sim_level,theta_marg1,theta_cond)
> xx1<-xx[,,1]
> xx2<-xx[,,2]
> xx<-rbind(xx1,xx2)
> xx<-data.frame(xx)
> colnames(xx)<-c("x1","x2","x3","x4")
> perm_sim<-sample(1:200,replace=FALSE)

# synthetic data matrix and true labels
> x_sim<-xx[perm_sim,]  # 200 four-dimensional synthetic data (200x4)
> sim_label<-xx_label[perm_sim] # true group label vector

4.1.2 Clustering the synthetic data into two groups with lambda = 0.6 

> result<-mnpcl(X=x_sim, g=2, n_level=sim_level, lambda=0.6, initmax=10, itmax=500, tol_logL=1.e-3, init_para = NULL)
initialization =  1 
initialization =  2 
initialization =  3 
initialization =  4 
initialization =  5 
initialization =  6 
initialization =  7 
initialization =  8 
initialization =  9 
initialization =  10 
n_g = 114 # number of observations allocated into group1
n_g = 86 # number of observations allocated into group2

> ACC <- 1-error.rate(result$cluster, sim_label) # ACC for MN-PCL
[1] 0.92 # The ACC may be different because the synthetic data is generated randomly.

4.2 Clustering Zoo data into seven groups with lambda = 0.2

> result<-mnpcl(xnew_zoo, 7, zoo_level, 0.2, 10, 500, 1.e-3, init_para = NULL)
initialization =  1 
initialization =  2 
initialization =  3 
initialization =  4 
initialization =  5 
initialization =  6 
initialization =  7 
initialization =  8 
initialization =  9 
initialization =  10 
n_g = 14 
n_g = 4 
n_g = 10 
n_g = 8 
n_g = 37 
n_g = 7 
n_g = 21 
> ACC <- 1-error.rate(result$cluster, zoonew_label)
[1] 0.8910891
> str(result)
List of 7
 $ pivec          : num [1:7] 0.1295 0.0447 0.0963 0.0875 0.3656 ...
 $ theta          : num [1:6, 1:6, 1:16, 1:16, 1:7] NA NA NA NA NA NA NA NA ...
 $ logL           : num -16344
 $ cluster        : int [1:101] 1 1 1 3 7 5 4 7 5 5 ...
 $ theta_tilde     : num [1:6, 1:6, 1:16, 1:16, 1:7] 1 NA NA NA NA NA NA 0 NA  ...
 $ thetacond_tilde: num [1:6, 1:6, 1:16, 1:16, 1:7] NA NA NA NA NA NA NAA NA ...
 $ call           : language mnpcl(X = xnew_zoo, g = 7, n_level = zoo_level, lambda = 0.2, initmax = 10, itmax = 500, tol_logL = 0.001, init_para = NULL)
 - attr(*, "class")= chr "mnpcl"
> 

