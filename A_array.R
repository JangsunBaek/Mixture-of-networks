A_array<-function(x,n_level){
  #
  # function for converting the input data to adjacency array
  #
  # input: x(nxp): input data matrix
  #       n_level(px1): number of response levels for categorical variables
  #
  # output: A(max(n_level),max(n_level),p,p,n): adjacency array A(k,l,r,s,i), 
  #               i in 1:n, r in 1: p-1, s in (r+1):p, k in 1:m_r, l in 1:m_s
  #
  n=nrow(x)
  p=length(n_level)
  m=max(n_level)
  A=array(rep(0,prod(m^2*p^2*n)),c(m,m,p,p,n))
  for (i in 1:n) {
    for (r in 1:(p-1)) {
      for (s in (r+1):p) {
        k=x[i,r]
        l=x[i,s]
        A[k,l,r,s,i]=1
      }
    }
  }
  A
}
