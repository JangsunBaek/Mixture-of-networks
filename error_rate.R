error.rate<-function(clust1,clust2)
{
  
  clust1 <- unclass(as.ordered(clust1))
  clust2 <- unclass(as.ordered(clust2))
  
  if((n=length(clust1))!=length(clust2))
  {warning("error: length not equal");return}
  
  if( (g=length(table(clust1)))!=length(table(clust2)))
  {warning("the number of clusters are not equal");return}
  
  permute<-function(a)
  {
    n<-length(a)
    if(n==1)
      f<-a
    else
    {
      nm<-gamma(n)
      f<-array(0,c(n,n*nm))
      j<-1
      
      for(i in a)
      {
        f[1, (j-1)*nm+1:nm]<-i
        f[-1,(j-1)*nm+1:nm]<-permute(setdiff(a,i))
        j<-j+1
      }
    }
    
    f
  }
  
  
  #
  id<-1:n
  
  cmb<-permute(1:g)
  
  nperm<-ncol(cmb)
  
  rate<-rep(0,nperm)
  
  #
  for(i in 1:nperm)
  {
    
    tmp<-rep(0,g)
    
    tc<-rep(0,n)
    
    for(j in 1:g)
      tc[clust2==j]=cmb[j,i]
    
    for(j in 1:g)
    {  
      tmp1<-0 
      
      for(k in (1:g)[-j])
        tmp1<-tmp1+length(intersect(id[clust1==j],id[tc==k]))
      
      tmp[j]<-tmp1
    }
    
    rate[i]<-sum(tmp)/n
  }
  
  min(rate)
}