next_data_generation<-function(x,level,transprob) {
  #
  # function for generating x_(r+1) data given x_r data using conditional prob. theta_(r,r+1)
  #
  # input: x(nx1): input data vector x_r
  #        level(m_(r+1),1): response level vector of x_(r+1)
  #        transprob: conditional prob. theta_(r,r+1)
  #
  # output: x_(r+1) (nx1)
  #
  n=length(x)
  y=rep(0,n)
  #
  r=rle(x)
  leng=r$lengths
  val=r$values
  k=0
  for ( i in 1:length(leng)){
    for (j in 1:leng[i]){
      y[k+j]=sample(level,1,prob=transprob[val[i],])
    }
  k=k+leng[i]
  }
  y
}

