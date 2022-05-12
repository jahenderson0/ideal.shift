############################################
# John Henderson
# - smoteFun: k-neighbor synthetic minority oversampling function
# - May 22, 2022
############################################

require(RcppHNSW)

############################################
# SMOTE (synthetic rare-class imputation) function
############################################
# - takes a matrix of covariates and vector of binary outcomes
# - augments data with synthetic cases using k nearest neighbor line projections
# - see, e.g., https://arxiv.org/abs/1106.1813
# - used to augment data when predicting rare class topics for bills

smoteFun=function(y,x,k=10){
  # which is the rare class?
  if(mean(y)<.5){
    mx=1
  } else {
    mx=0
  }

  # index x for the rare class side
  ix=which(y==mx)

  # collect indices for the k nearest neighborss for each index unit
  nn=hnsw_knn(x[ix,],k = k+1,distance = "cosine",
    M = 2000,ef_construction = 200,ef = 10,
    verbose = FALSE,progress = "bar",n_threads = 0,
    grain_size = 1)[[1]][,-c(1)]

  # imputed matrix object
  imat=matrix(NA,nrow(nn)*ncol(nn),ncol(x))
  set.seed(1005)
  tt=0
  # loop over rows of nn
  for(i in 1:nrow(nn)){
	   d=runif(0,1,n=k)
     for(j in 1:k){
       tt=tt+1
       imat[tt,]=x[ix[i],] + d[j]*(x[ix[nn[i,j]],]-x[ix[i],])
     }
   }

  xs=rbind(x,imat)
  ys=c(y,rep(mx,nrow(imat)))
  return(list("x"=xs,"y"=ys))
}

#END smoteFun
