
.distance.ethseq <- function(point,centroid)
{
  return(dist(rbind(centroid,point),method="euclidean"))
}

.centroid.ethseq <- function(k)
{
  t <- delaunayn(k,options="Qt");
  n <- dim(t)[1];
  w <- matrix(0,n,1) 
  c <- rep(0,dim(k)[2])
  for(m in 1:n)
  {
    sp = k[t[m,],]
    w[m] = convhulln(sp,options="FA")[[2]]
    c = c+w[m]*apply(sp,2,mean)
  }

  return(c/sum(w))
}

.inhull.ethseq <- function(testpts, 
                          calpts,
                          hull=convhulln(calpts))
{ 
  # R implementation of  http://www.mathworks.com/matlabcentral/fileexchange/10226-inhul
  # from https://github.com/cran/hypervolume/blob/master/R/inhull.R
  
  calpts <- as.matrix(calpts) 
  testpts <- as.matrix(testpts) 
  
  tol <- 1.e-13*mean(abs(calpts))
  
  p <- dim(calpts)[2] # columns in calpts
  cx <- dim(testpts)[1] # rows in testpts
  nt <- dim(hull)[1] # number of simplexes in hull 
  # find normal vectors to each simplex 
  nrmls <- matrix(NA, nt, p) # predefine each nrml as NA, degenerate
  
  degenflag <- matrix(TRUE, nt, 1) 
  for (i in 1:nt) { 
    nullsp <- t(Null(t(calpts[hull[i,-1],] - matrix(calpts[hull[i,1],],p-1,p, byrow=TRUE))))
    
    if (dim(nullsp)[1] == 1) { nrmls[i,] <- nullsp
    
    degenflag[i] <- FALSE}} 
  # Warn of degenerate faces, and remove corresponding normals 
  if(length(degenflag[degenflag]) > 0) warning(length(degenflag[degenflag])," degenerate faces in convex hull")
  
  nrmls <- nrmls[!degenflag,] 
  nt <- dim(nrmls)[1] 
  # find center point in hull, and any (1st) point in the plane of each simplex
  
  center = apply(calpts, 2, mean) 
  a <- calpts[hull[!degenflag,1],] 
  # scale normal vectors to unit length and ensure pointing inwards 
  nrmls <- nrmls/matrix(apply(nrmls, 1, function(x) sqrt(sum(x^2))), nt, p)
  dp <- sign(apply((matrix(center, nt, p, byrow=TRUE)-a) * nrmls, 1, sum))
  nrmls <- nrmls*matrix(dp, nt, p) 
  # if min across all faces of dot((x - a),nrml) is 
  # +ve then x is inside hull 
  # 0 then x is on hull 
  # -ve then x is outside hull 
  # Instead of dot((x - a),nrml) use dot(x,nrml) - dot(a, nrml) 
  aN <- diag(a %*% t(nrmls)) 
  val <- apply(testpts %*% t(nrmls) - matrix(aN, cx, nt, byrow=TRUE), 1,min) 
  # code values inside 'tol' to zero, return sign as integer 
  val[abs(val) < tol] <- 0 
  return(as.integer(sign(val))) 
}
