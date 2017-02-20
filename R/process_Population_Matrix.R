get.Position.Leafs <- function(m)
{
  leaf = rep(FALSE,nrow(m))
  leaf[length(leaf)] = TRUE
  if(nrow(m)>1)
  {
    for(i in 1:(nrow(m)-1))
    {
      if(which(m[i,]!="")>=which(m[i+1,]!=""))
        leaf[i] = TRUE
    }
  }
  return(leaf)
}

get.Position.Matrix <- function(m)
{
  n = matrix("",ncol=ncol(m),nrow=nrow(m))
  n[which(m[,1]!=""),1] = 1
  ind = 1
  prev = 1
  for(i in 1:nrow(n))
  {
    id = which(m[i,]!="")
    if(id!=1)
    {
      if(id!=prev)
        ind=ind+1
      n[i,id] = ind
      prev = id
    } else
    {
      prev = 1
    }
  }
  return(n)
}