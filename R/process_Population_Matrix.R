check.Matrix <- function(m,populations)
{
  pos = c()
  for(i in 1:nrow(m))
  {
    id = which(m[i,]!="")
    if(length(id)==0)
      return(FALSE)
    res = TRUE
    if(id>1)
      res = all(c(1:(id-1))%in%pos)
    pos = c(pos,id)
    pops = strsplit(m[i,id],"\\|")[[1]]
    if(!all(pops%in%populations)|!res)
      return(FALSE)
  }
  return(TRUE)
}

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