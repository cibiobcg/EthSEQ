.check.Sibling.Pops = function(pop.nodes)
{
  return(any(duplicated(unlist(strsplit(pop.nodes,'\\|')))))
}

.check.Parent.Pops = function(siblings, populations)
{
  !all(unlist(strsplit(siblings,"\\|"))%in%populations)
}

.check.Matrix <- function(m,populations)
{
  if(!is.null(nrow(m))&!is.null(ncol(m))){
    idx = which(m[,1]!="")
    
    olapsPops = .check.Sibling.Pops(m[idx,1])
    parentPops = .check.Parent.Pops(m[idx,1],populations)
    
    idx = c(idx,nrow(m)+1)
    
    r = lapply(1:(length(idx)-1),function(x)
    {
      .check.Matrix(m = m[idx[x]:(idx[x+1]-1),-1], populations = unlist(strsplit(m[idx[x],1],"\\|")))
    })
  } else {
    return(FALSE)
  }
  return(any(c(unlist(r),olapsPops,parentPops)))
}

.get.Position.Leafs <- function(m)
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

.get.Position.Matrix <- function(m)
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