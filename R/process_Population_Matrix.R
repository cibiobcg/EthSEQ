.check.Sibling.Pops = function(pop.nodes)
{
  return(any(duplicated(unlist(strsplit(pop.nodes,'\\|')))))
}

.check.Parent.Pops = function(siblings, populations)
{
  !all(unlist(strsplit(siblings,"\\|"))%in%populations)
}

.check.Matrix <- function(m,populations,col)
{
  if(!is.null(nrow(m))&!is.null(ncol(m))){
    idx = which(m[,col+1]!="")
    
    olapsPops = .check.Sibling.Pops(m[idx,col+1])
    parentPops = .check.Parent.Pops(m[idx,col+1],populations)
    
    idx = c(idx,nrow(m)+1)
    
    r = lapply(1:(length(idx)-1),function(x)
    {
      .check.Matrix(m = m[idx[x]:(idx[x+1]-1),], populations = unlist(strsplit(m[idx[x],col+1],"\\|")),col = col+1)
    })
  } else {
    return(FALSE)
  }
  any(c(unlist(r),olapsPops,parentPops))
}



.check.Matrix <- function(m,populations, col.m, row.start.m, row.stop.m)
{
  if(row.stop.m-row.start.m!=0&col.m<=ncol(m)){
    
    siblingNodes = m[row.start.m:row.stop.m,col.m]
    parentPops = .check.Parent.Pops(siblingNodes[siblingNodes!=""],populations)
    
    idx = which(siblingNodes!="")
    idx = (idx+row.start.m)-1
    olapsPops = .check.Sibling.Pops(m[idx,col.m])
    idx = c(idx,row.stop.m+1)
    r = lapply(1:(length(idx)-1),function(x)
    {
      .check.Matrix(m,populations = unlist(strsplit(siblingNodes[x],"\\|")),col.m = col.m+1, row.start.m = idx[x], row.stop.m = idx[x+1]-1)
    })
  } else {
    return(FALSE)
  }
  any(c(unlist(r),olapsPops,parentPops))
}

.check.Matrix <- function(m,populations)
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