get.Ethnicity <- function(tab,space="2D")
{
  ## get population polygons
  populations = as.character(unique(tab$pop))
  populations = populations[populations!="ND"]
  colors = rainbow(length(populations),alpha=0.8)
  polygons = list()
  polygons.area = list()
  
  pca.comp = 3:4
  if (space=="3D")
      pca.comp = 3:5
  
  cols = rep(rgb(0,0,0,alpha=0.8),nrow(tab))
  for(i in 1:length(populations))
  {
    ids = which(tab$sample.id%in%tab$sample.id[which(tab$pop==populations[i])])
    cols[ids] = colors[i]
    tmp = tab[ids,]
    polygons[[length(polygons)+1]] = tmp[unique(as.vector(convhulln(tmp[,pca.comp]))),pca.comp]
    polygons.area[[length(polygons.area)+1]] = tmp[as.vector(t(convhulln(tmp[,pca.comp]))),pca.comp]
  }
  names(polygons) = populations
  names(polygons.area) = populations
  
  if(sum(tab$sample.id=='ND')>0)
  {
    samples = tab[tab$sample.id=='ND',]
    samples$type = NA
    samples$pop = ""
    samples$contribution = ""
    
    centroids = list()
    for(j in 1:length(polygons))
    {
      centroids[[length(centroids)+1]] = centroid.ethseq(polygons[[j]])
    }
    
    for(i in 1:nrow(samples))
    {
      dist = c()
      pnts = data.frame(samples[i,pca.comp])
      
      for(j in 1:length(polygons))
      {
        res = inhull.ethseq(pnts,polygons[[j]])
        dist = c(dist,distance.ethseq(pnts,centroids[[j]]))
        
        if(res==1)
        {
          samples$pop[i] = paste(c(strsplit(as.character(samples$pop[i]),"\\|")[[1]],populations[j]),collapse="|")
          samples$type[i] = "INSIDE"
        }
      }
      
      if(as.character(samples$pop[i])=="")
      {
        samples$pop[i] = paste(populations[which(dist==min(dist))],collapse="|")
        samples$type[i] = "CLOSEST"
        ## contribution
        dist = min(dist)/dist
        idx = which(dist>0.33)
        if(length(idx)==1)
          idx = order(dist,decreasing = T)[1:2]
        dist = dist[idx]
        dist = dist/sum(dist)
        samples$contribution[i] = paste(paste(populations[idx],"(",round(dist*100,digits=2),"%)",sep=""),collapse="|")
      }
      
      if(length(grep("\\|",as.character(samples$pop[i])))>0)
      {
        pop.elems = strsplit(as.character(samples$pop[i]),"\\|")[[1]]
        ## contribution
        idx1 = which(populations%in%pop.elems)
        dist = dist[idx1]
        dist = min(dist)/dist
        dist = dist/sum(dist)
        samples$contribution[i] = paste(paste(populations[idx1],"(",round(dist*100,digits=2),"%)",sep=""),collapse="|")
      }
      
    }
    
    return(list(samples,polygons,polygons.area,TRUE))
  } else
  {
    return(list(NA,polygons,polygons.area,FALSE))
  }
}
