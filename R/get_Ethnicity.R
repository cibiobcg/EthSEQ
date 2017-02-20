get.Ethnicity <- function(tab)
{
  ## get population polygons
  populations = as.character(unique(tab$pop))
  populations = populations[populations!="ND"]
  colors = rainbow(length(populations),alpha=0.8)
  polygons = list()
  
  cols = rep(rgb(0,0,0,alpha=0.8),nrow(tab))
  for(i in 1:length(populations))
  {
    ids = which(tab$sample.id%in%tab$sample.id[which(tab$pop==populations[i])])
    cols[ids] = colors[i]
    tmp = tab[ids,]
    hpts = chull(tmp$EV1,tmp$EV2)
    hpts = c(hpts, hpts[1])
    polygons[[length(polygons)+1]] = tmp[hpts,3:4]
  }
  names(polygons) = populations
  
  samples = tab[grep("target.",tab$sample.id),]
  samples$type = NA
  samples$pop = ""
  samples$contribution = ""
  
  for(i in 1:nrow(samples))
  {
    dist = c()
    for(j in 1:length(polygons))
    {
      polypnts = polygons[[j]]
      pnts = data.frame(samples$EV1[i],samples$EV2[i])
      colnames(polypnts) = c("x","y")
      colnames(pnts) = c("x","y")
      
      polypnts = readWKT(paste("POLYGON((",paste(paste(polypnts$x,polypnts$y,sep=" "),collapse=","),"))",sep=""))
      pnts = readWKT(paste("POINT(",paste(paste(pnts$x,pnts$y,sep=" "),collapse=","),")",sep=""))
      res = gContains(polypnts,pnts)
      
      centroid = gCentroid(polypnts)
      dist = c(dist,gDistance(centroid,pnts))
      
      if(res)
      {
        samples$pop[i] = paste(c(strsplit(as.character(samples$pop[i]),"\\|")[[1]],populations[j]),collapse="|")
        samples$type[i] = "INSIDE"
      }
    }
    
    if(samples$pop[i]=="")
    {
      samples$pop[i] = paste(populations[which(dist==min(dist))],collapse="|")
      samples$type[i] = "CLOSEST"
      dist = min(dist)/dist
      idx = which(dist>0.33)
      if(length(idx)==1)
        idx = order(dist,decreasing = T)[1:2]
      dist = dist[idx]
      dist = dist/sum(dist)
      samples$contribution[i] = paste(paste(populations[idx],"(",round(dist*100,digits=2),"%)",sep=""),collapse="|")
    }
    
  }
  
  return(list(samples,polygons))
}