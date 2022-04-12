.plot.Report <-function(tab,samples,n.snps,polygons,polygons.area,out.dir,label="",space="2D")
{
  pdf(file.path(out.dir,paste("Report",label,".pdf",sep="")),11, 8)
  
  colors = c(rainbow(length(polygons),alpha=0.8),"black")
  colors[2] = "orange"
  
  idx = which(names(polygons)%in%samples$pop)
  colors.scatter = colors[idx]
  cols = rep(length(polygons[names(polygons)%in%samples$pop])+1,nrow(samples))
  
  for(i in 1:length(polygons[names(polygons)%in%samples$pop]))
  {
    cols[which(samples$pop==names(polygons[names(polygons)%in%samples$pop])[i])] = i
  }
  if(!all(samples$pop%in%names(polygons[names(polygons)%in%samples$pop])))
  {
    colors.scatter = c(colors.scatter,"black")
  }
  
  names = c(names(polygons),"Multiple")
  
  if(space=="3D")
  {
    scatter3D(samples$EV1,samples$EV2,samples$EV3,pch=19,col=colors.scatter,colvar=cols,cex=0.5,
              main=paste("Target samples with reference populations (#SNPs=",n.snps,")",sep=""),
              xlim=c(min(tab$EV1),max(tab$EV1)),ylim=c(min(tab$EV2),max(tab$EV2)),
              zlim=c(min(tab$EV3),max(tab$EV3)),xlab="PCA1",ylab="PCA2",zlab="PCA3",
              colkey = FALSE)
    
    x = 0
    unit = (max(tab$EV1)-min(tab$EV1))/10
    offset = (max(tab$EV2)-min(tab$EV2))/4
    for(i in 1:(length(polygons)+1))
    {
      text3D(min(tab$EV1)+x,min(tab$EV2)-offset,min(tab$EV3),labels = names[i],add=TRUE,theta = 90,cex=0.6,col=colors[i])
      x = x+unit
    }
    
    for(i in 1:length(polygons))
    {
      triangle3D(as.matrix(polygons.area[[i]]),colvar=rep(1,nrow(polygons.area[[i]])/3),col=colors[i],shade=0.1,add=TRUE,alpha=0.07)
    }
  }
  
  legend.space = (max(tab$EV1)-min(tab$EV1))/7
  a=samples[,3:4]
  scatter2D(a[,1],a[,2],xlim=c(min(tab$EV1)-legend.space,max(tab$EV1)),
            main=paste("Target samples with reference populations (#SNPs=",n.snps,")",sep=""),
            ylim=c(min(tab$EV2),max(tab$EV2)),col=colors.scatter,
            pch=19,colvar=cols,colkey=FALSE,xlab="PCA1",ylab="PCA2",cex=0.5)
  for(i in 1:length(polygons))
  {
    a = chull(as.matrix(polygons[[i]])[,1:2])
    a = as.matrix(polygons[[i]])[,1:2][c(a,a[1]),]
    lines2D(a[,1],a[,2],xlim=c(min(tab$EV1),max(tab$EV1)),ylim=c(min(tab$EV2),max(tab$EV2)),
            col=colors[i],pch=19,add=TRUE,colvar=rep(1,nrow(a)),lwd=2)
  }
  legend("bottomleft",legend=c(names(polygons),"Multiple"),pch=16, col=c(colors,"grey"),bty = "n",cex=0.7)
  
  if(space == "3D")
  {
    legend.space = (max(tab$EV2)-min(tab$EV2))/7
    a=samples[,4:5]
    scatter2D(a[,1],a[,2],xlim=c(min(tab$EV2)-legend.space,max(tab$EV2)),
              main=paste("Target samples with reference populations (#SNPs=",n.snps,")",sep=""),
              ylim=c(min(tab$EV3),max(tab$EV3)),col=colors.scatter,
              pch=19,colvar=cols,colkey=FALSE,xlab="PCA2",ylab="PCA3",cex=0.5)
    for(i in 1:length(polygons))
    {
      a = chull(as.matrix(polygons[[i]])[,2:3])
      a = as.matrix(polygons[[i]])[,2:3][c(a,a[1]),]
      lines2D(a[,1],a[,2],xlim=c(min(tab$EV2,max(tab$EV2))),ylim=c(min(tab$EV3),max(tab$EV3)),col=colors[i],
              pch=19,add=TRUE,colvar=rep(1,nrow(a)),lwd=2)
    }
    legend("bottomleft",legend=c(names(polygons),"Multiple"),pch=16, col=c(colors,"grey"),bty = "n",cex=0.7)
    
    legend.space = (max(tab$EV1)-min(tab$EV1))/7
    a=samples[,c(3,5)]
    scatter2D(a[,1],a[,2],xlim=c(min(tab$EV1)-legend.space,max(tab$EV1)),
              main=paste("Target samples with reference populations (#SNPs=",n.snps,")",sep=""),
              ylim=c(min(tab$EV3),max(tab$EV3)),col=colors.scatter,
              pch=19,colvar=cols,colkey=FALSE,xlab="PCA1",ylab="PCA3",cex=0.5)
    for(i in 1:length(polygons))
    {
      a = chull(as.matrix(polygons[[i]])[,c(1,3)])
      a = as.matrix(polygons[[i]])[,c(1,3)][c(a,a[1]),]
      lines2D(a[,1],a[,2],xlim=c(min(tab$EV1,max(tab$EV1))),ylim=c(min(tab$EV3),max(tab$EV3)),col=colors[i],
              pch=19,add=TRUE,colvar=rep(1,nrow(a)),lwd=2)
    }
    legend("bottomleft",legend=c(names(polygons),"Multiple"),pch=16, col=c(colors,"grey"),bty = "n",cex=0.7)
  }
  
  dev.off()
}