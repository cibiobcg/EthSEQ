plot.Report <-function(tab,samples,n.snps,polygons,out.dir,label="")
{
  colors = rainbow(length(polygons),alpha=0.8)
  pdf(file.path(out.dir,paste("Report",label,".pdf",sep="")),11, 8)
  
  legend.space = (max(tab$EV1)-min(tab$EV1))/7
  plot(NA,xlim=c(min(tab$EV1)-legend.space,max(tab$EV1)),ylim=c(min(tab$EV2),max(tab$EV2)),
       main=paste("Target samples with reference populations (#SNPs=",n.snps,")",sep=""),xlab="PC1",ylab="PC2")
  points(as.numeric(samples$EV1),as.numeric(samples$EV2),pch=16,cex=0.6)
  for(i in 1:length(polygons))
  {
    lines(polygons[[i]][,1],polygons[[i]][,2],lwd=3,col=colors[i])
  }
  legend("bottomleft",legend=names(polygons),pch=16, col=colors,bty = "n",cex=0.7)
  
  dev.off()
}