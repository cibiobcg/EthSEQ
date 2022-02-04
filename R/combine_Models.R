combine.Models <- function(reference.fn,target.fn,out.dir,composite.model.call.rate)
{
  reference = snpgdsOpen(reference.fn,allow.duplicate = T)
  reference.geno = snpgdsGetGeno(reference,verbose=F)
  reference.snps = read.gdsn(index.gdsn(reference,"snp.rs.id"))
  reference.samples = read.gdsn(index.gdsn(reference,"sample.id"))
  reference.alleles = read.gdsn(index.gdsn(reference,"snp.allele"))
  reference.annot = read.gdsn(index.gdsn(reference,"sample.annot"))
  reference.ref = read.gdsn(index.gdsn(reference,"snp.ref"))
  reference.alt = read.gdsn(index.gdsn(reference,"snp.alt"))
  reference.chr = read.gdsn(index.gdsn(reference,"snp.chromosome"))
  reference.pos = read.gdsn(index.gdsn(reference,"snp.position")) 
  reference.sign = paste(reference.chr,reference.pos,sep=":")
  
  target = snpgdsOpen(target.fn,allow.duplicate = T)
  target.geno = snpgdsGetGeno(target,verbose = F)
  target.snps = read.gdsn(index.gdsn(target,"snp.rs.id"))
  target.samples = read.gdsn(index.gdsn(target,"sample.id"))
  target.alleles = read.gdsn(index.gdsn(target,"snp.allele"))
  target.annot = read.gdsn(index.gdsn(target,"sample.annot"))
  target.ref = read.gdsn(index.gdsn(target,"snp.ref"))
  target.alt = read.gdsn(index.gdsn(target,"snp.alt"))
  target.chr = read.gdsn(index.gdsn(target,"snp.chromosome"))
  target.pos = read.gdsn(index.gdsn(target,"snp.position")) 
  target.sign = paste(target.chr,target.pos,sep=":")
  
  intersect.sign = intersect(target.sign,reference.sign)
  
  if(length(intersect.sign)>0)
  {
    ## Restrict reference model to target model SNPs
    idx = which(reference.sign%in%intersect.sign)
    reference.geno = reference.geno[idx,,drop=FALSE]
    reference.snps = reference.snps[idx]
    reference.alleles = reference.alleles[idx] 
    reference.ref = reference.ref[idx]
    reference.alt = reference.alt[idx]
    reference.chr = reference.chr[idx]
    reference.pos = reference.pos[idx]
    reference.sign = paste(reference.chr,reference.pos,sep=":")
    
    idx = which(target.sign%in%intersect.sign)
    target.geno = target.geno[idx,,drop=FALSE]
    target.snps = target.snps[idx]
    target.alleles = target.alleles[idx]
    target.ref = target.ref[idx]
    target.alt = target.alt[idx]
    target.chr = target.chr[idx]
    target.pos = target.pos[idx]
    target.sign = paste(target.chr,target.pos,sep=":")
    
    isort = match(reference.sign,target.sign)
    target.geno = target.geno[isort,,drop=FALSE]
    target.snps = target.snps[isort]
    target.alleles = target.alleles[isort] 
    target.ref = target.ref[isort]
    target.alt = target.alt[isort]
    target.chr = target.chr[isort]
    target.pos = target.pos[isort]
  } else
  {
    snpgdsClose(target)
    snpgdsClose(reference)
    return(FALSE)
  }
  
  if(length(target.ref)!=length(reference.ref))
  {
    snpgdsClose(target)
    snpgdsClose(reference)
    return(FALSE)
  }
  
  if(!all(target.alt==reference.alt))
  {
    #Exclude SNPs with different alternative  
    idx = which(target.alt==reference.alt)
    
    reference.geno = reference.geno[idx,,drop=FALSE]
    reference.snps = reference.snps[idx]
    reference.alleles = reference.alleles[idx] 
    reference.ref = reference.ref[idx]
    reference.alt = reference.alt[idx]
    reference.chr = reference.chr[idx]
    reference.pos = reference.pos[idx]
    reference.sign = paste(reference.chr,reference.pos,sep=":")
    
    target.geno = target.geno[idx,,drop=FALSE]
    target.snps = target.snps[idx]
    target.alleles = target.alleles[idx]
    target.ref = target.ref[idx]
    target.alt = target.alt[idx]
    target.chr = target.chr[idx]
    target.pos = target.pos[idx]
    target.sign = paste(target.chr,target.pos,sep=":")
  }
  
  if(length(target.ref)==0)
  {
    snpgdsClose(target)
    snpgdsClose(reference)
    return(FALSE)
  }
  
  if(!all(all(target.ref==reference.ref)&all(target.alt==reference.alt)&
          all(target.chr==reference.chr)&all(target.pos==reference.pos)))
  {
    snpgdsClose(target)
    snpgdsClose(reference)
    return(FALSE)
  }
  
  idx = which(target.alleles!=reference.alleles)
  if(length(idx)>0)
  {
    for(i in idx)
    {
      tmp = target.geno[i,]
      tmp[which(tmp==0)] = 4
      tmp[which(tmp==2)] = 0
      tmp[which(tmp==4)] = 2
      target.geno[i,] = tmp
    }
  }
  
  tmp.list = list(target.geno,reference.geno)
  genmat = do.call(cbind,tmp.list)
  
  #calls = apply(genmat,2,function(x) length(which(is.na(x)))/length(x))
  
  snpgdsCreateGeno(file.path(out.dir,"Aggregated.gds"),
                   genmat = genmat,
                   sample.id = c(paste("target.",target.samples,sep=""),paste("reference.",reference.samples,sep="")),
                   snp.id = 1:nrow(genmat),
                   snp.rs.id = target.snps,
                   snp.chromosome = reference.chr,
                   snp.position = reference.pos,
                   snp.allele = reference.alleles,
                   snpfirstdim=TRUE)
  
  snpgdsClose(target)
  snpgdsClose(reference)
  
  genofile <- snpgdsOpen(file.path(out.dir,"Aggregated.gds"),readonly = F)
  sample.annot <- rbind(target.annot,reference.annot)
  add.gdsn(genofile, "sample.annot", sample.annot)
  
  add.gdsn(genofile, "snp.ref", reference.ref)
  add.gdsn(genofile, "snp.alt", reference.alt)
  
  snpgdsClose(genofile)
  return(TRUE)
  
}
library(SNPRelate)
reference.fn = "/shares/CIBIO-Storage/BCGLAB/EthSEQ_salmon/modelCreationSnakemake/out/hg38/models/SS5_Regions_3.3.gds"
target.fn = "/home/ddalfovo/Desktop/Target.gds"
out.dir = "/home/ddalfovo/Desktop/"
combine.Models(reference.fn,target.fn,out.dir)

snpgdsCombineGeno(c(reference.fn,target.fn),out.fn = '/home/ddalfovo/Desktop/Agg.gds',
                  method = 'position',snpfirstdim = TRUE)
genofile <- snpgdsOpen('/home/ddalfovo/Desktop/Agg.gds',readonly = F,allow.duplicate = T)
one = read.gdsn(index.gdsn(snpgdsOpen(reference.fn),'sample.annot'))
two = read.gdsn(index.gdsn(snpgdsOpen(target.fn),'sample.annot'))
sample.annot <- rbind(one,two)
add.gdsn(genofile,"sample.annot",sample.annot)
snpgdsClose(genofile)

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
  
  if(sum(tab$pop=='ND')>0)
  {
    samples = tab[grep("ND",tab$pop),]
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

aggregated = snpgdsOpen("/home/ddalfovo/Desktop/Aggregated.gds")
model <- aggregated
pca <- snpgdsPCA(model,num.thread = 1,eigen.method = "DSPEVX",verbose = TRUE,missing.rate = 1-1,eigen.cnt=5)
sample.id <- read.gdsn(index.gdsn(model, "sample.id"))
pop_code <- read.gdsn(index.gdsn(model, "sample.annot/pop.group"))
tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],
                  EV2 = pca$eigenvect[,2],
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  EV5 = pca$eigenvect[,5],
                  stringsAsFactors = FALSE)
annotations = get.Ethnicity(tab,'3D')
out = annotations[[1]]
out = out[,c(1,2,8,9)]
out[,1] = gsub("target.","",out[,1])
out
snpgdsClose(aggregated)


agg = snpgdsOpen("/home/ddalfovo/Desktop/Agg.gds")
model <- agg
pca <- snpgdsPCA(model,num.thread = 1,eigen.method = "DSPEVX",verbose = TRUE,missing.rate = 1-1,eigen.cnt=5)
sample.id <- read.gdsn(index.gdsn(model, "sample.id"))
pop_code <- read.gdsn(index.gdsn(model, "sample.annot/pop.group"))
tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],
                  EV2 = pca$eigenvect[,2],
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  EV5 = pca$eigenvect[,5],
                  stringsAsFactors = FALSE)
annotations = get.Ethnicity(tab,'3D')
out = annotations[[1]]
out = out[,c(1,2,8,9)]
out[,1] = gsub("target.","",out[,1])
out
snpgdsClose(agg)

