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
