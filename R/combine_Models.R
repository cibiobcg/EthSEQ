combine.Models <- function(reference.fn,target.fn,out.dir,composite.model.call.rate)
{
  reference = snpgdsOpen(reference.fn)
  reference.geno = snpgdsGetGeno(reference,verbose=F)
  reference.snps = read.gdsn(index.gdsn(reference,"snp.rs.id"))
  reference.samples = read.gdsn(index.gdsn(reference,"sample.id"))
  reference.alleles = read.gdsn(index.gdsn(reference,"snp.allele"))
  reference.annot = read.gdsn(index.gdsn(reference,"sample.annot"))
  reference.ref = read.gdsn(index.gdsn(reference,"snp.ref"))
  reference.alt = read.gdsn(index.gdsn(reference,"snp.alt"))
  
  target = snpgdsOpen(target.fn)
  target.geno = snpgdsGetGeno(target,verbose = F)
  target.snps = read.gdsn(index.gdsn(target,"snp.rs.id"))
  target.samples = read.gdsn(index.gdsn(target,"sample.id"))
  target.alleles = read.gdsn(index.gdsn(target,"snp.allele"))
  target.annot = read.gdsn(index.gdsn(target,"sample.annot"))
  target.ref = read.gdsn(index.gdsn(target,"snp.ref"))
  target.alt = read.gdsn(index.gdsn(target,"snp.alt"))
  
  if(!all(all(target.snps==reference.snps)&all(target.ref==reference.ref)&all(target.alt==reference.alt)&
          all(read.gdsn(index.gdsn(target,"snp.chromosome"))==read.gdsn(index.gdsn(reference,"snp.chromosome")))&
          all(read.gdsn(index.gdsn(target,"snp.position"))==read.gdsn(index.gdsn(reference,"snp.position")))))
  {
    return(FALSE)
  }
  
  idx = which(target.alleles!=reference.alleles)
  if(length(idx)>0)
  {
    for(i in idx)
    {
      tmp = target.geno[,i]
      tmp[which(tmp==0)] = 4
      tmp[which(tmp==2)] = 0
      tmp[which(tmp==4)] = 2
      target.geno[,i] = tmp
    }
  }
  
  tmp.list = list(target.geno,reference.geno)
  genmat = as.numeric(rbind(tmp.list))
  #calls = apply(genmat,2,function(x) length(which(is.na(x)))/length(x))
  
  snpgdsCreateGeno(file.path(out.dir,"Aggregated.gds"),
                   genmat = genmat,
                   sample.id = c(paste("target.",target.samples,sep=""),paste("reference.",reference.samples,sep="")),
                   snp.id = 1:ncol(genmat),
                   snp.rs.id = reference.snps,
                   snp.chromosome = read.gdsn(index.gdsn(reference,"snp.chromosome")),
                   snp.position = read.gdsn(index.gdsn(reference,"snp.position")),
                   snp.allele = reference.alleles,
                   snpfirstdim=FALSE)
  
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
