.combine.Models <- function(reference.fn,target.fn,out.dir,composite.model.call.rate)
{
  
  snpgdsCombineGeno(c(target.fn,reference.fn),file.path(out.dir,"Aggregated.gds"),
                    method = 'position',snpfirstdim = TRUE)
  
  if(file.exists(file.path(out.dir,"Aggregated.gds")))
  {
    # ancestry DataFrame have to be manually combined into aggregated model
    genofile <- snpgdsOpen(file.path(out.dir,"Aggregated.gds"),readonly = F)
    target.model = snpgdsOpen(target.fn)
    target.sample.annot = read.gdsn(index.gdsn(target.model,'sample.annot'))
    target.snp.ref = read.gdsn(index.gdsn(target.model,'snp.ref'))
    target.snp.alt = read.gdsn(index.gdsn(target.model,'snp.alt'))
    reference.model = snpgdsOpen(reference.fn)
    reference.sample.annot = read.gdsn(index.gdsn(reference.model,'sample.annot'))
    reference.snp.ref = read.gdsn(index.gdsn(reference.model,'snp.ref'))
    reference.snp.alt = read.gdsn(index.gdsn(reference.model,'snp.alt'))
    sample.annot <- rbind(target.sample.annot,reference.sample.annot)
    snp.ref <- c(target.snp.ref,reference.snp.ref)
    snp.alt <- c(target.snp.alt,reference.snp.alt)
    add.gdsn(genofile,"sample.annot",sample.annot)
    add.gdsn(genofile,"snp.ref",snp.ref)
    add.gdsn(genofile,"snp.alt",snp.alt)
    snpgdsClose(genofile)
    snpgdsClose(target.model)
    snpgdsClose(reference.model)
    
    return(TRUE)
  } else 
  {
    return(FALSE)
  }
}
