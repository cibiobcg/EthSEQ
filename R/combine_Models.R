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
    # target.sign = paste(read.gdsn(index.gdsn(target.model,'snp.chromosome')),read.gdsn(index.gdsn(target.model,'snp.position')),sep=":")
    
    reference.model = snpgdsOpen(reference.fn)
    reference.sample.annot = read.gdsn(index.gdsn(reference.model,'sample.annot'))
    # reference.sign = paste(read.gdsn(index.gdsn(reference.model,'snp.chromosome')),read.gdsn(index.gdsn(reference.model,'snp.position')),sep=":")
    # ref.idx = reference.sign%in%target.sign
    # snp.ref = read.gdsn(index.gdsn(reference.model,'snp.ref'))[ref.idx]
    # snp.alt = read.gdsn(index.gdsn(reference.model,'snp.alt'))[ref.idx⎄]
    
    sample.annot <- rbind(target.sample.annot,reference.sample.annot)
    add.gdsn(genofile,"sample.annot",sample.annot)
    # add.gdsn(genofile,"snp.ref",snp.ref)
    # add.gdsn(genofile,"snp.alt",snp.alt)
    snpgdsClose(genofile)
    snpgdsClose(target.model)
    snpgdsClose(reference.model)
    
    return(TRUE)
  } else 
  {
    
    return(FALSE)
  }
}
