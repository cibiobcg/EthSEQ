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
    reference.model = snpgdsOpen(reference.fn)
    reference.sample.annot = read.gdsn(index.gdsn(reference.model,'sample.annot'))
    sample.annot <- rbind(target.sample.annot,reference.sample.annot)
    add.gdsn(genofile,"sample.annot",sample.annot)
    snpgdsClose(genofile)
    snpgdsClose(target.model)
    snpgdsClose(reference.model)
    
    return(TRUE)
  } else 
  {
    return(FALSE)
  }
}
