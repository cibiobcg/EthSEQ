#' Ethnicity analysis from whole-exome sequencing data
#'
#' This function performs ethnicity analysis of a set of samples ad reports the results.
#'
#' @param target.vcf Path to the sample's genotypes in VCF format
#' @param target.gds Path to the sample's genotypes in GDS format
#' @param bam.list Path to a file containing a list of BAM files paths
#' @param out.dir Path to the folder where the output of the analysis is saved
#' @param model.gds Path to a GDS file specifying the reference model
#' @param model.available String specifying the pre-computed reference model to use (SS2,SS4,NimblegenV3,HALO,Exonic)
#' @param model.folder Path to the folder where reference models are already present of downloaded when needed
#' @param run.genotype Logical values indicating wheter the ASEQ genotype should be run
#' @param aseq.path Path to the folder where ASEQ binary is available of is downloaded when needed
#' @param mbq Minmum base quality used in the pileup by ASEQ
#' @param mrq Minimum read quality used in the piluep by ASEQ
#' @param mdc Minimum read count accettable for genotype inference by ASEQ
#' @param cores Number of parallel cores used for the analysis
#' @param verbose Print detaild information
#' @param composite.model.call.rate SNP call rate used to run Principal Component Analysis (PCA)
#' @param refinement.analysis Matrix specyfing a tree of ethnicities
#' @param space Dimensions of PCA space used to infer ethnicity (2D or 3D)
#' @param bam.chr.encoding Logical value indicating whether input BAM files have chromosomes encoded with "chr" prefix
#' @return Logical value indicating the success of the analysis
#' @export
ethseq.Analysis <- function(
  target.vcf = NA,
  target.gds = NA,
  bam.list = NA,
  out.dir = "/tmp",
  model.gds = NA,
  model.available = NA,
  model.folder = "/tmp",
  run.genotype = FALSE,
  aseq.path = "/tmp",
  mbq=20,
  mrq=20,
  mdc=10,
  cores=1,
  verbose=TRUE,
  composite.model.call.rate = 1,
  refinement.analysis = NA,
  space = "2D",
  bam.chr.encoding = FALSE)
{
  
  ## Version
  message.Date("Running EthSEQ")
  message.Date(paste("Working directory: ",out.dir,sep=""))
  
  if(!file.exists(out.dir))
  {
    message.Date(paste("Create ",out.dir," folder",sep=""))
    dir.create(out.dir)
  }
  
  ## Select Model when model.gds is not specified
  if(is.na(model.gds))
  {
    if(is.na(model.available))
    {
      message.Date("ERROR: model.gds or model.available variables should be specified.")
      return(FALSE)
    }
    model.path = get.Model(model.available,model.folder)
  } else
  {
    model.path = model.gds
  }
  
  # Run ASEQ when target VCF and GDS are not provided
  if(is.na(target.vcf)&is.na(target.gds))
  {
    if(is.na(bam.list))
    {
      message.Date("ERROR: target.vcf or bam.list variables should be specified.")
      return(FALSE)
    }
    
    ### Check if all BAM files exists
    if(run.genotype&any(!file.exists(unique(readLines(bam.list)))))
    {
      message.Date("ERROR: one or more BAM files do not exist.")
      return(FALSE)
    }
    bam.files = unique(readLines(bam.list))
    sample.names = gsub(".bam","",basename(bam.files))
    
    ## Execute ASEQ to generate genotypes
    table = data.frame(sample=gsub(".bam","",basename(bam.files)),bam=bam.files)
    genotype.dir = file.path(out.dir,"ASEQGenotypes","")
    if(!file.exists(genotype.dir))
    {
      message.Date(paste("Create ",genotype.dir," folder",sep=""))
      dir.create(genotype.dir)
    }
    
    message.Date(paste("Run ASEQ to generate genotypes for ",length(bam.files)," samples",sep=""))
    if(run.genotype)
    {
      res = aseq.Run(bam.files,aseq.path,genotype.dir,out.dir,mbq,mrq,mdc,model.path,cores,bam.chr.encoding)
      if(!res)
      {
        message.Date("ERROR: Error while executing ASEQ.")
        return(FALSE)
      }
    }
    
    ## Create Target Model GDS file
    message.Date("Create target model")
    create.Target.Model(sample.names,genotype.dir,out.dir,cores,bam.chr.encoding)
    
    target.model = file.path(out.dir,"Target.gds")
  } else if(is.na(target.gds))
  {
    message.Date("Create target model from VCF")
    create.Target.Model.From.VCF(target.vcf,out.dir,cores)

    target.model = file.path(out.dir,"Target.gds")
  } else
  {
    target.model = target.gds
  }
  
  ### Create Composite model - DD Changed with SNPRelate function
  message.Date("Create aggregated model")
  res = combine.Models(model.path,target.model,out.dir,composite.model.call.rate)
  # snpgdsCombineGeno(c(target.model,model.path),file.path(out.dir,"Aggregated.gds"),
  #                   method = 'position',snpfirstdim = TRUE)
  # res = file.exists(file.path(out.dir,"Aggregated.gds"))
  # genofile <- snpgdsOpen(file.path(out.dir,"Aggregated.gds"),readonly = F)
  # target.sample.annot = read.gdsn(index.gdsn(snpgdsOpen(target.model),'sample.annot'))
  # reference.sample.annot = read.gdsn(index.gdsn(snpgdsOpen(model.path),'sample.annot'))
  # sample.annot <- rbind(target.sample.annot,reference.sample.annot)
  # add.gdsn(genofile,"sample.annot",sample.annot)
  # snpgdsClose(genofile)
  
  if(!res)
  {
    message.Date("ERROR: Target and reference models are not compatible.")
    return(FALSE)
  }
  
  ### Perform PCA
  message.Date("Perform PCA on aggregated model")

  model <- snpgdsOpen(file.path(out.dir,"Aggregated.gds"),readonly = F)
  pca <- snpgdsPCA(model,num.thread = cores,eigen.method = "DSPEVX",verbose = verbose,missing.rate = 1-composite.model.call.rate,eigen.cnt=5)
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
  snpgdsClose(model)
  
  ### Infer ethnicity
  message.Date("Infer ethnicities")
  annotations = get.Ethnicity(tab,space)
  
  if(!all(is.na(refinement.analysis)))
  {
    res = check.Matrix(refinement.analysis,names(annotations[[2]]))
    if(!res)
    {
      message.Date("ERROR: Refinement analysis matrix is wrongly formatted.")
      return(FALSE)
    }
    refinement.index = 1
    refinement.leafs = get.Position.Leafs(refinement.analysis)
    refinement.position = get.Position.Matrix(refinement.analysis)
    refinement.position = unlist(apply(refinement.position,1,function(x) x[which(x!="")]))
    refinement.position = as.numeric(refinement.position)
    refinement.analysis = unlist(apply(refinement.analysis,1,function(x) x[which(x!="")]))
    refinements = list()
    refinements[[1]] = annotations
    out = refinements[[1]][[1]][,c(1,2,8)]
    
    ### Plot visual PDF report
    message.Date("Plot visual report")
    coord = annotations[[1]]
    n.dim = as.numeric(gsub("D","",space))-1
    coord = coord[,c(1,3:(3+n.dim))]
    coord[,1] = gsub("target.","",coord[,1])
    write.table(coord,file.path(out.dir,"Report_Ref0.PCAcoord"),sep="\t",row.names=F,quote=F)
    plot.Report(tab,annotations[[1]],length(pca$snp.id),annotations[[2]],annotations[[3]],out.dir,label="_Ref0",space = space)
    
    for(s in 1:length(refinement.analysis))
    {
      message.Date(paste("Refinement step ",s," on populations (",paste(refinement.analysis[[s]],collapse=","),")",sep=""))
      model <- snpgdsOpen(file.path(out.dir,"Aggregated.gds"),readonly = F)
      samples <- refinements[[refinement.position[s]]][[1]]
      samples <- samples$sample.id[which(sapply(samples$pop,function(x) length(intersect(strsplit(x,"\\|")[[1]],strsplit(refinement.analysis[[s]],"\\|")[[1]]))>0))]
      sample.id <- read.gdsn(index.gdsn(model, "sample.id"))
      pop_code <- read.gdsn(index.gdsn(model, "sample.annot/pop.group"))
      samples = c(samples,sample.id[which(as.character(pop_code)%in%strsplit(refinement.analysis[[s]],"\\|")[[1]])])
      pca <- snpgdsPCA(model,num.thread = cores,eigen.method = "DSPEVX",verbose = verbose,sample.id = samples, missing.rate=1-composite.model.call.rate,eigen.cnt =5)
      
      idx = which(as.character(pop_code)%in%c("ND",strsplit(refinement.analysis[[s]],"\\|")[[1]]))
      tab <- data.frame(sample.id = pca$sample.id,
                        pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                        EV1 = pca$eigenvect[,1],
                        EV2 = pca$eigenvect[,2],
                        EV3 = pca$eigenvect[,3],
                        EV4 = pca$eigenvect[,4],
                        EV5 = pca$eigenvect[,5],
                        stringsAsFactors = FALSE)
      snpgdsClose(model)
      
      annotations = get.Ethnicity(tab,space)
      if(annotations[[4]])
      {
        coord = annotations[[1]]
        n.dim = as.numeric(gsub("D","",space))-1
        coord = coord[,c(1,3:(3+n.dim))]
        coord[,1] = gsub("target.","",coord[,1])
        write.table(coord,file.path(out.dir,paste("Report_Ref",s,".PCAcoord",sep="")),sep="\t",row.names=F,quote=F)
        plot.Report(tab,annotations[[1]],length(pca$snp.id),annotations[[2]],annotations[[3]],out.dir,label=paste("_Ref",s,sep=""),space = space)
      }
      if(!refinement.leafs[s])
      {
        refinements[[refinement.index+1]] = annotations
        refinement.index = refinement.index+1
      }
      
      if(annotations[[4]])
      {
        for(i in 1:nrow(out))
        {
          id = which(annotations[[1]]$sample.id==out$sample.id[i])
          if(length(id)>0)
          {
            out$pop[i] = annotations[[1]]$pop[id]
            out$type[i] = annotations[[1]]$type[id]
          }
        }
      }
    }
    
    message.Date("Print annotations")
    out[,1] = gsub("target.","",out[,1])
    write.table(out,file.path(out.dir,"Report.txt"),sep="\t",row.names=F,quote=F)
    
  } else
  {
    ### Print ethnicity annotations on tex tab-delimited file
    message.Date("Print annotations")
    out = annotations[[1]]
    out = out[,c(1,2,8,9)]
    out[,1] = gsub("target.","",out[,1])
    write.table(out,file.path(out.dir,"Report.txt"),sep="\t",row.names=F,quote=F)
    out = annotations[[1]]
    n.dim = as.numeric(gsub("D","",space))-1
    out = out[,c(1,3:(3+n.dim))]
    out[,1] = gsub("target.","",out[,1])
    write.table(out,file.path(out.dir,"Report.PCAcoord"),sep="\t",row.names=F,quote=F)
    
    ### Plot visual PDF report
    message.Date("Plot visual report")
    plot.Report(tab,annotations[[1]],length(pca$snp.id),annotations[[2]],annotations[[3]],out.dir,space=space)
  }
  
  message.Date("Computation end")
  return(TRUE)
}
