#' Create Reference Model for Ethnicity Analysis
#'
#' This function creates a GDS reference model that can be used to performe EthSEQ ethnicity analysis
#'
#' @param vcf.fn vector of paths to genotype files in VCF format
#' @param annotations data.frame with mapping of all samples names, ethnicities and gender
#' @param out.dir Path to output folder
#' @param model.name Name of the output model
#' @param bed.fn path to BED file with regions of interest
#' @param call.rate SNPs call rate cutoff for inclusion in the final reference model
#' @param cores How many parallel cores to use in the reference model generation
#' @return Logical value indicating the success of the analysis
#' @export


ethseq.RM <- function(
  vcf.fn,
  annotations,
  out.dir = "./",
  model.name = "Reference.Model",
  bed.fn = NA,
  call.rate = 1,
  cores = 1)
{
  if(!file.exists(out.dir))
  {
    .message.Date(paste("Create ",out.dir," folder",sep=""))
    dir.create(out.dir)
  }
  geno.list = list()
  common.snps = c()
  samples = c()
  for(geno in vcf.fn)
  {
    ### Load VCF file
    .message.Date(paste("Load genotype data in VCF file: ",geno,sep=""))
    vcf = fread(geno,data.table = FALSE,skip="#CHROM")
    ### Chromosomes without chr encoding
    vcf[,1] = gsub("chr","",as.character(vcf[,1]))
    vcf = vcf[which(vcf[,1]%in%1:22),]
    
    snp.allele = rep("A/B",nrow(vcf))
    .convertGeno(vcf,snp.allele)

    sample.id = colnames(vcf)[10:ncol(vcf)]
    ### Check names in VCF and annotations
    if(!all(sample.id%in%annotations$sample))
    {
      .message.Date("ERROR: Samples ID in VCF but not in annotations")
      return(FALSE)
    }

    
    ### Only SNPs in BED
    if(!is.na(bed.fn))
    {
      .message.Date("Load BED file")
      bed = fread(bed.fn,data.table = FALSE)
      bed[,1] = gsub("chr","",as.character(bed[,1]))
      .message.Date("Filter VCF file based on BED file coordinates")
      is.in = unlist(mclapply(1:nrow(vcf),function(i)
      {
        tmp = bed[which(bed[,1]==vcf[i,1]),]
        return(any(tmp[,2]<=vcf[i,2]&tmp[,3]>=vcf[i,2]))
      },mc.cores=cores))
      vcf = vcf[which(is.in),]
    }
    
    ### Select only variants with single reference and alternative bases
    .message.Date("Select variants with single reference and alternative bases")
    l1 = sapply(vcf[,4],function(x) nchar(x))
    l2 = sapply(vcf[,5],function(x) nchar(x))
    vcf = vcf[which(l1==1&l2==1),]
    
    ### Filter for call rate
    .message.Date("Filter VCF file based on SNPs call rates")
    cr = apply(vcf[,10:ncol(vcf)],1,function(x) length(which(!x%in%c("./.",".|.")))/length(x))
    vcf = vcf[which(cr>=call.rate),]
    
    
    
    vcf.info = vcf[,1:9]
    vcf = as.matrix(vcf[,10:ncol(vcf)])
    mode(vcf) = "integer"
    snpgdsCreateGeno(file.path(out.dir,paste0(basename(geno),".gds")),
                     genmat = vcf,
                     sample.id = sample.id,
                     snp.id = paste(vcf.info[,1],vcf.info[,2],vcf.info[,3],sep=":"),
                     snp.rs.id = vcf.info[,3],
                     snp.chromosome = vcf.info[,1],
                     snp.position = vcf.info[,2],
                     snp.allele = snp.allele,
                     snpfirstdim=TRUE)
    
    # 
    # 
    # 
    # 
    # 
    # 
    # geno.list[[length(geno.list)+1]] = vcf
    # 
    signature = paste(vcf.info[,1],
                      vcf.info[,2],
                      vcf.info[,3],
                      vcf.info[,4],
                      vcf.info[,5],sep="-")
    if(length(common.snps)==0)
    {
      common.snps = signature
    } else
    {
      common.snps = intersect(common.snps,signature)
    }

    if(any(colnames(vcf)%in%samples))
    {
      .message.Date("ERROR: duplicated samples in VCF files")
      return(FALSE)
    }
    samples = union(samples,colnames(vcf))
  }
  
  ### Merge all vcf files
  if(length(vcf.fn)>1)
  {
    .message.Date("Merge GDS files")
    snpgdsCombineGeno(file.path(out.dir,paste0(basename(vcf.fn),".gds")),file.path(out.dir,paste(model.name,".gds",sep="")), method="exact",snpfirstdim = TRUE)
    # vcf = geno.list[[1]]
    # signature = paste(vcf[,1],vcf[,2],vcf[,3],vcf[,4],vcf[,5],sep="-")
    # vcf = vcf[which(signature%in%common.snps),]
    # vcf = vcf[order(as.numeric(vcf[,1]),as.numeric(vcf[,2])),]
    # 
    # for(i in 2:length(geno.list))
    # {
    #   tmp = geno.list[[i]]
    #   signature = paste(tmp[,1],tmp[,2],tmp[,3],tmp[,4],tmp[,5],sep="-")
    #   tmp = tmp[which(signature%in%common.snps),]
    #   tmp = tmp[order(as.numeric(tmp[,1]),as.numeric(tmp[,2])),]
    #   vcf = cbind(vcf,tmp[,10:ncol(tmp)])
    # }
  } else
  {
    file.rename(file.path(out.dir,paste0(basename(geno),".gds")),
                file.path(out.dir,paste(model.name,".gds",sep="")))
    # vcf = geno.list[[1]]
    # vcf = vcf[order(as.numeric(vcf[,1]),as.numeric(vcf[,2])),]
  }
  
  # .message.Date("Create GDS reference model")
  # geno = t(vcf[,10:ncol(vcf)])
  # res = mclapply(1:ncol(geno),function(i)
  # {
  #   tmp = geno[,i]
  #   tmp[which(tmp=="./.")] = "3"
  #   tmp[which(tmp=="0/1")] = "1"
  #   tmp[which(tmp=="0/0")] = "0"
  #   tmp[which(tmp=="1/1")] = "2"
  #   return(as.numeric(tmp))
  # },mc.cores=cores)
  # 
  # geno = matrix(as.numeric(unlist(res)),nrow=length(res[[1]]),byrow=FALSE)
  # 
  # mafs = 1-apply(geno,2,function(x) (length(which(x==0))*2+length(which(x==1)))/(length(which(x!=3))*2))
  # idx = which(mafs>0.5)
  # for(i in idx)
  # {
  #   tmp = geno[,i]
  #   tmp[which(tmp==0)] = 4
  #   tmp[which(tmp==2)] = 0
  #   tmp[which(tmp==4)] = 2
  #   geno[,i] = tmp
  # }
  # 
  # snp.allele = rep("A/B",ncol(geno))
  # snp.allele[idx] = "B/A"
  # 
  # snpgdsCreateGeno(file.path(out.dir,paste(model.name,".gds",sep="")),
  #                  genmat = geno,
  #                  sample.id = colnames(vcf)[10:ncol(vcf)],
  #                  snp.id = 1:ncol(geno),
  #                  snp.rs.id = vcf[,3],
  #                  snp.chromosome = vcf[,1],
  #                  snp.position = vcf[,2],
  #                  snp.allele = snp.allele,
  #                  snpfirstdim=FALSE)
  # 
  
  genofile <- snpgdsOpen(file.path(out.dir,paste(model.name,".gds",sep="")),readonly = F)
  sample.id = read.gdsn(index.gdsn(genofile,'sample.id'))
  annotations = annotations[which(annotations$sample%in%sample.id),]
  isort = match(sample.id,annotations$sample)
  annotations = annotations[isort,]
  
  sample.annot <- data.frame(pop.group=annotations$pop,sex=annotations$gender)
  add.gdsn(genofile,"sample.annot",sample.annot)
  
  signature.aggregated = paste(read.gdsn(index.gdsn(genofile,'snp.chromosome')),
                               read.gdsn(index.gdsn(genofile,'snp.position')),
                               read.gdsn(index.gdsn(genofile,'snp.rs.id')),sep="-")
  common.snps.ref = sapply(strsplit(common.snps,"-"),'[[',4)
  common.snps.alt = sapply(strsplit(common.snps,"-"),'[[',5)
  common.snps = gsub("-[ACGT]-[ACGT]","",common.snps)
  common.snps.ref = common.snps.ref[which(common.snps%in%signature.aggregated)]
  common.snps.alt = common.snps.alt[which(common.snps%in%signature.aggregated)]
  common.snps = common.snps[which(common.snps%in%signature.aggregated)]
  isort = match(signature.aggregated,common.snps)

  add.gdsn(genofile,"snp.ref",common.snps.ref[isort])
  add.gdsn(genofile,"snp.alt",common.snps.alt[isort])
  # annotations
  snpgdsClose(genofile)
  
  # write.table(vcf[,1:8],paste(out.dir,"/Filtered_SNPs.vcf",sep=""),sep="\t",quote=F,row.names=F)
  return(TRUE)
}
