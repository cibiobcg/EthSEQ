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
  geno.list = list()
  common.snps = c()
  samples = c()
  for(geno in vcf.fn)
  {
    ### Load VCF file
    message.Date(paste("Load genotype data in VCF file: ",geno,sep=""))
    vcf = fread(geno,data.table = FALSE)
    if(ncol(vcf)<10)
      return(FALSE)
    
    vcf = vcf[which(vcf[,1]%in%1:22),]
    
    ### Check names in VCF and annotations
    if(!all(colnames(vcf)[10:ncol(vcf)]%in%annotations$sample))
    {
      message.Date("ERROR: Samples ID in VCF but not in annotations")
      return(FALSE)
    }
    
    ### Only SNPs in BED
    if(!is.na(bed.fn))
    {
      message.Date("Load BED file")
      bed = fread(bed.fn,data.table = FALSE)
      message.Date("Filter VCF file based on BED file coordinates")
      is.in = unlist(mclapply(1:nrow(vcf),function(i)
      {
        tmp = bed[which(bed[,1]==vcf[i,1]),]
        return(any(tmp[,2]<=vcf[i,2]&tmp[,3]>=vcf[i,2]))
      },mc.cores=cores))
      vcf = vcf[which(is.in),]
    }
    
    ### Select only variants with single reference and alternative bases
    message.Date("Select variants with single reference and alternative bases")
    l1 = sapply(vcf[,4],function(x) length(strsplit(x,"")[[1]]))
    l2 = sapply(vcf[,5],function(x) length(strsplit(x,"")[[1]]))
    vcf = vcf[which(l1==1&l2==1),]
    
    ### Filter for call rate
    message.Date("Filter VCF file based on SNPs call rates")
    cr = apply(vcf[,10:ncol(vcf)],1,function(x) length(which(!x%in%c("./.",".|.")))/length(x))
    vcf = vcf[which(cr>=call.rate),]
    
    geno.list[[length(geno.list)+1]] = vcf
    
    signature = paste(vcf[,1],vcf[,2],vcf[,3],vcf[,4],vcf[,5],sep="-")
    if(length(common.snps)==0)
    {
      common.snps = signature
    } else
    {
      common.snps = intersect(common.snps,signature)
    }
    
    if(any(colnames(vcf)[10:ncol(vcf)]%in%samples))
    {
      message.Date("ERROR: duplicated samples in VCF files")
      return(FALSE)
    }
    samples = union(samples,colnames(vcf)[10:ncol(vcf)])
  }
  
  ### Merge all vcf files
  if(length(geno.list)>1)
  {
    message.Date("Merge VCF files on intersecting SNPs")
    vcf = geno.list[[1]]
    signature = paste(vcf[,1],vcf[,2],vcf[,3],vcf[,4],vcf[,5],sep="-")
    vcf = vcf[which(signature%in%common.snps),]
    vcf = vcf[order(as.numeric(vcf[,1]),as.numeric(vcf[,2])),]
    
    for(i in 2:length(geno.list))
    {
      tmp = geno.list[[i]]
      signature = paste(tmp[,1],tmp[,2],tmp[,3],tmp[,4],tmp[,5],sep="-")
      tmp = tmp[which(signature%in%common.snps),]
      tmp = tmp[order(as.numeric(tmp[,1]),as.numeric(tmp[,2])),]
      vcf = cbind(vcf,tmp[,10:ncol(vcf)])
    }
  } else
  {
    vcf = geno.list[[1]]
    vcf = vcf[order(as.numeric(vcf[,1]),as.numeric(vcf[,2])),]
  }
  
  annotations = annotations[which(annotations$sample%in%colnames(vcf)[10:ncol(vcf)]),]
  isort = match(colnames(vcf)[10:ncol(vcf)],annotations$sample)
  annotations = annotations[isort,]
  
  message.Date("Create GDS reference model")
  geno = t(vcf[,10:ncol(vcf)])
  res = mclapply(1:ncol(geno),function(i)
  {
    tmp = geno[,i]
    tmp[which(tmp=="./.")] = "3"
    tmp[which(tmp=="0/1")] = "1"
    tmp[which(tmp=="0/0")] = "0"
    tmp[which(tmp=="1/1")] = "2"
    return(as.numeric(tmp))
  },mc.cores=cores)
  
  geno = matrix(as.numeric(unlist(res)),nrow=length(res[[1]]),byrow=FALSE)
  
  mafs = 1-apply(geno,2,function(x) (length(which(x==0))*2+length(which(x==1)))/(length(which(x!=3))*2))
  idx = which(mafs>0.5)
  for(i in idx)
  {
    tmp = geno[,i]
    tmp[which(tmp==0)] = 4
    tmp[which(tmp==2)] = 0
    tmp[which(tmp==4)] = 2
    geno[,i] = tmp
  }
  
  snp.allele = rep("A/B",ncol(geno))
  snp.allele[idx] = "B/A"
  
  snpgdsCreateGeno(file.path(out.dir,paste(model.name,".gds",sep="")),
                   genmat = geno,
                   sample.id = colnames(vcf)[10:ncol(vcf)],
                   snp.id = 1:ncol(geno),
                   snp.rs.id = vcf[,3],
                   snp.chromosome = vcf[,1],
                   snp.position = vcf[,2],
                   snp.allele = snp.allele,
                   snpfirstdim=FALSE)
  
  genofile <- snpgdsOpen(file.path(out.dir,paste(model.name,".gds",sep="")),readonly = F)
  sample.annot <- data.frame(pop.group=annotations$pop,sex=rep(annotations$gender,nrow(geno)))
  add.gdsn(genofile,"sample.annot",sample.annot)
  
  add.gdsn(genofile,"snp.ref",vcf[,4])
  add.gdsn(genofile,"snp.alt",vcf[,5])
  
  snpgdsClose(genofile)
  
  write.table(vcf,paste(out.dir,"/Filtered_",basename(vcf.fn),sep=""),
              sep="\t",quote=F,row.names=F)
  
  return(TRUE)
}
  
  
  
  
  
  
  
  
  