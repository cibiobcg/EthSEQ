
create.Target.Model <-function(sample.names,genotype.dir,out.dir,cores,bam.chr.encoding)
{
  files = list.files(genotype.dir,"genotype.vcf")
  files = files[which(files%in%paste(sample.names,".genotype.vcf",sep=""))]
  geno = fread(paste(genotype.dir,files[1],sep=""),sep="\t",header=T,data.table=FALSE,showProgress = F)
  
  res = mclapply(files,function(f)
  {
    geno = fread(paste(genotype.dir,f,sep=""),sep="\t",header=T,data.table=FALSE,showProgress=F)
    geno = geno[,10]
    geno[which(geno=="./.")] = "3"
    geno[which(geno=="0/1")] = "1"
    geno[which(geno=="0/0")] = "0"
    geno[which(geno=="1/1")] = "2"
    return(as.numeric(geno))
  },mc.cores=cores)
  
  genmat = matrix(unlist(res),ncol = length(res[[1]]), byrow = TRUE)
  rm(res)
  gc()
  mafs = 1-apply(genmat,2,function(x) (length(which(x==0))*2+length(which(x==1)))/(length(which(x!=3))*2))
  idx = which(mafs>0.5)
  for(i in idx)
  {
    tmp = genmat[,i]
    tmp[which(tmp==0)] = 4
    tmp[which(tmp==2)] = 0
    tmp[which(tmp==4)] = 2
    genmat[,i] = tmp
  }
  
  snp.allele = rep("A/B",ncol(genmat))
  snp.allele[idx] = "B/A"
  
  chrs = geno[,1]
  if(bam.chr.encoding)
    chrs = gsub("chr","",as.character(chrs))
  
  snpgdsCreateGeno(file.path(out.dir,"Target.gds"),
                   genmat = genmat,
                   sample.id = gsub(".genotype.vcf","",files),
                   snp.id = 1:nrow(geno),
                   snp.rs.id = geno[,3],
                   snp.chromosome = chrs,
                   snp.position = geno[,2],
                   snp.allele = snp.allele,
                   snpfirstdim=FALSE)
  
  genofile <- snpgdsOpen(file.path(out.dir,"Target.gds"),readonly = F)
  sample.annot <- data.frame(pop.group=rep("ND",nrow(genmat)),sex=rep("M",nrow(genmat)))
  add.gdsn(genofile,"sample.annot",sample.annot)
  
  add.gdsn(genofile,"snp.ref",geno[,4])
  add.gdsn(genofile,"snp.alt",geno[,5])
  
  snpgdsClose(genofile)
}

create.Target.Model.From.VCF <- function(vcf.fn,out.dir,cores)
{
  vcf = fread(vcf.fn,sep="\t",data.table=FALSE,showProgress=FALSE,skip="#CHR")
  ### Chromosomes without chr encoding
  vcf[,1] = gsub("chr","",as.character(vcf[,1]))
  
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
  
  snpgdsCreateGeno(file.path(out.dir,"Target.gds"),
                   genmat = geno,
                   sample.id = colnames(vcf)[10:ncol(vcf)],
                   snp.id = 1:ncol(geno),
                   snp.rs.id = vcf[,3],
                   snp.chromosome = vcf[,1],
                   snp.position = vcf[,2],
                   snp.allele = snp.allele,
                   snpfirstdim=FALSE)
  
  genofile <- snpgdsOpen(file.path(out.dir,"Target.gds"),readonly = F)
  sample.annot <- data.frame(pop.group=rep("ND",nrow(geno)),sex=rep("M",nrow(geno)))
  add.gdsn(genofile,"sample.annot",sample.annot)
  
  add.gdsn(genofile,"snp.ref",vcf[,4])
  add.gdsn(genofile,"snp.alt",vcf[,5])
  
  snpgdsClose(genofile)
}
