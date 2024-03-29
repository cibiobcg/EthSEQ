.create.Target.Model <-function(sample.names,genotype.dir,out.dir,cores,bam.chr.encoding)
{
  files = list.files(genotype.dir,"genotype.vcf")
  files = files[which(files%in%paste(sample.names,".genotype.vcf",sep=""))]
  geno = fread(paste(genotype.dir,files[1],sep=""),sep="\t",header=T,data.table=FALSE,showProgress = F,skip="#CHROM")
  
  res = mclapply(files,function(f)
  {
    geno = fread(paste(genotype.dir,f,sep=""),sep="\t",header=T,data.table=TRUE,showProgress=F,skip="#CHROM")
    geno[,10]
  },mc.cores=cores)
  
  vcf = cbind(geno[,1:9],data.frame(do.call(cbind,res)))
  
  vcf[,1] = gsub("chr","",as.character(vcf[,1]))
  snp.allele = rep("A/B",nrow(vcf))
  
  .convertGeno(vcf,snp.allele)

  vcf.info = vcf[,1:9]
  sample.id = colnames(vcf)[10:ncol(vcf)]
  vcf = as.matrix(vcf[,10:ncol(vcf)])
  mode(vcf) = "integer"
  
  snpgdsCreateGeno(file.path(out.dir,"Target.gds"),
                   genmat = vcf,
                   sample.id = sample.id,
                   snp.id = paste(vcf.info[,1],vcf.info[,2],vcf.info[,3],sep=":"),
                   snp.rs.id = vcf.info[,3],
                   snp.chromosome = vcf.info[,1],
                   snp.position = vcf.info[,2],
                   snp.allele = snp.allele,
                   snpfirstdim=TRUE)
  
  genofile <- snpgdsOpen(file.path(out.dir,"Target.gds"),readonly = F)
  sample.annot <- data.frame(pop.group=rep("ND",ncol(vcf)),sex=rep("M",ncol(vcf)))
  add.gdsn(genofile,"sample.annot",sample.annot)
  
  add.gdsn(genofile,"snp.ref",vcf.info[,4])
  add.gdsn(genofile,"snp.alt",vcf.info[,5])
  
  snpgdsClose(genofile)

}

.create.Target.Model.From.VCF <- function(vcf.fn,out.dir,cores)
{
  vcf = fread(vcf.fn,sep="\t",data.table=FALSE,showProgress=FALSE,skip="#CHROM")
  
  ### Chromosomes without chr encoding
  vcf[,1] = gsub("chr","",as.character(vcf[,1]))
  
  snp.allele = rep("A/B",nrow(vcf))
  .convertGeno(vcf,snp.allele)
  
  vcf.info = vcf[,1:9]
  sample.id = colnames(vcf)[10:ncol(vcf)]
  vcf = as.matrix(vcf[,10:ncol(vcf)])
  mode(vcf) = "integer"
  
  snpgdsCreateGeno(file.path(out.dir,"Target.gds"),
                   genmat = vcf,
                   sample.id = sample.id,
                   snp.id = paste(vcf.info[,1],vcf.info[,2],vcf.info[,3],sep=":"),
                   snp.rs.id = vcf.info[,3],
                   snp.chromosome = vcf.info[,1],
                   snp.position = vcf.info[,2],
                   snp.allele = snp.allele,
                   snpfirstdim=TRUE)
  
  genofile <- snpgdsOpen(file.path(out.dir,"Target.gds"),readonly = F)
  sample.annot <- data.frame(pop.group=rep("ND",ncol(vcf)),sex=rep("M",ncol(vcf)))
  add.gdsn(genofile,"sample.annot",sample.annot)
  
  add.gdsn(genofile,"snp.ref",vcf.info[,4])
  add.gdsn(genofile,"snp.alt",vcf.info[,5])
  
  snpgdsClose(genofile)
}
