aseq.Run <- function(bam.files,aseq.path,genotype.dir,out.dir,mbq,mrq,mdc,model.path,cores,bam.chr.encoding)
{
  tryCatch(
    {
      ## Create VCF file
      model = snpgdsOpen(model.path,readonly = F)
      snp.list = snpgdsSNPList(model)
      vcf = cbind(snp.list$chromosome,pos=snp.list$position,snp.list$snp.id,
                  as.character(read.gdsn(index.gdsn(model,"snp.ref"))),
                  as.character(read.gdsn(index.gdsn(model,"snp.alt"))),".",".",".")
      colnames(vcf)= c("CHR","POS","ID","REF","ALT","QUAL","FILTER","INFO")
      if(bam.chr.encoding)
        vcf[,1] = paste("chr",vcf[,1],sep="")
      write.table(vcf,file.path(out.dir,"ModelPositions.vcf"),sep="\t",quote=F,row.names=F)
      snpgdsClose(model)
      ## Check ASEQ path or download
      if(get.OS()=="linux")
      {
        aseq.exec = file.path(aseq.path,"ASEQ")
        if(!file.exists(aseq.exec))
        {
          download.file("https://github.com/cibiobcg/EthSEQ_Data/raw/master/ASEQ_binaries/linux64/ASEQ",file.path(aseq.path,"ASEQ"))
          #unzip(file.path(aseq.path,"ASEQ.zip"),exdir=aseq.path)
          Sys.chmod(aseq.exec, mode = "0755", use_umask = TRUE)
        }
        for (b in bam.files)
        {
          message.Date(paste("Computing pileup of BAM file ",b,sep=""))
          command = paste(aseq.exec," vcf=",file.path(out.dir,"ModelPositions.vcf")," bam=",b," mode=GENOTYPE threads=",cores," htperc=0.2 mbq=",mbq,
                          " mrq=",mrq," mdc=",mdc," out=",genotype.dir,sep="")
          system(command,ignore.stderr = T,ignore.stdout = T)
        }
      }
      if(get.OS()=="osx")
      {
        aseq.exec = file.path(aseq.path,"ASEQ")
        if(!file.exists(aseq.exec))
        {
          download.file("https://github.com/cibiobcg/EthSEQ_Data/raw/master/ASEQ_binaries/macosx/ASEQ",file.path(aseq.path,"ASEQ"))
          Sys.chmod(aseq.exec, mode = "0755", use_umask = TRUE)
        }
        for (b in bam.files)
        {
          command = paste(aseq.exec," vcf=",file.path(out.dir,"ModelPositions.vcf")," bam=",b," mode=GENOTYPE threads=",cores," htperc=0.2 mbq=",mbq,
                          " mrq=",mrq," mdc=",mdc," out=",genotype.dir,sep="")
          system(command,ignore.stderr = T,ignore.stdout = T)
        }
      }
      if(get.OS()=="windows")
      {
        aseq.exec = file.path(aseq.path,"ASEQ.exe")
        if(!file.exists(aseq.exec))
        {
          download.file("https://github.com/cibiobcg/EthSEQ_Data/raw/master/ASEQ_binaries/win32/ASEQ.exe",file.path(aseq.path,"ASEQ.exe"))
        }
        for (b in bam.files)
        {
          command = paste(aseq.exec," vcf=",file.path(out.dir,"ModelPositions.vcf")," bam=",b," mode=GENOTYPE threads=",cores," htperc=0.2 mbq=",mbq,
                          " mrq=",mrq," mdc=",mdc," out=",genotype.dir,sep="")
          system(command,ignore.stderr = T,ignore.stdout = T)
        }
      }
    }, error = function(e) {
      message.Date(e)
      return(FALSE)
    })
  return(TRUE)
}
