#' Perform Pileup with ASEQ
#'
#' This function performs pileup of set of BAM files using ASEQ on reference model positions.
#'
#' @param bam.files Vector of BAM files paths
#' @param aseq.path Path to ASEQ binary folder
#' @param genotype.dir Path to genotype output folder
#' @param out.dir Path to analysis output folder
#' @param mbq Minimum base quality for ASEQ pileup
#' @param mrq Minimum read quality for ASEQ pileup
#' @param mdc Minimum read count to call genotype
#' @param model.path Path to reference model GDS file
#' @param cores Number of cores used in the analysis
#' @return Logical value indicating the success of the analysis
aseq.Run <- function(bam.files,aseq.path,genotype.dir,out.dir,mbq,mrq,mdc,model.path,cores)
{
  tryCatch(
    {
      ## Create VCF file
      model = snpgdsOpen(model.path,readonly = F)
      snp.list = snpgdsSNPList(model)
      vcf = cbind(snp.list$chromosome,pos=snp.list$position,snp.list$rs.id,
                  as.character(read.gdsn(index.gdsn(model,"snp.ref"))),
                  as.character(read.gdsn(index.gdsn(model,"snp.alt"))),".",".",".")
      colnames(vcf)= c("CHR","POS","ID","REF","ALT","QUAL","FILTER","INFO")
      write.table(vcf,file.path(out.dir,"ModelPositions.vcf"),sep="\t",quote=F,row.names=F)
      snpgdsClose(model)
      ## Check ASEQ path or download
      if(get.OS()=="linux")
      {
        aseq.exec = file.path(aseq.path,"ASEQ")
        if(!file.exists(aseq.exec))
        {
          download.file("https://github.com/aromanel/EthSEQ_Data/raw/master/ASEQ_binaries/linux64/ASEQ",file.path(aseq.path,"ASEQ"))
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
          download.file("https://github.com/aromanel/EthSEQ_Data/raw/master/ASEQ_binaries/macosx/ASEQ",file.path(aseq.path,"ASEQ"))
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
          download.file("https://github.com/aromanel/EthSEQ_Data/raw/master/ASEQ_binaries/win32/ASEQ",file.path(aseq.path,"ASEQ.exe"))
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
      retunr(FALSE)
    })
  return(TRUE)
}
