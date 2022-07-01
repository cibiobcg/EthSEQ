#' List the models available
#'
#' This function prints the list of all available models.
#'
#' @param assembly Human genome assembly to use: hg19 or hg38
#' @param pop population to use: All, AFR, AMR, EAS, EUR or SAS
#' @return data.frame of all available models to use with specified assembly and population
#' @export

getModels <- function(assembly='hg38',pop='All')
{
  if(assembly%in%c('hg19','hg38')&pop%in%c('All','AFR','AMR','EAS','EUR','SAS')) {
    exKit = system.file("extdata", "exonic_kits_map.tsv",
                        package="EthSEQ")
    df = fread(exKit,data.table = F)
    df = df[sapply(strsplit(df$assembly,","),function(obj){any(obj==assembly)})&sapply(strsplit(df$pop,","),function(obj){any(obj==pop)}),1:2]
    return(df)
  } else {
    .message.Date(paste0("ERROR: No available models using assembly: ",assembly," and/or population: ",pop))
    return(FALSE)
  }
}

#' List the samples annotation
#'
#' This function prints the list of all available samples used to build the reference models.
#'
#' @export
getPops <- function()
{
  exKit = system.file("extdata", "samplesPop_map.tsv",
                      package="EthSEQ")
  fread(exKit)
}
