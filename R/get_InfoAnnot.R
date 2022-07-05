#' List the models available
#'
#' This function prints the list of all available models.
#'
#' @param assembly Human genome assembly to use: hg19 or hg38
#' @param pop population to use: All, AFR, AMR, EAS, EUR or SAS
#' @return data.frame of all available models to use with specified assembly and population
#' @export

getModelsList <- function()
{
  exKit = system.file("extdata", "exonic_kits_map.tsv",
                      package="EthSEQ")
  df = fread(exKit,data.table = F)
  df = df[sapply(strsplit(df$assembly,","),function(obj){any(obj==assembly)})&sapply(strsplit(df$pop,","),function(obj){any(obj==pop)}),]
  return(df)
}

#' List the samples annotation
#'
#' This function prints the list of all available samples used to build the reference models.
#'
#' @export
getSamplesInfo <- function()
{
  exKit = system.file("extdata", "samplesPop_map.tsv",
                      package="EthSEQ")
  fread(exKit)
}
