#' List the models available
#'
#' This function prints the list of all available models.
#'
#' @return data.frame of all available models to use with specified assembly and population
#' @export

getModelsList <- function()
{
  exKit = system.file("extdata", "exonic_kits_map.tsv",
                      package="EthSEQ")
  df = fread(exKit,data.table = F)
  return(df)
}

#' List the samples annotation
#'
#' This function prints the list of 1,000 Genomes Project samples used to build the reference models.
#'
#' @export
getSamplesInfo <- function()
{
  exKit = system.file("extdata", "samplesPop_map.tsv",
                      package="EthSEQ")
  fread(exKit)
}
