### Print a message with time info
.message.Date <- function(message)
{
  cat(paste("[",Sys.time() ,"] ",message,"\n",sep=""))
}

### Get operating system
.get.OS <- function()
{
  sysinf <- Sys.info()
  if (!is.null(sysinf))
  {
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else
  {
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  as.character(tolower(os))
}