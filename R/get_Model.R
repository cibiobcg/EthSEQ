get.Model <- function(model.available,model.folder)
{
  if(!model.available%in%c("SS2","SS4","HALO","NimblegenV3","Universal"))
  {
    message.Date(paste("Model ",model.available," not available",sep=""))
    return(NA)
  }
  if(!file.exists(file.path(model.folder,paste(model.available,".Model.gds",sep=""))))
  {
    download.file(paste("https://demichelislab.unitn.it/lib/exe/fetch.php?media=",model.available,".Model.gds",sep=""),
                  file.path(model.folder,paste(model.available,".Model.gds",sep="")))
  }
  return(file.path(model.folder,paste(model.available,".Model.gds",sep="")))
}