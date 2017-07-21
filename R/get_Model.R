get.Model <- function(model.available,model.folder)
{
  if(!model.available%in%c("SS2.All","SS2.Major","SS4.All","SS4.Major","HALO.All","HALO.Major","NimblegenV3.All","NimblegenV3.Major","Exonic.All","Exonic.Major"))
  {
    message.Date(paste("Model ",model.available," not available",sep=""))
    return(NA)
  }
  if(!file.exists(file.path(model.folder,paste(model.available,".Model.gds",sep=""))))
  {
    download.file(paste("https://github.com/aromanel/EthSEQ_Data/raw/master/EthSEQ_Models/",model.available,".Model.gds",sep=""),
                  file.path(model.folder,paste(model.available,".Model.gds",sep="")))
  }
  return(file.path(model.folder,paste(model.available,".Model.gds",sep="")))
}
