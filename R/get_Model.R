.get.Model <- function(model.available,model.folder, assembly)
{
  list.models = getModels()
  if(!any(list.models$name==model.available&list.models$assembly==assembly))
  {
    .message.Date(paste("Model ",model.available," (using assembly ",assembly,") not available",sep=""))
    return(NA)
  }
  if(!file.exists(file.path(model.folder,paste(model.available,".Model.gds",sep=""))))
  {
    download.file(paste("https://github.com/cibiobcg/EthSEQ_Data/raw/master/EthSEQ_Models_v3/",assembly,".",model.available,".Model.gds",sep=""),
                  file.path(model.folder,paste(model.available,".Model.gds",sep="")),mode='wb')
  }
  return(file.path(model.folder,paste(assembly,".",model.available,".Model.gds",sep="")))
}
