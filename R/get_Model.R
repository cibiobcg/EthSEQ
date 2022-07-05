.get.Model <- function(model.available, model.folder, assembly, pop)
{
  list.models = getModels(assembly, pop)
  if(is.data.frame(list.models)) {
    if(!any(list.models$name==model.available))
    {
      return(NA)
    }
    model.name = paste(model.available,".",assembly,".",pop,".Model.gds",sep="")
    if(!file.exists(file.path(model.folder,model.name)))
    {
      download.file(paste("https://github.com/cibiobcg/EthSEQ_Data/raw/master/EthSEQ_Models_v3/",model.name,sep=""),
                    file.path(model.folder,model.name),mode='wb')
    }
    return(file.path(model.folder,model.name))
  }
  else {
    return(NA)
  }
}
