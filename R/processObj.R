#' @title  processObj
#' @description process obj using Seurat workflow
#' @param obj obj data
#' @return obj processed obj data
#' @examples
#' obj=x
#' processObj(obj)
#'@export
################################################################################
processObj <- function(obj)
{
  print('processObj:...')


  # obj <- FindNeighbors(obj, dims = 1:dims)
  # obj<- FindClusters(obj, resolution = resolution)

  obj <- NormalizeData(obj,verbose=F)
  print('DONE Norm')

  obj <- FindVariableFeatures(obj,selection.method = "vst",nfeatures = min(3000,rownames(obj)),verbose=F) # or 2000 or 5000
  print('DONE VF')

  obj <- ScaleData(obj,verbose=FALSE) #features=rownames(obj) , default is variable genes
  print('DONE Scale')


  if(ncol(obj) < 50)
    obj <- RunPCA(obj, npcs=30, features = VariableFeatures(object = obj))
  else
    obj <- RunPCA(obj, features = VariableFeatures(object = obj))
  print('DONE PCA')

  obj <- RunUMAP(obj, dims = 1:30)

  obj <- UpdateSeuratObject(obj)

  return(obj)
}
