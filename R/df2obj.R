#' @title  df2obj
#' @description convert datframe to obj
#' @param df dataframe
#' @param annot metadata
#' @return obj object
#' @examples
#' df2obj(df,annot)
#' @export
################################################################################
df2obj <- function(df,annot)
{
  print('df2obj:...')
  n_pcs = 10
  no_feat=nrow(df)

  obj<- CreateSeuratObject(
    counts=df,
    project = "query",
    assay = "RNA",
    min.cells = 0,
    min.features =0,
    names.field = 1,
    names.delim = "_",
    meta.data = annot)

  obj = Seurat::ScaleData(obj)
  obj = Seurat::FindVariableFeatures(obj,selection.method = getVariableFeatures(obj,min(no_feat,3000)),verbose=F)
  print('DONE FindVariableFeatures')
  obj = Seurat::RunPCA(obj,verbose=F,npcs = n_pcs)
  print('DONE RunPCA')
  obj = Seurat::RunUMAP(obj,dims = 1:n_pcs, n.components = 2,verbose = F)
  print('DONE RunUMAP')

  VarFeat = obj@assays[["RNA"]]@var.features
  print(length(VarFeat))
  Idents(object = obj) <- 'type'
  print(obj)
  #saveRDS(obj, file=paste0(opfname,'_ObjVis.RDS'))
  return(obj)
}
################################################################################
