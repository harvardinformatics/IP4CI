#' @title  df2obj
#' @description convert datframe to obj
#' @param df dataframe
#' @param annot metadata
#' @param nVarF number variable features
#' @return obj object
#' @examples
#' df2obj(df,annot)
#' @export
################################################################################
df2obj <- function(df,annot,nVarF)
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
print(obj)
  obj = Seurat::ScaleData(obj)
  if(nVarF != 0) # all no featVar
  {
    obj = Seurat::FindVariableFeatures(obj,selection.method = getVariableFeatures(obj),nfeatures =min(no_feat,nVarF) ,verbose=F)
    print('DONE FindVariableFeatures')
    obj = Seurat::RunPCA(obj,verbose=F,npcs = n_pcs)
    print('DONE RunPCA')
  }
  if(nVarF == 0) # do featVar
  {
    obj = Seurat::RunPCA(obj,verbose=F,npcs = n_pcs,features = rownames(df))
    print('DONE RunPCA')
  }

  obj = Seurat::RunUMAP(obj,dims = 1:n_pcs, n.components = 2,verbose = F)
  print('DONE RunUMAP')

  # VarFeat = obj@assays[["RNA"]]@var.features
  # print(length(VarFeat))
  Idents(object = obj) <- 'type'
  print(obj)
  #saveRDS(obj, file=paste0(opfname,'_ObjVis.RDS'))
  return(obj)
}
################################################################################
