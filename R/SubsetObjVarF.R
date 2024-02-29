#' @title  df2objVarF
#' @description convert datframe to obj with selected variable features
#' @param obj object
#' @param VarF list variable features
#' @param cells list cells
#' @return obj object
#' @examples
#' df2objVarF(df,annot,VarF)
#' @export
################################################################################
SubsetObjVarF <- function(obj,VarF,cells)
{
  print('SubsetObjVarF:...')
  n_pcs = 10

  VarFeat = obj@assays[["RNA"]]@var.features
  print(length(VarFeat))
  Idents(object = obj) <- 'type'
  print(obj)
  #saveRDS(obj, file=paste0(opfname,'_ObjVis.RDS'))
  return(obj)
}
################################################################################
