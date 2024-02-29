#' @title  getVariableFeatures
#' @description get method of variable features from obj
#' @param sc_obj obj
#' @param no_feat number of features
#' @return name of method used to select fetaures
#' @examples
#' getVariableFeatures(sc_obj,no_feat)
#' @export
################################################################################
getVariableFeatures = function(sc,no_feat){
  tryCatch({
    sc = Seurat::FindVariableFeatures(sc,selection.method = 'vst',verbose=F,nfeatures=no_feat)
    return('vst')
  },error=function(e){
    tryCatch({
      sc = Seurat::FindVariableFeatures(sc,selection.method= 'disp',verbose=F,nfeatures=no_feat)
      return('disp')
    },error=function(e){
      sc = Seurat::FindVariableFeatures(sc,selection.method= 'mvp',verbose=F)
      return('mvp')
    })
  })
}
