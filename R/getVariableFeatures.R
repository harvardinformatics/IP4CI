#' @title  getVariableFeatures
#' @description get method of variable features from obj
#' @param sc obj
#' @return name of method used to select fetaures
#' @examples
#' getVariableFeatures(sc)
#' @export
################################################################################
getVariableFeatures = function(sc){
  tryCatch({
    sc = Seurat::FindVariableFeatures(sc,selection.method = 'vst',verbose=F)
    return('vst')
  },error=function(e){
    tryCatch({
      sc = Seurat::FindVariableFeatures(sc,selection.method= 'disp',verbose=F)
      return('disp')
    },error=function(e){
      sc = Seurat::FindVariableFeatures(sc,selection.method= 'mvp',verbose=F)
      return('mvp')
    })
  })
}
