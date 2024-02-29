#' @title  RnamingCTconv
#' @description Rnaming cell-type labels convienent for R
#' @param obj obj data
#' @return obj updated obj data
#' @examples
#' obj=x
#' RnamingCTconv(obj)
#' @export


RnamingCTconv <-function(obj)
{
  obj@meta.data$type<- gsub("/", ".",obj@meta.data$type)
  obj@meta.data$type <- gsub(" ", "_",obj@meta.data$type)
  obj@meta.data$type <- gsub("-", ".",obj@meta.data$type)
  print( unique(obj@meta.data$type))
  return(obj)
}
