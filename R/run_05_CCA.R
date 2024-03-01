#' @title  run_05_CCA
#' @description run Seurat CCA over two objectes
#' @param obj1 object of dataset 1
#' @param obj2 object of dataset 2
#' @param num.cc number of cc dimensions
#' @return ccaRes CCA results
#' @examples
#' run_05_CCA(obj1,obj2,num.cc)
#' @export
################################################################################
run_05_CCA <- function(obj1,obj2,num.cc)
{
  print('run_05_CCA:...')
  cca_res = RunCCA(obj1, obj2, num.cc = num.cc   ,renormalize = FALSE,rescale = FALSE,verbose = TRUE) #rownames(obj1@assays$RNA@counts)
 return(cca_res)

}
