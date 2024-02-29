#' @title  run_05_CCA
#' @description run Seurat CCA over two objectes
#' @param obj1 object of dataset 1
#' @param obj2 object of dataset 2
#' @return CCA.obj CCA results
#' @examples
#' run_05_CCA(obj1,obj2)
#' @export
################################################################################
run_05_CCA <- function(obj1,obj2)
{
  print('run_05_CCA:...')
  cc_dims = 30
  cca_res = RunCCA(obj1, obj2, num.cc = cc_dims,features =rownames(x = obj1)   ,renormalize = FALSE,rescale = FALSE,verbose = TRUE) #rownames(obj1@assays$RNA@counts)
  return(cca_res)

}
