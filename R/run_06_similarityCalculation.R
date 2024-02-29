#' @title  run_06_similarityCalculation
#' @description calculate similarity between two datasets based on CCA info
#' @param cca_res cca results
#' @param annot metadata
#' @param id.list ids for the first and second dataset
#' @param ct rcell-type
#' @param res_dir results directory
#' @examples
#' run_06_similarityCalculation(cca_res,annot,id.list,ct,res_dir)
#' @export

###############################################################################

run_06_similarityCalculation <- function(cca_res,annot,id.list,ct,res_dir)
{

  proj_p_c = rbind((cca_res$u %*% t(cca_res$xscores)),(cca_res$v %*% t(cca_res$yscores)))
  calcSim_avgCT(df=t(proj_p_c),annot, id.list, ct ,res_dir)

}
