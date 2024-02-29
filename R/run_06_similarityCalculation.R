#' @title  run_06_similarityCalculation
#' @description calculate similarity between two datasets based on CCA info
#' @param merged_res merged_res results either from cca or integration
#' @param annot metadata
#' @param id.list ids for the first and second dataset
#' @param ct rcell-type
#' @param res_dir results directory
#' @examples
#' run_06_similarityCalculation(merged_res,annot,id.list,ct,res_dir)
#' @export

###############################################################################

run_06_similarityCalculation <- function(merged_res,annot,id.list,ct,res_dir)
{
  proj_p_c = rbind((merged_res$u %*% t(merged_res$xscores)),(merged_res$v %*% t(merged_res$yscores))) # cellxp
  calcSim_avgCT(df=t(proj_p_c),annot, id.list,ct ,paste0(res_dir,'1_'))

  proj_p_c = rbind((merged_res$u %*% t(merged_res$pembd)),(merged_res$v %*% t(merged_res$pembd))) # cellxp
  calcSim_avgCT(df=t(proj_p_c),annot, id.list,ct ,paste0(res_dir,'2_'))


  proj_cell_c = rbind(merged_res$ul,merged_res$vl) # cellxcc using calculated p embd for each x/y
  calcSim_avgCT(df=t(proj_cell_c),annot, id.list,ct ,paste0(res_dir,'3_'))

  proj_cell_c = rbind(merged_res$ul2,merged_res$vl2)# cellxcc using fixed p embd from cca
  calcSim_avgCT(df=t(proj_cell_c),annot, id.list,ct ,paste0(res_dir,'4_'))

  proj_cell_c = rbind(merged_res$u ,merged_res$v) # cellxp
  calcSim_avgCT(df=t(proj_cell_c),annot, id.list,ct ,paste0(res_dir,'5_'))

  proj_cell_c = rbind(merged_res$ulW ,merged_res$vlW) # cellxp
  calcSim_avgCT(df=t(proj_cell_c),annot, id.list,ct ,paste0(res_dir,'6_'))
}
