#' @title  run_06_similarityCalculationIntegration
#' @description calculate similarity between two datasets based on CCA info
#' @param merged_res merged_res results either from cca or integration
#' @param annot metadata
#' @param id.list ids for the first and second dataset
#' @param ct rcell-type
#' @param res_dir results directory
#' @examples
#' run_06_similarityCalculationIntegration(merged_res,annot,id.list,ct,res_dir)
#' @export

###############################################################################

run_06_similarityCalculationIntegration <- function(merged_res,annot,id.list,ct,res_dir)
{
  print('run_06_similarityCalculationIntegration:...')
print(names(merged_res))
    proj_p_c = merged_res$integrated
    print(head(colnames(proj_p_c)))
    print(head(rownames(proj_p_c)))
    print(head((annot)))

    calcSim_avgCT(df=(proj_p_c),annot, id.list,ct ,res_dir)

}
