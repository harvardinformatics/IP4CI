#' @title  run_07_rankPathways
#' @description rank pathways
#' @param cca_res CCA res
#' @param p.info pathway info
#' @param w_opt to use weighted pathways scores
#' @return rankP ranked pathways with info
#' @examples
#' run_07_rankPathways(cca_res,p.info,w_opt)
#' @export
################################################################################
run_07_rankPathways <- function(cca_res,p.info,w_opt)
{
  rank.reg= rankPathwaysRegression(cca_res,p.info,w_opt)
  print('rank.reg')
  print(rank.reg)
  rank.cor= rankPathwaysCorrelation(cca_res,p.info,w_opt)
  print('rank.cor')
  print(rank.cor)
  rankP =  merge(rank.reg,rank.cor, by="p_id")
  rownames(rankP)= rankP$Row.names
  rankP$Row.names=NULL

print(rankP)
  return(rankP)
}
