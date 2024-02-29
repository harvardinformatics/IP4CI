#' @title  run_07_rankPathways
#' @description rank pathways
#' @param cca_res CCA res
#' @param p.info pathway info
#' @return rankP ranked pathways with info
#' @examples
#' run_07_rankPathways()
#' @export
################################################################################
run_07_rankPathways <- function(cca_res,p.info)
{
  rank.reg= rankPathwaysRegression(cca_res,p.info)
  print('rank.reg')
  print(rank.reg)
  rank.cor= rankPathwaysCorrelation(cca_res,p.info)
  print('rank.cor')
  print(rank.cor)
  rankP =  merge(rank.reg,rank.cor, by="p_id")
  rownames(rankP)= rankP$Row.names
  rankP$Row.names=NULL

print(rankP)
  return(rankP)
}
