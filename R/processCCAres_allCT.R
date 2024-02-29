#' @title  processCCAres_allCT
#' @description process CCA res to incluide extra info from all cell-types
#' @param cca_res CCA res
#' @param annot metadata
#' @return cca_res_updated updated cca_res with info
#' @examples
#' processCCAres_allCT(df,annot)
#' @export
################################################################################

processCCAres_allCT<-function(cca_res, x, obj1,y,obj2)
{
  print('processCCAres_allCT:..')

  x.aux.calc = obj1@assays$RNA@scale.data #(pxcell)
  y.aux.calc = obj2@assays$RNA@scale.data
  print(dim(x.aux.calc))
  print(dim(y.aux.calc))
  cca_call= cca_res[['cca']]
  u = cca_call@cell.embeddings[colnames(x),] #colnames(x = obj1)
  v = cca_call@cell.embeddings[colnames(y),] #colnames(y = obj2)
  print(dim(u))
  print(dim(v))
  xscores <- x.aux.calc %*% u# MASS::ginv(sm[[1]]) %*% res$u  pathwaysXcell cellXcc
  yscores <- y.aux.calc %*% v
  cor=diag(cor(xscores,yscores))
  print(cor)
  no_cc_selected = length(cor)
  ul <- t(x.aux.calc) %*% xscores# sampleXpathway . pathwaysXcc == sampleXcc
  vl <- t(y.aux.calc) %*% yscores

  cca_res_updated = list()
  cca_res_updated$u =u
  cca_res_updated$v =v
  cca_res_updated$ul =ul
  cca_res_updated$vl =vl
  cca_res_updated$xscores =xscores
  cca_res_updated$yscores =yscores
  cca_res_updated$cor =cor
  cca_res_updated$x =x
  cca_res_updated$y =y
  cca_res_updated$x.aux.calc = x.aux.calc
  cca_res_updated$y.aux.calc = y.aux.calc
  colnames(cca_res_updated$u)= seq(1:no_cc_selected)
  colnames(cca_res_updated$v)= seq(1:no_cc_selected)
  colnames(cca_res_updated$ul)= seq(1:no_cc_selected)
  colnames(cca_res_updated$vl)= seq(1:no_cc_selected)
  colnames(cca_res_updated$xscores)= seq(1:no_cc_selected)
  colnames(cca_res_updated$yscores)= seq(1:no_cc_selected)

  return(cca_res_updated)
}
