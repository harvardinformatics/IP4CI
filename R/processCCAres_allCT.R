#' @title  processCCAres_allCT
#' @description process CCA res to incluide extra info from all cell-types
#' @param cca_res CCA res
#' @param x df of obj1
#' @param obj1 obj1
#' @param y df of obj2
#' @param obj2 obj2
#' @return cca_res_updated updated cca_res with info
#' @examples
#' processCCAres_allCT(cca_res, x, obj1,y,obj2)
#' @export
################################################################################

processCCAres_allCT<-function(cca_res, x, obj1,y,obj2)
{
  print('processCCAres_allCT:..')

  cca_call= cca_res[['cca']]
  varF = rownames(cca_call@feature.loadings)
print(dim(obj1@assays$RNA@scale.data))
  # in cca: xscores = p x cc , u= cell x cc
  x.aux.calc = obj1@assays$RNA@scale.data[varF,] #(pxcell)
  y.aux.calc = obj2@assays$RNA@scale.data[varF,]
  print(dim(x.aux.calc))
  print(dim(y.aux.calc))

  u = cca_call@cell.embeddings[colnames(x.aux.calc),] #colnames(x = obj1)
  v = cca_call@cell.embeddings[colnames(y.aux.calc),] #colnames(y = obj2)
  print(dim(u))
  print(dim(v))
  xscores <- x.aux.calc %*% u# MASS::ginv(sm[[1]]) %*% res$u  pathwaysXcell cellXcc
  yscores <- y.aux.calc %*% v
  xscoresW <- apply(xscores,2,function(x) x/sqrt(sum(x*x)))
  yscoresW <- apply(yscores,2,function(x) x/sqrt(sum(x*x)))
  cor=diag(cor(xscores,yscores))
  print(cor)
  corW=diag(cor(xscoresW,yscoresW))
  no_cc_selected = length(cor)
  ul <- t(x.aux.calc) %*% xscores# sampleXpathway . pathwaysXcc == sampleXcc
  vl <- t(y.aux.calc) %*% yscores
  ulW <- t(x.aux.calc) %*% xscoresW# sampleXpathway . pathwaysXcc == sampleXcc
  vlW <- t(y.aux.calc) %*% yscoresW
  # using gene.embed
  ul2 <- t(x.aux.calc) %*% cca_call@feature.loadings #cellxcc
  vl2 <- t(y.aux.calc) %*% cca_call@feature.loadings


  no_cc_selected = length(cor)
  cca_res_updated = list()
  cca_res_updated$pembd = cca_call@feature.loadings
  cca_res_updated$u =u
  cca_res_updated$v =v
  cca_res_updated$ul =ul
  cca_res_updated$vl =vl
  cca_res_updated$ulW =ulW
  cca_res_updated$vlW =vlW
  cca_res_updated$ul2 =ul2
  cca_res_updated$vl2 =vl2
  cca_res_updated$xscores =xscores
  cca_res_updated$yscores =yscores
  cca_res_updated$cor =cor
  cca_res_updated$xscoresW =xscoresW
  cca_res_updated$yscoresW =yscoresW
  cca_res_updated$corW =corW
  cca_res_updated$x =x
  cca_res_updated$y =y
  cca_res_updated$x.aux.calc = x.aux.calc
  cca_res_updated$y.aux.calc = y.aux.calc
  no_cc_selected = length(cor)
  colnames(cca_res_updated$u)= seq(1:no_cc_selected)
  colnames(cca_res_updated$v)= seq(1:no_cc_selected)
  colnames(cca_res_updated$ul)= seq(1:no_cc_selected)
  colnames(cca_res_updated$vl)= seq(1:no_cc_selected)
  colnames(cca_res_updated$ul2)= seq(1:no_cc_selected)
  colnames(cca_res_updated$vl2)= seq(1:no_cc_selected)
  colnames(cca_res_updated$xscores)= seq(1:no_cc_selected)
  colnames(cca_res_updated$yscores)= seq(1:no_cc_selected)
  colnames(cca_res_updated$pembd)= seq(1:no_cc_selected)
  return(cca_res_updated)
}
