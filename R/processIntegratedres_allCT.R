#' @title  processIntegratedres_allCT
#' @description process integration res to incluide extra info from all cell-types
#' @param integrated_res CCA res
#' @param x metadata
#' @return integrated_res_updated updated integrated_res with info
#' @examples
#' processIntegratedres_allCT(df,annot)
#' @export
################################################################################

processIntegratedres_allCT<-function(integrated_res, x, obj1,y,obj2)
{
  print('processIntegratedres_allCT:..')

  varF = integrated_res@assays$integrated@var.features
  x.aux.calc = obj1@assays$RNA@scale.data[varF,] #(pxcell)
  y.aux.calc = obj2@assays$RNA@scale.data[varF,]
  print(dim(x.aux.calc))
  print(dim(y.aux.calc))
  cell.embeddings= Embeddings(integrated_res) # cells X pcs
  u = cell.embeddings[colnames(x),] #colnames(x = obj1)
  v = cell.embeddings[colnames(y),] #colnames(y = obj2)
  print(dim(u))
  print(dim(v))
  xscores <- x.aux.calc %*% u# MASS::ginv(sm[[1]]) %*% res$u  pathwaysXcell cellXcc
  yscores <- y.aux.calc %*% v
  cor=diag(cor(xscores,yscores))
  print(cor)
  no_cc_selected = length(cor)
  ul <- t(x.aux.calc) %*% xscores# sampleXpathway . pathwaysXcc == sampleXcc
  vl <- t(y.aux.calc) %*% yscores

  integrated_res_updated = list()
  integrated_res_updated$u =u
  integrated_res_updated$v =v
  integrated_res_updated$ul =ul
  integrated_res_updated$vl =vl
  integrated_res_updated$xscores =xscores
  integrated_res_updated$yscores =yscores
  integrated_res_updated$cor =cor
  integrated_res_updated$x =x
  integrated_res_updated$y =y
  integrated_res_updated$x.aux.calc = x.aux.calc
  integrated_res_updated$y.aux.calc = y.aux.calc
  integrated_res_updated$integrated = integrated_res[['integrated']]@data
  integrated_res_updated$meta.data = integrated_res@meta.data
  colnames(integrated_res_updated$u)= seq(1:no_cc_selected)
  colnames(integrated_res_updated$v)= seq(1:no_cc_selected)
  colnames(integrated_res_updated$ul)= seq(1:no_cc_selected)
  colnames(integrated_res_updated$vl)= seq(1:no_cc_selected)
  colnames(integrated_res_updated$xscores)= seq(1:no_cc_selected)
  colnames(integrated_res_updated$yscores)= seq(1:no_cc_selected)

  return(integrated_res_updated)
}
