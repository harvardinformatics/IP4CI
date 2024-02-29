#' @title  processCCAres_perCT
#' @description get ct-specific info from CCA res
#' @param cca_res CCA res
#' @param x dataframe score of dataset 1
#' @param y dataframe score of dataset 2
#' @return cca_res cca results updated for this x & y ct-based datasets
#' @examples
#' processCCAres_perCT(cca_res,x,y)
#' @export
################################################################################
processCCAres_perCT <-function(cca_res,x,y)
{
  print('processCCAres_perCT:..')

  x.aux.calc = cca_res$x.aux.calc[,colnames(x)]
  y.aux.calc = cca_res$y.aux.calc[,colnames(y)]
  print(dim(x.aux.calc))
  print(dim(y.aux.calc))
  u = cca_res$u[colnames(x.aux.calc),]
  v = cca_res$v[colnames(y.aux.calc),]
  print(dim(u))
  print(dim(v))
  xscores <- x.aux.calc %*% u# MASS::ginv(sm[[1]]) %*% res$u  pathwaysXcell cellXcc
  yscores <- y.aux.calc %*% v
  cor=diag(cor(xscores,yscores))
  print(cor)

  xscoresW <- apply(xscores,2,function(x) x/sqrt(sum(x*x)))
  yscoresW <- apply(yscores,2,function(x) x/sqrt(sum(x*x)))
  corW=diag(cor(xscoresW,yscoresW))


  no_cc_selected = length(cor)
  ul <- t(x.aux.calc) %*% xscores# sampleXpathway . pathwaysXcc == sampleXcc
  vl <- t(y.aux.calc) %*% yscores

  cca_res = list()
  cca_res$u =u
  cca_res$v =v
  cca_res$ul =ul
  cca_res$vl =vl
  cca_res$xscores =xscores
  cca_res$yscores =yscores
  cca_res$cor =cor

  cca_res$xscoresW =xscoresW
  cca_res$yscoresW =yscoresW
  cca_res$corW =corW

  cca_res$x =x
  cca_res$y =y
  colnames(cca_res$u)= seq(1:no_cc_selected)
  colnames(cca_res$v)= seq(1:no_cc_selected)
  colnames(cca_res$ul)= seq(1:no_cc_selected)
  colnames(cca_res$vl)= seq(1:no_cc_selected)
  colnames(cca_res$xscores)= seq(1:no_cc_selected)
  colnames(cca_res$yscores)= seq(1:no_cc_selected)
  colnames(cca_res$xscoresW)= seq(1:no_cc_selected)
  colnames(cca_res$yscoresW)= seq(1:no_cc_selected)

  return(cca_res)
}

################################################################################
