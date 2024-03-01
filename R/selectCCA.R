#' @title  selectCCA
#' @description select important CCA dims
#' @param cca_res CCA results
#' @param select_opt option to select dims
#' @param select_opt_cutoff cutoff to select dims
#' @return selected_dims selected dims to use
#' @examples
#' selectCCA(cca_res,select_opt,select_opt_cutoff)
#' @export
################################################################################
selectCCA <- function(cca_res,select_opt,select_opt_cutoff)
{
  print('selectCCA:...')
  selected_dims = c(1:30)
  print(cca_res$cor)
  names(cca_res$cor) = seq(1:length(cca_res$cor))

  # select based on cor
  if(select_opt == 'cor')
  {
    idx = names(cca_res$cor[which(cca_res$cor >= select_opt_cutoff) ])
    print(idx)
    selected_dims = idx
  }

  # select baaed on var
  if(select_opt == 'var')
  {
    cor_sq = (cca_res$cor^2)
    var = cor_sq #/sum(cor_sq)
    print(var)
    idx = names(var[which(var >= select_opt_cutoff) ])
    print(idx)
    selected_dims = idx
  }

  # select baaed on var
  if(select_opt == 'CC')
  {
    idx = sort(cca_res$cor, index.return=TRUE,decreasing=T)$ix
    print(idx)
    selected_dims = idx[1:select_opt_cutoff]
  }
  # select based on both
print(selected_dims)

  #How do you explain the variance of a correlation?
  #The strength of the relationship between X and Y is sometimes expressed by squaring the correlation coefficient and multiplying by 100. The resulting statistic is known as variance explained (or R2). Example: a correlation of 0.5 means 0.52x100 = 25% of the variance in Y is "explained" or predicted by the X variable.

  # update
  cca_res_updated = list()
  cca_res_updated$pembd = cca_res$pembd[,selected_dims]
  cca_res_updated$u =cca_res$u[,selected_dims]
  cca_res_updated$v =cca_res$v[,selected_dims]
  cca_res_updated$ul =cca_res$ul[,selected_dims]
  cca_res_updated$vl =cca_res$vl[,selected_dims]
  cca_res_updated$ul2 =cca_res$ul2[,selected_dims]
  cca_res_updated$vl2 =cca_res$vl2[,selected_dims]
  cca_res_updated$xscores =cca_res$xscores[,selected_dims]
  cca_res_updated$yscores =cca_res$yscores[,selected_dims]
  cca_res_updated$cor =cca_res$cor[selected_dims]
  cca_res_updated$xscoresW =cca_res$xscoresW[,selected_dims]
  cca_res_updated$yscoresW =cca_res$yscoresW[,selected_dims]
  cca_res_updated$corW =cca_res$corW[selected_dims]
  cca_res_updated$x =cca_res$x
  cca_res_updated$y =cca_res$y
  cca_res_updated$x.aux.calc = cca_res$x.aux.calc
  cca_res_updated$y.aux.calc = cca_res$y.aux.calc
  no_cc_selected = length(selected_dims)
  colnames(cca_res_updated$u)= seq(1:no_cc_selected)
  colnames(cca_res_updated$v)= seq(1:no_cc_selected)
  colnames(cca_res_updated$ul)= seq(1:no_cc_selected)
  colnames(cca_res_updated$vl)= seq(1:no_cc_selected)
  colnames(cca_res_updated$xscores)= seq(1:no_cc_selected)
  colnames(cca_res_updated$yscores)= seq(1:no_cc_selected)

  return(cca_res_updated)

}
