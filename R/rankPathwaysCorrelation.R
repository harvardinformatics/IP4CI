#' @title  rankPathwaysCorrelation
#' @description rank pathways based on correlation
#' @param cca_res CCA res
#' @param id.list id names for the objects
#' @param p.info pathways info
#' @return rankP.cor  ranked pathways based on correlation
#' @examples
#' rankPathwaysRegression()
#' @export
################################################################################
rankPathwaysCorrelation <- function(cca_res,p.info)
{
  print('rankPathwaysCorrelation:...')
  x_p = cca_res$xscores
  y_p = cca_res$yscores
  print(dim(x_p))
  print(dim(y_p))
  cnames=c('p_id','cor')
  rankP.cor =  setNames(data.frame(matrix(ncol = length(cnames), nrow = nrow(x_p))),cnames)
  rownames(rankP.cor) = rownames(x_p)
  print((rankP.cor))

  for(p1 in rownames(x_p))
  {
    print(p1)

    px=x_p[p1,]
    py=y_p[p1,]
    cor_res = cor.test(px,py)
    message(cor_res$estimate,'\t',cor_res$p.value)
    rankP.cor[p1,'p_id']=p1
    rankP.cor[p1,'cor']=cor_res$estimate
  }
  print(rankP.cor)
  rankP.cor$norm = (rankP.cor$cor-min(rankP.cor$cor))/(max(rankP.cor$cor)-min(rankP.cor$cor))
  print(rankP.cor)
  # topP = (subset(rankP.cor,rankP.cor$cor >=0.20))$X
  # print(topP)
  # downP = (subset(rankP.cor,rankP.cor$cor < 0.20))$X
  # print(downP)
  # others = rankP.cor$X[!  rankP.cor$X %in% topP & !  rankP.cor$X %in% downP]
  # rankP.cor$status = ifelse(rankP.cor$X %in% downP , 'Diverged',ifelse(rankP.cor$X %in% topP , 'Conserved','Others') )

  rankP.cor =  rankP.cor[order(rankP.cor$cor, decreasing = T),]
  print(head(rankP.cor))
  print(tail(rankP.cor))

  rankP.cor = merge(rankP.cor,p.info,by="p_id")


  return(rankP.cor)

}
################################################################################
