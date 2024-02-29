#' @title  run_08_rankGenes
#' @description rank genes in rankedPathways based on the expression of genes between the two datasets
#' @param expr1 gene expression or logFC of dataset 1
#' @param expr2 gene expression or logFC of dataset 2
#' @param commonG common genes per pathway between the two datasets
#' @param rankedPathways ranked pathways with info
#' @return rankG ranked genes with info
#' @examples
#' run_08_rankGenes(expr1,expr2,commonG,rankedrankedPathways)
#' @export
################################################################################
run_08_rankGenes <- function(expr1,expr2,commonG,rankedPathways)
{
  print('run_08_rankGenes:...')
  print(dim(expr1))
  print(dim(expr2))
  expr1.avg<- rowMeans(expr1, na.rm = TRUE)
  expr2.avg<- rowMeans(expr2, na.rm = TRUE)
message('avg')
  print(length(expr1.avg))
print(length(expr2.avg))
  cnames = c('g_id','g_lfc1','g_lfc2','g_diff','g_sign','g_rank','g_expr_rank',
             'p_id')
  rankG = setNames(data.frame(matrix(ncol = length(cnames), nrow = 0)), cnames)
  rankG.idx=1

  for(p in commonG$p_id) #rownames(rankedrankedPathways)
  {
    print(p)
    common_g = commonG[commonG$p_id == p,'g_common']
    #print(common_g)
    common_g =  as.list(strsplit(common_g, ', ')[[1]])
    print(length(common_g))
    expr1.avg.s = expr1.avg[names(expr1.avg) %in% common_g]
    print(length(expr1.avg.s))
    expr2.avg.s = expr2.avg[names(expr2.avg) %in% common_g]
    print(length(expr2.avg.s))
    expr.avg <- merge(as.data.frame(expr1.avg.s), as.data.frame(expr2.avg.s), by=0, all=TRUE)
    rownames(expr.avg)= expr.avg$Row.names
    expr.avg$Row.names=NULL
    expr.avg = expr.avg %>% replace(is.na(.), 0)
    colnames(expr.avg) =c('X','Y')
    print(head(expr.avg))
    print(tail(expr.avg))
    print(dim(expr.avg))

    expr.avg.fltr =  expr.avg[abs(expr.avg$X) > 0.1 & abs(expr.avg$Y) > 0.1,]# half lfc or 0.25
message('expr.avg.fltr')
print(head(expr.avg.fltr))
print(dim(expr.avg.fltr))
    # loop over genes
    rankG.s = setNames(data.frame(matrix(ncol = length(cnames), nrow = 0)), cnames)

    for(idx in 1:nrow(expr.avg.fltr))
    {
      g= rownames(expr.avg.fltr)[idx]
      print(g)

      rankG.s[idx,'g_id']=g
      rankG.s[idx,'g_lfc1']=expr.avg.fltr[g,'X']
      rankG.s[idx,'g_lfc2']=expr.avg.fltr[g,'Y']
      rankG.s[idx,'p_id']=p
      rankG.s[idx,'g_diff']= ifelse((sign(expr.avg.fltr[g,'X']) * (sign(expr.avg.fltr[g,'Y']))) ==1, abs(abs(expr.avg.fltr[g,'X'])-abs(expr.avg.fltr[g,'Y'])),(abs(expr.avg.fltr[g,'X'])+abs(expr.avg.fltr[g,'Y'])))
      rankG.s[idx,'g_sign']=ifelse((sign(expr.avg.fltr[g,'X']) * (sign(expr.avg.fltr[g,'Y']))) ==1,sign(expr.avg.fltr[g,'X']) ,0)
      rankG.s[idx,'g_expr_rank']= abs(expr.avg.fltr[g,'X']) + abs(expr.avg.fltr[g,'Y'])
    }
    # rank genes
    print('rankG.s')
    print('rankG.s')
    g_diff_sameDirection = rankG.s[rankG.s$g_sign != 0,]
    print(g_diff_sameDirection)
    g_diff_OppositeDirection = rankG.s[rankG.s$g_sign == 0,]
    print(g_diff_OppositeDirection)
    rank_idx = 1

    if(nrow(g_diff_sameDirection) >0)
    {
      g_diff_sameDirection = g_diff_sameDirection[order(g_diff_sameDirection$g_diff, decreasing = F),]
      # rank first the higher lfc
      for(g in g_diff_sameDirection$g_id)
      {
        rankG.s[rankG.s$g_id == g,'g_rank']=rank_idx
        rank_idx=rank_idx+1
      }
    }

    if(nrow(g_diff_OppositeDirection) >0)
    {
      g_diff_OppositeDirection = g_diff_OppositeDirection[order(g_diff_OppositeDirection$g_diff, decreasing = F),]

      for(g in g_diff_OppositeDirection$g_id)
      {
        rankG.s[rankG.s$g_id == g,'g_rank']=rank_idx
        rank_idx=rank_idx+1
      }
    }
    # rank by total expr in both dataset for each gene
    rank_idx = 1
    g_expr_rank = rankG.s[order(rankG.s$g_expr_rank, decreasing = T),]
    # rank first the higher lfc
    for(g in g_expr_rank$g_id)
    {
      rankG.s[rankG.s$g_id == g,'g_expr_rank']=rank_idx
      rank_idx=rank_idx+1
    }

    rankG=rbind(rankG,rankG.s)
  }

  rankG =  merge(rankG,commonG, by="p_id")
  drop=c('p_sim.x')
  rankG=rankG[,!(names(rankG) %in% drop)]
  names(rankG)[names(rankG) == "p_sim.y"] <- "p_sim"
  print(head(rankG))

  # write.table(rankG, file = paste0(opfname,'_rankingG.txt'), sep = "\t", quote = TRUE, row.names = TRUE, col.names = T)
  return(rankG)
}

################################################################################
