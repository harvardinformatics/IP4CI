#' @title  getCommonGenesPerPathway
#' @description get gene members common per pathway between two datasets
#' @param enrichedG_1 genes enriched per pathway for dataset 1
#' @param enrichedG_2 genes enriched per pathway for dataset 2
#' @return p_g_common common genes per pathway
#' @examples
#' getCommonGenesPerPathway(rankedPathways,enrichedG_1,enrichedG_2)
#' @export
################################################################################
getCommonGenesPerPathway <- function(rankedPathways,enrichedG_1,enrichedG_2)
{
  print('getCommonGenesPerPathway:...')

  cn=c("p_id","g_common","g_common_no")
  commonP = intersect(enrichedG_1$p_id,enrichedG_2$p_id)
  print(length(commonP))
  p_g_common =  setNames(data.frame(matrix(ncol = length(cn), nrow = length(commonP))), cn)
  rownames(p_g_common) =  commonP
  print(dim(p_g_common))


  for(p in commonP) #rownames(rankedPathways)
  {
    print(p)
    p_g_common[p,'p_id']=p
    # cnames=c("p_simScore","g_count","g_rank_EntrezID","g_rank_Symbol","g_activity","g_activityCount")
    g1 = enrichedG_1[enrichedG_1$p_id %in% p,'g_id']
    g2= enrichedG_2[enrichedG_2$p_id %in% p,'g_id']
    common_g = intersect(g1,g2)
    print(common_g)
    p_g_common[p,'g_common']=toString(common_g)
    p_g_common[p,'g_common_no']=length(common_g)
  }

  p_g_common = merge(p_g_common,rankedPathways,by = 'p_id')
  print(head(p_g_common))
  print(dim(p_g_common))
  return(p_g_common)
}
