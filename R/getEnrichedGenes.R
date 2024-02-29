#' @title  getEnrichedGenes
#' @description get gene members in pathways enriched in samples of interest based on gsea_res results
#' @param gsea_res gsea_res res
#' @param rankedPathways ranked pathways data.frame
#' @return enrichedG ranked pathways with info
#' @examples
#' getEnrichedGenes(gsea_res,rankedPathways)
#' @export
################################################################################
getEnrichedGenes <- function(gsea_res,rankedPathways)
{
  print('getEnrichedGenes:...')

  # ids<-bitr(genes_i, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
  cnames=c("g_EntrezID","g_id","p_activity_gsea","p_id","s_id")
  enrichedG = setNames(data.frame(matrix(ncol = length(cnames), nrow = 0)), cnames)
  idx=1

  for(i in 1:length(gsea_res)) # for s  length(gsea_res)
  {
    df = gsea_res[[i]]
    s=names(gsea_res)[i]

if(df != 0)
{
  print(class(df))
  print(df)
  df_p = subset(df,df$ID %in% rankedPathways$p_id)
  if(nrow(df_p)>0)
  {
    for(p in rownames(df_p))#for(j in 1:nrow(df)) # for p
    {
      print(s)
      #p=df[j,'ID']
      #print(dim(df))
      df_s = subset(df_p,df_p$ID == p)
      #print(head(df_s))
      if(nrow(df_s) > 0)
      {
        g_lst= df_s$core_enrichment#df[j,'core_enrichment']
        g_lst = as.list(strsplit(g_lst[1], '/')[[1]])

        p_act= df_s$NES
        #fill-in enrichedG
        #cnames=c("g_EntrezID","g_id","g_activity","p_id","s_id")
        for(g in g_lst)
        {
          enrichedG[idx,'g_EntrezID']=g
          enrichedG[idx,'p_activity_gsea']=p_act
          enrichedG[idx,'p_id']=p
          enrichedG[idx,'s_id']=s

          idx=idx+1
        }

      }

}

        }
      }

  }


  ids<-bitr(unique(enrichedG$g_EntrezID), fromType = "ENTREZID", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
  dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
  for(i in 1:nrow(enrichedG))
  {
    g =enrichedG[i,'g_EntrezID']
    enrichedG[i,'g_id']=dedup_ids[dedup_ids$ENTREZID == g,'SYMBOL']
  }


  enrichedG = merge(enrichedG,rankedPathways,by="p_id")
  #save(enrichedG,file = paste0(opfname,'_enrichedG.RData'))

  return(enrichedG)
}
################################################################################
