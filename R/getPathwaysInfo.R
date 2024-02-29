#' @title  getPathwaysInfo
#' @description merge sample scores from two datasets
#' @param pathwaysToAnnotate set of pathways to annotate
#' @return  P_annnot data.frame of pathways info
#' @examples
#' getPathwaysInfo(pathwaysToAnnotate)
#' @export
################################################################################

getPathwaysInfo <-function(pathwaysToAnnotate)
{
  print('getPathwaysInfo:..')
  p.df = data.frame()
  idx=1
  for(p1 in pathwaysToAnnotate)
  {
    topLevel_df='notFound'
    tname='notFound'
    topLevel_df_s='notFound'
    p_level_top='notFound'
    tryCatch({
      topLevel_df = rba_reactome_event_ancestors(p1)[[1]]
      #print(topLevel_df)
      tname=rba_reactome_query(ids = p1, attribute_name = "displayName")
      #print(tname)
      topLevel_df_s = subset(topLevel_df, topLevel_df$schemaClass == 'TopLevelPathway')
      #print(topLevel_df_s)
      p_level_top=toString(topLevel_df_s$displayName)
      print(p_level_top)
    },
    error=function(e){})

    p.df[idx,'p_name']=tname
    p.df[idx,'p_level']=p_level_top
    p.df[idx,'p_id']=p1
    idx=idx+1
  }

  rownames(p.df)=p.df$p_id

  return(p.df)

}
##############################################################################################

