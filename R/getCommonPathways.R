#' @title  getCommonPathways
#' @description get common  pathways between two lists of pathways
#' @param p1 CCA res
#' @param p2 metadata
#' @return commpnP common set of pathways
#' @examples
#' run_07_rankPathways()
#' @export
################################################################################
getCommonPathways <- function(p1,p2)
{
  print('getCommonSelectedPathways')

  p_ct_selected_1 = filter(p_ct_selected_1, grepl("R-HSA-*", p_ct_selected_1[,1]))
  print(dim(p_ct_selected_1))
  p_ct_selected_2 = filter(p_ct_selected_2, grepl("R-HSA-*", p_ct_selected_2[,1]))
  print(dim(p_ct_selected_2))
  message('common:..')
  if(opt == 'common')
    common_selected_pathways = intersect(p_ct_selected_1[,1],p_ct_selected_2[,1])
  if(opt == 'union')
    common_selected_pathways = union(p_ct_selected_1[,1],p_ct_selected_2[,1])
  print(length(common_selected_pathways))

  return(common_selected_pathways)
}
################################################################################################
