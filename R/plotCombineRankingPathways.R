#' @title  plotCombineRankingPathways
#' @description combine rankingP across all cell-types
#' @param fname.lst list of ranking pathways files
#' @param ct.lst list of cell-types according to the files
#' @param res_dir result directory to save the plot
#' @examples
#' fname.lst = paste0(data_dir,c('rankP_quiescent_stellate.txt','rankP_macrophage.txt'))
#' ct.lst = c("quiescent_stellate","macrophage")
#' plotCombineRankingPathways(fname.lst,ct.lst,res_dir)
#' @export
################################################################################
plotCombineRankingPathways <- function(fname.lst,ct.lst,res_dir)
{
  cnames=c("p_id","p_level","p_name","cor")
  res.df =setNames(data.frame(matrix(NA,ncol = length(cnames), nrow = 0)),cnames )
  df.lst = list()
  idx=1
  col = 'cor'
  for(fname in fname.lst)
  {
    ct =ct.lst[[idx]]
    print(ct)
    print(fname)
    data =read.table(file = fname, header=T, sep="\t", row.names = 1)
    data = data
    print(head(data))
    print(dim(data))
    data = data[,cnames]
    colnames(data)[which(names(data) == col)] <- ct
    print(head(data))
    print(dim(data))

    df.lst[[idx]]=data
    df.lst[[idx]]$ROWNAMES  <- rownames(df.lst[[idx]])
    idx=idx+1
  }
  message('df.lst:...')
  print(length(df.lst))

  merged_df =   join_all( df.lst, by='ROWNAMES', type="full" )
  rownames(merged_df) <- merged_df$ROWNAMES; merged_df$ROWNAMES <- NULL
  # merged_df = merged_df[rowSums(is.na(merged_df)) != ncol(merged_df),]
  #merged_df[is.na(merged_df)] <- 0

  #res.df = merge(res.df,data, by =c("p_id","p_level","p_name"), all = TRUE)
  drop=c('cor')
  merged_df=merged_df[,!(names(merged_df) %in% drop)]
  merged_df= merged_df %>% distinct(p_id, .keep_all = TRUE)
  print(head(merged_df))
  print(dim(merged_df))

  set.seed(123)

  print(dim(merged_df))
  rownames(merged_df) = merged_df$p_id

  ct_names = sort(ct.lst)
  percentage.df = setNames(data.frame(matrix(ncol = 2, nrow = length(ct_names))), c('cons','div'))
  rownames(percentage.df)=ct_names
  for(ct in ct_names)
  {
    print(ct)
    tmp=as.numeric(merged_df[,ct])
    tmp = tmp[! is.na(tmp)]
    print(tmp)
    percentage.df[ct,'cons'] =(100)* length(tmp[tmp > 0])/length(tmp)
    percentage.df[ct,'div'] = (100)*length(tmp[tmp < 0])/length(tmp)
  }

  #merged_df=merged_df %>% replace(is.na(.), 0)
  column_ha = HeatmapAnnotation(name = merged_df$p_name, pathway_level = merged_df$p_level,
                                annotation_legend_param = list(pathway_level = list(title = "pathway_level", title_gp = gpar(face="bold",col = "darkgreen",fontsize =8),title_position = "leftcenter-rot",labels_gp = gpar(fontsize = 6))))
  drop=c('p_id',	'p_level',	'p_name')
  merged_df=merged_df[,!(names(merged_df) %in% drop)]
  percentage.df = percentage.df[colnames(merged_df),]

  row_ha = rowAnnotation(cellType = colnames(merged_df),cons = anno_barplot(percentage.df$cons),div = anno_barplot(percentage.df$div),show_legend = FALSE)
  pdf(paste0(res_dir,'rankP_combined.pdf'), width=10, height = 6)
  p1=Heatmap(t(merged_df), name = "rankP_combined", top_annotation = column_ha, right_annotation = row_ha, cluster_rows = TRUE,clustering_distance_rows= 'pearson',cluster_columns=FALSE,
             row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 2),
             heatmap_legend_param = list(title = paste0("Pathway signatures",'\n',"across cell types"), title_gp = gpar(face="bold",col = "darkgreen",fontsize =8),title_position = "leftcenter-rot",labels_gp = gpar(fontsize = 6),
                                         at = c(-1, 1),labels = c("diverged", "conserved"),legend_height = unit(2, "cm")))
  print(p1)
  p1=Heatmap(t(merged_df), name = "rankP_combined", top_annotation = column_ha, right_annotation = row_ha, cluster_rows = F,cluster_columns=FALSE,
             row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 2),
             heatmap_legend_param = list(title = paste0("Pathway signatures",'\n',"across cell types"), title_gp = gpar(face="bold",col = "darkgreen",fontsize =8),title_position = "leftcenter-rot",labels_gp = gpar(fontsize = 6),
                                         at = c(-1, 1),labels = c("diverged", "conserved"),legend_height = unit(2, "cm")))
  print(p1)
  dev.off()
}
