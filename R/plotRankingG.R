#' @title  plotRankingG
#' @description plot ranking genes per pathway
#' @param df dataframe to plot
#' @param col_annot1 first column annotation
#' @param col_annot2 second column annotation
#' @param title title of plot
#' @param opfname file to print plot
#' @examples
#' plotRankingG(df,col_annot1,col_annot2,title,opfname)
#' @export
################################################################################
plotRankingG <- function(df,col_annot1,col_annot2,title,opfname)
{
  print('plotRankingG:...')

  pdf(opfname, width = 12, height = 12)

  column_ha = HeatmapAnnotation(g_rank = (col_annot1),g_diff =  anno_barplot(col_annot2))
  row_ha = rowAnnotation(dataset_name =c('mouse','human'),show_legend = TRUE)
  ph=Heatmap(
    df, name =title, top_annotation = column_ha,right_annotation = row_ha ,cluster_rows = F,cluster_columns=FALSE,show_row_names=T,show_column_names=T,
    heatmap_legend_param = list( title_gp = gpar(face="bold",fontsize =12))
  )
  print(ph)

  dev.off()

}
