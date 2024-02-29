#' @title  plotRankingP
#' @description plot pathways ranking
#' @param df dataframe to plot
#' @param col_annot1 first column annotation
#' @param col_annot2 second column annotation
#' @param row_annot1 first row annotation
#' @param title title of plot
#' @param opfname file to print plot
#' @examples
#' plotRankingP(df,col_annot1,col_annot2,row_annot1,title,opfname)
#' @export
################################################################################
plotRankingP <- function(df,col_annot1,col_annot2,row_annot1,title,opfname)
{
  print('plotRankingP:...')

  pdf(opfname, width = 20, height = 12)

  column_ha = HeatmapAnnotation(p_level = (col_annot1),p_sim =  anno_barplot(col_annot2))
  row_ha = rowAnnotation(dataset_name =row_annot1,show_legend = TRUE)
  ph=Heatmap(
    df, name =title, top_annotation = column_ha,right_annotation = row_ha ,cluster_rows = F,cluster_columns=FALSE,show_row_names=F,show_column_names=T,
    heatmap_legend_param = list( title_gp = gpar(face="bold",fontsize =12))
  )
  print(ph)

  dev.off()
}
