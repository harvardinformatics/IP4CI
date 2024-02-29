#' @title  visData
#' @description visualize dataframe
#' @examples
#' obj=x
#' processObj(obj)
#'@export
################################################################################

visData<-function(df,annot,p.df,opfname)
{
  print('visData_df:...')
  df = t(df)
  print(dim(df))

  pdf(opfname,width = 20, height =10)

  library(uwot)
  for(k in c(7,10,15))
  {
    ump <- uwot::umap(df,n_neighbors = k, metric =  'correlation',verbose = F, ret_model = F,n_components = 2)# umap
    message('ump')
    print(dim(ump))
    print(head(ump))
    ump_plot <- data.frame(x =ump[,1], y = ump[,2])#data.frame(x = ump[,1], y = ump[,2])
    message('ump_plot')

    if(length(unique(annot$dataset_name)) >1 ) # 2 datasets
    {
      p1=ggplot(ump_plot,aes(ump_plot$x, ump_plot$y , color = factor(annot$dataset_name))) + #, label = annot$Sex
        geom_point(aes(x=ump_plot[,1], y=ump_plot[,2] , shape = annot$type), size=3)+
        # geom_text_repel(aes(x=ump_plot[,1], y=ump_plot[,2]), max.overlaps = Inf,min.segment.length = Inf, box.padding = 0.5)+
        theme(text = element_text(size = 10),axis.text.x = element_text(size=10,angle = 90, hjust = 1),axis.text.y = element_text(size=10,vjust = 1),legend.position="right")+
        xlab("UMAP 1") + ylab("UMAP 2")+ggtitle(paste0('gsea with n_neighbors = ',k))

    }
    else
    {
      p1=ggplot(ump_plot,aes(ump_plot$x, ump_plot$y , color = factor(annot$type))) + #, label = annot$Sex
        geom_point(aes(x=ump_plot[,1], y=ump_plot[,2] , shape = annot$dataset_name), size=3)+
        # geom_text_repel(aes(x=ump_plot[,1], y=ump_plot[,2]), max.overlaps = Inf,min.segment.length = Inf, box.padding = 0.5)+
        theme(text = element_text(size = 10),axis.text.x = element_text(size=10,angle = 90, hjust = 1),axis.text.y = element_text(size=10,vjust = 1),legend.position="right")+
        xlab("UMAP 1") + ylab("UMAP 2")+ggtitle(paste0('gsea with n_neighbors = ',k))

    }
    print(p1)

  }
  column_ha = HeatmapAnnotation(p_level = p.df$p_level)
  row_ha = rowAnnotation(cell_type = annot$type,dataset_name = annot$dataset_name,show_legend = TRUE)
  ph=Heatmap(df, name = "heatmap of gsea", top_annotation = column_ha,right_annotation = row_ha ,cluster_rows = T,cluster_columns=FALSE,show_row_names=FALSE,show_column_names=FALSE,
             heatmap_legend_param = list( title_gp = gpar(face="bold",col = "darkgreen",fontsize =20),title_position = "leftcenter-rot",labels_gp = gpar(fontsize = 18),legend_height = unit(5, "cm"))
  ) #
  print(ph)
  ph=Heatmap(df, name = "heatmap of gsea", top_annotation = column_ha,right_annotation = row_ha ,cluster_rows =F,cluster_columns=FALSE,show_row_names=FALSE,show_column_names=FALSE,
             heatmap_legend_param = list( title_gp = gpar(face="bold",col = "darkgreen",fontsize =20),title_position = "leftcenter-rot",labels_gp = gpar(fontsize = 18),legend_height = unit(5, "cm"))
  ) #
  print(ph)
  # pos_df = data.frame("SamplesLabels" = annot$type)
  # rownames(pos_df) = rownames(df)# name matching
  # p5=pheatmap(df,
  #             clustering_distance_cols = "correlation",
  #             cluster_cols = F,
  #             cluster_rows = F,
  #             show_rownames = F,
  #             show_colnames = F,
  #             legend = T,
  #             fontsize =12,
  #             fontsize_row = 14,
  #             fontsize_col=14,
  #             treeheight_row = 5,
  #             treeheight_col = 5,
  #             #color = colorRampPalette(brewer.pal(9,"Reds"))(400),
  #             annotation_row = pos_df,
  #             main = 'heatmaps of GSEA scores using pheatmap')
  #
  dev.off()

}
