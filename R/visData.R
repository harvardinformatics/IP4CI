#' @title  visData
#' @description visualize dataframe
#' @param df dataframe
#' @param annot metadata for df
#' @param p.df dataframe holding pathways info
#' @param opfname file name to save plot
#' @examples
#' obj=x
#' visData(df,annot,p.df,opfname)
#'@export
################################################################################

visData<-function(df,annot,p.df,opfname)
{
  print('visData_df:...')
  df = t(df) # cell x p
  print(dim(df))

  pdf(opfname,width = 20, height =10)

# avg
    id1 = unique(annot$dataset_name)[1]
    print(id1)
    s1 = rownames(subset(annot, annot$dataset_name == id1))
    df_1 = as.matrix(df[rownames(df) %in% s1, ] )# pathway x cell
    annot_1 = annot[rownames(annot) %in% rownames(df_1),]
    print(dim(df_1))
    d1_avg =  setNames(data.frame(matrix(ncol = ncol(df_1), nrow = length(unique(annot_1$type)))),colnames(df_1)) # avg_ct x p
    row.names(d1_avg) = unique(sort(annot_1$type))
    for(l in unique(annot_1$type))
    {
      l_s = rownames(subset(annot_1, annot_1$type == l))
      df_l = df_1[l_s,]
      for(p in colnames(df_1))
      {
        d1_avg[l,p]=mean(df_l[,p],na.rm = T)
      }
    }
    rownames(d1_avg) = paste0(rownames(d1_avg),':',id1)
    id2 = unique(annot$dataset_name)[2]
    s2 = rownames(subset(annot, annot$dataset_name == id2))
    df_2 = as.matrix(df[rownames(df) %in% s2,] )# pathway x cell
    annot_2 = annot[rownames(annot) %in% rownames(df_2),]
    d2_avg =  setNames(data.frame(matrix(ncol = ncol(df_2), nrow = length(unique(annot_2$type)))),colnames(df_2)) # avg_ct x p
    row.names(d2_avg) = unique(sort(annot_2$type))
    for(l in unique(annot_2$type))
    {
      l_s = rownames(subset(annot_2, annot_2$type == l))
      df_l = df_2[l_s,]
      for(p in colnames(df_2))
      {
        d2_avg[l,p]=mean(df_l[,p],na.rm = T)
      }
    }
    rownames(d2_avg) = paste0(rownames(d2_avg),':',id2)
    df_avg = rbind(d1_avg,d2_avg)

    df_avg_annot_col = c('type','dataset')
    df_avg_annot =  setNames(data.frame(matrix(ncol = length(df_avg_annot_col), nrow =nrow(df_avg))),df_avg_annot_col) # avg_ct x p
    rownames(df_avg_annot)=rownames(df_avg)
    message('df_avg_annot')
    print(head(df_avg_annot))
    print(dim(df_avg_annot))

    df_avg_annot[c('type', 'dataset')] <- str_split_fixed(rownames(df_avg), ':', 2)
    df_avg_annot$both = paste0(df_avg_annot$type,':',df_avg_annot$dataset)
    print(head(df_avg_annot))
    print(dim(df_avg_annot))

    #d.umap = umap(df)
    for(k in c(3,5,7,10,20, 30))
    {
      if(k < nrow(df_avg))
      {
   d.umap <- uwot::umap(df_avg,n_neighbors = k, metric =  'correlation',verbose = F, ret_model = F,n_components = 2)# umap

   # avg.embedding = umap_transform(as.matrix(df_avg),d.umap)
   # print(head(avg.embedding))
   # print(dim(avg.embedding))
   # ump_plot <- data.frame(x =avg.embedding[,1], y = avg.embedding[,2])#
   ump_plot <- data.frame(x =d.umap[,1], y = d.umap[,2])
   p1= ump_plot%>%ggplot(aes(x, y , color = factor(df_avg_annot$type), label =df_avg_annot$both, shape =df_avg_annot$dataset) ) + #, label = annot$Sex
     geom_point(aes(x, y), size=4)+
     geom_text_repel(aes(x, y), max.overlaps = Inf,min.segment.length = Inf, box.padding = 0.5)+
     theme(text = element_text(size = 10),axis.text.x = element_text(size=10,angle = 90, hjust = 1),axis.text.y = element_text(size=10,vjust = 1),legend.position="right")+
     xlab("UMAP 1") + ylab("UMAP 2")+ggtitle(paste0('avg of gsea with n_neighbors = ',k))
print(p1)
      }
    }
    column_ha = HeatmapAnnotation(p_level = p.df$p_level)
    row_ha = rowAnnotation(cell_type = df_avg_annot$type,dataset_name = df_avg_annot$dataset,show_legend = TRUE)
    ph=Heatmap(df_avg, name = "heatmap of avg gsea", top_annotation = column_ha,right_annotation = row_ha ,cluster_rows =F,cluster_columns=FALSE,show_row_names=FALSE,show_column_names=FALSE,
               heatmap_legend_param = list( title_gp = gpar(face="bold",col = "darkgreen",fontsize =20),title_position = "leftcenter-rot",labels_gp = gpar(fontsize = 18),legend_height = unit(5, "cm"))
    ) #
    print(ph)


  if(length(unique(annot$type)) <10)
  {
  for(k in c(7,10,15))
  {
    ump <- uwot::umap(df,n_neighbors = k, metric =  'correlation',verbose = F, ret_model = F,n_components = 2)# umap
    message('ump')
    print(dim(ump))
    print(head(ump))
    ump_plot <- data.frame(x =ump[,1], y = ump[,2])#data.frame(x = ump[,1], y = ump[,2])
    message('ump_plot')
    ump_plot_annot = cbind(ump_plot,annot)
    if(length(unique(annot$dataset_name)) >1 ) # 2 datasets
    {
      p1= ump_plot_annot%>%ggplot(aes(x, y , color = factor(type))) + #, label = annot$Sex
        geom_point(aes(x, y))+
        #geom_text_repel(aes(x, y), max.overlaps = Inf,min.segment.length = Inf, box.padding = 0.5)+
        theme(text = element_text(size = 10),axis.text.x = element_text(size=10,angle = 90, hjust = 1),axis.text.y = element_text(size=10,vjust = 1),legend.position="right")+
        xlab("UMAP 1") + ylab("UMAP 2")+ggtitle(paste0('gsea with n_neighbors = ',k))
      p2= ump_plot_annot%>%ggplot(aes(x, y , color = dataset_name)) + #, label = annot$Sex
        geom_point(aes(x, y))+
        #geom_text_repel(aes(x, y), max.overlaps = Inf,min.segment.length = Inf, box.padding = 0.5)+
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
    print(p1 + p2)

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
}
  dev.off()

}
