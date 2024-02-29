#' @title  getOverlapMarkers
#' @description get Overlap Markers between two objs
#' @param Markers1 Markers of obj1
#' @param Markers2 Markers of obj2
#' @param id1 id of markers1
#' @param id2 id of markers2
#' @param res_dir results_directory to save results and plot
#' @examples
#'getOverlapMarkers(markers1,markers2,res_dir)
#'@export
###############################################################################
getOverlapMarkers <- function(markers1,id1,markers2,id2,res_dir)
{
  print('getOverlapMarkers:...')

  print(markers1)
  print(markers2)
  # 1. visualize markers
  top1 = markers1 %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
  top2 = markers2 %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)

  # overlp markers bw obj1 & obj2
  df =  setNames(data.frame(matrix(ncol = length(unique(top1$cluster)), nrow = length(unique(top2$cluster)))),sort(unique(top1$cluster)))
  row.names(df) = sort(unique(top2$cluster))
  print(df)
  df.c =  setNames(data.frame(matrix(ncol = length(unique(top1$cluster)), nrow = length(unique(top2$cluster)))),sort(unique(top1$cluster)))
  row.names(df.c) = sort(unique(top2$cluster))
  for(i in unique(top1$cluster))
  {
    g1= (subset(top1,top1$cluster == i))$gene
    #i= paste0('mouse:',i)
    #message('ct: ',i)

    for(j in unique(top2$cluster))
    {
      g2= (subset(top2,top2$cluster == j))$gene
      #j= paste0('human:',j)
      #message('ct: ',j)
      overlapG = jaccard(g1,g2)#intersect(g1,g2)
      df[j,i]=toString(intersect(g1,g2))
      df.c[j,i]=overlapG#length(overlapG)
    }
  }
  df <- df[ order(rownames(df)),order(colnames(df)) ]
  df.c <- df.c[ order(rownames(df.c)),order(colnames(df.c)) ]
  write.table(df, paste0(res_dir,'markersOverlapG.txt'),sep='\t',row.names = T, col.names = TRUE,quote = FALSE)
  write.table(df.c, paste0(res_dir,'markersOverlapGcount.txt'),sep='\t',row.names = T, col.names = TRUE,quote = FALSE)
  df.mlt = melt(as.matrix(df.c)) #     Var1                Var2 value
  print(df.mlt)
  pdf(paste0(res_dir,'markersOverlapG.pdf'), width=15, height = 15)
  heatmap_plot= ggplot(data = df.mlt, aes(x = X1, y = X2)) +
    theme_classic(base_size = 25)+
    theme(legend.position = "none",axis.text.x=element_text(angle=45,hjust=1))+
    geom_point(aes(size = value))+
    coord_fixed()+
    xlab(id2) + ylab(id1)
  print(heatmap_plot)
  dev.off()

  # overlp markers bw cell-type in obj1
  df =  setNames(data.frame(matrix(ncol = length(unique(top1$cluster)), nrow = length(unique(top1$cluster)))),sort(unique(top1$cluster)))
  row.names(df) = sort(unique(top1$cluster))
  print(df)
  df.c =  setNames(data.frame(matrix(ncol = length(unique(top1$cluster)), nrow = length(unique(top1$cluster)))),sort(unique(top1$cluster)))
  row.names(df.c) = sort(unique(top1$cluster))
  for(i in unique(top1$cluster))
  {
    g1= (subset(top1,top1$cluster == i))$gene
    #i= paste0('mouse:',i)
    #message('ct: ',i)

    for(j in unique(top1$cluster))
    {
      g2= (subset(top1,top1$cluster == j))$gene
      #j= paste0('human:',j)
      #message('ct: ',j)
      overlapG = jaccard(g1,g2)#intersect(g1,g2)
      df[j,i]=toString(intersect(g1,g2))
      df.c[j,i]=overlapG#length(overlapG)
    }
  }
  df <- df[ order(rownames(df)),order(colnames(df)) ]
  df.c <- df.c[ order(rownames(df.c)),order(colnames(df.c)) ]
  write.table(df, paste0(res_dir,'markersOverlapG_1.txt'),sep='\t',row.names = T, col.names = TRUE,quote = FALSE)
  write.table(df.c, paste0(res_dir,'markersOverlapGcount_1.txt'),sep='\t',row.names = T, col.names = TRUE,quote = FALSE)
  df.mlt = melt(as.matrix(df.c)) #     Var1                Var2 value
  print(df.mlt)
  pdf(paste0(res_dir,'markersOverlapG_1.pdf'), width=15, height = 15)
  heatmap_plot= ggplot(data = df.mlt, aes(x = X1, y = X2)) +
    theme_classic(base_size = 25)+
    theme(legend.position = "none",axis.text.x=element_text(angle=45,hjust=1))+
    geom_point(aes(size = value))+
    coord_fixed()+
    xlab(id1) + ylab(id1)
  print(heatmap_plot)
  dev.off()

  # overlp markers bw cell-type in obj2
  df =  setNames(data.frame(matrix(ncol = length(unique(top2$cluster)), nrow = length(unique(top2$cluster)))),sort(unique(top2$cluster)))
  row.names(df) = sort(unique(top2$cluster))
  print(df)
  df.c =  setNames(data.frame(matrix(ncol = length(unique(top2$cluster)), nrow = length(unique(top2$cluster)))),sort(unique(top2$cluster)))
  row.names(df.c) = sort(unique(top2$cluster))
  for(i in unique(top2$cluster))
  {
    g1= (subset(top2,top2$cluster == i))$gene
    #i= paste0('mouse:',i)
    #message('ct: ',i)

    for(j in unique(top2$cluster))
    {
      g2= (subset(top2,top2$cluster == j))$gene
      #j= paste0('human:',j)
      #message('ct: ',j)
      overlapG = jaccard(g1,g2)#intersect(g1,g2)
      df[j,i]=toString(intersect(g1,g2))
      df.c[j,i]=overlapG#length(overlapG)
    }
  }
  df <- df[ order(rownames(df)),order(colnames(df)) ]
  df.c <- df.c[ order(rownames(df.c)),order(colnames(df.c)) ]
  write.table(df, paste0(res_dir,'markersOverlapG_2.txt'),sep='\t',row.names = T, col.names = TRUE,quote = FALSE)
  write.table(df.c, paste0(res_dir,'markersOverlapGcount_2.txt'),sep='\t',row.names = T, col.names = TRUE,quote = FALSE)
  df.mlt = melt(as.matrix(df.c)) #     Var1                Var2 value
  print(df.mlt)
  pdf(paste0(res_dir,'markersOverlapG_2.pdf'), width=15, height = 15)
  heatmap_plot= ggplot(data = df.mlt, aes(x = X1, y = X2)) +
    theme_classic(base_size = 25)+
    theme(legend.position = "none",axis.text.x=element_text(angle=45,hjust=1))+
    geom_point(aes(size = value))+
    coord_fixed()+
    xlab(id2) + ylab(id2)
  print(heatmap_plot)
  dev.off()

}
