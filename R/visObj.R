#' @title  visObj
#' @description visualize obj using Seurat
#' @param obj object
#' @param dim_name dim name to use 'umap' or 'cca'
#' @param opfname file name to save plot
#' @examples
#'visObj(obj,dim_name,opfname)
#'@export

################################################################################
visObj <- function(obj,dim_name,opfname)
{
  print('visObj:...')
  print(obj)
  print(Idents(object = obj))

  message('obj dataset name : ',length(unique(obj@meta.data$dataset_name)))
  pdf(opfname, width = 12, height = 6)
  if(length(unique(obj@meta.data$dataset_name))>1) # 2 datasets
  {
    p1=DimPlot(obj, group.by   ="dataset_name",pt.size =  0.9,reduction = dim_name, label = TRUE,repel = TRUE) + ggtitle(paste0(unique(obj@meta.data$dataset_name)[1],'&',unique(obj@meta.data$dataset_name)[2]))+
      theme(plot.title = element_text(size = 25, face = "bold"))+ NoLegend()

    obj@meta.data$type_dataset = paste0(obj@meta.data$dataset_name,':',obj@meta.data$type)
    p2=DimPlot(obj, group.by  ="type_dataset",pt.size =  0.9,reduction = dim_name, label = TRUE,repel = TRUE) + ggtitle(paste0(unique(obj@meta.data$dataset_name)[1],'&',unique(obj@meta.data$dataset_name)[2]))+
      theme(plot.title = element_text(size = 25, face = "bold")) +NoLegend()

    p3=DimPlot(obj, split.by   ="dataset_name",pt.size =  0.9,reduction = dim_name, label = TRUE,repel = TRUE) + ggtitle(paste0(unique(obj@meta.data$dataset_name)[1],'&',unique(obj@meta.data$dataset_name)[2]))+
      theme(plot.title = element_text(size = 25, face = "bold"))+ NoLegend()

    p4=DimPlot(obj, group.by   ="type",pt.size =  0.9,reduction = dim_name, label = TRUE,repel = TRUE) + ggtitle(paste0(unique(obj@meta.data$dataset_name)[1],'&',unique(obj@meta.data$dataset_name)[2]))+
      theme(plot.title = element_text(size = 25, face = "bold"))+ NoLegend()

    # p3= DimPlot(obj, reduction = "umap", group.by = "type", label = TRUE, repel = TRUE) +
    #   NoLegend()

    # my_colors <- scales::hue_pal()(length(unique(obj@meta.data$type)))
    # p2=DoHeatmap(object = obj)+ ggtitle(paste0(unique(obj@meta.data$dataset_name)[1],'&',unique(obj@meta.data$dataset_name)[2]))+
    #   scale_color_manual(values = my_colors,limits = unique(obj@meta.data$type))#features =
    print(p1)
    print(p2)
    print(p3)
    print(p4)
  }
  else # one dataset
  {
    p1=DimPlot(obj, group.by='type' , shape.by ="dataset_name",pt.size =  0.9,reduction =dim_name, label = TRUE,repel = TRUE) + ggtitle(unique(obj@meta.data$dataset_name)[1])+
      theme(plot.title = element_text(size = 25, face = "bold")) #oLegend()+
    my_colors <- scales::hue_pal()(length(unique(obj@meta.data$type)))
    p2=DoHeatmap(object = obj)+ ggtitle(unique(obj@meta.data$dataset_name))+
      scale_color_manual(values = my_colors,limits = unique(obj@meta.data$type))#features =
    print(p1)
    print(p2)

  }

  # plot varFeat
  # top10 <- head(VariableFeatures(obj), 10)
  # p3 <- VariableFeaturePlot(obj)
  # p4 <- LabelPoints(plot = p3, points = top10, repel = TRUE)
  # print(p3 + p4)

  # plot markers
 #  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
 #  markers %>%
 #    group_by(cluster) %>%
 #    slice_max(n = 2, order_by = avg_log2FC)
 #  markers %>%
 #    group_by(cluster) %>%
 #    top_n(n = 10, wt = avg_log2FC) -> top10
 # p5= DoHeatmap(obj, features = top10$gene) + NoLegend()
 #  print(p5)

  dev.off()
}
