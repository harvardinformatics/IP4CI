
#' @title  run_05_dataIntegration
#' @description run Seurat data integarte over two objectes
#' @param obj1 object of dataset 1
#' @param obj2 object of dataset 2
#' @return integrated_res  integration obj
#' @examples
#' run_05_dataIntegration(obj1,obj2)
#' @export
################################################################################
run_05_dataIntegration <- function(obj1,obj2)
{
  print('run_05_dataIntegration:...')

  anchors <- FindIntegrationAnchors(object.list = list(obj1,obj2),reduction = "cca",dims = 1:30,anchor.features = 500) #union(rownames(obj1),rownames(obj2))   ,anchor.features = Fno
  integrated_res <- IntegrateData(anchorset = anchors,dims = 1:30)
  integrated_res <- ScaleData(integrated_res, assay = "integrated")
  integrated_res <- RunPCA(integrated_res, assay = "integrated")
  integrated_res <- RunUMAP(integrated_res,  dims = 1:30) #reduction = "pca",
  integrated_res <- RunTSNE(integrated_res, dims = 1:30)
  return(integrated_res)

}
