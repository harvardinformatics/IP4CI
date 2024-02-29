#' @title  run_01_processData
#' @description process data
#' @param data_dir data directory folder name
#' @param res_dir result directory folder name
#' @param id.list id names for the objects
#' @param obj.fname.list file names for the objects
#' @param genes2keep which genes to keep between the two objects
#' @return processedData list of processedData: obj1,obj2, and annot
#' @examples
#' data_name = 'panc'
#' expr_dir =  paste0('/Users/rao198/OneDrive\ -\ Harvard\ University/Informatics/doc/paper/DA4LT/V1/analysis/',data_name)
#' data_dir =  paste0(expr_dir,'/data/')
#' res_dir =  paste0(expr_dir,'/01_processData/')
#' createDir(res_dir)
#' id1='mouse'
#' id2='human'
#' obj1_fname='mouse_sc.RDS'
#' obj2_fname='human_sc.RDS'
#' genes2keep='commonVarG'
#' run_01_processData(data_dir,res_dir, id.list=c(id1,id2),obj.fname.list =c(obj1_fname,obj2_fname),genes2keep)
#' @export

################################################################################
run_01_processData <- function(data_dir,res_dir, id.list,obj.fname.list,genes2keep)
{
  print('run_01_processData:...')

  # 1. read obj
  # obj1 = readRDS(paste0(data_dir,obj.fname.list[[1]]))
  # print(obj1)
  # obj2 = readRDS(paste0(data_dir,obj.fname.list[[2]]))
  # print(obj2)

  # 2. add dataset_name, type col to metadata
  # obj1@meta.data$dataset_name=id.list[[1]]
  # obj2@meta.data$dataset_name=id.list[[2]]
  # Idents(object = obj1) = obj1@meta.data$type
  # Idents(object = obj2) = obj2@meta.data$type


  # 4. keep genes

if(genes2keep == 'unionVarG')
{
  print('consider each var. feat')
 # vf12 = union(vf1,vf2)
}

#   if(genes2keep == 'commonVarG')
  {
    #   vf1 = obj1@assays$DA4CA.1@var.features
    #   print(length(vf1))
    #   vf2 = obj2@assays$DA4CA.2@var.features
    #   print(length(vf2))


#     vf12 = intersect(vf1,vf2)
# message('vf12:  ',length(vf12))

  # obj1
  # obj1@assays$DA4CA.1@var.features= vf12
  # print('DONE VAR.FEAT')
  # print(dim( obj1@assays$DA4CA.1@counts ))
  # obj1@assays$DA4CA.1@counts =  obj1@assays$DA4CA.1@counts[vf12,] # OR GetAssayData(object = obj1, assay = "DA4CA.1", slot = "counts")
  # print('DONE counts')
  # obj1@assays$DA4CA.1@scale.data = obj1@assays$DA4CA.1@scale.data[vf12,]
  # print('DONE scale.data')
  # obj1@assays$DA4CA.1@data = obj1@assays$DA4CA.1@data[vf12,]
  # print('DONE data')


  # obj2
  # obj2@assays$DA4CA.2@var.features= vf12
  # obj2@assays$DA4CA.2@counts =  obj2@assays$DA4CA.2@counts[vf12,]
  # obj2@assays$DA4CA.2@scale.data = obj2@assays$DA4CA.2@scale.data[vf12,]
  # obj2@assays$DA4CA.2@data = obj2@assays$DA4CA.2@data[vf12,]
# print('DONE obj2')

}
  # 5. make sure cell-type names are R valid
  # obj1=RnamingCTconv(obj1)
  # obj2=RnamingCTconv(obj2)

  # saveRDS(obj1, file=paste0(res_dir,genes2keep,obj1_fname))
  # saveRDS(obj2, file=paste0(res_dir,genes2keep,obj2_fname))

  obj1 =readRDS(file=paste0(res_dir,genes2keep,obj1_fname))
  print(obj1)
  obj2 =readRDS(file=paste0(res_dir,genes2keep,obj2_fname))
  print(obj2)

  # 6.annot stat
  # print(colnames(obj1@meta.data))
  # print(colnames(obj2@meta.data))
  # commonCol = intersect(colnames(obj1@meta.data),colnames(obj2@meta.data))
  # obj1@meta.data = obj1@meta.data[,commonCol]
  # obj2@meta.data = obj2@meta.data[,commonCol]
  annot =rbind(obj1@meta.data,obj2@meta.data)
  # plotMetadataStat(annot, data_dir)
  # annot_stat = getMetadataStat(annot)
  # write.table(annot_stat, file = paste0(data_dir,'annot_stat.txt'), sep = "\t", quote = FALSE, row.names = TRUE)

  # 7.vis
  dim_name = 'umap'
  # obj1 = visObj(obj1,dim_name,paste0(res_dir,id.list[[1]],'_vis.pdf'))
  # obj2 = visObj(obj2,dim_name,paste0(res_dir,id.list[[2]],'_vis.pdf'))

  # 8. general analysis on obj list including: marker genes, gsea, Go & pathway enrichment analysis
  # markers1 = FindAllMarkers(obj1, logfc.threshold = 0.25, min.pct = 0.1, test.use="MAST")
  # saveRDS(markers1, file=paste0(res_dir,id.list[[1]],'_markers.RDS'))
  # markers1 = readRDS(file=paste0(res_dir,id.list[[1]],'_markers.RDS'))
  # markers2 = FindAllMarkers(obj2, logfc.threshold = 0.25, min.pct = 0.1, test.use="MAST")
  # saveRDS(markers2, file=paste0(res_dir,id.list[[2]],'_markers.RDS'))
  # markers2 = readRDS(file=paste0(res_dir,id.list[[2]],'_markers.RDS'))
  # performMarkersAnalysis(markers1,id.list[[1]],res_dir)
  # performMarkersAnalysis(markers2,id.list[[2]],res_dir)
  # getOverlapMarkers(markers1,id.list[[1]],markers2,id.list[[2]],res_dir)

  processedData = list("obj1" = obj1,"obj2" =obj2, "annot"= annot)
  return(processedData)
}
################################################################################

