#' @title  run_00_preprocessData
#' @description preprocess data
#' @param data_dir data directory folder name
#' @param res_dir result directory folder name
#' @param id.list id names for the objects
#' @param obj.fname.list file names for the objects
#' @param typeColumnToUse.list list of column name holds the cell-type labels
#' @param assaySlotToUse.list list of name of assay slot
#' @examples
#' data_name = 'panc'
#' expr_dir =  paste0('/Users/rao198/OneDrive\ -\ Harvard\ University/Informatics/doc/paper/DA4LT/V1/analysis/',data_name)
#' data_dir =  paste0(expr_dir,'/data/')
#' res_dir =  paste0(expr_dir,'/01_processData/')
#' createDir(res_dir)
#' id.list=c('mouse','human')
#' obj.fname.list=c('mouse_sc.RDS','human_sc.RDS')
#' typeColumnToUse.list=c('cell.names','cell.names') # later replaced by 'type'
#' assaySlotToUse.list=list('integrated','SCT')
#' run_00_preprocessData(data_dir,res_dir, id.list,obj.fname.list,typeColumnToUse.list,assaySlotToUse.list)
#' @export

################################################################################
run_00_preprocessData <- function(data_dir,res_dir,
                                  id.list,obj.fname.list,
                                  typeColumnToUse.list,assaySlotToUse.list)
{
print('run_00_preprocessData:...')

  # 0. read
  obj1 = readRDS(paste0(data_dir,obj.fname.list[[1]]))
  print(obj1)
  obj2 = readRDS(paste0(data_dir,obj.fname.list[[2]]))
  print(obj2 )

  # 1. subset
  # df1 = obj1.tmp@assays[[assaySlotToUse.list[[1]]]]@counts
  # annot1 = obj1.tmp@meta.data
  # obj1<- CreateSeuratObject(
  #   counts=df1,
  #   project = "DA4CA.1",
  #   assay = "DA4CA.1",
  #   min.cells = 0,
  #   min.features =0,
  #   names.field = 1,
  #   names.delim = "_",
  #   meta.data = annot1)
  # df2 = obj2.tmp@assays[[assaySlotToUse.list[[2]]]]@counts
  # annot2 = obj2.tmp@meta.data
  # obj2<- CreateSeuratObject(
  #   counts=df2,
  #   project = "DA4CA.2",
  #   assay = "DA4CA.2",
  #   min.cells = 0,
  #   min.features =0,
  #   names.field = 1,
  #   names.delim = "_",
  #   meta.data = annot2)

  # OR 1. as recommended by Seurat
  assay.name.1 = assaySlotToUse.list[[1]]
  obj1[['DA4CA.1']] = obj1[[assay.name.1]]
  obj1[[assay.name.1]] = NULL
  assay.name.2 = assaySlotToUse.list[[2]]
  obj2[['DA4CA.2']] = obj2[[assay.name.2]]
  obj2[[assay.name.2]] = NULL

  # 2. deals with commonCol for each annot
  commonCol = intersect(colnames(obj1@meta.data),colnames(obj2@meta.data))
  obj1@meta.data = obj1@meta.data[,commonCol]
  obj2@meta.data = obj2@meta.data[,commonCol]
  if((! 'type' %in% typeColumnToUse.list) & ('type' %in% commonCol))
  {
    colnames(obj1@meta.data)[which(names(obj1@meta.data) == "type")] <- "type.org"
    colnames(obj2@meta.data)[which(names(obj2@meta.data) == "type")] <- "type.org"
  }
  if('type' != typeColumnToUse.list[[1]])
    obj1@meta.data$type=obj1@meta.data[,typeColumnToUse.list[[1]]]
  if('type' != typeColumnToUse.list[[2]])
    obj2@meta.data$type=obj2@meta.data[,typeColumnToUse.list[[2]]]

  saveRDS(obj1, file=paste0(res_dir,obj.fname.list[[1]]))
  saveRDS(obj2, file=paste0(res_dir,obj.fname.list[[2]]))

  print(obj1)
  print(obj2)

}
