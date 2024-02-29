#' @title  run_02_gsea_fast
#' @description score samples with gsea for a single cell-type at a time
#' @param data_dir data directory folder name
#' @param res_dir result directory folder name
#' @param dataset_name id name for the object
#' @param obj object
#' @param pathwaydatabase pathway database
#' @param dataType which obj data slot to use
#' @param geneType which genes to use
#' @param ct name of cell-type
#' @examples
#' data_dir =  paste0(expr_dir,'/data/')
#' res_dir =  paste0(expr_dir,'/02_runGSEA/')
#' createDir(res_dir)
#' dataset_name='mouse'
#' obj='mouse_sc.RDS'
#' dataType = 'norm' # use normalized obj data slot
#' geneType = 'var' # use variable genes
#' pathwaydatabase = 'reactome.db'
#' run_02_gsea_fast(data_dir,res_dir, dataset_name,obj,dataType,geneType,pathwaydatabase,ct)
#' @export

###############################################################################
run_02_gsea_fast <- function(data_dir,res_dir, dataset_name,obj,dataType,geneType,pathwaydatabase,ct)
{
  print('run_02_gsea_fast:...')

  set.seed(1234)


  # select obj
  varfeat <- obj@assays[["RNA"]]@var.features

  if(geneType =='var' & dataType == 'norm')
    expr = obj@assays$RNA@data[varfeat,]
  if(geneType =='none' & dataType == 'norm')
    expr = obj@assays$RNA@data
  if(geneType =='var' & dataType == 'norm_scale')
    expr = obj@assays$RNA@scale.data[varfeat,]
  if(geneType =='none' & dataType == 'counts')
    expr = obj@assays$RNA@counts
  if(geneType =='var' & dataType == 'counts')
    expr = obj@assays$RNA@counts[varfeat,]
  annot =  obj@meta.data

  expr= as.matrix(expr)
  print(dim(expr))
  print(expr[1:3,1:3])
  print(class(expr))
  print(dim(annot))

  # serial run over cell-type of obj
    print(ct)

    calcLogFC_res = calcLogFC(expr,annot,ct)
    logFCs=calcLogFC_res[[1]]
    reference_samples = calcLogFC_res[[2]]
    expr_logFC = calcLogFC_res[[3]]
    GSEA_res = runGSEA(pathwaydatabase,logFCs,reference_samples)
    print(names(GSEA_res))
    opfname = paste(pathwaydatabase,dataset_name,dataType,geneType,ct,sep='_')
    save(GSEA_res,logFCs,expr_logFC,file = paste0(res_dir,opfname,'.RData'));

}

################################################################################
# data_name = 'panc'
# expr_dir =  paste0('/Users/rao198/OneDrive\ -\ Harvard\ University/Informatics/doc/paper/DA4LT/V1/analysis/',data_name,'/')
# data_dir =  paste0(expr_dir,'/01_processData/')
# res_dir =  paste0(expr_dir,'/02_runGSEA/')
# createDir(res_dir)
# id='mouse'
# obj='commonVarGmouse_sc.RDS'
# dataType = 'norm' # use normalized obj data slot
# geneType = 'var' # use variable genes
# pathwaydatabase = 'reactome.db'
# run_02_gsea(data_dir,res_dir, id,obj,dataType,geneType,pathwaydatabase)



