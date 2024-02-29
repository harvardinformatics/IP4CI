#' @title  run_02_gseaFast
#' @description score samples with gsea for a single cell-type at a time
#' @param data_dir data directory folder name
#' @param res_dir result directory folder name
#' @param dataset_name id name for the object
#' @param expr data from object
#' @param annot annot from object
#' @param dataType which obj data slot to use
#' @param geneType which genes to use
#' @param pathwaydatabase pathway database
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
#' run_02_gseaFast(data_dir,res_dir, dataset_name,expr,annot,dataType,geneType,pathwaydatabase,ct)
#' @export

###############################################################################
run_02_gseaFast <- function(data_dir,res_dir, dataset_name,expr,annot,dataType,geneType,pathwaydatabase,ct)
{
  print('run_02_gsea_fast:...')

  set.seed(1234)

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
    save(GSEA_res,expr_logFC,file = paste0(res_dir,opfname,'.RData')); #logFCs,

}

################################################################################



