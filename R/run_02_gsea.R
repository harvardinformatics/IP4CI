#' @title  run_02_gsea
#' @description score samples with gsea for serial run on cell-types of an obj
#' @param data_dir data directory folder name
#' @param res_dir result directory folder name
#' @param dataset_name id name for the object
#' @param expr data from object
#' @param annot annot from object
#' @param dataType which obj data slot to use
#' @param geneType which genes to use
#' @param pathwaydatabase pathway database
#' @examples
#' run_02_gsea(data_dir,res_dir, dataset_name,expr,annot,dataType,geneType,pathwaydatabase)
#' @export

###############################################################################
run_02_gsea <- function(data_dir,res_dir, dataset_name,expr,annot,dataType,geneType,pathwaydatabase)
{
  print('run_02_gsea:...')

  set.seed(1234)

  expr= as.matrix(expr)
  print(dim(expr))
  print(expr[1:3,1:3])
  print(class(expr))
  print(dim(annot))

  # serial run over cell-type of obj
  ct_obj_lst =  unique(annot$type)
      for(ct in ct_obj_lst ) # ct_obj_lst or for faster run use ct_obj_lst[ct_i] or the name of var c('macrophage')
      {
        print(ct)

      calcLogFC_res = calcLogFC(expr,annot,ct)
     logFCs=calcLogFC_res[[1]]
      reference_samples = calcLogFC_res[[2]]
      expr_logFC = calcLogFC_res[[3]]
    GSEA_res = runGSEA(pathwaydatabase,logFCs,reference_samples)
    print(names(GSEA_res))
    opfname = paste(pathwaydatabase,dataset_name,dataType,geneType,ct,sep='_')
    save(GSEA_res,logFCs,expr_logFC,file = paste0(res_dir,opfname,'.RData'));

   } # ct


}

################################################################################


