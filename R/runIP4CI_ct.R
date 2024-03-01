#' @title  runIP4CI_ct
#' @description run runIP4CI_ct main function to run the analysis using two objects for step 03

#' @param expr_dir experiment directory folder name
#' @param id.list id names for the objects
#' @param obj.fname.list file names for the objects
#' @param convertG.list whether to convert mouse gene to human genes T/F
#' @param typeColumnToUse.list list of column name holds the cell-type labels
#' @param assaySlotToUse.list list of name of assay slot
#' @param genes2keep which genes to keep between the two objects
#' @param processObjOpt T or F to process object with visualization
#' @param pathwaydatabase pathway database
#' @param dataType which obj data slot to use
#' @param geneType which genes to use
#' @param filter_opt column of GSEA to filter based on i.e., pvalue
#' @param filter_cutoff 0.05
#' @param sample_cutoff percentage of samples meets the cutoff of filter i.e. 0.1
#' @param topRanked 10 top ranked biological findings to investigate, used for pathways & genes
#' @param ct cell-type lable


#' @export
################################################################################

runIP4CI_ct <- function(expr_dir, id.list,obj.fname.list,
                     convertG.list,
                     typeColumnToUse.list,assaySlotToUse.list,
                     genes2keep,
                     processObjOpt,
                     pathwaydatabase,dataType,geneType,
                     filter_opt,filter_cutoff,sample_cutoff,
                     topRanked,
                     ct)
{
  print('start running the main function of IP4CI :...')


  # 02.run gsea
  data_dir =  paste0(expr_dir, '01_processData/')
  res_dir =  paste0(expr_dir, '02_runGSEA/')
  dir.create(res_dir)
  if(dataType == 'norm' & geneType == 'var')
  {
    df1=obj1@assays$IP4CI.1@data
    annot1 = obj1@meta.data
    df2=obj2@assays$IP4CI.2@data
    annot2 = obj2@meta.data
  }
  #run_02_gsea(data_dir,res_dir, id.list[[1]] ,df1,annot1,dataType,geneType,pathwaydatabase)# slow
  run_02_gseaFast(data_dir,res_dir,id.list[[1]] ,df1,annot1,dataType,geneType,pathwaydatabase,ct) # fast
  #run_02_gsea(data_dir,res_dir, id.list[[2]],df2,annot2,dataType,geneType,pathwaydatabase)
  run_02_gseaFast(data_dir,res_dir,id.list[[2]],df2,annot2,dataType,geneType,pathwaydatabase,ct) # fast
  ################################################################################
}
