#' @title  run_09_detectUnknownCT
#' @description detect unknown cell-type between two datasets (Reference and Query)
#' @param obj1 obj of the first dataset
#' @param obj2 obj of the second dataset
#' @param annot metadata
#' @param id.list name list of the two datasets
#' @param ct cell-type to detect
#' @param res_dir results directory
#' @example
#' run_09_detectUnknownCT(obj1,obj2,annot,id.list,ct, res_dir)
#' @export
###############################################################################
run_09_detectUnknownCT <- function(obj1,obj2,annot,id.list, ct,res_dir)
{
  print('run_09_detectUnknownCT:...')


  #01 run CCA
  print(head(rownames(x = obj1)))
  print(length(rownames(x = obj1)))

  print(length(colnames(x = obj1)))
  print(length(colnames(x = obj2)))

  print(length(unique(obj1@meta.data$type)))
  print(length(unique(obj2@meta.data$type)))

  # run CCA
  no_cc=30
  cca_res = RunCCA(obj1, obj2, num.cc = no_cc,features = rownames(x = obj1),renormalize = FALSE,rescale = FALSE,verbose = TRUE) #features = rownames(df)
  cca_call= cca_res[['cca']]

# process CCA results
  cca_res = processCCAres_allCT(cca_res, obj1@counts, obj1,obj2@counts,obj2)

  # #B. avg sim
  proj_p_c = rbind((cca_res$u %*% t(cca_res$xscores)),(cca_res$v %*% t(cca_res$yscores)))
  calcSim_avgCT(df=t(proj_p_c),annot,id.list, ct,res_dir)
}
################################################################################
