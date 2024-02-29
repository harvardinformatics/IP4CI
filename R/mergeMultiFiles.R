#' @title  mergeMultiFiles
#' @description merge sample scores from one dataset
#' @param my.files sample scores files
#' @return mergedFile merged file of all files
#' @examples
#' data_dir =  paste0(expr_dir, '03_sampleScoring/')
#'res_dir =  paste0(expr_dir, '04_mergeTwoDatasets/')
#' dir.create(res_dir)
#'pattern_s = paste0(filter_opt,filter_cutoff,sample_cutoff,'*',id1,'*.RData')
#'print(pattern_s)
#'mergingFiles_lst =  list.files(data_dir, pattern=glob2rx(pattern_s))
#' mergeMultiFiles(my.files)
#' @export
################################################################################

mergeMultiFiles <-function(my.files)
{
  df_lst = list()
  for(f in 1:length(my.files))
  {
    message('f')
    print(my.files[[f]])


    fname = my.files[[f]]
    load(fname)

    message('df_p_s_score:  ')
    print(dim(sampleScore$df_p_s_score))

    df_lst[[f]]=sampleScore$df_p_s_score
    df_lst[[f]]$ROWNAMES  <- rownames(df_lst[[f]])
  }
  message('df_lst:...')
  print(length(df_lst))

  merged_df =   join_all( df_lst, by='ROWNAMES', type="full" )
  rownames(merged_df) <- merged_df$ROWNAMES; merged_df$ROWNAMES <- NULL
  # merged_df = merged_df[rowSums(is.na(merged_df)) != ncol(merged_df),]
  merged_df[is.na(merged_df)] <- 0
  message('merged_df')
  print(dim(merged_df))

  return(merged_df)
}
################################################################################
