#' @title  run_04_mergeTwoDatasets
#' @description merge sample scores from two datasets
#' @param df1 merged scoring files for dataset 1
#' @param df2 merged scoring files for dataset 2
#' @return  mergedFile merged file from two datasets
#' @examples
#' run_04_mergeTwoDatasets(merged_score1,merged_score2)
#' @export
################################################################################
run_04_mergeTwoDatasets <- function(df1,df2)
{

    print('mergeTwoDatasets:...')
    df_lst = list()
    idx = 1
    FeatureMergingType = 'common'

    if(FeatureMergingType == 'union')
      features = union(rownames(df1),rownames(df2))
    if(FeatureMergingType == 'common')
      features = intersect(rownames(df1),rownames(df2))

    df1 = df1[rownames(df1) %in% features, ]
    df2 = df2[rownames(df2) %in% features, ]
    #print(dim(df1))
    #print(dim(df2))

    df_lst[[1]]=as.data.frame(df1)
    df_lst[[1]]$ROWNAMES  <- rownames(df_lst[[1]])

    df_lst[[2]]=as.data.frame(df2)
    df_lst[[2]]$ROWNAMES  <- rownames(df_lst[[2]])
    #message('features:...',length(features))

    merged_df =   join_all( df_lst, by='ROWNAMES', type="full" )
    rownames(merged_df) <- merged_df$ROWNAMES; merged_df$ROWNAMES <- NULL
    # merged_df = merged_df[rowSums(is.na(merged_df)) != ncol(merged_df),]
    merged_df[is.na(merged_df)] <- 0


    return(merged_df)

  }

################################################################################
