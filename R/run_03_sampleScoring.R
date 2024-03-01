#' @title  run_03_sampleScoring
#' @description score and filter samples based on gsea results
#' @param GSEA_res gsea results
#' @param filter_opt column of GSEA to filter based on i.e., pvalue
#' @param filter_cutoff 0.05
#' @param sample_cutoff percentage of samples meets the cutoff of filter i.e. 0.1
#' @examples
#' run_03_sampleScoring(GSEA_res,filter_opt,filter_cutoff,sample_cutoff)
#' @export
###############################################################################
run_03_sampleScoring <- function(GSEA_res,filter_opt,filter_cutoff,sample_cutoff)
{

  print('run_03_sampleScoring:...')

  gsea_res = GSEA_res$gsea_res
  pathways = GSEA_res$pathways_lst
    samples = names(gsea_res)
      df_p_s_score = setNames(data.frame(matrix(0,ncol = length(samples), nrow = length(pathways))),samples )
      rownames(df_p_s_score) = pathways
      df_p_s_sig = setNames(data.frame(matrix(0,ncol = length(samples), nrow = length(pathways))),samples )
      rownames(df_p_s_sig) = pathways
      df_p_s_sig_c = setNames(data.frame(matrix(0,ncol = length(samples), nrow = length(pathways))),samples )
      rownames(df_p_s_sig_c) = pathways

print(length(gsea_res))
      for(i in 1:length(gsea_res))
      {
        message('i:  ',i)
        s_name = names(gsea_res)[i]
        s_gsea = gsea_res[[i]]
#print(s_gsea)
        if(s_gsea != 0)
        {
          for(p in s_gsea$ID)
          {
            message('p:..',p,'    i:   ',i,'     s_name:   ',s_name)
            tmp = subset(s_gsea, s_gsea$ID == p)

            es = tmp$NES
sig = tmp[,filter_opt]

            df_p_s_score[p,s_name]=as.numeric(es)
            df_p_s_sig[p,s_name]=as.numeric(sig)

          }
        }
        else
        {
          df_p_s_score[,s_name]=0
          df_p_s_sig[,s_name]=100
        }
        message('DONE WITH SAMPLE# ',s_name)
      }

      print(dim(df_p_s_score))
      message('min:  ',min(df_p_s_score),max(df_p_s_score))

      print(dim(df_p_s_sig))
      message('min:  ',min(df_p_s_sig),max(df_p_s_sig))

      selected_pathways_sig = list()
      for(r in rownames(df_p_s_sig))
      {
        #message('r:..',r)
        tmp = df_p_s_sig[r,]
        df_sig_c = length(tmp[tmp < filter_cutoff])
        t =  round(length(tmp) * sample_cutoff)
        if(df_sig_c >= t) # df_sig_c > round(length(tmp) * 0.1) #10% samples > 3
          selected_pathways_sig = c(selected_pathways_sig,r)
      }
      selected_pathways_sig = unique(unlist(selected_pathways_sig))
      df_p_s_score = df_p_s_score[selected_pathways_sig ,]
      message('selected_pathways_sig')
      print(length(selected_pathways_sig))
      message('df_p_s_score')
      print(dim(df_p_s_score))

      return(sampleScore = list('selected_pathways_sig'=selected_pathways_sig,'df_p_s_score'=df_p_s_score,'df_p_s_sig'=df_p_s_sig))

}
###############################################################################
# expr_dir = '/Users/rao198/OneDrive\ -\ Harvard\ University/Informatics/doc/paper/DA4LT/V1/analysis/panc'
# data_dir =  paste0(expr_dir,'/02_runGSEA/')
# res_dir =  paste0(expr_dir,'/03_sampleScoring/')
# createDir(res_dir)
# d1
# scoringFiles_lst = c('reactome.db_mouse_norm_var_macrophage.RData')
# for(f in scoringFiles_lst)
# {
#   load(paste0(data_dir,f))
#   sampleScore=run_03_sampleScoring(data_dir,res_dir,GSEA_res,filter_opt,filter_cutoff,sample_cutoff)
#   print(names(sampleScore))
#   save(sampleScore,file =paste0(res_dir,filter_opt,filter_cutoff,sample_cutoff,'_',f))
# }
