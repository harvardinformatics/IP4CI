rankMetric <- function(fname,tool,metric_name)
{
  print('rankMetric:... ')

  stat_tool =  read.table(fname, sep='\t', header=T,row.names = 1)
  print(head(stat_tool))
  #Sensitivity	Specificity	Pos Pred Value	Neg Pred Value	Precision	Recall	F1	Prevalence	Detection Rate	Detection Prevalence	Balanced Accuracy
  stat_tool_metric = as.data.frame(stat_tool[,metric_name])
  colnames(stat_tool_metric) = metric_name
  rownames(stat_tool_metric)=rownames(stat_tool)
  stat_tool_metric$tool = tool
  stat_tool_metric$cellType =rownames(stat_tool_metric)
  stat_tool_metric$cellType =gsub("Class: ", "", stat_tool_metric$cellType )
  stat_tool_metric$metric=round(as.numeric(stat_tool_metric[,metric_name]),2)
  print(head(stat_tool_metric))
  print(dim(stat_tool_metric))


  return(stat_tool_metric)
}


###############################################################################

rankMetric_acrossTools <-function(data_dir,ref_fldr,data_name,fname.pattern,tools.list,metric_names,T_lab)
{
  for(metric_name in metric_names)
  {
    print(metric_name)
    #Sensitivity	Specificity	Pos Pred Value	Neg Pred Value	Precision	Recall	F1	Prevalence	Detection Rate	Detection Prevalence	Balanced Accuracy
    #metric_name = 'Sensitivity'
    opfname = paste0(data_dir,'LT_metrics_Ranked_by_',data_name,"_",metric_name,"_",ref_fldr)
    df_stat_tools = list()
    for(i in 1:length(tools.list))
    {
      tool = tools.list[i]
      print(tool)

      fname = paste0(data_dir,tool,'/',fname.pattern)
      print(fname)
      df_tool = rankMetric(fname,tool,metric_name)
      df_stat_tools[[i]]= df_tool
    }
    print(df_stat_tools)
    df_merged = reshape::merge_recurse(df_stat_tools)
    df_merged[is.na(df_merged)] <- 0
    df_merged = df_merged[df_merged$cellType %in% T_lab ,]
    print(head(df_merged))
    print(dim(df_merged))
    write.table(df_merged, file = paste0(opfname,".txt"), sep = "\t", quote=FALSE, row.names=TRUE,col.names=TRUE)
    pdf(paste0(opfname,".pdf"), width= 20, height = 20)
    p=ggplot(data=df_merged, aes(x=tool, y=metric, fill=tool)) +
      geom_bar(position="dodge", stat="identity") +
      ggtitle(metric_name) +
      geom_text(aes(x=tool, y=metric/2,label =  metric), vjust = 0,  size = 10) +
      facet_wrap(~cellType) +
      theme_bw()+ theme(text = element_text(size = 35,face = "bold"),axis.text.x = element_text(angle = 90))



    print(p)
    dev.off()
  }
}
###############################################################################

#rank using all tools
# ref_fldr = 'MH'
# data_dir = paste0('/Users/rao198/OneDrive\ -\ Harvard\ University/Informatics/doc/paper/DA4LT/V1/scripts/DA4CA/relatedWork/brain/',ref_fldr,'/')
#
# data_name = 'brain'
# #T_lab = c('beta', 'ductal', 'delta', 'quiescent_stellate', 'endothelial', 'gamma', 'alpha', 'macrophage')
# T_lab=c('Vip_2','Lamp5_Rosehip','Vip_5','Lamp5_1','Sst_5','Vip_Sncg','Pax6','Vip_1','Vip_4','Vip_3','Pvalb_2','Lamp5_2','Sst_1','Sst_4','Chandelier','Sst_3','Sst_Chodl','Pvalb_1','Exc_L6_CT','Exc_L6b','Exc_L6_IT_1','Exc_L5.6_IT_1','Sst_2','Exc_L2.3_IT','Lamp5_Lhx6','Exc_L5_PT','Exc_L4.5_IT','Exc_L5.6_NP','Exc_L3.5_IT','Exc_L5.6_IT_2','Exc_L6_IT_2','Oligo','Exc_L5.6_IT_3','Astrocyte')
# fname.pattern = 'labelTransfer_Metrics.txt'
# tools.list = c('S3','scGCN','HieRFIT')
# metric_names=c('Sensitivity','Precision','Recall','F1','Balanced.Accuracy')
# rankMetric_acrossTools(data_dir,ref_fldr,data_name,fname.pattern,tools.list,metric_names,T_lab)
