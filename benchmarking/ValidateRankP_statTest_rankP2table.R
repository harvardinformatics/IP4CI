library(stringr)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)


################################################################################


################################################################################

visData<-function(df,annot,p.df,opfname)
{
  print('visData_df:...')
  df = t(df) # cell x p
  print(dim(df))
  
  pdf(opfname,width = 20, height =10)
  
  if(length(unique(annot$dataset_name))>1)
  {
    print('if(length(unique(annot$dataset_name))>1)')
    # avg
    id1 = unique(annot$dataset_name)[1]
    print(id1)
    s1 = rownames(subset(annot, annot$dataset_name == id1))
    df_1 = as.matrix(df[rownames(df) %in% s1, ] )# pathway x cell
    print(dim(df_1))
    annot_1 = annot[rownames(annot) %in% rownames(df_1),]
    print(dim(annot_1))
    d1_avg =  setNames(data.frame(matrix(ncol = ncol(df_1), nrow = length(unique(annot_1$type)))),colnames(df_1)) # avg_ct x p
    row.names(d1_avg) = unique(sort(annot_1$type))
    for(l in unique(annot_1$type))
    {
      l_s = rownames(subset(annot_1, annot_1$type == l))
      df_l = df_1[l_s,]
      
      for(p in colnames(df_1))
      {
        
        d1_avg[l,p]=mean(df_l[,p],na.rm = T)
      }
    }
    rownames(d1_avg) = paste0(rownames(d1_avg),':',id1)
    print('DONE d1_avg : ...')
    
    id2 = unique(annot$dataset_name)[2]
    s2 = rownames(subset(annot, annot$dataset_name == id2))
    df_2 = as.matrix(df[rownames(df) %in% s2,] )# pathway x cell
    annot_2 = annot[rownames(annot) %in% rownames(df_2),]
    d2_avg =  setNames(data.frame(matrix(ncol = ncol(df_2), nrow = length(unique(annot_2$type)))),colnames(df_2)) # avg_ct x p
    row.names(d2_avg) = unique(sort(annot_2$type))
    for(l in unique(annot_2$type))
    {
      l_s = rownames(subset(annot_2, annot_2$type == l))
      df_l = df_2[l_s,]
      print(dim(df_l))
      for(p in colnames(df_2))
      {
        d2_avg[l,p]=mean(df_l[,p],na.rm = T)
      }
    }
    rownames(d2_avg) = paste0(rownames(d2_avg),':',id2)
    df_avg = rbind(d1_avg,d2_avg)
    
    df_avg_annot_col = c('type','dataset')
    df_avg_annot =  setNames(data.frame(matrix(ncol = length(df_avg_annot_col), nrow =nrow(df_avg))),df_avg_annot_col) # avg_ct x p
    rownames(df_avg_annot)=rownames(df_avg)
    message('df_avg_annot')
    print(head(df_avg_annot))
    print(dim(df_avg_annot))
    
    df_avg_annot[c('type', 'dataset')] <- str_split_fixed(rownames(df_avg), ':', 2)
    df_avg_annot$both = paste0(df_avg_annot$type,':',df_avg_annot$dataset)
    print(head(df_avg_annot))
    print(dim(df_avg_annot))
    
    #d.umap = umap(df)
    for(k in c(3,5,7,10,20, 30))
    {
      if(k < nrow(df_avg))
      {
        d.umap <- uwot::umap(df_avg,n_neighbors = k, metric =  'correlation',verbose = F, ret_model = F,n_components = 2)# umap
        
        # avg.embedding = umap_transform(as.matrix(df_avg),d.umap)
        # print(head(avg.embedding))
        # print(dim(avg.embedding))
        # ump_plot <- data.frame(x =avg.embedding[,1], y = avg.embedding[,2])#
        ump_plot <- data.frame(x =d.umap[,1], y = d.umap[,2])
        p1= ump_plot%>%ggplot(aes(x, y , color = factor(df_avg_annot$type), label =df_avg_annot$both, shape =df_avg_annot$dataset) ) + #, label = annot$Sex
          geom_point(aes(x, y), size=4)+
          geom_text_repel(aes(x, y), max.overlaps = Inf,min.segment.length = Inf, box.padding = 0.5)+
          theme(text = element_text(size = 10),axis.text.x = element_text(size=10,angle = 90, hjust = 1),axis.text.y = element_text(size=10,vjust = 1),legend.position="right")+
          xlab("UMAP 1") + ylab("UMAP 2")+ggtitle(paste0('avg of gsea with n_neighbors = ',k))
        print(p1)
      }
    }
    column_ha = HeatmapAnnotation(p_level = p.df$p_level)
    row_ha = rowAnnotation(cell_type = df_avg_annot$type,dataset_name = df_avg_annot$dataset,show_legend = TRUE)
    ph=Heatmap(df_avg, name = "heatmap of avg gsea", top_annotation = column_ha,right_annotation = row_ha ,cluster_rows =F,cluster_columns=FALSE,show_row_names=FALSE,show_column_names=FALSE,
               heatmap_legend_param = list( title_gp = gpar(face="bold",col = "darkgreen",fontsize =20),title_position = "leftcenter-rot",labels_gp = gpar(fontsize = 18),legend_height = unit(5, "cm"))
    ) #
    print(ph)
  }
  if(length(unique(annot$dataset_name))==1)
  {
    # avg
    id1 = unique(annot$dataset_name)[1]
    print(id1)
    s1 = rownames(df)#rownames(subset(annot, annot$dataset_name == id1))
    df_1 = as.matrix(df)#as.matrix(df[rownames(df) %in% s1, ] )# pathway x cell
    print(dim(df_1))
    annot_1 = annot[rownames(annot) %in% rownames(df_1),]
    print(dim(annot_1))
    d1_avg =  setNames(data.frame(matrix(ncol = ncol(df_1), nrow = length(unique(annot_1$type)))),colnames(df_1)) # avg_ct x p
    row.names(d1_avg) = unique(sort(annot_1$type))
    for(l in unique(annot_1$type))
    {
      l_s = rownames(subset(annot_1, annot_1$type == l))
      df_l = df_1[l_s,]
      
      for(p in colnames(df_1))
      {
        
        d1_avg[l,p]=mean(df_l[,p],na.rm = T)
      }
    }
    rownames(d1_avg) = paste0(rownames(d1_avg),':',id1)
    print('DONE d1_avg : ...')
    
    
    df_avg = d1_avg
    df_avg_annot_col = c('type','dataset')
    df_avg_annot =  setNames(data.frame(matrix(ncol = length(df_avg_annot_col), nrow =nrow(df_avg))),df_avg_annot_col) # avg_ct x p
    rownames(df_avg_annot)=rownames(df_avg)
    message('df_avg_annot')
    print(head(df_avg_annot))
    print(dim(df_avg_annot))
    
    df_avg_annot[c('type', 'dataset')] <- str_split_fixed(rownames(df_avg), ':', 2)
    df_avg_annot$both = paste0(df_avg_annot$type,':',df_avg_annot$dataset)
    print(head(df_avg_annot))
    print(dim(df_avg_annot))
    
    #d.umap = umap(df)
    for(k in c(3,5,7,10,20, 30))
    {
      if(k < nrow(df_avg))
      {
        d.umap <- uwot::umap(df_avg,n_neighbors = k, metric =  'correlation',verbose = F, ret_model = F,n_components = 2)# umap
        
        # avg.embedding = umap_transform(as.matrix(df_avg),d.umap)
        # print(head(avg.embedding))
        # print(dim(avg.embedding))
        # ump_plot <- data.frame(x =avg.embedding[,1], y = avg.embedding[,2])#
        ump_plot <- data.frame(x =d.umap[,1], y = d.umap[,2])
        p1= ump_plot%>%ggplot(aes(x, y , color = factor(df_avg_annot$type), label =df_avg_annot$both, shape =df_avg_annot$dataset) ) + #, label = annot$Sex
          geom_point(aes(x, y), size=4)+
          geom_text_repel(aes(x, y), max.overlaps = Inf,min.segment.length = Inf, box.padding = 0.5)+
          theme(text = element_text(size = 10),axis.text.x = element_text(size=10,angle = 90, hjust = 1),axis.text.y = element_text(size=10,vjust = 1),legend.position="right")+
          xlab("UMAP 1") + ylab("UMAP 2")+ggtitle(paste0('avg of gsea with n_neighbors = ',k))
        print(p1)
      }
    }
    column_ha = HeatmapAnnotation(p_level = p.df$p_level)
    row_ha = rowAnnotation(cell_type = df_avg_annot$type,dataset_name = df_avg_annot$dataset,show_legend = TRUE)
    ph=Heatmap(df_avg, name = "heatmap of avg gsea", top_annotation = column_ha,right_annotation = row_ha ,cluster_rows =F,cluster_columns=FALSE,show_row_names=FALSE,show_column_names=FALSE,
               heatmap_legend_param = list( title_gp = gpar(face="bold",col = "darkgreen",fontsize =20),title_position = "leftcenter-rot",labels_gp = gpar(fontsize = 18),legend_height = unit(5, "cm"))
    ) #
    print(ph)
  }
  
  if(length(unique(annot$type)) <10 | length(unique(annot$dataset_name))==1 )
  {
    for(k in c(7,10,15))
    {
      ump <- uwot::umap(df,n_neighbors = k, metric =  'correlation',verbose = F, ret_model = F,n_components = 2)# umap
      message('ump')
      print(dim(ump))
      print(head(ump))
      ump_plot <- data.frame(x =ump[,1], y = ump[,2])#data.frame(x = ump[,1], y = ump[,2])
      message('ump_plot')
      ump_plot_annot = cbind(ump_plot,annot)
      if(length(unique(annot$dataset_name)) >1 ) # 2 datasets
      {
        p1= ump_plot_annot%>%ggplot(aes(x, y , color = factor(type))) + #, label = annot$Sex
          geom_point(aes(x, y))+
          #geom_text_repel(aes(x, y), max.overlaps = Inf,min.segment.length = Inf, box.padding = 0.5)+
          theme(text = element_text(size = 10),axis.text.x = element_text(size=10,angle = 90, hjust = 1),axis.text.y = element_text(size=10,vjust = 1),legend.position="right")+
          xlab("UMAP 1") + ylab("UMAP 2")+ggtitle(paste0('gsea with n_neighbors = ',k))
        p2= ump_plot_annot%>%ggplot(aes(x, y , color = dataset_name)) + #, label = annot$Sex
          geom_point(aes(x, y))+
          #geom_text_repel(aes(x, y), max.overlaps = Inf,min.segment.length = Inf, box.padding = 0.5)+
          theme(text = element_text(size = 10),axis.text.x = element_text(size=10,angle = 90, hjust = 1),axis.text.y = element_text(size=10,vjust = 1),legend.position="right")+
          xlab("UMAP 1") + ylab("UMAP 2")+ggtitle(paste0('gsea with n_neighbors = ',k))
        
        print(p1 + p2)
      }
      else
      {
        p1=ggplot(ump_plot,aes(ump_plot$x, ump_plot$y , color = factor(annot$type))) + #, label = annot$Sex
          geom_point(aes(x=ump_plot[,1], y=ump_plot[,2] , shape = annot$dataset_name), size=3)+
          # geom_text_repel(aes(x=ump_plot[,1], y=ump_plot[,2]), max.overlaps = Inf,min.segment.length = Inf, box.padding = 0.5)+
          theme(text = element_text(size = 10),axis.text.x = element_text(size=10,angle = 90, hjust = 1),axis.text.y = element_text(size=10,vjust = 1),legend.position="right")+
          xlab("UMAP 1") + ylab("UMAP 2")+ggtitle(paste0('gsea with n_neighbors = ',k))
        print(p1)
      }
      
      
    }
    column_ha = HeatmapAnnotation(p_level = p.df$p_level)
    row_ha = rowAnnotation(cell_type = annot$type,dataset_name = annot$dataset_name,show_legend = TRUE)
    ph=Heatmap(df, name = "heatmap of gsea", top_annotation = column_ha,right_annotation = row_ha ,cluster_rows = T,cluster_columns=FALSE,show_row_names=FALSE,show_column_names=FALSE,
               heatmap_legend_param = list( title_gp = gpar(face="bold",col = "darkgreen",fontsize =20),title_position = "leftcenter-rot",labels_gp = gpar(fontsize = 18),legend_height = unit(5, "cm"))
    ) #
    print(ph)
    ph=Heatmap(df, name = "heatmap of gsea", top_annotation = column_ha,right_annotation = row_ha ,cluster_rows =F,cluster_columns=FALSE,show_row_names=FALSE,show_column_names=FALSE,
               heatmap_legend_param = list( title_gp = gpar(face="bold",col = "darkgreen",fontsize =20),title_position = "leftcenter-rot",labels_gp = gpar(fontsize = 18),legend_height = unit(5, "cm"))
    ) #
    print(ph)
  }
  dev.off()
  
}
################################################################################

expr_dir ='/rolayan/newBrain/'
region = 'AMY/'

################################################################################
# var needed
id.list = c('mouse','human')
FeatureMergingType ='union'
filter_opt = 'pvalue'
filter_cutoff = 0.05
sample_cutoff = 0.1
num.cc = 30
selected_dims_opt = 'none'
selected_dims_opt_cutoff =0
nVarF='all'
ObjType='single'
################################################################################

# files needed

data_dir =  paste0(expr_dir,region, '04_mergeTwoDatasets/')
annot = read.table(file = paste0(data_dir,'annot_Merged_',FeatureMergingType,'.txt'), header=T, sep="\t", row.names = 1)
print(dim(annot))
merged_score12.Pinfo = read.table(file = paste0(data_dir,id.list[[1]],'_',id.list[[2]],'_Merged_',FeatureMergingType,'_PathwaysInfo.txt'), header=T, sep="\t", row.names = 1)
print(dim(merged_score12.Pinfo))
merged_score12 = read.table(file = paste0(data_dir,id.list[[1]],'_',id.list[[2]],'_Merged_',FeatureMergingType,'.txt'), header=T, sep="\t", row.names = 1,check.names = FALSE)
print(dim(merged_score12))
data_dir =  paste0(expr_dir,region, '07_rankPathways/')
rankP = read.table(file = paste0(data_dir,'rankP_allCT.txt'), header=T, sep="\t", row.names = 1)
print(dim(rankP))
merged_score12.p = merged_score12[rankP$p_id,]
print(dim(merged_score12.p))
data_dir =  paste0(expr_dir,region, '05_CCA/')
cca_res = readRDS(paste0(data_dir,ObjType,'_',num.cc,'_',nVarF,'_',selected_dims_opt,selected_dims_opt_cutoff,'.RDS'))
print(dim(cca_res$pembd))
print(dim(cca_res$u))
print(dim(cca_res$v))

################################################################################
pathwayCC = cca_res$pembd
cellCC= rbind(cca_res$u,cca_res$v)
pcCC= (pathwayCC) %*% t(cellCC)
# rownames(pcCC) = rownames(pathwayCC)
# colnames(pcCC) = rownames(cellCC)
# print(dim(pcCC)) # pxcell
data_dir =  paste0(expr_dir,region, '07_rankPathways/')
visData(pcCC,annot[colnames(pcCC),],merged_score12.Pinfo[rownames(pcCC),],paste0(data_dir,'projectedCC_visData.pdf'))


annot.1=subset(annot, annot$dataset_name == id.list[[1]])
annot.2=subset(annot, annot$dataset_name == id.list[[2]])
ct1.lst = unique(annot.1$type)
ct2.lst = unique(annot.2$type)

# t.test based on GSEA
# t.test based on projection of cellxCC & pathwaysxCC


for(ct1 in ct1.lst)
{
  cnames=c("ct_id","p_id","p_level","p_name","p_all_cor","p_ct_cor",
           "p_selected_1","p_selected_2","p_selected_common",
           'p_ttest_gsea','p_ttpval_gsea','p_ftest_gsea','p_ftpval_gsea',
           'p_avg_gsea_1','p_avg_gsea_2',
           'p_med_gsea_1','p_med_gsea_2',
           'p_diff_gsea','p_diffSign_gsea',
           'p_ttest_CPproj','p_ttpval_CPproj','p_ftest_CPproj','p_ftpval_CPproj')
  c_p_org_proj_info =  setNames(data.frame(matrix(ncol = length(cnames), nrow = 0)), cnames)
  idx=1
  print(ct1)
  s1 =rownames(subset(annot.1, annot.1$type == ct1))
  print(length(s1))
  load(paste0(paste0(expr_dir, region,'03_sampleScoring/'),filter_opt,filter_cutoff,sample_cutoff,'_reactome.db_',id.list[[1]] ,'_norm_var_',ct1,'.RData')) #pvalue0.050.07_reactome.db_human_norm_var_Astrocyte.RData
  selected_p1 = sampleScore$selected_pathways_sig
  message('selected_p1:...',length(selected_p1))
  for(ct2 in ct2.lst)
  {
    print(ct2)
  s2 = rownames(subset(annot.2, annot.2$type == ct2))
  print(length(s2))
  load(paste0(paste0(expr_dir, region,'03_sampleScoring/'),filter_opt,filter_cutoff,sample_cutoff,'_reactome.db_',id.list[[2]] ,'_norm_var_',ct2,'.RData')) #pvalue0.050.07_reactome.db_human_norm_var_Astrocyte.RData
  selected_p2 = sampleScore$selected_pathways_sig
  message('selected_p2:...',length(selected_p2))
  rankP_ct =read.table(file = paste0(paste0(expr_dir, region,'07_rankPathways/'),'rankP_',paste0(ct1,'_',ct2),'_all.txt'), header=T, sep="\t", row.names = 1)
  rankP_ct = rankP_ct[order(rankP_ct$cor,decreasing = TRUE),]
  for(p in rankP_ct$p_id)   #rownames(merged_score12)
  {
  s1_p=as.numeric(merged_score12[p,s1])
  s2_p=as.numeric(merged_score12[p,s2])
  s1_p_avg=mean(s1_p, na.rm = TRUE)
  s2_p_avg=mean(s2_p, na.rm = TRUE)
  s1_p_med=median(s1_p, na.rm = TRUE)
  s2_p_med=median(s2_p, na.rm = TRUE)
  s1_p_proj=as.numeric(pcCC[p,s1])
  s2_p_proj=as.numeric(pcCC[p,s2])
  #A large t-score, or t-value, indicates that the groups are different while a small t-score indicates that the groups are similar.
  tt =t.test(s1_p, s2_p, var.equal=T)
  tf =t.test(s1_p, s2_p, var.equal=F)
  tt_proj =t.test(s1_p_proj, s2_p_proj, var.equal=T)
  tf_proj =t.test(s1_p_proj, s2_p_proj, var.equal=F)
  c_p_org_proj_info[idx,'ct_id']=paste0(ct1,':',ct2)
  c_p_org_proj_info[idx,'p_id']=p
  c_p_org_proj_info[idx,'p_level']=ifelse(p %in% merged_score12.Pinfo$p_id,merged_score12.Pinfo[merged_score12.Pinfo$p_id == p,'p_level'],'none')
  c_p_org_proj_info[idx,'p_name']=ifelse(p %in% merged_score12.Pinfo$p_id,merged_score12.Pinfo[merged_score12.Pinfo$p_id == p,'p_name'],'none')
  c_p_org_proj_info[idx,'p_all_cor']=rankP[rankP$p_id == p,'cor']
  c_p_org_proj_info[idx,'p_ct_cor']=rankP_ct[rankP_ct$p_id == p,'cor']
  c_p_org_proj_info[idx,'p_selected_1']=ifelse(p %in% selected_p1, 'yes', 'no')
  c_p_org_proj_info[idx,'p_selected_2']=ifelse(p %in% selected_p2, 'yes', 'no')
  c_p_org_proj_info[idx,'p_selected_common']=ifelse(p %in% intersect(selected_p1,selected_p2) , 'yes', 'no')
  c_p_org_proj_info[idx,'p_ttest_gsea']=toString(tt$statistic)
  c_p_org_proj_info[idx,'p_ttpval_gsea']=round(tt$p.value,2)
  c_p_org_proj_info[idx,'p_ftest_gsea']=toString(tf$statistic)
  c_p_org_proj_info[idx,'p_ftpval_gsea']=round(tf$p.value,2)
  c_p_org_proj_info[idx,'p_ttest_CPproj']=toString(tt_proj$statistic)
  c_p_org_proj_info[idx,'p_ttpval_CPproj']=round(tt_proj$p.value,2)
  c_p_org_proj_info[idx,'p_ftest_CPproj']=toString(tf_proj$statistic)
  c_p_org_proj_info[idx,'p_ftpval_CPproj']=round(tf_proj$p.value,2)
  c_p_org_proj_info[idx,'p_avg_gsea_1']=s1_p_avg
  c_p_org_proj_info[idx,'p_avg_gsea_2']=s2_p_avg
  c_p_org_proj_info[idx,'p_med_gsea_1']=s1_p_med
  c_p_org_proj_info[idx,'p_med_gsea_2']=s2_p_med
  c_p_org_proj_info[idx,'p_diff_gsea']=ifelse((sign(s1_p_avg) * (sign(s2_p_avg))) ==1, abs(abs(s1_p_avg)-abs(s2_p_avg)),(abs(s1_p_avg)+abs(s2_p_avg)))
  c_p_org_proj_info[idx,'p_diffSign_gsea']=ifelse((sign(s1_p_avg) * (sign(s2_p_avg))) ==1,sign(s1_p_avg) ,0)

  
  idx=idx+1
  }
  }
  # print('c_p_org_proj_info')
  # print(dim(c_p_org_proj_info))
  res_dir =  paste0(expr_dir,region, '07_rankPathways/')
  write.table(c_p_org_proj_info,file = paste0(res_dir,'c_p_GSEAorg_GSEAproj_info_',ct1,'.txt'), sep = "\t", quote = TRUE, row.names = TRUE, col.names = T)

 }




