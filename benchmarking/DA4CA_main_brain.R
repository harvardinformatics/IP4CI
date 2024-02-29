# ---
#   title: "DA4CA-vignette"
# output: rmarkdown::html_vignette
# vignette: >
#   %\VignetteIndexEntry{my-vignette}
# %\VignetteEngine{knitr::rmarkdown}
# %\VignetteEncoding{UTF-8}
# ---
#
#   #```{r, include = FALSE}
#   knitr::opts_chunk$set(
#     collapse = TRUE ,
#     comment = "#>"
#   )
# #```

################################  main run #####################################
#```{r clean_env}
rm(list = ls())
#```


# run all required pkgs
#```{r run_packages}
run_packages()
#```
################################################################################
# folders & variables needed by the main scripts
#```{r var}
expr_dir = '/Users/rao198/OneDrive\ -\ Harvard\ University/Informatics/doc/paper/DA4LT/V1/analysis/brain/'
data_dir = '/Users/rao198/OneDrive\ -\ Harvard\ University/Informatics/doc/paper/DA4LT/V1/analysis/brain/data/'

# variables needed by the main scripts

id1='mouse'
id2='human'
obj1_fname='mouse_sc.RDS'
obj2_fname='human_sc.RDS'
typeColumnToUse.list=c('cell.names','cell.names') # later replaced by 'type'
assaySlotToUse.list=list('RNA','RNA')#list('integrated','SCT')
processObjOpt = F
genes2keep='commonVarG'

dataType = 'norm' # use normalized obj data slot
geneType = 'var' # use variable genes
pathwaydatabase = 'reactome.db'

filter_opt = 'pvalue'
filter_cutoff = 0.05
sample_cutoff = 0.07

topRanked = 10

print('hello end of var')
#```
################################################################################
# 00. preprocess data
#```{r 00.preprocess data}
data_dir = data_dir
res_dir =  paste0(expr_dir,'00_preprocessData/')
createDir(res_dir)
# run_00_preprocessData(data_dir,res_dir, id.list=c(id1,id2),obj.fname.list =c(obj1_fname,obj2_fname),typeColumnToUse.list,assaySlotToUse.list)
#```
################################################################################
# 01.process data
#```{r 01.process data}
data_dir = paste0(expr_dir,'00_preprocessData/') # or paste0(expr_dir,'/00_preprocessData/')
res_dir =  paste0(expr_dir,'01_processData/')
createDir(res_dir)
processObjOpt = F
# processedData = run_01_processData(data_dir,res_dir, id.list=c(id1,id2),obj.fname.list =c(obj1_fname,obj2_fname),processObjOpt,genes2keep)
# obj1 = processedData$obj1
# obj2 = processedData$obj2
# annot = processedData$annot
#```
################################################################################

# 02.run gsea
#```{r 02.run gsea}
data_dir =  paste0(expr_dir, '01_processData/')
res_dir =  paste0(expr_dir, '02_runGSEA/')
dir.create(res_dir)

# obj1
# obj1 = readRDS(paste0(data_dir,genes2keep,obj1_fname)) #human_sc_commonVarG.RDS(res_dir,genes2keep,obj1_fname)
# print(obj1)
# obj2 = readRDS(paste0(data_dir,genes2keep,obj2_fname))
# print(obj2)
# if(dataType == 'norm' & geneType == 'var')
# {
#   df1=obj1@assays$DA4CA.1@data
#   annot1 = obj1@meta.data
#   df2=obj2@assays$DA4CA.2@data
#   annot2 = obj2@meta.data
# }
#run_02_gsea(data_dir,res_dir, id1,df1,annot1,dataType,geneType,pathwaydatabase)# slow
# for(ct in 'macrophage') #unique(obj1@meta.data$type)
#   run_02_gseaFast(data_dir,res_dir,id1,df1,annot1,dataType,geneType,pathwaydatabase,ct) # fast
# #run_02_gsea(data_dir,res_dir, id2,df2,annot2,dataType,geneType,pathwaydatabase)
# for(ct in 'macrophage') #unique(obj2@meta.data$type)
#   run_02_gseaFast(data_dir,res_dir,id2,df2,annot2,dataType,geneType,pathwaydatabase,ct) # fast
#```
################################################################################

# 03 score samples
#```{r 03.score samples}
data_dir =  paste0(expr_dir, '02_runGSEA/')
res_dir =  paste0(expr_dir, '03_sampleScoring/')
dir.create(res_dir)
#d1
# pattern_s = paste0(pathwaydatabase,'_',id1,'*.RData')
# print(pattern_s)
# scoringFiles_lst =  list.files(data_dir, pattern=glob2rx(pattern_s))
# print(scoringFiles_lst)
# scoringFiles_lst = c('reactome.db_mouse_norm_var_Astrocyte.RData','reactome.db_mouse_norm_var_Chandelier.RData')
# for(f in scoringFiles_lst)
# {
#   load(paste0(data_dir,f))
#   sampleScore=run_03_sampleScoring(GSEA_res,filter_opt,filter_cutoff,sample_cutoff)
#   print(names(sampleScore))
#   save(sampleScore,file =paste0(res_dir,filter_opt,filter_cutoff,sample_cutoff,'_',f))
# }
# d2
# pattern_s = paste0(pathwaydatabase,'_',id2,'*.RData')
# print(pattern_s)
# scoringFiles_lst =  list.files(data_dir, pattern=glob2rx(pattern_s))
# print(scoringFiles_lst)
# scoringFiles_lst = c('reactome.db_human_norm_var_Astrocyte.RData','reactome.db_human_norm_var_Chandelier.RData')
# for(f in scoringFiles_lst)
# {
#   load(paste0(data_dir,f))
#   print(GSEA_res)
#   sampleScore=run_03_sampleScoring(GSEA_res,filter_opt,filter_cutoff,sample_cutoff)
#   print(names(sampleScore))
#   save(sampleScore,file =paste0(res_dir,filter_opt,filter_cutoff,sample_cutoff,'_',f))
# }
#```
################################################################################

# 04 merge datasets
#```{r 04.merge datasets}
data_dir =  paste0(expr_dir, '03_sampleScoring/')
res_dir =  paste0(expr_dir, '04_mergeTwoDatasets/')
dir.create(res_dir)
annot = read.table(file = '/Users/rao198/OneDrive\ -\ Harvard\ University/Informatics/Analysis/DA4LT/brain/04_mergeTwoDatasets/reactome.db_mouse_human_norm_var_GSEAS_pvalue_FALSE_0.05_0.07_common_Merged_annot.txt', header=T, sep="\t", row.names = 1)
print(dim(annot))

# id1
# pattern_s = paste0(filter_opt,filter_cutoff,sample_cutoff,'*',id1,'*.RData')
# print(pattern_s)
# mergingFiles_lst =  list.files(data_dir, pattern=glob2rx(pattern_s))
# print(mergingFiles_lst)
# merged_score1 = mergeMultiFiles(paste0(data_dir,mergingFiles_lst))
# write.table(merged_score1, file = paste0(res_dir,id1,'_Merged.txt'), sep = "\t", quote = FALSE, row.names = TRUE)
# merged_score1 = read.table(file = '/Users/rao198/OneDrive\ -\ Harvard\ University/Informatics/Analysis/DA4LT/brain/04_mergeTwoDatasets/reactome.db_mouse_norm_var_GSEAS_pvalue_FALSE_0.05_0.07_Merged.txt', header=T, sep="\t", row.names = 1)
# print(dim(merged_score1))
# merged_score1.Pinfo = getPathwaysInfo(rownames(merged_score1))
# write.table(merged_score1.Pinfo, file =paste0(res_dir,id1,'_PathwaysInfo.txt'), sep = "\t", quote = FALSE, row.names = TRUE)
# merged_score1.Pinfo = read.table(file = paste0(res_dir,id1,'_PathwaysInfo.txt'), header=T, sep="\t", row.names = 1)
# print(dim(merged_score1.Pinfo))
# annot.s = annot[colnames(merged_score1),]
# p.df.s= merged_score1.Pinfo[merged_score1.Pinfo$p_id %in% rownames(merged_score1),]
# print(dim(p.df.s))
# df.s = merged_score1[rownames(merged_score1) %in% rownames(p.df.s),]
# print(dim(df.s))
# visData(df.s,annot.s,p.df.s,paste0(res_dir,id1,'_visData.pdf'))
# nVarF=nrow(merged_score12)
# df.s.obj = df2obj(df.s,annot.s,nVarF)
# print(df.s.obj@meta.data)
# dim_name = 'umap'
# visObj(df.s.obj,dim_name,paste0(res_dir,id1,'_vis_varF_',nVarF,'.pdf'))
# # id2
# pattern_s = paste0(filter_opt,filter_cutoff,sample_cutoff,'*',id2,'*.RData')
# print(pattern_s)
# mergingFiles_lst =  list.files(data_dir, pattern=glob2rx(pattern_s))
# print(mergingFiles_lst)
# merged_score2 = mergeMultiFiles(paste0(data_dir,mergingFiles_lst))
# write.table(merged_score2, file = paste0(res_dir,id2,'_Merged.txt'), sep = "\t", quote = FALSE, row.names = TRUE)
# merged_score2 = read.table(file = '/Users/rao198/OneDrive\ -\ Harvard\ University/Informatics/Analysis/DA4LT/brain/04_mergeTwoDatasets/reactome.db_human_norm_var_GSEAS_pvalue_FALSE_0.05_0.07_Merged.txt', header=T, sep="\t", row.names = 1)
# print(dim(merged_score2))
# merged_score2.Pinfo = getPathwaysInfo(rownames(merged_score2))
# write.table(merged_score2.Pinfo, file =paste0(res_dir,id2,'_PathwaysInfo.txt'), sep = "\t", quote = FALSE, row.names = TRUE)
# annot.s = annot[colnames(merged_score2),]
# p.df.s= merged_score2.Pinfo[merged_score2.Pinfo$p_id %in% rownames(merged_score2),]
# print(dim(p.df.s))
# df.s = merged_score2[rownames(merged_score2) %in% rownames(p.df.s),]
# print(dim(df.s))
# visData(df.s,annot.s,p.df.s,paste0(res_dir,id2,'_visData.pdf'))
# df.s.obj = df2obj(df.s,annot.s)
# print(df.s.obj@meta.data)
dim_name = 'umap'
# visObj(df.s.obj,dim_name,paste0(res_dir,id2,'_vis.pdf'))
# # id1 & id2
# merged_score12=run_04_mergeTwoDatasets(merged_score1,merged_score2)
merged_score12 = read.table(file = '/Users/rao198/OneDrive\ -\ Harvard\ University/Informatics/Analysis/DA4LT/brain/04_mergeTwoDatasets/reactome.db_mouse_human_norm_var_GSEAS_pvalue_FALSE_0.05_0.07_common_Merged.txt', header=T, sep="\t", row.names = 1)
print(dim(merged_score12))
# merged_score12.Pinfo = getPathwaysInfo(rownames(merged_score12))
# write.table(merged_score12.Pinfo, file =paste0(res_dir,id1,'_',id2,'_PathwaysInfo.txt'), sep = "\t", quote = FALSE, row.names = TRUE)
merged_score12.Pinfo = read.table(file = paste0(res_dir,id1,'_',id2,'_PathwaysInfo.txt'), header=T, sep="\t", row.names = 1)
print(dim(merged_score12.Pinfo))
annot.s = annot[colnames(merged_score12),]
p.df.s= merged_score12.Pinfo[merged_score12.Pinfo$p_id %in% rownames(merged_score12),]
print(dim(p.df.s))
df.s = merged_score12[rownames(merged_score12) %in% rownames(p.df.s),]
print(dim(df.s))
visData(df.s,annot.s,p.df.s,paste0(res_dir,id1,'_',id2,'_visDataAvg.pdf'))
# nVarF=nrow(merged_score12)
# df.s.obj = df2obj(df.s,annot.s,nVarF)
# # print(nrow(df.s.obj))
# saveRDS(df.s.obj,paste0(res_dir,id1,'_',id2,'_obj_varF_',nVarF,'.RDS'))
# # df.s.obj = readRDS(paste0(res_dir,id1,'_',id2,'_obj_varF_',nVarF,'.RDS'))
# # print(df.s.obj@assays[["RNA"]]@var.features)
# dim_name = 'umap'
# visObj(df.s.obj,dim_name,paste0(res_dir,id1,'_',id2,'_vis_varF_',nVarF,'.pdf'))
# annot = annot[colnames(merged_score12),]
# write.table(annot, file = paste0(res_dir,'annot_Merged.txt'), sep = "\t", quote = FALSE, row.names = TRUE)
# Feature selection

#```
################################################################################

# 05 CCA
#```{r  05.CCA}
data_dir =  paste0(expr_dir, '04_mergeTwoDatasets/')
res_dir =  paste0(expr_dir, '05_CCA/')
dir.create(res_dir)

# 1. VF
ObjType = 'merged' # 'merged = from gsea , take x,y' 'single = each dataset x/y has FV'

# s1 = rownames(subset(annot, annot$dataset_name == id1))
# s2 = rownames(subset(annot, annot$dataset_name == id2))
# if(ObjType == 'merged')
# {
#   x=merged_score12[,s1]
#   y=merged_score12[,s2]
# }
# if(ObjType == 'single')
# {
#   x=merged_score1[,s1]
#   y=merged_score2[,s2]
# }

# x=merged_score1[,s1]
# y=merged_score2[,s2]
nVarF ='all'
if(nVarF == 'all')
{
    # obj.x =df2obj(x,annot[s1,],nrow(x))
    # print(dim(obj.x))
    # saveRDS(obj.x,paste0(res_dir,id1,'_',ObjType,'obj_varF_',nVarF,'.RDS'))
    # obj.y =df2obj(y,annot[s2,],nrow(y))
    # print(dim(obj.y))
    # saveRDS(obj.y,paste0(res_dir,id2,'_',ObjType,'obj_varF_',nVarF,'.RDS'))

  # obj.x = readRDS(paste0(res_dir,id1,'_',ObjType,'obj_varF_',nVarF,'.RDS'))
  # print(obj.x)
  # dim_name = 'umap'
  # visObj(obj.x,dim_name,paste0(res_dir,id1,'_vis_varF_',nVarF,'.pdf'))

  # obj.y = readRDS(paste0(res_dir,id2,'_',ObjType,'obj_varF_',nVarF,'.RDS'))
  # print(obj.y)
  # dim_name = 'umap'
  # visObj(obj.y,dim_name,paste0(res_dir,id2,'_vis_varF_',nVarF,'.pdf'))


}
if(nVarF != 'all')
{
  # obj.x =df2obj(x,annot[s1,],nVarF)
  # print(dim(obj.x))
  # saveRDS(obj.x,paste0(res_dir,id1,'_',ObjType,'obj_varF_',nVarF,'.RDS'))
  # obj.y =df2obj(y,annot[s2,],nVarF)
  # print(dim(obj.y))
  # saveRDS(obj.y,paste0(res_dir,id2,'_',ObjType,'obj_varF_',nVarF,'.RDS'))

  obj.x = readRDS(paste0(res_dir,id1,'_',ObjType,'obj_varF_',nVarF,'.RDS'))
  print(obj.x)
  dim_name = 'umap'
  visObj(obj.x,dim_name,paste0(res_dir,id1,'_vis_varF_',nVarF,'.pdf'))

  obj.y = readRDS(paste0(res_dir,id2,'_',ObjType,'obj_varF_',nVarF,'.RDS'))
  print(obj.y)
  dim_name = 'umap'
  visObj(obj.y,dim_name,paste0(res_dir,id2,'_vis_varF_',nVarF,'.pdf'))

}

# 2 CCA
num.cc = 30
selected_dims_opt = 'none' #cor'
selected_dims_opt_cutoff =0
# cca_res =run_05_CCA(obj.x,obj.y,num.cc)
# saveRDS(cca_res,paste0(res_dir,ObjType,'_',num.cc,'_',nVarF,'_1.RDS') )
# cca_res= readRDS(paste0(res_dir,ObjType,'_',num.cc,'_',nVarF,'_1.RDS') )
# print(unique(cca_res@meta.data$dataset_name))
# dim_name = 'umap_cca'
# print(length(rownames(cca_res)))
# #cca_res = Seurat::ScaleData(cca_res)
# cca_res = RunPCA(cca_res,verbose=F,npcs = 10,features = rownames(cca_res))
# print(cca_res)
# cca_res = RunUMAP(cca_res,dims = 1:10, n.components = 2,verbose = F,reduction.name = dim_name)
# print(cca_res)
# visObj(cca_res,dim_name,paste0(res_dir,ObjType,'_',num.cc,'_',nVarF,'_',selected_dims_opt,selected_dims_opt_cutoff,'_obj_',nVarF,'_',dim_name,'.pdf'))
# dim_name = 'cca'
# visObj(cca_res,dim_name,paste0(res_dir,ObjType,'_',num.cc,'_',nVarF,'_',selected_dims_opt,selected_dims_opt_cutoff,'_obj_',nVarF,'_',dim_name,'.pdf'))

# cca_res =processCCAres_allCT(cca_res, x, obj.x,y,obj.y) # add stability over cc
# saveRDS(cca_res,paste0(res_dir,ObjType,'_',num.cc,'_',nVarF,'_','none','0','.RDS'))
# cca_res = readRDS(paste0(res_dir,ObjType,'_',num.cc,'_',nVarF,'_','none','0','.RDS'))
# if(selected_dims_opt != 'none')
# {
#   cca_res = selectCCA(cca_res,selected_dims_opt,selected_dims_opt_cutoff)
#   saveRDS(cca_res,paste0(res_dir,ObjType,'_',num.cc,'_',nVarF,'_',selected_dims_opt,selected_dims_opt_cutoff,'.RDS') )
# }

# 3 visualize after cca results
# cca_res= readRDS(paste0(res_dir,ObjType,'_',num.cc,'_',nVarF,'_',selected_dims_opt,selected_dims_opt_cutoff,'.RDS') )
# embed2use = 'uembedvembed'
# if(embed2use == 'uv')
# {
#   xy = t(rbind(cca_res$u,cca_res$v))
#   rownames(xy)=paste0('cc',rownames(xy))
#   xy.annot = rbind(annot[rownames(cca_res$u),],annot[rownames(cca_res$v),])
# }
# if(embed2use == 'ulvl')
# {
#   xy = t(rbind(cca_res$ul,cca_res$vl))
#   rownames(xy)=paste0('cc',rownames(xy))
#   xy.annot = rbind(annot[rownames(cca_res$ul),],annot[rownames(cca_res$vl),])
# }
# if(embed2use == 'uembedvembed')
# {
#   xy = t(rbind((cca_res$u %*% t(cca_res$pembd)),(cca_res$v %*% t(cca_res$pembd))) )
#   rownames(xy)=paste0('cc',rownames(xy))
#   xy.annot = rbind(annot[rownames(cca_res$u),],annot[rownames(cca_res$v),])
# }
# obj.xy =df2obj(xy,xy.annot,0)
# saveRDS(obj.xy,paste0(res_dir,ObjType,'_',num.cc,'_',nVarF,'_',selected_dims_opt,selected_dims_opt_cutoff,'_obj_',nVarF,'2.RDS' ))
# # obj.xy =  readRDS(paste0(res_dir,ObjType,'_',num.cc,'_',nVarF,'_',selected_dims_opt,selected_dims_opt_cutoff,'_obj.RDS' ))
dim_name = 'umap'
# visObj(obj.xy,dim_name,paste0(res_dir,ObjType,'_',num.cc,'_',nVarF,'_',selected_dims_opt,selected_dims_opt_cutoff,'_obj_',nVarF,'_',dim_name,'_',embed2use,'.pdf'))
# nVarF=500
# obj.xy =df2obj(xy,xy.annot,nVarF)
# saveRDS(obj.xy,paste0(res_dir,ObjType,'_',num.cc,'_',nVarF,'_',selected_dims_opt,selected_dims_opt_cutoff,'_obj_',nVarF,'.RDS' ))
# # obj.xy =  readRDS(paste0(res_dir,ObjType,'_',num.cc,'_',nVarF,'_',selected_dims_opt,selected_dims_opt_cutoff,'_obj.RDS' ))
dim_name = 'umap'
# visObj(obj.xy,dim_name,paste0(res_dir,ObjType,'_',num.cc,'_',nVarF,'_',selected_dims_opt,selected_dims_opt_cutoff,'_obj_',nVarF,'.pdf'))

# 1. run Seurat DataIntegrate
# integrated_res =run_05_dataIntegration(obj.x,obj.y)
# saveRDS(integrated_res,paste0(res_dir,ObjType,'_integrated','_varF_',nVarF,'_1.RDS') )
# integrated_res = readRDS(paste0(res_dir,ObjType,'_integrated','_varF_',nVarF,'_1.RDS'))
# visObj(integrated_res,dim_name,paste0(res_dir,ObjType,'_integrated_',nVarF,'_obj.pdf'))
dim_name = 'umap'
# pdf(paste0(res_dir,'merged2_integrated','_obj_varF_',nVarF,'_1_vis_updated.pdf'), width = 12, height = 12)
#           print(DimPlot(integrated_res, reduction = "pca", group.by = "type"))
#           print(DimPlot(integrated_res, reduction = "TSNE_on_CCA", group.by = "type"))
#           print(DimPlot(integrated_res, reduction = "UMAP_on_CCA", group.by = "type", label = TRUE,repel = TRUE))
#           print(DimPlot(integrated_res, reduction = "UMAP_on_CCA", split.by   ="dataset_name", label = TRUE,repel = TRUE))
#           print(DimPlot(integrated_res, reduction = "UMAP_on_CCA", group.by   ="dataset_name", label = TRUE,repel = TRUE))
# dev.off()
# integrated_res =processIntegratedres_allCT(integrated_res, x, obj.x,y,obj.y)
# saveRDS(integrated_res,paste0(res_dir,ObjType,'_integrated','_obj_varF_',nVarF,'.RDS') )
# integrated_res = readRDS(paste0(res_dir,ObjType,'_integrated','_obj_varF_',nVarF,'.RDS'))


# old
# old cca res file
# cca_res = readRDS("/Users/rao198/Documents/DataTransfer/data/Intg_cca_30_ccaObj.RDS")
# cca_res =processCCAres_allCT(cca_res, x, obj.x,y,obj.y)
# saveRDS(cca_res,paste0(res_dir,'merged_cca.RDS') )
# cca_res = readRDS(paste0(res_dir,'merged_cca.RDS'))
# old
################################################################################

# 06 CCA similarity
#```{r  06.CCA similarity}
data_dir =  paste0(expr_dir, '05_CCA/')
res_dir =  paste0(expr_dir, '06_similarityCalculation/')
dir.create(res_dir)
# run_06_similarityCalculation(cca_res,annot,id.list=c(id1,id2),ct='all',res_dir)

# in case of Dataintegrate
# run_06_similarityCalculation(integrated_res,integrated_res$meta.data,id.list=c(id1,id2),ct='all',paste0(res_dir,ObjType,'_integrated_',nVarF))
# run_06_similarityCalculationIntegration(integrated_res,integrated_res$meta.data,id.list=c(id1,id2),ct='all',paste0(res_dir,ObjType,'_integratedCorrected_',nVarF))

# cca_res = readRDS(paste0(data_dir,'merged_cca_',nVarF,'_selected.RDS') )
# run_06_similarityCalculation(cca_res,annot,id.list=c(id1,id2),ct='all',paste0(res_dir,nVarF,'_selected_'))

# in case of runCCA
# cca_res = readRDS(paste0(data_dir,ObjType,'_',num.cc,'_',nVarF,'_',selected_dims_opt,selected_dims_opt_cutoff,'.RDS'))
# run_06_similarityCalculation(cca_res,annot,id.list=c(id1,id2),ct='all',paste0(res_dir,ObjType,'_',num.cc,'_',nVarF,'_',selected_dims_opt,selected_dims_opt_cutoff))
#```
################################################################################

# 07 rank pathways
#```{r  07.rank pathways}
data_dir =  paste0(expr_dir, '05_CCA/')
res_dir =  paste0(expr_dir, '07_rankPathways/')
dir.create(res_dir)

w_opt =F  # no need to weight the pathways scores
#a. rank all pathways
# rankP = run_07_rankPathways(cca_res,merged_score12.Pinfo,w_opt)
# write.table(rankP,file = paste0(res_dir,'rankP_allCT',.txt'), sep = "\t", quote = TRUE, row.names = TRUE, col.names = T)
# rankP = read.table(file = paste0(res_dir,'rankP_allCT_new.txt'), header=T, sep="\t", row.names = 1)
#plot
# for(sim_col in c('cor')) #'reg',
# {
#   rankP.tmp = rankP
#   print(head(rankP.tmp))
#   colnames(rankP.tmp)[which(names(rankP.tmp) == sim_col)] <- "p_sim"
#   rankP.tmp= rankP.tmp[order(rankP.tmp$p_sim, decreasing = T),]
#   rankP.tmp = rankP.tmp[1:topRanked,c('p_name','p_level','p_id','p_sim')]
#   rownames(rankP.tmp) = rankP.tmp$p_id
#   # print(head(rankP.tmp))
#   # print(rownames(cca_res$x)[1:10])
#   # print(cca_res$x[(rankP.tmp$p_id),1:10])
#   # print(cca_res$y[(rankP.tmp$p_id),1:10])
#   df.s=as.data.frame(rbind(t(cca_res$x[(rankP.tmp$p_id),]),t(cca_res$y[(rankP.tmp$p_id),])))
#   print(head(df.s))
#   annot.s = annot[rownames(df.s),]
#   title = paste0('heatmap of gsea for pathways','\n','across all cell-types between two datasets:')
#   print(title)
#   plotRankingP(df.s ,rankP.tmp$p_level,rankP.tmp$p_sim,annot.s$dataset_name,title,paste0(res_dir,'rankP_',sim_col,'_all_new.pdf'))
# }
#b. rank per cell-type pathways
# CTs = c("Astrocyte","Chandelier")#unique(annot$type)
# for(ct in CTs)
# {
#   print(ct)
#   # get samples
#   s1 = rownames(subset(annot, annot$dataset_name == id1 & annot$type == ct))
#   s2 = rownames(subset(annot, annot$dataset_name == id2 & annot$type == ct))
#   print(length(s1))
#   print(length(s2))
#   x.s = x[,s1]
#   y.s = y[,s2]
#   cc_res_ct= processCCAres_perCT(cca_res,x.s,y.s)
# saveRDS(cc_res_ct,paste0(res_dir,ct,'_merged_cca.RDS') )
#   cc_res_ct = readRDS(paste0(res_dir,ct,'_merged_cca.RDS') )
#   w_opt = F
# rankP_ct = run_07_rankPathways(cc_res_ct,merged_score12.Pinfo,w_opt)
# write.table(rankP_ct,file = paste0(res_dir,'rankP_',ct,'.txt'), sep = "\t", quote = TRUE, row.names = TRUE, col.names = T)
# #plot
# for(sim_col in c('cor'))
# {
#   rankP.tmp = rankP_ct
#   colnames(rankP.tmp)[which(names(rankP.tmp) == sim_col)] <- "p_sim"
#   rankP.tmp= rankP.tmp[order(rankP.tmp$p_sim, decreasing = T),]
#   rankP.tmp = rankP.tmp[1:topRanked,c('p_name','p_level','p_id','p_sim')]
#   rownames(rankP.tmp) = rankP.tmp$p_id
#   df.s=as.data.frame(rbind(t(cc_res_ct$x[(rankP.tmp$p_id),]),t(cc_res_ct$y[(rankP.tmp$p_id),])))
#   print(head(df.s))
#   annot.s = annot[rownames(df.s),]
#   title = paste0('heatmap of gsea for pathways','\n','across all cell-types between two datasets:')
#   print(title)
#   plotRankingP(df.s ,rankP.tmp$p_level,rankP.tmp$p_sim,annot.s$dataset_name,title,paste0(res_dir,'rankP_',sim_col,'_',ct,'.pdf'))
# }
# }

# c. combine rankingP
# fname.lst = paste0(res_dir,c('rankP_Astrocyte.txt','rankP_Chandelier.txt'))
# ct.lst = c("Astrocyte","Chandelier")
# plotCombineRankingPathways(fname.lst,ct.lst,res_dir)
################################################################################


# 07 rank pathways
#```{r  07.rank pathways}
data_dir =  paste0(expr_dir, '05_CCA/')
res_dir =  paste0(expr_dir, '07_rankPathways/')
dir.create(res_dir)
w_opt=F
#a. rank all pathways
# rankP = run_07_rankPathways(cca_res,merged_score12.Pinfo,w_opt)
# write.table(rankP,file = paste0(res_dir,'rankP_allCT.txt'), sep = "\t", quote = TRUE, row.names = TRUE, col.names = T)
#plot
# for(sim_col in c('reg','cor'))
# {
#   rankP.tmp = rankP
#   colnames(rankP.tmp)[which(names(rankP.tmp) == sim_col)] <- "p_sim"
#   rankP.tmp= rankP.tmp[order(rankP.tmp$p_sim, decreasing = T),]
#   rankP.tmp = rankP.tmp[1:topRanked,c('p_name','p_level','p_id','p_sim')]
#   rownames(rankP.tmp) = rankP.tmp$p_id
#   df.s=as.data.frame(rbind(t(cca_res$x[(rankP.tmp$p_id),]),t(cca_res$y[(rankP.tmp$p_id),])))
#   print(head(df.s))
#   annot.s = annot[rownames(df.s),]
#   title = 'heatmap of gsea for pathways across all cell-types between two datasets:'
#   print(title)
#   plotRankingP(df.s ,rankP.tmp$p_level,rankP.tmp$p_sim,annot.s$dataset_name,title,paste0(res_dir,'rankP_',sim_col,'_all.pdf'))
# }
#b. rank per cell-type pathways
# CTs = c('Chandelier')#unique(annot$type)
# for(ct in CTs)
# {
#   print(ct)
# #get samples
# s1 = rownames(subset(annot, annot$dataset_name == id1 & annot$type == ct))
# s2 = rownames(subset(annot, annot$dataset_name == id2 & annot$type == ct))
# print(length(s1))
# print(length(s2))
# x.s = cca_res$x[,s1]
# y.s = cca_res$y[,s2]
# cc_res_ct= processCCAres_perCT(cca_res,x.s,y.s)
# saveRDS(cc_res_ct,paste0(res_dir,ct,'_merged_cca.RDS') )
# cc_res_ct = readRDS(paste0(res_dir,ct,'_merged_cca.RDS') )
# print(names(cc_res_ct))
# load(paste0(paste0(expr_dir, '03_sampleScoring/'),filter_opt,filter_cutoff,sample_cutoff,'_reactome.db_',id1,'_norm_var_',ct,'.RData')) #pvalue0.050.07_reactome.db_human_norm_var_Astrocyte.RData
# print(length(sampleScore$selected_pathways_sig))
# selected_p1 = sampleScore$selected_pathways_sig
# load(paste0(paste0(expr_dir, '03_sampleScoring/'),filter_opt,filter_cutoff,sample_cutoff,'_reactome.db_',id2,'_norm_var_',ct,'.RData'))
# print(length(sampleScore$selected_pathways_sig))
# selected_p2 = sampleScore$selected_pathways_sig
# rankP_ct = run_07_rankPathways(cc_res_ct,merged_score12.Pinfo,w_opt)
# write.table(rankP_ct,file = paste0(res_dir,'rankP_',ct,'_all.txt'), sep = "\t", quote = TRUE, row.names = TRUE, col.names = T)
# rankP_ct = rankP_ct[rankP_ct$p_id %in% intersect(selected_p1,selected_p2),]
# write.table(rankP_ct,file = paste0(res_dir,'rankP_',ct,'_selected.txt'), sep = "\t", quote = TRUE, row.names = TRUE, col.names = T)
# rankP_ct = read.table(file = paste0(res_dir,'rankP_',ct,'_selected.txt'), header=T, sep="\t", row.names = 1)
#
#     # #plot
#   for(sim_col in c('reg','cor'))
#   {
#     rankP.tmp = rankP_ct
#     colnames(rankP.tmp)[which(names(rankP.tmp) == sim_col)] <- "p_sim"
#     rankP.tmp= rankP.tmp[order(rankP.tmp$p_sim, decreasing = T),]
#     rankP.tmp = rankP.tmp[1:topRanked,c('p_name','p_level','p_id','p_sim')]
#     rownames(rankP.tmp) = rankP.tmp$p_id
#     df.s=as.data.frame(rbind(t(cc_res_ct$x[(rankP.tmp$p_id),]),t(cc_res_ct$y[(rankP.tmp$p_id),])))
#     print(head(df.s))
#     annot.s = annot[rownames(df.s),]
#     title = paste0('heatmap of gsea for pathways', '\n','across all cell-types', '\n','between two datasets:')
#     print(title)
#     plotRankingP(df.s ,rankP.tmp$p_level,rankP.tmp$p_sim,annot.s$dataset_name,title,paste0(res_dir,'rankP_',sim_col,'_',ct,'.pdf'))
#     colnames(cc_res_ct$xscores)= paste0('X:',colnames(cc_res_ct$xscores))
#     colnames(cc_res_ct$yscores)= paste0('Y:',colnames(cc_res_ct$yscores))
#     df.s=as.data.frame(rbind(t(cc_res_ct$xscores[(rankP.tmp$p_id),]),t(cc_res_ct$yscores[(rankP.tmp$p_id),])))
#     print(head(df.s))
#     dataset_name = c(rep(id1,each=ncol(cc_res_ct$xscores)),rep(id2,each=ncol(cc_res_ct$yscores)))
#     annot.s = data.frame(dataset_name)
#     rownames(annot.s) = rownames(df.s)
#     title = paste0('heatmap of CCA scores for pathways', '\n','across all cell-types', '\n','between two datasets:')
#     print(title)
#     plotRankingP(df.s ,rankP.tmp$p_level,rankP.tmp$p_sim,annot.s$dataset_name,title,paste0(res_dir,'rankP_',sim_col,'_',ct,'_XYscores.pdf'))
#   }
# }
# c. combine rankingP
# fname.lst = paste0(res_dir,c('rankP_Chandelier_selected.txt','rankP_Astrocyte_selected.txt'))
# ct.lst = c("Chandelier","Astrocyte")
# plotCombineRankingPathways(fname.lst,ct.lst,res_dir)


################################################################################
########################      NOT TO RUN YET      ##############################
# using loadings from integration not cca ---- NEED TESTING ---
#a. rank all pathways
# re-run old brain
# cca_res = readRDS("/Users/rao198/Documents/DataTransfer/data/Intg_cca_30_ccaObj.RDS")
# cca_res =processCCAres_allCT(cca_res, x, obj.x,y,obj.y)
# saveRDS(cca_res,paste0(res_dir,'merged_cca.RDS') )
# rankP = run_07_rankPathways(integrated_res,merged_score12.Pinfo)
# write.table(rankP,file = paste0(res_dir,'rankP_allCT_integrated.txt'), sep = "\t", quote = TRUE, row.names = TRUE, col.names = T)
#plot
# for(sim_col in c('reg','cor'))
# {
#   rankP.tmp = rankP
#   colnames(rankP.tmp)[which(names(rankP.tmp) == sim_col)] <- "p_sim"
#   rankP.tmp= rankP.tmp[order(rankP.tmp$p_sim, decreasing = T),]
#   rankP.tmp = rankP.tmp[1:topRanked,c('p_name','p_level','p_id','p_sim')]
#   rownames(rankP.tmp) = rankP.tmp$p_id
#   df.s=as.data.frame(rbind(t(integrated_res$x[(rankP.tmp$p_id),]),t(integrated_res$y[(rankP.tmp$p_id),])))
#   print(head(df.s))
#   annot.s = annot[rownames(df.s),]
#   title = 'heatmap of gsea for pathways across all cell-types between two datasets:'
#   print(title)
#   plotRankingP(df.s ,rankP.tmp$p_level,rankP.tmp$p_sim,annot.s$dataset_name,title,paste0(res_dir,'rankP_',sim_col,'_all_integrated.pdf'))
# }
#b. rank per cell-type pathways
# CTs = c('macrophage','quiescent_stellate')#unique(annot$type)
# for(ct in CTs)
# {
#   print(ct)
#   # get samples
#   s1 = rownames(subset(annot, annot$dataset_name == id1 & annot$type == ct))
#   s2 = rownames(subset(annot, annot$dataset_name == id2 & annot$type == ct))
#   print(length(s1))
#   print(length(s2))
#   x.s = x[,s1]
#   y.s = y[,s2]
#   integrated_res_ct= processCCAres_perCT(integrated_res,x.s,y.s)
# saveRDS(integrated_res_ct,paste0(res_dir,ct,'merged_integrated.RDS') )
#   rankP_ct = run_07_rankPathways(integrated_res_ct,merged_score12.Pinfo)
#   write.table(rankP_ct,file = paste0(res_dir,'rankP_',ct,'_integrated.txt'), sep = "\t", quote = TRUE, row.names = TRUE, col.names = T)
#   #plot
#   for(sim_col in c('reg','cor'))
#   {
#     rankP.tmp = rankP_ct
#     colnames(rankP.tmp)[which(names(rankP.tmp) == sim_col)] <- "p_sim"
#     rankP.tmp= rankP.tmp[order(rankP.tmp$p_sim, decreasing = T),]
#     rankP.tmp = rankP.tmp[1:topRanked,c('p_name','p_level','p_id','p_sim')]
#     rownames(rankP.tmp) = rankP.tmp$p_id
#     df.s=as.data.frame(rbind(t(integrated_res_ct$x[(rankP.tmp$p_id),]),t(integrated_res_ct$y[(rankP.tmp$p_id),])))
#     print(head(df.s))
#     annot.s = annot[rownames(df.s),]
#     title = 'heatmap of gsea for pathways across all cell-types between two datasets:'
#     print(title)
#     plotRankingP(df.s ,rankP.tmp$p_level,rankP.tmp$p_sim,annot.s$dataset_name,title,paste0(res_dir,'rankP_',sim_col,'_',ct,'.pdf'))
#   }
# }
# c. combine rankingP
# fname.lst = paste0(res_dir,c('rankP_quiescent_stellate_integrated.txt','rankP_macrophage_integrated.txt'))
# ct.lst = c("quiescent_stellate","macrophage")
# plotCombineRankingPathways(fname.lst,ct.lst,res_dir)
#```
########################      NOT TO RUN YET      ##############################
################################################################################

# 08 rank genes
#```{r  08. rank genes}
data_dir =  paste0(expr_dir, '02_runGSEA/')
res_dir =  paste0(expr_dir, '08_rankGenes/')
if (!dir.exists(res_dir))
  dir.create(res_dir)
print(res_dir)

# CTs = c('Chandelier')#unique(annot$type)
# for(ct in CTs)
#  {
#a. get pathways ranked by 'cor' or 'reg' for topRanked for ct
# rankP.fname =  paste0( paste0(expr_dir, '07_rankPathways/'),'rankP_',ct,'.txt')
# rankP = read.table(file = rankP.fname, header=T, sep="\t", row.names = 1)
#b. get enriched pathway genes members from GSEA results, for human, for mouse
# id1
# gsea.fname = paste0(paste0(expr_dir, '02_runGSEA/'),'reactome.db_',id1,'_norm_var_',ct,'.RData')
# load(gsea.fname)
# expr_logFC1 = expr_logFC
# enrichedG_1 = getEnrichedGenes(GSEA_res$gsea_res,rankP)
# write.table(enrichedG_1,file = paste0(res_dir,id1,'_enrichedG_',ct,'.txt'), sep = "\t", quote = TRUE, row.names = TRUE, col.names = T)
# enrichedG_1 = read.table(file = paste0(res_dir,id1,'_enrichedG_',ct,'.txt'), header=T, sep="\t", row.names = 1)
# id2
# gsea.fname = paste0(paste0(expr_dir, '02_runGSEA/'),'reactome.db_',id2,'_norm_var_',ct,'.RData')
# load(gsea.fname)
# expr_logFC2 = expr_logFC
# enrichedG_2 = getEnrichedGenes(GSEA_res$gsea_res,rankP)
# write.table(enrichedG_2,file = paste0(res_dir,id2,'_enrichedG_',ct,'.txt'), sep = "\t", quote = TRUE, row.names = TRUE, col.names = T)
# enrichedG_2 = read.table(file = paste0(res_dir,id2,'_enrichedG_',ct,'.txt'), header=T, sep="\t", row.names = 1)
# c. get common genes between human and mouse
# commonG = getCommonGenesPerPathway(rankP,enrichedG_1,enrichedG_2)
# print(dim(commonG))
# d. rank genes
# rankG = run_08_rankGenes(expr_logFC1,expr_logFC2,commonG,rankP)
# write.table(rankG,file = paste0(res_dir,'rankG_',ct,'.txt'), sep = "\t", quote = TRUE, row.names = TRUE, col.names = T)
# rankG = read.table(file = paste0(res_dir,'rankG_',ct,'.txt'), header=T, sep="\t", row.names = 1)
# e. plot ranking
# for(p in unique(rankG$p_id))
# {
#   print(p)
#     df.s = rankG[rankG$p_id == p,]
#     name = unique(df.s[df.s$p_id == p,'p_name'])
#     print(name)
#     rownames(df.s)=df.s$g_id
#     df.s = df.s[order(df.s$g_expr_rank, decreasing = F),]
#     title = paste0('heatmap of genes in pathway:\n',name, '\nfor cell-type:',ct)
#     print(title)
#     plotRankingG(t(df.s[, c('g_lfc1','g_lfc2')]) ,df.s$g_rank,df.s$g_diff,title,paste0(res_dir,'rankG_',p,'_',ct,'.pdf'))
# }
# }
#```
################################################################################
# 09 detect Unknown CellType
#```{r  09. detect Unknown CellType}
data_dir =  paste0(expr_dir, '01_processData/')
res_dir =  paste0(expr_dir, '09_detectUnknownCT/')
if (!dir.exists(res_dir))
  dir.create(res_dir)
print(res_dir)
# CTs = c('Astrocyte','Chandelier','Oligo')#unique(annot$type)
# for(ct in CTs)
# {
# print(ct)
#  obj1.s =subset(obj.x,cells=colnames(obj.x)[Idents(obj.x)!=ct])
#  print(obj1.s)
#  run_09_detectUnknownCT(obj1.s,obj.y,annot,id.list=c(id1,id2),ct,res_dir)
# }
#```
#############################################################



