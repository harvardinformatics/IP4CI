#' @title  runDA4CA
#' @description run DA4CA main function to run the analysis using two objects

#' @param expr_dir experiment directory folder name
#' @param id.list id names for the objects
#' @param obj.fname.list file names for the objects
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

#' @examples
#' expr_dir = ''
#' createDir(expr_dir)
#' id.list=c('mouse','human')
#' obj.fname.list=c('mouse_sc.RDS','human_sc.RDS')
#' typeColumnToUse.list=c('cell.names','cell.names')
#' assaySlotToUse.list=list('integrated','SCT')
#' convertG.list = c(T,F)
#' processObjOpt = F
#' convertG.list = c(FALSE,FALSE)
#' normObj.list = c(FALSE,FALSE)
#' genes2keep = 'commonVarG'
#' pathwaydatabase = 'reactome.db'
#' dataType = 'norm'
#' geneType = 'var'
#' filter_opt = 'pvalue'
#' filter_cutoff = 0.05
#' sample_cutoff = 0.07
#' topRanked = 10

#' runDA4CA(expr_dir, id.list,obj.fname.list,convertG.list,normObj.list,processObjOpt,typeColumnToUse.list,assaySlotToUse.list,
#' genes2keep,
#' pathwaydatabase,dataType,geneType,
#' filter_opt,filter_cutoff,sample_cutoff,
#' topRanked)

#' @export
################################################################################

runDA4CA <- function()
{
  print('start running the main function of DA4CA :...')


  ################################  main run #####################################
  # clean_env
  rm(list = ls())
  args = commandArgs(trailingOnly=TRUE)

  # run_packages
  run_packages()
  library(DA4CA)
  ################################################################################
  # 00. preprocess data
  data_dir = paste0(expr_dir,'data/')
  res_dir =  paste0(expr_dir,'00_preprocessData/')
  createDir(res_dir)
  run_00_preprocessData(data_dir,res_dir, id.list,obj.fname.list,convertG.list,normObj.list,processObjOpt,typeColumnToUse.list,assaySlotToUse.list)
  ################################################################################

   # 01.process data
  data_dir = paste0(expr_dir,'00_preprocessData/') # or paste0(expr_dir,'/00_preprocessData/')
  res_dir =  paste0(expr_dir,'01_processData/')
  createDir(res_dir)
  processedData = run_01_processData(data_dir,res_dir, id.list,obj.fname.list,genes2keep)
  obj1 = processedData$obj1
  obj2 = processedData$obj2
  annot = processedData$annot
  ################################################################################

  # 02.run gsea
  data_dir =  paste0(expr_dir, '01_processData/')
  res_dir =  paste0(expr_dir, '02_runGSEA/')
  dir.create(res_dir)
  if(dataType == 'norm' & geneType == 'var')
  {
    df1=obj1@assays$DA4CA.1@data
    annot1 = obj1@meta.data
    df2=obj2@assays$DA4CA.2@data
    annot2 = obj2@meta.data
  }
  run_02_gsea(data_dir,res_dir, id.list[[1]] ,df1,annot1,dataType,geneType,pathwaydatabase)# slow
  #run_02_gseaFast(data_dir,res_dir,id.list[[1]] ,df1,annot1,dataType,geneType,pathwaydatabase,ct) # fast
  run_02_gsea(data_dir,res_dir, id.list[[2]],df2,annot2,dataType,geneType,pathwaydatabase)
  #run_02_gseaFast(data_dir,res_dir,id.list[[2]],df2,annot2,dataType,geneType,pathwaydatabase,ct) # fast
  ################################################################################

   # 03 score samples
  data_dir =  paste0(expr_dir, '02_runGSEA/')
  res_dir =  paste0(expr_dir, '03_sampleScoring/')
  dir.create(res_dir)
  pattern_s = paste0(pathwaydatabase,'_',id.list[[1]] ,'*.RData')
  scoringFiles_lst =  list.files(data_dir, pattern=glob2rx(pattern_s))
  for(f in scoringFiles_lst)
  {
    load(paste0(data_dir,f))
    sampleScore=run_03_sampleScoring(GSEA_res,filter_opt,filter_cutoff,sample_cutoff)
    save(sampleScore,file =paste0(res_dir,filter_opt,filter_cutoff,sample_cutoff,'_',f))
  }
  pattern_s = paste0(pathwaydatabase,'_',id.list[[2]],'*.RData')
  scoringFiles_lst =  list.files(data_dir, pattern=glob2rx(pattern_s))
  for(f in scoringFiles_lst)
  {
    load(paste0(data_dir,f))
    sampleScore=run_03_sampleScoring(GSEA_res,filter_opt,filter_cutoff,sample_cutoff)
    save(sampleScore,file =paste0(res_dir,filter_opt,filter_cutoff,sample_cutoff,'_',f))
  }
  ################################################################################

  # 04 merge datasets
 data_dir =  paste0(expr_dir, '03_sampleScoring/')
  res_dir =  paste0(expr_dir, '04_mergeTwoDatasets/')
  dir.create(res_dir)

  commonCol = intersect(colnames(obj1@meta.data),colnames(obj2@meta.data))
  obj1@meta.data = obj1@meta.data[,commonCol]
  obj2@meta.data = obj2@meta.data[,commonCol]
  annot =rbind(obj1@meta.data,obj2@meta.data)
  write.table(annot, file = paste0(res_dir,'annot_Merged.txt'), sep = "\t", quote = FALSE, row.names = TRUE)
  # id1
  pattern_s = paste0(filter_opt,filter_cutoff,sample_cutoff,'*',id1,'*.RData')
  mergingFiles_lst =  list.files(data_dir, pattern=glob2rx(pattern_s))
  merged_score1 = mergeMultiFiles(paste0(data_dir,mergingFiles_lst))
  write.table(merged_score1, file = paste0(res_dir,id1,'_Merged.txt'), sep = "\t", quote = FALSE, row.names = TRUE)
  merged_score1.Pinfo = getPathwaysInfo(rownames(merged_score1))
  write.table(merged_score1.Pinfo, file =paste0(res_dir,id1,'_PathwaysInfo.txt'), sep = "\t", quote = FALSE, row.names = TRUE)
  annot.s = annot[colnames(merged_score1),]
  p.df.s= merged_score1.Pinfo[merged_score1.Pinfo$p_id %in% rownames(merged_score1),]
  df.s = merged_score1[rownames(merged_score1) %in% rownames(p.df.s),]
  visData(df.s,annot.s,p.df.s,paste0(res_dir,id1,'_visData.pdf'))
  nVarF=nrow(df.s)
  df.s.obj = df2obj(df.s,annot.s,nVarF)
  dim_name = 'umap'
  visObj(df.s.obj,dim_name,paste0(res_dir,id1,'_vis_varF_',nVarF,'.pdf'))
  # # id2
  pattern_s = paste0(filter_opt,filter_cutoff,sample_cutoff,'*',id2,'*.RData')
  mergingFiles_lst =  list.files(data_dir, pattern=glob2rx(pattern_s))
  merged_score2 = mergeMultiFiles(paste0(data_dir,mergingFiles_lst))
  write.table(merged_score2, file = paste0(res_dir,id2,'_Merged.txt'), sep = "\t", quote = FALSE, row.names = TRUE)
  merged_score2.Pinfo = getPathwaysInfo(rownames(merged_score2))
  write.table(merged_score2.Pinfo, file =paste0(res_dir,id2,'_PathwaysInfo.txt'), sep = "\t", quote = FALSE, row.names = TRUE)
  annot.s = annot[colnames(merged_score2),]
  p.df.s= merged_score2.Pinfo[merged_score2.Pinfo$p_id %in% rownames(merged_score2),]
  df.s = merged_score2[rownames(merged_score2) %in% rownames(p.df.s),]
  visData(df.s,annot.s,p.df.s,paste0(res_dir,id2,'_visData.pdf'))
  nVarF=nrow(df.s)
  df.s.obj = df2obj(df.s,annot.s,nVarF)
  dim_name = 'umap'
  visObj(df.s.obj,dim_name,paste0(res_dir,id2,'_vis_varF_',nVarF,'.pdf'))
  # id1 & id2
  merged_score12=run_04_mergeTwoDatasets(merged_score1,merged_score2)
  merged_score12.Pinfo = getPathwaysInfo(rownames(merged_score12))
  write.table(merged_score12.Pinfo, file =paste0(res_dir,id1,'_',id2,'_PathwaysInfo.txt'), sep = "\t", quote = FALSE, row.names = TRUE)
  annot.s = annot[colnames(merged_score12),]
  p.df.s= merged_score12.Pinfo[merged_score12.Pinfo$p_id %in% rownames(merged_score12),]
  df.s = merged_score12[rownames(merged_score12) %in% rownames(p.df.s),]
  visData(df.s,annot.s,p.df.s,paste0(res_dir,id1,'_',id2,'_visData.pdf'))
  nVarF=nrow(df.s)
  df.s.obj = df2obj(df.s,annot.s,nVarF)
  dim_name = 'umap'
  visObj(df.s.obj,dim_name,paste0(res_dir,id1,'_',id2,'_vis_varF_',nVarF,'.pdf'))
  annot = annot[colnames(merged_score12),]
  write.table(annot, file = paste0(res_dir,'annot_Merged.txt'), sep = "\t", quote = FALSE, row.names = TRUE)
  ################################################################################

  # 05 CCA
  data_dir =  paste0(expr_dir, '04_mergeTwoDatasets/')
  res_dir =  paste0(expr_dir, '05_CCA/')
  dir.create(res_dir)

  # 1. VF
  VF_opt = T
  nVarF=nrow(merged_score12)
  ObjType = 'merged' # 'merged = from gsea , take x,y' 'single = each dataset x/y has FV'
  s1 = rownames(subset(annot, annot$dataset_name == id1))
  s2 = rownames(subset(annot, annot$dataset_name == id2))
  if(ObjType == 'merged')
  {
    x=merged_score12[,s1]
    y=merged_score12[,s2]
  if(VF_opt)
  {
  annot.s = annot[colnames(merged_score12),]
  p.df.s= merged_score12.Pinfo[merged_score12.Pinfo$p_id %in% rownames(merged_score12),]
  df.s = merged_score12[rownames(merged_score12) %in% rownames(p.df.s),]
  visData(df.s,annot.s,p.df.s,paste0(res_dir,id1,'_',id2,'_visDataAvg.pdf'))
  df.s.obj = df2obj(df.s,annot.s,nVarF)
  saveRDS(df.s.obj,paste0(res_dir,id1,'_',id2,'_obj_varF_',nVarF,'.RDS'))
   vf = df.s.obj@assays$RNA@var.features
   x=x[vf,]
   y=y[vf,]
  }
  }
  if(ObjType == 'single')
  {
    x=merged_score1[,s1]
    y=merged_score2[,s2]
    # more to add
  }
  nVarFx = nrow(x)
  obj.x =df2obj(x,annot[s1,],nVarFx)
  saveRDS(obj.x,paste0(res_dir,id1,'_',ObjType,'obj_varF_',nVarFx,'.RDS'))
  nVarFy = nrow(y)
  obj.y =df2obj(y,annot[s2,],nVarFy)
  saveRDS(obj.y,paste0(res_dir,id2,'_',ObjType,'obj_varF_',nVarFy,'.RDS'))
  dim_name = 'umap'
  visObj(obj.x,dim_name,paste0(res_dir,id1,'_vis_varF_',nVarFx,'.pdf'))
  dim_name = 'umap'
  visObj(obj.y,dim_name,paste0(res_dir,id2,'_vis_varF_',nVarFy,'.pdf'))

  # 2 CCA
  num.cc = 30
  selected_dims_opt = 'none' #cor'
  selected_dims_opt_cutoff =0
  cca_res =run_05_CCA(obj.x,obj.y,num.cc)
  saveRDS(cca_res,paste0(res_dir,ObjType,'_',num.cc,'_',nVarF,'_1.RDS') )
  dim_name = 'umap_cca'
  #cca_res = Seurat::ScaleData(cca_res)
  cca_res = RunPCA(cca_res,verbose=F,npcs = 10,features = rownames(cca_res))
  cca_res = RunUMAP(cca_res,dims = 1:10, n.components = 2,verbose = F,reduction.name = dim_name)
  visObj(cca_res,dim_name,paste0(res_dir,ObjType,'_',num.cc,'_',nVarF,'_',selected_dims_opt,selected_dims_opt_cutoff,'_obj_',nVarF,'_',dim_name,'.pdf'))
  dim_name = 'cca'
  visObj(cca_res,dim_name,paste0(res_dir,ObjType,'_',num.cc,'_',nVarF,'_',selected_dims_opt,selected_dims_opt_cutoff,'_obj_',nVarF,'_',dim_name,'.pdf'))
  cca_res =processCCAres_allCT(cca_res, x, obj.x,y,obj.y) # add stability over cc
  saveRDS(cca_res,paste0(res_dir,ObjType,'_',num.cc,'_',nVarF,'_','none','0','.RDS'))
  if(selected_dims_opt != 'none')
  {
    cca_res = selectCCA(cca_res,selected_dims_opt,selected_dims_opt_cutoff)
    saveRDS(cca_res,paste0(res_dir,ObjType,'_',num.cc,'_',nVarF,'_',selected_dims_opt,selected_dims_opt_cutoff,'.RDS') )
  }
  # 3 visualize after cca results
 embed2use = 'uv'
  if(embed2use == 'uv')
  {
    xy = t(rbind(cca_res$u,cca_res$v))
    rownames(xy)=paste0('cc',rownames(xy))
    xy.annot = rbind(annot[rownames(cca_res$u),],annot[rownames(cca_res$v),])
  }
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
  obj.xy =df2obj(xy,xy.annot,0)
  saveRDS(obj.xy,paste0(res_dir,ObjType,'_',num.cc,'_',nVarF,'_',selected_dims_opt,selected_dims_opt_cutoff,'_obj_',nVarF,'2.RDS' ))
  dim_name = 'umap'
  visObj(obj.xy,dim_name,paste0(res_dir,ObjType,'_',num.cc,'_',nVarF,'_',selected_dims_opt,selected_dims_opt_cutoff,'_obj_',nVarF,'_',dim_name,'_',embed2use,'.pdf'))

  # 1. run Seurat DataIntegrate
  # integrated_res =run_05_dataIntegration(obj.x,obj.y)
  # saveRDS(integrated_res,paste0(res_dir,ObjType,'_integrated','_varF_',nVarF,'_1.RDS') )
  # dim_name = 'umap'
  #visObj(integrated_res,dim_name,paste0(res_dir,ObjType,'_integrated_',nVarF,'_',dim_name,'.pdf'))
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

  ################################################################################

  # 06 CCA similarity
  data_dir =  paste0(expr_dir, '05_CCA/')
  res_dir =  paste0(expr_dir, '06_similarityCalculation/')
  dir.create(res_dir)
  # cca_res = readRDS(paste0(data_dir,ObjType,'_',num.cc,'_',nVarF,'_',selected_dims_opt,selected_dims_opt_cutoff,'.RDS'))
  run_06_similarityCalculation(cca_res,annot,id.list=c(id1,id2),ct='all',paste0(res_dir,ObjType,'_',num.cc,'_',nVarF,'_',selected_dims_opt,selected_dims_opt_cutoff))

  # in case of Dataintegrate
  # run_06_similarityCalculation(integrated_res,integrated_res$meta.data,id.list=c(id1,id2),ct='all',paste0(res_dir,ObjType,'_integrated_',nVarF))
  # run_06_similarityCalculationIntegration(integrated_res,integrated_res$meta.data,id.list=c(id1,id2),ct='all',paste0(res_dir,ObjType,'_integratedCorrected_',nVarF))

  ################################################################################

  # 07 rank pathways
  data_dir =  paste0(expr_dir, '05_CCA/')
  res_dir =  paste0(expr_dir, '07_rankPathways/')
  dir.create(res_dir)
  w_opt =F
  #a. rank all pathways
  rankP = run_07_rankPathways(cca_res,merged_score12.Pinfo,w_opt)
  write.table(rankP,file = paste0(res_dir,'rankP_allCT.txt'), sep = "\t", quote = TRUE, row.names = TRUE, col.names = T)
  #plot
  for(sim_col in c('cor')) #'reg',
  {
    rankP.tmp = rankP
    colnames(rankP.tmp)[which(names(rankP.tmp) == sim_col)] <- "p_sim"
    rankP.tmp= rankP.tmp[order(rankP.tmp$p_sim, decreasing = T),]
    rankP.tmp = rankP.tmp[1:topRanked,c('p_name','p_level','p_id','p_sim')]
    rownames(rankP.tmp) = rankP.tmp$p_id
    df.s=as.data.frame(rbind(t(cca_res$x[(rankP.tmp$p_id),]),t(cca_res$y[(rankP.tmp$p_id),])))
    annot.s = annot[rownames(df.s),]
    title = paste0('heatmap of gsea for pathways','\n','across all cell-types between two datasets:')
    plotRankingP(df.s ,rankP.tmp$p_level,rankP.tmp$p_sim,annot.s$dataset_name,title,paste0(res_dir,'rankP_',sim_col,'_all.pdf'))
  }
  #b. rank per cell-type pathways
  CTs = unique(annot$type)
  fname.lst = list()
  ct.lst = list()
  for(ct in CTs)
  {
    print(ct)
  s1 = rownames(subset(annot, annot$dataset_name == id.list[[1]]  & annot$type == ct))
  s2 = rownames(subset(annot, annot$dataset_name == id.list[[2]] & annot$type == ct))
  x.s = cca_res$x[,s1]
  y.s = cca_res$y[,s2]
  cc_res_ct= processCCAres_perCT(cca_res,x.s,y.s)
  saveRDS(cc_res_ct,paste0(res_dir,ct,'_merged_cca.RDS') )
  load(paste0(paste0(expr_dir, '03_sampleScoring/'),filter_opt,filter_cutoff,sample_cutoff,'_reactome.db_',id.list[[1]] ,'_norm_var_',ct,'.RData')) #pvalue0.050.07_reactome.db_human_norm_var_Astrocyte.RData
  selected_p1 = sampleScore$selected_pathways_sig
  load(paste0(paste0(expr_dir, '03_sampleScoring/'),filter_opt,filter_cutoff,sample_cutoff,'_reactome.db_',id.list[[2]],'_norm_var_',ct,'.RData'))
  selected_p2 = sampleScore$selected_pathways_sig
  rankP_ct = run_07_rankPathways(cc_res_ct,merged_score12.Pinfo,w_opt)
  write.table(rankP_ct,file = paste0(res_dir,'rankP_',ct,'_all.txt'), sep = "\t", quote = TRUE, row.names = TRUE, col.names = T)
  rankP_ct = rankP_ct[rankP_ct$p_id %in% intersect(selected_p1,selected_p2),]
  write.table(rankP_ct,file = paste0(res_dir,'rankP_',ct,'_selected.txt'), sep = "\t", quote = TRUE, row.names = TRUE, col.names = T)
  fname.lst = append(fname.lst,paste0(res_dir,'rankP_',ct,'_selected.txt'))
  ct.lst = append(ct.lst,ct)
  #plot
    for(sim_col in c('cor')) #'reg',
    {
      rankP.tmp = rankP_ct
      colnames(rankP.tmp)[which(names(rankP.tmp) == sim_col)] <- "p_sim"
      rankP.tmp= rankP.tmp[order(rankP.tmp$p_sim, decreasing = T),]
      rankP.tmp = rankP.tmp[1:topRanked,c('p_name','p_level','p_id','p_sim')]
      rownames(rankP.tmp) = rankP.tmp$p_id
      df.s=as.data.frame(rbind(t(cc_res_ct$x[(rankP.tmp$p_id),]),t(cc_res_ct$y[(rankP.tmp$p_id),])))
      annot.s = annot[rownames(df.s),]
      title = paste0('heatmap of gsea for pathways', '\n','across all cell-types', '\n','between two datasets:')
      plotRankingP(df.s ,rankP.tmp$p_level,rankP.tmp$p_sim,annot.s$dataset_name,title,paste0(res_dir,'rankP_',sim_col,'_',ct,'.pdf'))
      colnames(cc_res_ct$xscores)= paste0('X:',colnames(cc_res_ct$xscores))
      colnames(cc_res_ct$yscores)= paste0('Y:',colnames(cc_res_ct$yscores))
      df.s=as.data.frame(rbind(t(cc_res_ct$xscores[(rankP.tmp$p_id),]),t(cc_res_ct$yscores[(rankP.tmp$p_id),])))
      dataset_name = c(rep(id.list[[1]] ,each=ncol(cc_res_ct$xscores)),rep(id.list[[2]],each=ncol(cc_res_ct$yscores)))
      annot.s = data.frame(dataset_name)
      rownames(annot.s) = rownames(df.s)
      title = paste0('heatmap of CCA scores for pathways', '\n','across all cell-types', '\n','between two datasets:')
      plotRankingP(df.s ,rankP.tmp$p_level,rankP.tmp$p_sim,annot.s$dataset_name,title,paste0(res_dir,'rankP_',sim_col,'_',ct,'_XYscores.pdf'))
    }
  }
 # c. combine rankingP
  plotCombineRankingPathways(fname.lst,ct.lst,res_dir)
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
  # CTs = unique(annot$type)
  # for(ct in CTs)
  # {
  #   print(ct)
  #   # get samples
  #   s1 = rownames(subset(annot, annot$dataset_name == id.list[[1]]  & annot$type == ct))
  #   s2 = rownames(subset(annot, annot$dataset_name == id.list[[2]] & annot$type == ct))
  #   print(length(s1))
  #   print(length(s2))
  #   x.s = x[,s1]
  #   y.s = y[,s2]
  #   integrated_res_ct= processCCAres_perCT(integrated_res,x.s,y.s)
  # saveRDS(integrated_res_ct,paste0(res_dir,ct,'merged_integrated.RDS') )
  #   rankP_ct = run_07_rankPathways(integrated_res_ct,merged_score12.Pinfo)
  #   write.table(rankP_ct,file = paste0(res_dir,'rankP_',ct,'_integrated.txt'), sep = "\t", quote = TRUE, row.names = TRUE, col.names = T)
  #fname.lst = append(fname.lst,paste0(res_dir,'rankP_',ct,'_integrated.txt'))
  #ct.lst = append(ct.lst,ct)
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
  # plotCombineRankingPathways(fname.lst,ct.lst,res_dir)
  ########################      NOT TO RUN YET      ##############################

  ################################################################################

  # 08 rank genes
  #```{r  08. rank genes}
  data_dir =  paste0(expr_dir, '02_runGSEA/')
  res_dir =  paste0(expr_dir, '08_rankGenes/')
  if (!dir.exists(res_dir))
    dir.create(res_dir)
  print(res_dir)
  CTs = unique(annot$type)
  for(ct in CTs)
   {
  # a. get pathways ranked by 'cor' or 'reg' for topRanked for ct
  rankP.fname =  paste0( paste0(expr_dir, '07_rankPathways/'),'rankP_',ct,'.txt')
  # b. get enriched pathway genes members from GSEA results, for human, for mouse
  gsea.fname = paste0(paste0(expr_dir, '02_runGSEA/'),'reactome.db_',id.list[[1]] ,'_norm_var_',ct,'.RData')
  load(gsea.fname)
  expr_logFC1 = expr_logFC
  enrichedG_1 = getEnrichedGenes(GSEA_res$gsea_res,rankP)
  write.table(enrichedG_1,file = paste0(res_dir,id.list[[1]] ,'_enrichedG_',ct,'.txt'), sep = "\t", quote = TRUE, row.names = TRUE, col.names = T)
  gsea.fname = paste0(paste0(expr_dir, '02_runGSEA/'),'reactome.db_',id.list[[2]],'_norm_var_',ct,'.RData')
  load(gsea.fname)
  expr_logFC2 = expr_logFC
  enrichedG_2 = getEnrichedGenes(GSEA_res$gsea_res,rankP)
  write.table(enrichedG_2,file = paste0(res_dir,id.list[[2]],'_enrichedG_',ct,'.txt'), sep = "\t", quote = TRUE, row.names = TRUE, col.names = T)
  # c. get common genes between human and mouse
  commonG = getCommonGenesPerPathway(rankP,enrichedG_1,enrichedG_2)
  print(dim(commonG))
  # d. rank genes
  rankG = run_08_rankGenes(expr_logFC1,expr_logFC2,commonG,rankP)
  write.table(rankG,file = paste0(res_dir,'rankG_',ct,'.txt'), sep = "\t", quote = TRUE, row.names = TRUE, col.names = T)
  # e. plot ranking
  for(p in unique(rankG$p_id))
  {
    print(p)
      df.s = rankG[rankG$p_id == p,]
      name = unique(df.s[df.s$p_id == p,'p_name'])
      rownames(df.s)=df.s$g_id
      df.s = df.s[order(df.s$g_expr_rank, decreasing = F),]
      title = paste0('heatmap of genes in pathway:\n',name, '\nfor cell-type:',ct)
      plotRankingG(t(df.s[, c('g_lfc1','g_lfc2')]) ,df.s$g_rank,df.s$g_diff,title,paste0(res_dir,'rankG_',p,'_',ct,'.pdf'))
  }
  }
  ################################################################################

  # 09 detect Unknown CellType
  data_dir =  paste0(expr_dir, '01_processData/')
  res_dir =  paste0(expr_dir, '09_detectUnknownCT/')
  if (!dir.exists(res_dir))
    dir.create(res_dir)
  print(res_dir)
  CTs = unique(annot$type)
  for(ct in CTs)
  {
  obj1.s =subset(obj.x,cells=colnames(obj.x)[Idents(obj.x)!=ct])
   run_09_detectUnknownCT(obj1.s,obj.y,annot,id.list,ct,res_dir)
  }
  #############################################################

}
