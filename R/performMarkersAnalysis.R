#' @title  performMarkersAnalysis
#' @description perform Analysis using Markers of obj
#' @param markers Markers of obj
#' @param dataset_id dataset_id of obj
#' @param res_dir results_directory to save results and plot
#' @examples
#'performMarkersAnalysis(Markers,dataset_id,opfname)
#'@export
###############################################################################
performMarkersAnalysis <- function(markers,dataset_id,res_dir)
{

  print(dim(markers))
  markers.top = markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
  print(dim(markers.top))


uq_cellTypes = unique(markers$cluster)
gene_lst_all = list()
gene_lst_all_logFC = list()
for(i in 1:length(uq_cellTypes))
{
  cellType = uq_cellTypes[i]
  message('cellType:..,', cellType)
  df = subset(markers,markers$cluster ==cellType )

  # id-maping for gsea and enrichment anlysis
  ids<-bitr(df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
  dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
  df2 = df[df$gene %in% dedup_ids$SYMBOL,]
  df2$Y = dedup_ids$ENTREZID
  original_gene_list <- df2$avg_log2FC
  names(original_gene_list) <- df2$Y
  gene_list<-na.omit(original_gene_list)
  gene_list = sort(gene_list, decreasing = TRUE) # or (e.g.-log10p-valuemultiplied by the sign of log-transformedfold-change)
  print(gene_list)
  message('cellType: before:  ',cellType)
  cellType <- gsub("/", ".",cellType)
  cellType <- gsub(" ", "_",cellType)
  cellType<- gsub("-", ".",cellType)
  message('cellType: after:  ',cellType)
  gene_lst_all[i]=as.list(strsplit(toString(names(gene_list)), ", "))
  names(gene_lst_all)[i]=cellType
  gene_lst_all = gene_lst_all[sort(names(gene_lst_all))]

  #gseGO
  # ids<-bitr(dfL$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
  # dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
  # df2 = dfL[dfL$gene %in% dedup_ids$SYMBOL,]
  # df2$Y = dedup_ids$ENTREZID
  # original_gene_list <- df2$avg_log2FC
  # names(original_gene_list) <- df2$Y
  # gene_list<-na.omit(original_gene_list)
  # gene_list = sort(gene_list, decreasing = TRUE) # or (e.g.-log10p-valuemultiplied by the sign of log-transformedfold-change)
  # gene_lst_all_logFC[[i]] = gene_list
  # names(gene_lst_all_logFC)[i]=cellType
  # gene_lst_all_logFC = gene_lst_all_logFC[sort(names(gene_lst_all_logFC))]

  # if(length(gene_list)>=3)
  # {
  #
  #   kk2 <- gsePathway(geneList     = gene_list,
  #                     organism     = 'human',
  #                     nPerm        = 10000,
  #                     minGSSize    = 3,
  #                     maxGSSize    = 800,
  #                     pvalueCutoff = 0.25,
  #                     pAdjustMethod = "BH",
  #                     seed=TRUE)
  #   if(nrow(as.data.frame(kk2))>0)
  #   {
  #     pdf(paste0(res_dir,dataset_id,'_',cellType,'_GSEA_Reactome.pdf'), width=6, height = 6)
  #     print(dotplot(kk2, showCategory=50 ,split = '.sign')+facet_grid(.~.sign))
  #     kk2 = as.data.frame(kk2)
  #     print(dim(kk2))
  #     write.table(kk2, file = paste0(res_dir,dataset_id,'_',cellType,'_GSEA_Reactome.txt'), sep = "\t", quote=FALSE, row.names=TRUE,col.names=TRUE)
  #     dev.off()
  #   }
  # }

}

# analysis over all cell-types using enrichment
# print(length(gene_lst_all_logFC))
# print(gene_lst_all_logFC)
# for(i in 1:length(gene_lst_all_logFC))
# {
#   gene_lst_all_logFC[[i]]= as.numeric(gene_lst_all_logFC[[i]])
# }
# print('gene_lst_all_logFC')
# print(gene_lst_all_logFC)
# pdf(paste0(opfname,'_markers_gseGO.pdf'), width=70, height = 55)
# ck <- compareCluster(geneCluster = gene_lst_all_logFC, fun = gseGO,OrgDb = org.Hs.eg.db,keyType="ENTREZID")
# # ck <- setReadable(ck, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
#  write.table(ck, file = paste0(opfname,'_markers_gseGO.txt'), sep = "\t", quote=FALSE, row.names=TRUE,col.names=TRUE)
#  dotplot(ck,font.size = 25)
#  dev.off()
#  print('done gseGO')
# print(gene_lst_all)


# pdf(paste0(res_dir,dataset_id,'_enrichReactome.pdf'), width=70, height = 55)
#  ck <- compareCluster(geneCluster = gene_lst_all, fun = enrichPathway,organism = 'human',pvalueCutoff = 0.25, pAdjustMethod = "BH")
#  #ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
#  write.table(ck, file = paste0(res_dir,dataset_id,'_enrichReactome.txt'), sep = "\t", quote=FALSE, row.names=TRUE,col.names=TRUE)
#  print(dotplot(ck,font.size = 25)) #showCategory=30
#  dev.off()
#
# pdf(paste0(res_dir,dataset_id,'_enrichKEGG.pdf'), width=70, height = 55)
#  ck <- compareCluster(geneCluster = gene_lst_all, fun = enrichKEGG,organism = 'hsa',pvalueCutoff = 0.25, pAdjustMethod = "BH")
#  ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
#  write.table(ck, file = paste0(res_dir,dataset_id,'_enrichKEGG.txt'), sep = "\t", quote=FALSE, row.names=TRUE,col.names=TRUE)
#  dotplot(ck,font.size = 25)
#  dev.off()
# #
pdf(paste0(res_dir,dataset_id,'_enrichDO.pdf'), width=25, height = 20)
ck <- compareCluster(geneCluster = gene_lst_all, fun = enrichDO)
#ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(ck, file = paste0(res_dir,dataset_id,'_enrichDO.txt'), sep = "\t", quote=FALSE, row.names=TRUE,col.names=TRUE)
dotplot(ck,font.size = 25)
dev.off()
#
for(ont_go in c('BP','CC','MF'))
{
  pdf(paste0(res_dir,dataset_id,'_enrichGO',ont_go,'.pdf'), width=70, height = 55)
  ck <- compareCluster(geneCluster = gene_lst_all, fun = enrichGO,OrgDb = 'org.Hs.eg.db',ont =ont_go )
  #ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  write.table(ck, file = paste0(res_dir,dataset_id,'_enrichGO',ont_go,'.txt'), sep = "\t", quote=FALSE, row.names=TRUE,col.names=TRUE)
  dotplot(ck,font.size = 25)
  dev.off()
}

}
