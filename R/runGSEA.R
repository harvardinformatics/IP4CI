#' @title  run_02_gsea
#' @description score samples with gsea
#' @param pathwaydatabase pathway database
#' @param DE_lst list of logFCs
#' @param reference_samples control samples
#' @examples
#' run_02_gsea(pathwaydatabase,DE_lst,reference_samples)
#' @export

###############################################################################
runGSEA <- function(pathwaydatabase,DE_lst,reference_samples)
{
  print('runGSEA:...')
  pathways_lst = list()
  gsea_res_tmp = list()

  idx = 1
  for(i in 1:length(DE_lst)) #length(DE_lst)
  {
    message('sample: ',i)
    DE = DE_lst[[i]]
    genes_i = rownames(DE)
    #print(genes_i)
    s_i = names(DE_lst)[i]
    #print(s_i)
    dataset_name ='human'
    if(dataset_name != 'mouse')
    {
      organism_name = 'human'
      message('organism_name:..',organism_name)
      ids<-bitr(genes_i, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
    }

    if(dataset_name == 'mouse')
    {
      organism_name = 'mouse'
      ids<-bitr(genes_i, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Mm.eg.db)
    }
    #print(head(ids))
    dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
    df2 = DE[rownames(DE) %in% dedup_ids$SYMBOL,]
    df2$Y = dedup_ids$ENTREZID
    original_gene_list <- df2$logFC
    names(original_gene_list) <- df2$Y
    gene_list<-na.omit(original_gene_list)
    gene_list = sort(gene_list, decreasing = TRUE) # or (e.g.-log10p-valuemultiplied by the sign of log-transformedfold-change)


    if(pathwaydatabase == 'reactome.db')
    {

      kk2 <- gsePathway(geneList     = gene_list,
                        organism     = organism_name,
                        nPerm        = 10000,
                        minGSSize    = 7, # before to 3
                        maxGSSize    = 500, # before 800
                        pvalueCutoff = 0.25, # before 2
                        pAdjustMethod = "BH",
                        seed=TRUE,BPPARAM = MulticoreParam(workers=3))
    }
    if(pathwaydatabase == 'KEGG.db')
    {
      kk2 <- gseKEGG(geneList     = gene_list,
                     organism     = organism_name,
                     nPerm        = 10000,
                     minGSSize    = 3,
                     maxGSSize    = 500,
                     pvalueCutoff = 0.25,
                     pAdjustMethod = "BH",
                     seed=TRUE)
    }
    if(pathwaydatabase %in% c('BP','MF','CC'))
    {
      if(organism_name == 'mouse')
        organism_db = 'org.Mm.eg.db'
      if(organism_name != 'mouse')
        organism_db = 'org.Hs.eg.db'
      kk2 <- gseGO(geneList     = gene_list,
                   OrgDb     = organism_db,
                   ont = pathwaydatabase,
                   nPerm        = 10000,
                   minGSSize    = 7,
                   maxGSSize    = 500,
                   pvalueCutoff = 0.25,
                   pAdjustMethod = "BH",
                   seed=TRUE)
      #print(head(kk2))
      print(dim(kk2))
      kk2 = simplify(kk2,cutoff = 0.7,by = "pvalue",select_fun = min,measure = "Wang",semData = NULL)
      #print(head(kk2))
      print(dim(kk2))
    }

    kk2 = as.data.frame(kk2)
    print(dim(kk2))


    if(nrow(kk2)== 0)
    {
      #write.table(kk2, paste0(data_dir,var,'_','s',s_i,'_',pathwaydatabase,'_gsea.txt'),sep='\t',row.names = T, col.names = TRUE,quote = FALSE)

      gsea_res_tmp[[i]]=0
      names(gsea_res_tmp)[[i]] = s_i
    }

    if(nrow(kk2)>0)
    {
      # write.table(kk2, paste0(data_dir,var,'_','s',s_i,'_',pathwaydatabase,'_gsea.txt'),sep='\t',row.names = T, col.names = TRUE,quote = FALSE)

      pathways_lst = c(pathways_lst,as.character(kk2$ID)) #Description
      gsea_res_tmp[[i]]=kk2
      names(gsea_res_tmp)[[i]] = s_i
    }
  }
  message('gsea_res_tmp')
  print(length(gsea_res_tmp))

  pathways_lst =  unique(as.character(pathways_lst))
  print(length(pathways_lst))

return(GSEA_res = list('gsea_res' = gsea_res_tmp,'pathways_lst' = pathways_lst,'reference_samples'=reference_samples))
}
