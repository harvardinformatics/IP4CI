#' @title  calcSim_avgCT
#' @description calculate average similarity between two datasets based on CCA info for each cell-type
#' @param df dataframe of cell embeddings
#' @param annot metadata
#' @param id.list ids for the first and second datasets
#' @param ct cell-type
#' @param res_dir results directory
#' @examples
#' calcSim_avgCT(df, annot, id.list,df_type, ct, res_dir)
#' @export
###############################################################################
calcSim_avgCT <- function(df, annot, id.list,ct,res_dir)
{
  print('calcSim_avgCT:...')
  id1=id.list[[1]]
  id2=id.list[[2]]

  df = as.matrix(df)
  # print(dim(df))
  s1 = rownames(subset(annot, annot$dataset_name == id1))
  df_1 = as.matrix(df[,colnames(df) %in% s1] )# pathway x cell
  # print(dim(df_1))
  annot_1 = annot[rownames(annot) %in% colnames(df_1),]
  # print(dim(annot_1))
  df_1 = t(df_1) # cell x pathways
  # message('df1')
  # print(dim(df_1))
  coefx_avg =  setNames(data.frame(matrix(ncol = ncol(df_1), nrow = length(unique(annot_1$type)))),colnames(df_1))
  row.names(coefx_avg) = unique(sort(annot_1$type))
  # print(dim(coefx_avg))
  for(l in unique(annot_1$type))
  {
    l_s = rownames(subset(annot_1, annot_1$type == l))
    df_l = df_1[l_s,]
    # message('df_l')
    # print(dim(df_l))
    for(p in colnames(df_1))
    {
      coefx_avg[l,p]=mean(df_l[,p],na.rm = T)

    }
  }
  rownames(coefx_avg)=paste0(id1,": ",rownames(coefx_avg))
  # print(head(coefx_avg))
  # print(dim(coefx_avg))

  s2 = rownames(subset(annot, annot$dataset_name == id2))
  df_2 = as.matrix(df[,colnames(df) %in% s2] )# pathway x cell
  print(dim(df_2))
  annot_2 = annot[rownames(annot) %in% colnames(df_2),]
  df_2 = t(df_2) # cell x pathways
  # message('df2')
  # print(dim(df_2))
  coefy_avg =  setNames(data.frame(matrix(ncol = ncol(df_2), nrow = length(unique(annot_2$type)))),colnames(df_2))
  row.names(coefy_avg) =  unique(sort(annot_2$type))
  for(l in unique(annot_2$type))
  {
    l_s = rownames(subset(annot_2, annot_2$type == l))
    df_l = df_2[l_s,]
    # message('df_l')
    # print(dim(df_l))
    for(p in colnames(df_2))
    {
      coefy_avg[l,p]=mean(df_l[,p],na.rm = T)

    }
  }
  rownames(coefy_avg)=paste0(id2,": ",rownames(coefy_avg))
  # print(head(coefy_avg))
  # print(dim(coefy_avg))

  # cor
  df = rbind(coefx_avg,coefy_avg)
  dist_res = corr.test(t(df), adjust = 'none')
  dist_sim = dist_res$r[rownames(coefy_avg),rownames(coefx_avg)]
  p=dist_res$p[rownames(coefy_avg),rownames(coefx_avg)]
  dist_sim_norm = (dist_sim-min(dist_sim))/(max(dist_sim)-min(dist_sim))

  dist_res11 = corr.test(t(coefx_avg), adjust = 'none')
  dist_sim1 = dist_res11$r
  p11=dist_res11$p
  dist_sim1_norm = (dist_sim1-min(dist_sim1))/(max(dist_sim1)-min(dist_sim1))

  dist_res12 = corr.test(t(coefy_avg), adjust = 'none')
  dist_sim2 = dist_res12$r
  p12=dist_res12$p
  dist_sim2_norm = (dist_sim2-min(dist_sim2))/(max(dist_sim2)-min(dist_sim2))

  dist_name = 'pearson_cor'
  pdf(paste0(res_dir,ct,'_corSim.pdf'), width = 12, height = 12)
  corrplot(as.matrix(dist_sim), type = "full",is.corr = FALSE, tl.col="black",
                 insig = "blank",cl.ratio = .2, cl.align = "c",tl.cex = 1.2,col=brewer.pal(n=8, name="RdBu"),
                 cl.cex = 1/par("cex"),tl.srt = 90,
                 title=paste0(dist_name,'\n',paste0(id2,'&',id1)) ,mar=c(0,0,1,0)) #, method="number"

  corrplot(as.matrix(dist_sim1), type = "full",is.corr = FALSE, tl.col="black",
                 insig = "blank",cl.ratio = .2, cl.align = "c",tl.cex = 1.2,col=brewer.pal(n=8, name="RdBu"),
                 cl.cex = 1/par("cex"),tl.srt = 90,
                 title=paste0(dist_name,'\n',id1) ,mar=c(0,0,1,0))#, method="number"

  corrplot(as.matrix(dist_sim2), type = "full",is.corr = FALSE, tl.col="black",
                 insig = "blank",cl.ratio = .2, cl.align = "c",tl.cex = 1.2,col=brewer.pal(n=8, name="RdBu"),
                 cl.cex = 1/par("cex"),tl.srt = 90,
                 title=paste0(dist_name,'\n',id2) ,mar=c(0,0,1,0))

  # norm sim
  corrplot(as.matrix(dist_sim_norm), type = "full",is.corr = FALSE, tl.col="black",
                 insig = "blank",cl.ratio = .2, cl.align = "c",tl.cex = 1.2,col=brewer.pal(n=8, name="Blues"),
                 method = 'color',
                 cl.cex = 1/par("cex"),tl.srt = 90,
                 title=paste0(dist_name,'\n',paste0(id2,'&',id1)) ,mar=c(0,0,1,0)) #, method="number"

  corrplot(as.matrix(dist_sim1_norm), type = "full",is.corr = FALSE, tl.col="black",
                 insig = "blank",cl.ratio = .2, cl.align = "c",tl.cex = 1.2,col=brewer.pal(n=8, name="Blues"),
                 method = 'color',
                 cl.cex = 1/par("cex"),tl.srt = 90,
                 title=paste0(dist_name,'\n',id1) ,mar=c(0,0,1,0)) #, method="number"

  corrplot(as.matrix(dist_sim2_norm), type = "full",is.corr = FALSE, tl.col="black",
                 insig = "blank",cl.ratio = .2, cl.align = "c",tl.cex = 1.2,col=brewer.pal(n=8, name="Blues"),
                 method = 'color', #addCoef.col = 'black',
                 cl.cex = 1/par("cex"),tl.srt = 90,
                 title=paste0(dist_name,'\n',id2) ,mar=c(0,0,1,0))

  dev.off()
}
################################################################################




