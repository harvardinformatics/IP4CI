#' @title  rankPathwaysRegression
#' @description rank pathways based on regression
#' @param cca_res CCA res
#' @param p.info pathways info
#' @return rankP.reg  ranked pathways based on regression
#' @examples
#' rankPathwaysRegression()
#' @export
################################################################################
rankPathwaysRegression <- function(cca_res, p.info)
{
  print('rankPathwaysRegression:...')

  x_scores = cca_res$xscores
  y_scores = cca_res$yscores
  cor = cca_res$cor
  print(dim(x_scores))
  print(dim(y_scores))

  # plot
  cd.df=data.frame()
  resd.df=data.frame()
  rankP.reg =  setNames(data.frame(matrix(ncol = ncol(x_scores), nrow = nrow(x_scores))), seq(1:ncol(x_scores)))
  rownames(rankP.reg)= rownames(x_scores)

  print(cor)
  cc_r2_total = sum(cor^2)
  for(cc in 1:min(ncol(x_scores),30)) #no_cc
  {
    message("cc:  ",cc)
    data_cc =  setNames(data.frame(matrix(ncol = 3, nrow = length(x_scores[,cc]))), c("X","Y",'label'))
    rownames(data_cc)= rownames(x_scores)
    data_cc$X = x_scores[,cc]
    data_cc$Y = y_scores[,cc]
    cor_r = round(cor.test(data_cc$X,data_cc$Y)$estimate,2)
    cor_pval=round(cor.test(data_cc$X,data_cc$Y)$p.value,2)
    message('cor_r:..',cor_r ,' cor_pval:...',cor_pval,' cc_r2_total:...',cc_r2_total)

    # print(ggplot(data_cc, aes(x=X, y=Y,  color=rownames(data_cc)))
    #       +geom_smooth(method=lm, se=TRUE, linetype="dashed",color="darkred",fullrange=TRUE, formula = y ~ 0 + x,conf.int = TRUE,
    #                    cor.coef = TRUE)
    #       +geom_point(size = 2)
    #       +geom_label_repel(data=data_cc,
    #                         aes(label=rownames(data_cc)),
    #                         colour='black',segment.color = NA,size = 5,label.size = NA,max.overlaps = Inf,min.segment.length = Inf, box.padding = 0.5)
    #       +ggtitle(paste0("pathways scores for CC",cc,'\n','cor estimate= ',cor_r,'\t','    cor pvalue= ',cor_pval))
    #       +theme(legend.position="none",plot.title = element_text(size = 25, face = "bold"),
    #              legend.title=element_text(size=12),
    #              legend.text=element_text(size=12))+
    #         xlab(paste0(id1," xscores for CC",cc)) + ylab(paste0(id2," yscores for CC",cc))+geom_jitter()
    # )
    # print(ggplot(data_cc, aes(x=X, y=Y,  color=rownames(data_cc)))
    #       +geom_smooth()
    #       +geom_point(size = 2)
    #       +geom_label_repel(data=data_cc,
    #                         aes(label=rownames(data_cc)),
    #                         colour='black',segment.color = NA,size = 5,label.size = NA,max.overlaps = Inf,min.segment.length = Inf, box.padding = 0.5)
    #       +ggtitle(paste0("pathways scores for CC",cc,'\n','cor estimate= ',cor_r,'\t','    cor pvalue= ',cor_pval))
    #       +theme(legend.position="none",plot.title = element_text(size = 25, face = "bold"),
    #              legend.title=element_text(size=12),
    #              legend.text=element_text(size=12))+
    #         xlab(paste0(id1," xscores for CC",cc)) + ylab(paste0(id2," yscores for CC",cc, ' smooth()'))+geom_jitter()
    # )

    print('done reg line')
    # study reg
    print(names(data_cc))
    reg.fit = lm(Y ~ 0 + X , data= data_cc)
    res = ols_plot_cooksd_chart(reg.fit)
    # print(names(res))
    # print(res$data)
    cd.cc = res$data
    reg.fit$resid = reg.fit$resid[match(rownames(cd.cc), names(reg.fit$resid))]
    cd.cc$cc = cc
    cd.cc$Pname = rownames(cd.cc)
    cd.cc$Cor = cor_r
    cd.cc$CorP = cor_pval
    cd.cc$resd= as.numeric(abs(reg.fit$resid))
    cd.cc$resdName = names(reg.fit$resid)
    cd.cc$resdNorm = ( cd.cc$resd- min( cd.cc$resd)) /(max(cd.cc$resd)-min( cd.cc$resd))
    cd.cc$resdNormSim = 1-cd.cc$resdNorm
    if(cor_r < 0)
      cd.cc$resdNormSim = -1 * cd.cc$resdNormSim
    #https://stats.stackexchange.com/questions/146069/canonical-correlation-analysis-taking-into-account-the-negatieve-correlations
    cc_r2 = (cor_r^2)/cc_r2_total
    message('cc_r2:...',cc_r2,'\t',cor[cc],'\t',cor_r,'\t',cc_r2_total)
    cd.cc$importance = cc_r2 *cd.cc$resdNormSim


    cd.df=rbind(cd.df,cd.cc)
    rankP.reg[,cc]= cd.cc$importance

    outliers = rownames(subset(cd.cc, cd.cc$fct_color == 'outlier'))
    data_cc.wo.Outliers = data_cc[!rownames(data_cc)%in% outliers ,]
    cor_r = round(cor.test(data_cc.wo.Outliers$X,data_cc.wo.Outliers$Y)$estimate,2)
    cor_pval=round(cor.test(data_cc.wo.Outliers$X,data_cc.wo.Outliers$Y)$p.value,2)
    # print(ggplot(data_cc.wo.Outliers, aes(x=X, y=Y,  color=rownames(data_cc.wo.Outliers)))
    #       +geom_smooth(method=lm, se=TRUE, linetype="dashed",color="darkred",fullrange=TRUE, formula = y ~ 0 + x,conf.int = TRUE,
    #                    cor.coef = TRUE)
    #       +geom_point(size = 2)
    #       +geom_label_repel(data=data_cc.wo.Outliers,
    #                         aes(label=rownames(data_cc.wo.Outliers)),
    #                         colour='black',segment.color = NA,size = 5,label.size = NA,max.overlaps = Inf,min.segment.length = Inf, box.padding = 0.5)
    #       +ggtitle(paste0("pathways scores for CC",cc,'\n','cor estimate= ',cor_r,'\t','    cor pvalue= ',cor_pval))
    #       +theme(legend.position="none",plot.title = element_text(size = 25, face = "bold"),
    #              legend.title=element_text(size=12),
    #              legend.text=element_text(size=12))+
    #         xlab(paste0(id1," xscores for CC",cc)) + ylab(paste0(id2," yscores for CC",cc))+geom_jitter()
    # )

    # ols_plot_dfbetas(reg.fit)
    # ols_plot_dffits(reg.fit)
    # ols_plot_resid_stud(reg.fit)
    # ols_plot_resid_stand(reg.fit)
    # ols_plot_resid_lev(reg.fit)
    #ols_plot_resid_stud_fit(reg.fit)


  }# CC
  #dev.off()

  cd.df$p_id = rownames(cd.df)
  # cd.df = merge(cd.df,p.info,by="p_id")
  rankP.reg$total = rowSums(rankP.reg)
  rankP.reg$p_id =rownames(rankP.reg)
  rankP.reg = merge(rankP.reg,p.info,by="p_id")
  rankP.reg = rankP.reg[order(rankP.reg$total, decreasing = T),]
  colnames(rankP.reg)[which(names(rankP.reg) == "total")] <- "reg"
  rankP.reg = rankP.reg[,c('p_id','reg')]
  return(rankP.reg)

}
################################################################################
