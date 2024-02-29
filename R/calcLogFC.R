#' @title  calcLogFC
#' @description calculate LogFC of expr
#' @param expr expr data from obj slot
#' @param annot metadata
#' @param var names of variable/cell-type to consider
#' @examples
#' expr= obj@assays$RNA@data[varfeat,]
#' annot= obj@meta.data
#' var='macrophage'
#' calcLogFC(expr,annot,var)
#' @export
################################################################################
calcLogFC <- function(expr,annot,var)
{

  print('calcLogFC:...')
print(class(expr))
  DE_lst = list()

  c = rownames(subset(annot, annot$type != var))
  n = rownames(annot)[! rownames(annot) %in% c] #non-control
  message('c:n')
  print(length(c))
  print(length(n))

  expr_avg_control = expr[,colnames(expr) %in% c]
  message('expr_avg_control:  ')
  print(dim(expr_avg_control))

  rowmeansmean<- rowMeans(expr_avg_control, na.rm = TRUE)
  expr_avg_control<- cbind(expr_avg_control, rowmeansmean)
  print(dim(expr_avg_control))
  #write.table(expr_avg_control,  paste0(opfname,'_expr_avg_control_manual.txt'),sep='\t',row.names = T, col.names = T,quote = FALSE)


  #https://mkempenaar.github.io/gene_expression_analysis/chapter-4.html
  # # # Get top differentially expressed genes
  expr_logFC = (expr - expr_avg_control[,'rowmeansmean'])
  # rownames(expr_logFC) = rownames(expr)
  print(dim(expr_logFC))

  # in case scoring only n samples
  n_idx=1
  message('n_s:..', length(n))
  for(n_s in n) #ncol(expr_logFC)expr_n_logFC
  {

    df_s = as.data.frame(expr_logFC[,n_s])
    colnames(df_s) = 'logFC'
    df_s$genes  = rownames(df_s)

    DE_lst[[n_idx]] = df_s
    names(DE_lst)[[n_idx]] = n_s
    n_idx=n_idx+1

  }
  message('DE_lst:..',length(DE_lst))
  #print(DE_lst)
  return(list(DE_lst, c,as.data.frame(expr_logFC[,n])))

}
################################################################################

