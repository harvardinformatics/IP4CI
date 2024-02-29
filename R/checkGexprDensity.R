#' @title  checkGexprDensity
#' @description check Gene expression Density
#' @param obj object
#' @param dataset_name dataset name
#' @param ct cell-type
#' @param res_dir results directory to save plot
#' @examples
#' expr_dir = '/Users/rao198/OneDrive\ -\ Harvard\ University/Informatics/doc/paper/DA4LT/V1/analysis/brain/'
#' res_dir =  paste0(expr_dir, '01_processData/')
#' obj = readRDS(paste0(res_dir,'mouse_sc_commonVarG.RDS'))  #commonVarGmouse_sc.RDS
#' ct='Chandelier'
#' dataset_name = 'mouse'
#' checkGexprDensity(obj,dataset_name,ct,res_dir)
#' @export
################################################################################


checkGexprDensity <- function(obj,dataset_name,ct,res_dir)
{

Idents(object = obj) = 'type' #obj1@meta.data$type
obj.s =subset(obj,cells=colnames(obj)[Idents(obj)==ct])

all_markers <- FindAllMarkers(object = obj)
top5_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 5, named_vector = FALSE,
                                    make_unique = TRUE)
cluster.s.markers <- FindMarkers(obj, ident.1 = ct, min.pct = 0.25)
gene2consider = intersect(rownames(cluster.s.markers),VariableFeatures(object = obj))
print(length(gene2consider))

pdf(paste0(res_dir,ct,'_',dataset_name,'.pdf'),width=12,height = 12)
print(DotPlot(obj,features = gene2consider )+scale_colour_gradient2(low="blue",high="red",mid="white"))
dev.off()

}
