#' @title  getMetadataStat
#' @description get Stat of metadata
#' @param annot metadata
#' @return metadata_updated updated metadata
#' @examples
#' getMetadataStat(annot)
#' @export
################################################################################
getMetadataStat <- function(annot)
{
  print('getMetadataStat:...')
  print(head(annot))

  dbs = unique(annot$dataset_name)
  cts1= (subset(annot, annot$dataset_name == dbs[1]))$type
  cts2 =  (subset(annot, annot$dataset_name == dbs[2]))$type
  cts = union(cts1,cts2)
  annot_stat =  setNames(data.frame(matrix(0,ncol = length(cts), nrow = length(dbs))), cts)
  rownames(annot_stat) = dbs

  for(ct in colnames(annot_stat))
  {
    annot.ct = subset(annot, annot$type == ct)
    ct1.no= nrow(annot.ct[annot.ct$dataset_name ==dbs[1] ,])
    annot_stat[dbs[1],ct]=ct1.no
    ct2.no= nrow(annot.ct[annot.ct$dataset_name ==dbs[2] ,])
    annot_stat[dbs[2],ct]=ct2.no
  }
  annot_stat=as.data.frame(t(annot_stat))
  print(annot_stat)
  t1= sum(annot_stat[,1])
  t2= sum(annot_stat[,2])
  annot_stat['total',2] = t2
  annot_stat['total',1]  = t1
  print(annot_stat)

  return(annot_stat)
}
################################################################################
