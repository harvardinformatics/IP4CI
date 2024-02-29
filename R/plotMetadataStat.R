#' @title  plotMetadataStat
#' @description plot Stat of metadata
#' @param metadata metadata
#' @param res_dir results directory to save plot
#' @examples
#' plotMetadataStat(metadata,opfname)
#' @export
################################################################################

plotMetadataStat <- function(annot,res_dir)
{
  print('plotMetadataStat:...')

  for(db in unique(annot$dataset_name))
  {
    print(db)
    pdf(paste0(res_dir,db,'.pdf'), width=12, height = 6)

    annot.s = subset(annot, annot$dataset_name ==db)
    tbl = table(annot.s$type)
    print(dim(tbl))
    df = as.data.frame(tbl)
    names(df)=c('celltype','count')

    colourCount = length(unique(df$celltype))
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))
    p1=df %>%
      ggplot(aes(celltype,count,fill=celltype))+
      geom_col() +
      labs(title=paste0("counts of celltypes : ", db),
           x="celltype", y= "count")+
      geom_label(aes(label = count), position = position_dodge(width = 0.9), angle = 45, color = 'darkblue', size=3,fontface="bold" )+
      scale_fill_manual(values = getPalette(colourCount)) +theme_set(theme_classic(base_size = 12))+theme(legend.position = "none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

    print(p1)
    dev.off()
  }
}

