% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runIP4CI.R
\name{runIP4CI}
\alias{runIP4CI}
\title{runIP4CI}
\usage{
runIP4CI(
  expr_dir,
  id.list,
  obj.fname.list,
  convertG.list,
  typeColumnToUse.list,
  assaySlotToUse.list,
  genes2keep,
  processObjOpt,
  pathwaydatabase,
  dataType,
  geneType,
  filter_opt,
  filter_cutoff,
  sample_cutoff,
  topRanked
)
}
\arguments{
\item{expr_dir}{experiment directory folder name}

\item{id.list}{id names for the objects}

\item{obj.fname.list}{file names for the objects}

\item{typeColumnToUse.list}{list of column name holds the cell-type labels}

\item{assaySlotToUse.list}{list of name of assay slot}

\item{genes2keep}{which genes to keep between the two objects}

\item{processObjOpt}{T or F to process object with visualization}

\item{pathwaydatabase}{pathway database}

\item{dataType}{which obj data slot to use}

\item{geneType}{which genes to use}

\item{filter_opt}{column of GSEA to filter based on i.e., pvalue}

\item{filter_cutoff}{0.05}

\item{sample_cutoff}{percentage of samples meets the cutoff of filter i.e. 0.1}

\item{topRanked}{10 top ranked biological findings to investigate, used for pathways & genes}
}
\description{
run IP4CI main function to run the analysis using two objects
}
\examples{
expr_dir = ''
createDir(expr_dir)
id.list=c('mouse','human')
obj.fname.list=c('mouse_sc.RDS','human_sc.RDS')
typeColumnToUse.list=c('cell.names','cell.names')
assaySlotToUse.list=list('integrated','SCT')
convertG.list = c(T,F)
processObjOpt = F
convertG.list = c(FALSE,FALSE)
normObj.list = c(FALSE,FALSE)
genes2keep = 'commonVarG'
pathwaydatabase = 'reactome.db'
dataType = 'norm'
geneType = 'var'
filter_opt = 'pvalue'
filter_cutoff = 0.05
sample_cutoff = 0.07
topRanked = 10
runIP4CI(expr_dir, id.list,obj.fname.list,convertG.list,normObj.list,processObjOpt,typeColumnToUse.list,assaySlotToUse.list,
genes2keep,
pathwaydatabase,dataType,geneType,
filter_opt,filter_cutoff,sample_cutoff,
topRanked)
}