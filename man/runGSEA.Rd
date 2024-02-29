% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runGSEA.R
\name{runGSEA}
\alias{runGSEA}
\title{run_02_gsea}
\usage{
runGSEA(pathwaydatabase, DE_lst, reference_samples)
}
\arguments{
\item{pathwaydatabase}{pathway database}

\item{DE_lst}{list of logFCs}

\item{reference_samples}{control samples}

\item{dataset_name}{datasetname}
}
\description{
score samples with gsea
}
\examples{
data_dir =  paste0(expr_dir,'/data/')
res_dir =  paste0(expr_dir,'/02_runGSEA/')
createDir(res_dir)
dataset_name='mouse'
obj_fname='mouse_sc.RDS'
dataType = 'norm' # use normalized obj data slot
geneType = 'var' # use variable genes
pathwaydatabase = 'reactome.db'
run_02_gsea(pathwaydatabase,DE_lst,reference_samples,dataset_name)
}