% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_06_similarityCalculation.R
\name{run_06_similarityCalculation}
\alias{run_06_similarityCalculation}
\title{run_06_similarityCalculation}
\usage{
run_06_similarityCalculation(merged_res, annot, id.list, ct, res_dir)
}
\arguments{
\item{merged_res}{merged_res results either from cca or integration}

\item{annot}{metadata}

\item{id.list}{ids for the first and second dataset}

\item{ct}{rcell-type}

\item{res_dir}{results directory}
}
\description{
calculate similarity between two datasets based on CCA info
}
\examples{
run_06_similarityCalculation(merged_res,annot,id.list,ct,res_dir)
}
