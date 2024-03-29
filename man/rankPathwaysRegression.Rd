% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rankPathwaysRegression.R
\name{rankPathwaysRegression}
\alias{rankPathwaysRegression}
\title{rankPathwaysRegression}
\usage{
rankPathwaysRegression(cca_res, p.info, w_opt)
}
\arguments{
\item{cca_res}{CCA res}

\item{p.info}{pathways info}

\item{w_opt}{to use weighted pathways scores}
}
\value{
rankP.reg  ranked pathways based on regression
}
\description{
rank pathways based on regression
}
\examples{
rankPathwaysRegression(cca_res, p.info, w_opt)
}
