% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functionsConventionalAccumulations.R
\name{conventionalAssessmentSummary}
\alias{conventionalAssessmentSummary}
\title{Print USGS Conventional Assessment Summary}
\usage{
conventionalAssessmentSummary(dataList)
}
\arguments{
\item{dataList}{The \code{risked} or \code{unrisked} list returned from \code{\link{conventionalAssessment}}}
}
\value{
Printed table in the console
}
\description{
Returns a matrix table of the accumulation mean and fractiles typically published in USGS Factsheets.
}
\note{
Edited by CDMartinez  29 Dec 15
}
\examples{
OGasmt <- conventionalAssessment(auMC = 5000, auType = "Oil", auCharge = 1, auRocks = 1,
auTiming = 1, auFieldNum = c(1, 2, 10), auFieldSize = c(.5, .8, 10), auGOR = c(200, 400, 600),
auNGLGR = c(35, 85, 115), year = 2013)
round(conventionalAssessmentSummary(OGasmt$risked))
round(conventionalAssessmentSummary(OGasmt$unrisked))
}
\author{
Created by CDMartinez  1 Nov 15
}
