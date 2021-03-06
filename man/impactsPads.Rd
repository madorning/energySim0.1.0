% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/impactFunctions.R
\name{impactsPads}
\alias{impactsPads}
\title{Calculate changes to raster dataset based on pad placement}
\usage{
impactsPads(ecoInputs, padsIn)
}
\arguments{
\item{ecoInputs}{List returned from \code{\link{prepareChangeRaster}} or \code{\link{prepareChangeShape}}}

\item{padsIn}{List returned from a single call to \code{\link{placePads}}}
}
\value{
Vector of cumulative change for each unique classification specified in \code{ecoInputs$ecoItems}
}
\description{
Calculates the changes to input raster dataset based on pad placement.
}
\note{
Edited by CDMartinez 8 Feb 17
}
\examples{
library(raster)
OGasmt <- continuousAssessment(auMC = 50,
                               auType = 'Gas',
                               auProbability = 1,
                               auAreaProductive = c(100,400,800),
                               auAreaDrainage = c(10,20,40),
                               auPercAreaUntested = c(93,96,99),
                               auPercAreaSweet = c(100,100,100),
                               auPercFutureSS = c(20,40,50),
                               auEURss = c(0.15,0.4,0.65),
                               auLGR = c(.08,.5,1),
                               year = 2016)

OGasmt <- convertAcre2sqMeter(OGasmt)

rBase <- raster(resolution = c(10,10), xmn = 0, xmx = 2000, ymn = 0, ymx = 2000)
values(rBase) <- sample(1:10, 40000, replace = TRUE)

points <- rbind(c(250,250),c(250,1750),c(1750,1750),c(1750,250),c(250,250))
shape <- SpatialPolygons(list(Polygons(list(Polygon(points)), 'auOutline')))

spatialPrep <- prepareSimSpatial(surfaceRaster = rBase, shape, OGasmt)
distributionPrep <- prepareSimDistributions(spatialPrep,wellsPerPad = 1,
padArea = 300, EA = OGasmt, numIterations=5)

pads <- placePads(distributionPrep, 5)

rHabitat <- rBase
prepHab <- prepareChangeRaster(rRaster = rHabitat, spatialList = spatialPrep)

padConversion <- impactsPads(ecoInputs = prepHab, padsIn = pads)
}
\author{
Created by CDMartinez 28 June 16
}
