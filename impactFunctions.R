
#' Assign as base raster
#' @description The input raster defines the resolution for all outputs modeled
#' at the land surface. Conversion of the raster to the matrix format used by other
#' simulation functions also occurs.
#' @param rBase The input raster with desired resolution and spatial extent.
#' @return A list with the raster data as a matrix and spatial extent object.
#' @details Do not change the names of items included in the list; these will be used
#' for other functionality.
#' @author Created by CDMartinez 28 June 16
#' @note Edited by CDMartinez 28 June 16
#' @importFrom raster as.matrix extent raster
#' @examples library(raster)
#' rLogo = raster(system.file("external/rlogo.grd", package="raster"))
#' plot(rLogo)
#' logoList = surfaceBase(rLogo)
#' summary(logoList)
#' logoList$rectBB
#' image(logoList$mBase)
#' @export
surfaceBase = function(rBase){
mBase = as.matrix(rBase)
rectBB = extent(rBase)
return(list(mBase = mBase, rectBB = rectBB))
}

#' Convert raster to matrix
#' @description Conversion of a raster to the matrix format used by other
#' simulation functions.
#' @param rSurface The input raster with resolution and spatial extent associated
#' with the base raster returned by \code{\link{surfaceBase}}.
#' @return A list with the raster data as a matrix and spatial extent object.
#' @author Created by CDMartinez 28 June 16
#' @note Edited by CDMartinez 28 June 16
#' @importFrom raster as.matrix raster
#' @examples library(raster)
#' rLogo = raster(system.file("external/rlogo.grd", package="raster"))
#' plot(rLogo)
#' mLogo = surfaceConvert(rLogo)
#' image(mLogo)
#' @export
surfaceConvert = function(rSurface){
  as.matrix(rSurface)
#  mSurface = as.matrix(rSurface)
#  rectBB = extent(rSurface)
#  return(list(mSurface = mSurface, rectBB = rectBB))
}

#' Resample raster
#' @description Resamples one raster to the resolution of the second raster.
#' @param rSurface Raster to be resampled.
#' @param rBase Raster with desired resolution and extent.
#' @return Raster sampled from data in \code{rSurface} with resolution and extent of \code{rBase}
#' @author Created by CDMartinez 28 June 16
#' @note Edited by CDMartinez 28 June 16
#' @importFrom raster resample raster
#' @examples library(raster)
#' rLogo = raster(system.file("external/rlogo.grd", package="raster"))
#' rBase = raster(resolution = c(0.5,0.5), xmn = 0, xmx = 101, ymn = 0, ymx = 77)
#' rLogoFine = surfaceResample(rLogo, rBase)
#' plot(rLogoFine)
#' plot(rLogo)
#' @export
surfaceResample = function(rSurface, rBase){
  rOut = resample(rSurface, rBase)
}

#' Prepare Raster for simulation
#' @description Prepares an input raster for change simulation. The raster values represent surface
#' classifications. The simulation will aggregate change according to the classifications.
#' @param rRaster The raster to prepare for simulation. Raster must have integer values representing
#' the desired classifications.
#' @param spatialList The list returned from \code{\link{prepareSimSpatial}}
#' @param shapeData Optional list used in conjunction with call to \code{\link{prepareChangeShape}}
#' @return List with the matrix version of raster with required spatial extents,
#' number of classifications present, and integer classification values.
#' @author Created by CDMartinez 21 Oct 16
#' @note Edited by CDMartinez 8 Feb 17
#' @importFrom raster resample extent as.matrix values crop res
#' @examples library(raster)
#' OGasmt = continuousAssessment(auMC = 50000,
#' auType = 'Gas',
#' auProbability = 1,
#' auAreaProductive = c(2220000,2800000,3742000),
#' auAreaDrainage = c(10,20,40),
#' auPercAreaUntested = c(93,96,99),
#' auPercAreaSweet = c(100,100,100),
#' auPercFutureSS = c(80,92,94),
#' auEURss = c(0.15,0.4,0.65),
#' auLGR = c(.08,.5,1),
#' year = 2012)
#'
#' rBase = raster(resolution = c(10,10), xmn = 0, xmx = 100, ymn = 0, ymx = 100)
#' values(rBase) = seq(1:100)
#' points = rbind(c(30, 40), c(30, 50), c(40, 50), c(40, 60), c(50, 60), c(60, 60),
#' c(60, 50), c(50, 50), c(50, 40), c(40, 40), c(30, 40))
#' shape = spPolygons(points)
#'
#' plot(rBase)
#' lines(shape)
#'
#' prepSpatial = prepareSimSpatial(rBase,shape, OGasmt)
#' rOwnerTypes <- raster(matrix(sample(1:4, 225, replace = TRUE), nrow = 15),
#' xmn = -500, xmx = 2500, ymn = -500, ymx = 2500)
#' prepOwnerTypes <- prepareChangeRaster(rRaster = rOwnerTypes, spatialList = prepSpatial)
#' @export
prepareChangeRaster = function(rRaster, spatialList, shapeData = NULL){
  rBase = spatialList$rGrid
  rastBB = spatialList$rastBoundBox
  rChange = rRaster
  resCompare = res(rBase) == res(rRaster)
  if(rastBB == extent(rRaster) && resCompare[1] == TRUE && resCompare[2] == TRUE){ # if input raster already has desired extent and resolution
    mEco = as.matrix(rRaster)    # then get raster numeric values as a matrix (for faster processing),
    ecoItems = unique(values(rRaster)) # and the number of unique values,
    nEco = length(ecoItems) # stored in nEco to prevent redundant calls later
  }else{ # if input raster does not have desired extent
    rTemp = crop(rRaster, rBase) # first crop the input raster to base raster extent
      if(rastBB == extent(rTemp) && resCompare[1] == TRUE && resCompare[2] == TRUE){ # and check if the extent is now the same
        rChange = rTemp
      }else{ # if extent is not the same, raster resolutions may vary
        rChange = resample(rRaster, rBase, method = "ngb") # and input raster is resampled to base raster resolution
      }
    mEco = as.matrix(rChange)
    ecoItems = unique(values(rChange))
    nEco = length(ecoItems)
  }

  ecoItems = sort(setdiff(ecoItems, 0))

  if(missing(shapeData) == TRUE){
    rastReturn = list(mEco = mEco, ecoItems = ecoItems, nEco = nEco, rastRes = res(rChange),
                      ecoNames = ecoItems)
  }else{
    rastReturn = list(mEco = mEco, ecoItems = ecoItems, nEco = nEco, rastRes = res(rChange),
                      ecoNames = shapeData$classesChar[ecoItems])
  }

  return(rastReturn)
}

#' Prepare Shapefile for Change Comparison
#' @description Prepares a shapefile values and parameters for change comparison
#' @param shape The shapefile to prepare for simulation. Shapefile must have an attribute
#' representing the surface classifications.
#' @param shapeAttribute String providing the field name of the attribute to be used
#' for classification.
#' @param spatialList The list returned from \code{\link{prepareSimSpatial}} containing the base
#' raster resolution to be used for conversion.
#' @return List returned from \code{\link{prepareChangeRaster}} with the matrix version of raster with required spatial extents,
#' number of classifications present, integer classification values
#' @author Created by CDMartinez 21 Oct 16
#' @note Edited by CDMartinez 8 Feb 17
#' @importFrom raster rasterize
#' @examples library(raster)
#' OGasmt = continuousAssessment(auMC = 50000,
#' auType = 'Gas',
#' auProbability = 1,
#' auAreaProductive = c(2220000,2800000,3742000),
#' auAreaDrainage = c(10,20,40),
#' auPercAreaUntested = c(93,96,99),
#' auPercAreaSweet = c(100,100,100),
#' auPercFutureSS = c(80,92,94),
#' auEURss = c(0.15,0.4,0.65),
#' auLGR = c(.08,.5,1),
#' year = 2012)
#'
#' rBase = raster(resolution = c(10,10), xmn = 0, xmx = 100, ymn = 0, ymx = 100)
#' values(rBase) = seq(1:100)
#' points = rbind(c(30, 40), c(30, 50), c(40, 50), c(40, 60), c(50, 60), c(60, 60),
#' c(60, 50), c(50, 50), c(50, 40), c(40, 40), c(30, 40))
#' shape = spPolygons(points)
#'
#' plot(rBase)
#' lines(shape)
#'
#' prepSpatial = prepareSimSpatial(rBase,shape, OGasmt)
#'
#' criticalArea <- SpatialPolygons(list(
#' Polygons(list(Polygon(rbind(c(0,500),c(1000,2000),c(1500,2000),c(0,0),c(0,500)))), 'A'),
#' Polygons(list(Polygon(rbind(c(1300,700),c(1100,950),c(900,700),c(1100,450),c(1300,700)))),'B')))
#' criticalAreaShape <- SpatialPolygonsDataFrame(criticalArea, data.frame(
#' Years = rbind(3,5), Designation = rbind('Wildlife Corridor', 'Nest Site'),
#' row.names = c('A','B')), match.ID = TRUE)
#' prepCriticalArea <- prepareChangeShape(shape = criticalAreaShape,
#' shapeAttribute = 'Designation', spatialList = prepSpatial)
#' @export
prepareChangeShape = function(shape, shapeAttribute = NULL, spatialList){
  rBase = spatialList$rGrid
  rastBB = spatialList$rastBoundBox

  if(missing(shapeAttribute) == TRUE){

    message('Defaulting to binary classification: present (1) or not present (NA)')
    rRaster = rasterize(shape, rBase, field = 1)
    return(prepareChangeRaster(rRaster, spatialList = spatialList))

  }else{

    message(paste('Using', eval(parse(text = "shapeAttribute")), 'as integer classification.'))
    classesChar = eval(parse(text = paste("shape$", eval(parse(text = "shapeAttribute")), sep = "")))
    classesInt = seq(1:length(classesChar))
    rRaster = rasterize(shape, rBase, field =  classesInt)
    # rRaster = rasterize(shape,rBase,field = eval(parse(text=paste("shape$",eval(parse(text = "shapeAttribute")),sep=""))))
    return(prepareChangeRaster(rRaster, spatialList = spatialList,
                               shapeData = list(classesChar = classesChar, classesInt = classesInt)))
  }

}

# Calculate total surface impacts
# @param ecoInputs List returned from \code{\link{prepareChangeRaster}}
# @param mPads Binary change matrix applied to \code{mEco}
# @return Vector of cumulative change for each unique classification
# @author Created by CDMartinez 28 June 16
# @note Edited by CDMartinez 28 Nov 16
# @export
# surfaceImpact = function(ecoInputs, mPads){
#   mEco = ecoInputs$mEco
#   nEco = ecoInputs$nEco
#   cellArea = ecoInputs$rastRes[1]*ecoInputs$rastRes[2]
#   table(mPads*mEco)*cellArea
# }

#' Calculate changes to raster dataset based on pad placement
#' @description Calculates the changes to input raster dataset based on pad placement.
#' @param ecoInputs List returned from \code{\link{prepareChangeRaster}} or \code{\link{prepareChangeShape}}
#' @param padsIn List returned from a single call to \code{\link{placePads}}
#' @return Vector of cumulative change for each unique classification specified in \code{ecoInputs$ecoItems}
#' @author Created by CDMartinez 28 June 16
#' @note Edited by CDMartinez 8 Feb 17
#' @examples library(raster)
#' OGasmt <- continuousAssessment(auMC = 50,
#'                                auType = 'Gas',
#'                                auProbability = 1,
#'                                auAreaProductive = c(100,400,800),
#'                                auAreaDrainage = c(10,20,40),
#'                                auPercAreaUntested = c(93,96,99),
#'                                auPercAreaSweet = c(100,100,100),
#'                                auPercFutureSS = c(20,40,50),
#'                                auEURss = c(0.15,0.4,0.65),
#'                                auLGR = c(.08,.5,1),
#'                                year = 2016)
#'
#' OGasmt <- convertAcre2sqMeter(OGasmt)
#'
#' rBase <- raster(resolution = c(10,10), xmn = 0, xmx = 2000, ymn = 0, ymx = 2000)
#' values(rBase) <- sample(1:10, 40000, replace = TRUE)
#'
#' points <- rbind(c(250,250),c(250,1750),c(1750,1750),c(1750,250),c(250,250))
#' shape <- SpatialPolygons(list(Polygons(list(Polygon(points)), 'auOutline')))
#'
#' spatialPrep <- prepareSimSpatial(surfaceRaster = rBase, shape, OGasmt)
#' distributionPrep <- prepareSimDistributions(spatialPrep,wellsPerPad = 1,
#' padArea = 300, EA = OGasmt, numIterations=5)
#'
#' pads <- placePads(distributionPrep, 5)
#'
#' rHabitat <- rBase
#' prepHab <- prepareChangeRaster(rRaster = rHabitat, spatialList = spatialPrep)
#'
#' padConversion <- impactsPads(ecoInputs = prepHab, padsIn = pads)
#' @export
impactsPads = function(ecoInputs, padsIn){
  mPads = padsIn$mPads
  mEco = ecoInputs$mEco
  nEco = ecoInputs$nEco
  ecoItems = ecoInputs$ecoItems
  cellArea = ecoInputs$rastRes[1]*ecoInputs$rastRes[2]
  impact = mPads*mEco
  tallyMatrix = tabulate(impact)
  tallyMatrix[ecoItems]*cellArea
}

#' Calculate changes to raster dataset based on road placement
#' @description Calculates changes to raster dataset based on random road placement.
#' @param ecoInputs List returned from \code{\link{prepareChangeRaster}} or \code{\link{prepareChangeShape}}
#' @param roadsIn List returned from \code{\link{placeRoads}}
#' @return Vector of cumulative change for each unique classification specified in \code{ecoInputs$ecoItems}
#' @author Created by CDMartinez 28 June 16
#' @note Edited by CDMartinez 8 Feb 17
#' @examples library(raster)
#' OGasmt <- continuousAssessment(auMC = 50,
#'                                auType = 'Gas',
#'                                auProbability = 1,
#'                                auAreaProductive = c(100,400,800),
#'                                auAreaDrainage = c(10,20,40),
#'                                auPercAreaUntested = c(93,96,99),
#'                                auPercAreaSweet = c(100,100,100),
#'                                auPercFutureSS = c(20,40,50),
#'                                auEURss = c(0.15,0.4,0.65),
#'                                auLGR = c(.08,.5,1),
#'                                year = 2016)
#'
#' OGasmt <- convertAcre2sqMeter(OGasmt)
#'
#' rBase <- raster(resolution = c(10,10), xmn = 0, xmx = 2000, ymn = 0, ymx = 2000)
#' values(rBase) <- sample(1:10, 40000, replace = TRUE)
#'
#' points <- rbind(c(250,250),c(250,1750),c(1750,1750),c(1750,250),c(250,250))
#' shape <- SpatialPolygons(list(Polygons(list(Polygon(points)), 'auOutline')))
#'
#' spatialPrep <- prepareSimSpatial(surfaceRaster = rBase, shape, OGasmt)
#' distributionPrep <- prepareSimDistributions(spatialPrep,wellsPerPad = 1,
#' padArea = 300, EA = OGasmt, numIterations=5)
#'
#' nVertices <- 500
#' road1 <- cbind(seq(0, 2000, length.out = nVertices), seq(0, 100,
#' length.out = nVertices)*sin(seq(-pi, 1.5*pi, length.out = nVertices)) + 600)
#' road2 <- cbind(200*cos(seq(-pi, 1.5*pi, length.out = nVertices)) +
#' seq(200, 1800, length.out = nVertices), seq(0, 2000, length.out = nVertices))
#'
#' prepRoads <- rbind(road1, road2, cbind(road1[,1],rev(road1[,2]) + 700))
#'
#' pads <- placePads(distributionPrep, 5)
#' roadLength <- makeRoads(xyStarts = pads$xyPadCenter, roadNodes = prepRoads)
#' roads <- placeRoads(padsIn = pads, simList = distributionPrep, totalRoadLength = roadLength)
#'
#' rHabitat <- rBase
#' prepHab <- prepareChangeRaster(rRaster = rHabitat, spatialList = spatialPrep)
#'
#' padConversion <- impactsPads(ecoInputs = prepHab, padsIn = pads)
#' roadConversion <- impactsRoads(ecoInputs = prepHab, roadsIn = roads)
#' @export
impactsRoads = function(ecoInputs, roadsIn){
  impact = ecoInputs$mEco[roadsIn$roadCells]
  ecoItems = ecoInputs$ecoItems
  tallyMatrix = tabulate(impact)
  tallyMatrix[ecoItems]*roadsIn$roadArea
}
