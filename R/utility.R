
#' Prepare simulation inputs related to the spatial extents
#' @description Generate spatial variables needed for simulation. The \code{surfaceRaster} resolution
#' will be used for all calculated impacts and defines the smallest resolution of placing well pads.
#' @param surfaceRaster Surface raster with desired resolution for modeling impacts
#' @param auOutline The AU outline shapefile.
#' @param EA Oil or gas continuous assessment from \code{\link{continuousAssessment}} and converted
#' to square meters using \code{\link{convertAcre2sqMeter}}.
#' @param sweetOutline Optional. The Sweet Spot shapefile. If not provided, \code{\link{continuousAssessment}}
#' should likewise be an all Sweet Spot assessment.
#' @param affineMatrix Optional. The affine matrix to  use in \code{\link{scaleShapeTransformation}}.
#' @details If \code{affineMatrix} is not provided, default scaling of the AU extent will occur using
#' the shapefile centroid and \code{\link{scaleShape}}.
#' @return List for use in spatial simulation. Do not alter output from this function, as the list
#' is required for use of other functions such as \code{\link{prepareSimDistributions}} or \code{\link{prepareRusle}}.
#' @import raster
#' @note Edited by CDMartinez  6 April 17
#' @author Created by CDMartinez 15 Sept 16
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
#' spatialPrep <- prepareSimSpatial(rBase,shape, OGasmt)
#' @export
prepareSimSpatial = function(surfaceRaster,auOutline,EA,sweetOutline = NULL, affineMatrix = NULL){

  maxArea = max(EA$risked$ProductiveArea)
  if(missing(affineMatrix)){
    message('Using centroid scaling of AU outline for simulation.')
    maxOutline = scaleShape(auOutline,maxArea)
  }else{
    message('Using specified affine scaling of AU outline for simulation.')
    maxOutline = scaleShapeTransformation(inShape=auOutline,newArea=maxArea,affine=affineMatrix)
  }

  if(sum(EA$risked$NonSweetArea) > 0){
    if(missing(sweetOutline)){stop('Sweet spot outline missing.')}
  }

  rBase = raster(ext = extent(surfaceRaster),
                 nrows = nrow(surfaceRaster), ncols = ncol(surfaceRaster))
  # rectBB=extent(maxOutline)
  # rCrop = crop(surfaceRaster,rectBB)
  # rBase = rCrop
  values(rBase) = 0L
  rGrid = crop(rBase,maxOutline)
  rastBB = extent(rGrid)
  names(rGrid) = 'AU maximum AOI'

  if(!missing(sweetOutline)){
    return(list(rGrid = rGrid, rastBoundBox = rastBB, auOutline = auOutline, sweetOutline = sweetOutline))
  }else{
    return(list(rGrid = rGrid, rastBoundBox = rastBB, auOutline = auOutline, sweetOutline = NULL))
  }
}


#' Prepare simulation inputs related to the oil and gas variables
#' @description Calculates constant variables for each iteration prior to simulation.
#' @param simSpatial List returned from \code{\link{prepareSimSpatial}}
#' @param wellsPerPad Number of wells per pad. May be either a constant integer or a vector of integers
#' the same length as the number of simulation iterations (\code{numIterations}) to perform.
#' @param padArea Area of each pad. Units should be consistent with \code{EA} units.
#' May be either a constant integer or a vector of integers the same length as the number of
#' simulation iterations (\code{numIterations}) to perform.
#' @param EA Oil or gas assessment list returned from \code{\link{continuousAssessment}}. Suggest
#' converting to square meters using \code{\link{convertAcre2sqMeter}} prior to using this function.
#' @param numIterations Number of MC iterations to perform
#' @details The \code{numIterations} to perform should be less than or equal to the number of
#' Monte Carlo iterations used to generate the energy
#' assessment list \code{EA}. If less than the number of energy assessments iterations available,
#' random sampling will determine the
#' energy assessment combinations to use.
#' @return List for use in spatial simulation [area units must be square meters]
#' @import raster
#' @note Edited by CDMartinez  15 Nov 16
#' @author Created by CDMartinez 15 Sept 16
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
#' spatialPrep = prepareSimSpatial(rBase,shape, OGasmt)
#' distributionPrep = prepareSimDistributions(spatialPrep,wellsPerPad = 3,
#' padArea = 100, EA = OGasmt, numIterations=5)
#' @export
prepareSimDistributions = function(simSpatial,wellsPerPad,padArea,EA,numIterations=NULL){

  auMC = length(EA$risked$DrainageArea)

  if(missing(numIterations)){
    numIterations = auMC
    message(cat('Note: Defaulting iterations to',auMC))}
  if(numIterations > auMC){
    numIterations = auMC
    message(cat('Note: Specified iterations larger than available oil and gas iterations. Defaulting iterations to',auMC))}

  if(length(wellsPerPad) == 1){
    nWellsPerPad = rep.int(wellsPerPad,times = numIterations)
    message(cat('Note: Constant number of',wellsPerPad,'wells per pad will be used for simulation.'))
  }else{
    nWellsPerPad = wellsPerPad
    message(cat('Note: Variable number of wells per pad will be used for simulation.'))
  }

  if(length(padArea) == 1){
    nPadArea = rep(padArea,times = numIterations)
    message(cat('Note: Constant pad area of',padArea,'square meters will be used for simulation.'))
  }else{
    nPadArea = padArea
    message(cat('Note: Variable pad area will be used for simulation.'))
  }

  if(length(which(nPadArea <= 0)) > 0){stop('Error: Variable pad area given contains zero or negative area.')}
  rGrid = simSpatial$rGrid

  minx = xmin(rGrid)
  maxx = xmax(rGrid)
  miny = ymin(rGrid)
  maxy = ymax(rGrid)
  mPads = as.matrix(rGrid)

  deltax = xres(rGrid)
  deltay = yres(rGrid)
  cellArea = deltay*deltax
  cellBuffer = ceiling(sqrt(nPadArea/(cellArea)))
  halfBuffer = ceiling(cellBuffer/2)
  cellBufferDist = cellBuffer*deltax

  constVars = list(mPads=mPads,minx=minx,maxx=maxx,miny=miny,maxy=maxy,
                   deltax=deltax,deltay=deltay,
                   cellArea=cellArea)#mEco=mEco,nEco=nEco,ecoItems=ecoItems,

  EA$risked$mcNumber = seq(1,auMC,by=1)
  useIters = sample.int(auMC,size=numIterations,replace=FALSE)
  simVars = EA$risked[useIters,]

  simVars$nWellsPerPad = nWellsPerPad
  simVars$nPadArea = nPadArea
  simVars$halfBuffer = halfBuffer
  simVars$cellBuffer = cellBuffer
  simVars$cellBufferDist = cellBufferDist

  simVars$nPadsSweet = ceiling(simVars$NumSweetWells/nWellsPerPad)
  simVars$nPadsNonSweet = ceiling(simVars$NumNonSweetWells/nWellsPerPad)
  simVars$nPads = simVars$nPadsSweet + simVars$nPadsNonSweet
  simVars$padDrainage = simVars$DrainageArea*nWellsPerPad
  simVars$drainCellSize = sqrt(simVars$padDrainage)

  simVars$xextent = maxx - minx
  simVars$yextent = maxy - miny

  return(list(rGrid = rGrid, assessmentVars = simVars, constantVars = constVars,
              rastBoundBox = simSpatial$rastBoundBox, rastRes = res(rGrid),
              anOutline = simSpatial$auOutline, inOutline = simSpatial$sweetOutline))
}


#' Prepare RUSLE inputs
#' @description Calculates base soil loss given initial conditions for the factors.
#' @param spatialList Output generated by \code{\link{prepareSimSpatial}}
#' @param R Rainfall erosivity factor raster
#' @param K Soil erodibility factor raster
#' @param LS Topographic slope length factor raster
#' @param C Vegetative cover factor raster
#' @param shapes Optional shapefile. If provided, simulation will return total soil loss relative to baseline as
#' aggregated to input polygons specified by shapefile. \code{shapes} will invoke call to \code{\link{prepareChangeShape}}
#' @param shapeAttribute Optional attribute. If provided along with \code{shapes} shapefile,
#' will invoke call to \code{\link{prepareChangeShape}}
#' @return List containing inputs used for soil loss simulation, \code{\link{rusle}}. List contains
#' baseline input (R*K*LS*C) and C as matrices and optional input shapes for aggregated output.
#' @author Created by CDMartinez 21 Oct 16
#' @note Edited by CDMartinez 2 Feb 17
#' @importFrom raster crop extent resample as.matrix rasterize values
#' @examples library(raster)
#' set.seed(46)
#' OGasmt <- continuousAssessment(auMC = 5,
#' auType = 'Gas',
#' auProbability = 1,
#' auAreaProductive = c(100,400,800),
#' auAreaDrainage = c(10,20,40),
#' auPercAreaUntested = c(93,96,99),
#' auPercAreaSweet = c(100,100,100),
#' auPercFutureSS = c(20,40,50),
#' auEURss = c(0.15,0.4,0.65),
#' auLGR = c(.08,.5,1),
#' year = 2016)
#'
#' OGasmt <- convertAcre2sqMeter(OGasmt)
#'
#' rBase <- raster(resolution = c(10,10), xmn = 0, xmx = 2000, ymn = 0, ymx = 2000)
#' values(rBase) <- sample(1:10, 40000, replace = TRUE)
#'
#' points <- rbind(c(250,250),c(250,1750),c(1750,1750),c(1750,250),c(250,250))
#' shape <- SpatialPolygons(list(Polygons(list(Polygon(points)), 'auOutline')))
#'
#' plot(rBase, xlim = c(0,2000), ylim = c(0,2000))
#' lines(shape)
#'
#' prepSpatial <- prepareSimSpatial(surfaceRaster = rBase, shape, OGasmt)
#' dxdy <- 400
#'
#' # Slope-length factor
#' tempI <- matrix(complex( real = rep(seq(0.4, .47, length.out = dxdy), each = dxdy ),
#' imag = rep(seq(.3, .42, length.out = dxdy), dxdy)), ncol= dxdy, nrow = dxdy)
#' tempZ <- 0
#' for(k in 1:20){tempZ <- tempZ^2+tempI}
#' rLS <- raster(exp(-abs(tempZ))*20)
#' extent(rLS) = extent(prepSpatial$rGrid)
#'
#' # Vegetative cover factor
#' rC <- raster(matrix(sample(seq(.001, .05, length.out = 7),
#' size = dxdy, replace = TRUE), nrow = sqrt(dxdy)), xmn = -500,
#' xmx = 2500, ymn = -500, ymx = 2500)
#'
#' # Soil erodibility factor
#' rK <- raster(matrix(sample(seq(.01, .7, length.out = 15),
#' size = dxdy, replace = TRUE), nrow = sqrt(dxdy)), xmn = -500,
#' xmx = 2500, ymn = -500, ymx = 2500)
#'
#' # Rainfall erosivity
#' rR <- raster(matrix(rep(10, dxdy), nrow = sqrt(dxdy)),
#' xmn = -500, xmx = 2500, ymn = -500, ymx = 2500)
#'
#' # Soil Loss will be aggregated across entire area
#' prepRusle <- prepareRusle(prepSpatial, R = rR, K = rK, LS = rLS, C = rC)
#' @export
prepareRusle = function(spatialList, R, K, LS, C, shapes = NULL, shapeAttribute = NULL){

  if(missing(shapes) == TRUE){
    message('Simulation will return total soil loss for region.')
  }else{
    message('Simulation will return total soil loss aggregated to input polygons specified by shapefile.')
  }

  rBase = spatialList$rGrid
  rastBB = spatialList$rastBoundBox

  # isoerodent raster check
  resCompare = res(rBase) == res(R)
  if(rastBB == extent(R) && resCompare[1] == TRUE && resCompare[2] == TRUE){
    message('R factor raster same as base rGrid!')
  }else{
    message('Cropping raster: R')
    rTemp = crop(R, rBase) # first crop the input raster to base raster extent
    resCompare = res(rBase) == res(rTemp)
    if(resCompare[1] == TRUE && resCompare[2] == TRUE){ # and check if the resolution is the same
      rR = rTemp
    }else{ # if raster resolutions vary
      message('Raster resolution differs: resampling R to base raster resolution')
      rR = resample(rTemp, rBase, method = "ngb") # input raster resampled to base raster resolution
    }
  }

  # soils raster check
  resCompare = res(rBase) == res(K)
  if(rastBB == extent(K) && resCompare[1] == TRUE && resCompare[2] == TRUE){
    message('K factor raster same as base rGrid!')
  }else{
    message('Cropping raster: K')
    rTemp = crop(K, rBase) # first crop the input raster to base raster extent
    resCompare = res(rBase) == res(rTemp)
    if(resCompare[1] == TRUE && resCompare[2] == TRUE){ # and check if the resolution is the same
      rK = rTemp
    }else{ # if raster resolutions vary
      message('Raster resolution differs: resampling K to base raster resolution')
      rK = resample(rTemp, rBase, method = "ngb") # input raster resampled to base raster resolution
    }
  }

  # slope-length raster check
  resCompare = res(rBase) == res(LS)
  if(rastBB == extent(LS) && resCompare[1] == TRUE && resCompare[2] == TRUE){
    message('LS factor raster same as base rGrid!')
  }else{
    message('Cropping raster: LS')
    rTemp = crop(LS, rBase) # first crop the input raster to base raster extent
    resCompare = res(rBase) == res(rTemp)
    if(resCompare[1] == TRUE && resCompare[2] == TRUE){ # and check if the resolution is the same
      rLS = rTemp
    }else{ # if raster resolutions vary
      message('Raster resolution differs: resampling LS to base raster resolution')
      rLS = resample(rTemp, rBase, method = "ngb") # input raster resampled to base raster resolution
    }
  }

  # vege cover check
  resCompare = res(rBase) == res(C)
  if(rastBB == extent(C) && resCompare[1] == TRUE && resCompare[2] == TRUE){
    message('C factor raster same as base rGrid!')
  }else{
    message('Cropping raster: C')
    rTemp = crop(C, rBase) # first crop the input raster to base raster extent
    resCompare = res(rBase) == res(rTemp)
    if(resCompare[1] == TRUE && resCompare[2] == TRUE){ # and check if the resolution is the same
      rC = rTemp
    }else{ # if raster resolutions vary
      message('Raster resolution differs: resampling C to base raster resolution')
      rC = resample(rTemp, rBase, method = "ngb") # input raster resampled to base raster resolution
    }
  }

  # base raster
  mR = as.matrix(rR)
  mK = as.matrix(rK)
  mLS = as.matrix(rLS)
  mC = as.matrix(rC)

  mRKLSC = mR*mK*mLS*mC

  if(missing(shapes) == TRUE){
    return(list(mRKLSC = mRKLSC, mC = mC, shapeList = NULL))
  }else{
    pShapes = prepareChangeShape(shapes, shapeAttribute = shapeAttribute, spatialList = spatialList)
    return(list(mRKLSC = mRKLSC, mC = mC, shapeList = pShapes))
  }
}


#' Commands for a single iteration
#' @description Optional utility function for a single iteration: computes surface disturbance,
#' road disturbance, and soil loss.
#' @param X Integer specifying iteration index in \code{simDistList}
#' @param simDistList The list returned by \code{\link{prepareSimDistributions}}
#' @param rasterSurfaceInputs Optional. List containing all surface rasters prepared using
#' \code{\link{prepareChangeShape}} or \code{\link{prepareChangeRaster}}
#' @param roads Optional matrix or vector. If vector is provided, road length is calculated by 
#' drawing a road length for each simulated pad from the provided vector of values. If matrix is provided, 
#' road length is calculated using distance between simulated pad and existing road network. 
#' The matrix format should be two column matrix with easting, northing vertices of known road network.
#' @param soilLoss Optional. The list returned by \code{\link{prepareRusle}}
#' @param totalLength Optional. Argument passed to \code{\link{makeRoads}}. \code{TRUE}
#' or \code{FALSE} boolean for returning total road length.
#' For use with \code{placeRoads}, leave as default value of \code{TRUE}.
#' @param roadWidth Optional. Argument passed to \code{\link{placeRoads}}.
#' @param cellProportion Optional. Argument passed to \code{\link{placeRoads}}.
#' @return List for each iteration containing a outputs for each of the three inputs.
#' @author Created by CDMartinez 15 February 17
#' @note Edited by CDMartinez 23 Feb 17
#' @export
#'
simIteration = function(X, simDistList, rasterSurfaceInputs, roads = NULL, soilLoss = NULL,
                        totalLength = TRUE, roadWidth = NULL, cellProportion = NULL){
  simPads = placePads(simDistList, X)
  nSurfaces = length(rasterSurfaceInputs)
  surfacePadStats = lapply(X = rasterSurfaceInputs, FUN = impactsPads, padsIn = simPads)

  simOutput = list(surfacePadStats = surfacePadStats)

  if(is.null(roads) == FALSE){
    if(is.matrix(roads) == TRUE){
      roadStats = makeRoads(simPads$xyPadCenter, roads, totalLength = totalLength)
    }  
    if(is.vector(roads) == TRUE){
      roadStats = makeRoadsD(simPads$xyPadCenter, roads, totalLength = totalLength)
    }
    simRoads = placeRoads(simPads, simDistList, totalRoadLength = roadStats, roadWidth = roadWidth, cellProportion = cellProportion)
    surfaceRoadStats = lapply(X = rasterSurfaceInputs, FUN = impactsRoads, roadsIn = simRoads)

    simOutput$surfaceRoadStats = surfaceRoadStats
    simOutput$roadStats = roadStats

    if(is.null(soilLoss) == FALSE){
      rusleRoadStats = rusle(rusleIn = soilLoss, padsIn = simPads, roadsIn = simRoads)
      simOutput$rusleRoadStats = rusleRoadStats
    }
  }

  if(is.null(roads) == FALSE){
    ruslePadStats = rusle(rusleIn = soilLoss, padsIn = simPads)
    simOutput$ruslePadStats = ruslePadStats
  }

  simOutput
}

#' Aggregate output distributions
#' @description Extracts alike distribution data from list simulation output from
#' \code{\link{simIteration}} that was used in conjuction with a list apply function
#' such as \code{\link{lapply}}, \code{\link{mclapply}}, or \code{\link{parLapply}}.
#' @param simOutput List returned from apply call to \code{\link{simIteration}}.
#' @param rasterSurfaceInputs Optional. List of surface inputs to \code{\link{simIteration}}.
#' If provided, names will be automatically ported to the output.
#' @param soilLoss Optional. List returned from \code{\link{prepareRusle}}
#' If provided, names will be automatically ported to the output.
#' @return List of aggregated impact distributions.
#' @author Created by CDMartinez 23 March 17
#' @note Edited by CDMartinez 23 March 17
#' @export
simAggregate = function(simOutput, rasterSurfaceInputs = NULL, soilLoss = NULL){

  nMC = length(simOutput)
  simDistributions = list()

  ## Surface
  if(is.null(simOutput[[1]]$surfacePadStats) == FALSE){
    simDistributions$statsSurfacePads = list()
    for(i in 1:length(simOutput[[1]]$surfacePadStats)){
      simDistributions$statsSurfacePads[[i]] = matrix(NA,nrow = nMC, ncol = length(simOutput[[1]]$surfacePadStats[[i]]))
      if(is.null(rasterSurfaceInputs) == FALSE){
      colnames(simDistributions$statsSurfacePads[[i]]) = rasterSurfaceInputs[[i]]$ecoNames
      }
    }
    names(simDistributions$statsSurfacePads) =  names(simOutput[[1]]$surfacePadStats)
  }

  ## Roads
  if(is.null(simOutput[[1]]$roadStats) == FALSE){
    rCols = length(simOutput[[1]]$roadStats)
    if(rCols == 1){
      simDistributions$statsRoads = vector(mode = "numeric", length = nMC)
    }else{
      simDistributions$statsRoads = vector(mode = "numeric", length = nMC)
      simDistributions$roadNodes = vector(mode = "list", length = nMC)
    }
  }

  if(is.null(simOutput[[1]]$surfaceRoadStats) == FALSE){
    simDistributions$statsSurfaceRoads = list()
    for(i in 1:length(simOutput[[1]]$surfaceRoadStats)){
      simDistributions$statsSurfaceRoads[[i]] = matrix(NA,nrow = nMC, ncol = length(simOutput[[1]]$surfaceRoadStats[[i]]))
      if(is.null(rasterSurfaceInputs) == FALSE){
      colnames(simDistributions$statsSurfaceRoads[[i]]) = rasterSurfaceInputs[[i]]$ecoNames
      }
    }
    names(simDistributions$statsSurfaceRoads) =  names(simOutput[[1]]$surfaceRoadStats)

    if(is.null(simOutput[[1]]$rusleRoadStats) == FALSE){
      simDistributions$statsRusleRoads = list()
      for(i in 1:3){
        simDistributions$statsRusleRoads[[i]] = matrix(NA, nrow = nMC, ncol = dim(simOutput[[1]]$rusleRoadStats)[1])
        if(is.null(rasterSurfaceInputs) == FALSE){
        colnames(simDistributions$statsRusleRoads[[i]]) = soilLoss$shapeList$ecoNames
        }
      }
      names(simDistributions$statsRusleRoads) = c('Pond', 'Rough', 'Smooth')
    }
  }

  ## RUSLE
  if(is.null(simOutput[[1]]$ruslePadStats) == FALSE){
    simDistributions$statsRuslePads = list()
    for(i in 1:3){
      simDistributions$statsRuslePads[[i]] = matrix(NA, nrow = nMC, ncol = dim(simOutput[[1]]$ruslePadStats)[1])
      if(is.null(rasterSurfaceInputs) == FALSE){
        colnames(simDistributions$statsRuslePads[[i]]) = soilLoss$shapeList$ecoNames
      }
    }
    names(simDistributions$statsRuslePads) = c('Pond', 'Rough', 'Smooth')
  }

  ## Extract and aggregate results
  for(i in 1:nMC){
    temp = simOutput[[i]]

    # If present, extract road length stats
    if(is.null(simDistributions$statsRoads) == FALSE){
      if(rCols == 1){
        simDistributions$statsRoads[i] = temp$roadStats
      }else{
        simDistributions$statsRoads[i] = temp$roadStats$roadLength
        simDistributions$roadNodes[[i]] = temp$roadStats$xyEnds
      }
    }

    # If present, extract rusle stats for pads
    if(is.null(simDistributions$statsRuslePads) == FALSE){
      simDistributions$statsRuslePads[[1]][i,] = t(temp$ruslePadStats[,1])
      simDistributions$statsRuslePads[[2]][i,] = t(temp$ruslePadStats[,2])
      simDistributions$statsRuslePads[[3]][i,] = t(temp$ruslePadStats[,3])
    }

    # If present, extract rusle stats for roads
    if(is.null(simDistributions$statsRusleRoads) == FALSE){
      simDistributions$statsRusleRoads[[1]][i,] = t(temp$rusleRoadStats[,1])
      simDistributions$statsRusleRoads[[2]][i,] = t(temp$rusleRoadStats[,2])
      simDistributions$statsRusleRoads[[3]][i,] = t(temp$rusleRoadStats[,3])
    }

    for(j in 1:length(simDistributions$statsSurfacePads)){
      # Extract surface stats for pads
      simDistributions$statsSurfacePads[[j]][i,] = temp$surfacePadStats[[j]]

      # If present, extract surface stats for roads
      if(is.null(simDistributions$statsSurfaceRoads) == FALSE){
        simDistributions$statsSurfaceRoads[[j]][i,] = temp$surfaceRoadStats[[j]]
      }
    }# end for, surface

  }# end for, nMC
  simDistributions
}

