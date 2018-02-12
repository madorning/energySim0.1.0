# placePads

##Order of ops
# continuousAssessment()
# convertAcre2sqMeter()
# prepareSimSpatial()
# prepareSimDistributions()
# for(i in MC){
# placePads()
# roads()
# ecoCompare1()
# ecoCompare2()
# ecoCompare3()

########################## placePads ################################
#' Place pads for one set of sampled oil and gas parameters
#' @description Places pads for a single Monte Carlo iteration.
#' @param simList List returned from \code{\link{prepareSimDistributions}}.
#' @param entry The index of constant values to extract from \code{simList}.
#' @details Optimal use of this single iteration function is
#' utilizing a preset function \code{\link{simIteration}}
#' @return Returns a list with the following items: rDrained,mPads,xyPadCenter,mcOutline,mcSweet
#' @note Edited by CDMartinez  14 Feb 17. Replaces runSingleIterationSweet
#' which is original implementation.
#' @author Created by CDMartinez 15 Mar 16
#' @import raster
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
#' values(rBase) <- seq(1:40000)
#'
#' points <- rbind(c(250,250),c(250,1750),c(1750,1750),c(1750,250),c(250,250))
#' shape <- SpatialPolygons(list(Polygons(list(Polygon(points)), 'auOutline')))
#'
#' plot(rBase, xlim = c(0,2000), ylim = c(0,2000))
#' lines(shape)
#'
#' spatialPrep <- prepareSimSpatial(surfaceRaster = rBase, shape, OGasmt)
#' distributionPrep <- prepareSimDistributions(spatialPrep,wellsPerPad = 3,
#' padArea = 500, EA = OGasmt, numIterations=5)
#'
#' disturbance <- placePads(distributionPrep, 5)
#' rPads = raster(disturbance$mPads)
#' extent(rPads) = extent(spatialPrep$rGrid)
#' plot(rPads)
#' points(disturbance$xyPadCenter, col = 'white', pch = '.')
#' legend('topleft', legend = 'Gas Pad', col = terrain.colors(12), pch = 15)
#' @export
placePads = function(simList, entry){

  # Extract needed variables from list

  ecoGridInput = simList$rGrid
  assessmentVars = simList$assessmentVars[entry,]
  constantVars = simList$constantVars
  rastBoundBox = simList$rastBoundBox
  anOutline = simList$anOutline
  inOutline = simList$inOutline

  # If present in USGS NOGA Assessment
  if(is.null(inOutline) == TRUE){
    sweetFlag = FALSE
  }else{
    sweetFlag = TRUE
  }

  minx = constantVars$minx
  maxx = constantVars$maxx
  miny = constantVars$miny
  maxy = constantVars$maxy
  deltax = constantVars$deltax
  deltay = constantVars$deltay
  mPads = constantVars$mPads
  cellArea = constantVars$cellArea

  #   NumSweetWells = unlist(assessmentVars[['NumSweetWells']])
  #   NumNonSweetWells = unlist(assessmentVars[['NumNonSweetWells']])
  ProductiveArea = unlist(assessmentVars[['ProductiveArea']])
  DrainageArea = unlist(assessmentVars[['DrainageArea']])
  drainCellSize = unlist(assessmentVars[['drainCellSize']])
  nWells = unlist(assessmentVars[['nWells']])
  nPads = unlist(assessmentVars[['nPads']])
  xextent = unlist(assessmentVars[['xextent']])
  yextent = unlist(assessmentVars[['yextent']])
  tracking = unlist(assessmentVars[['mcNumber']])
  # nWellsPerPad = unlist(assessmentVars[['nWellsPerPad']])
  # nPadArea = unlist(assessmentVars[['nPadArea']])
  halfBuffer = unlist(assessmentVars[['halfBuffer']])
  cellBuffer = unlist(assessmentVars[['cellBuffer']])
  cellBufferDist = unlist(assessmentVars[['cellBufferDist']])

  ## Scale AU Outline area to productive area; centroid hard-wired for now
  mcOutline = scaleShape(anOutline, ProductiveArea)

  ## If sweet spot is available, get constant variables
  if(sweetFlag == TRUE){
    SweetArea = unlist(assessmentVars[['SweetArea']])
    nPadsSweet = unlist(assessmentVars[['nPadsSweet']])
    nPadsNonSweet = unlist(assessmentVars[['nPadsNonSweet']])

    mcSweet = scaleShape2Shape(inOutline, mcOutline, SweetArea)
  }

  ## Generate grid of square drainage areas
  origin = runif(n = 2, min = (-drainCellSize/2), max = (drainCellSize)/2)
  xorigin = origin[1] + minx
  yorigin = origin[2] + miny
  xextent = xextent + xorigin
  yextent = yextent + yorigin
  rTemp = raster(resolution = drainCellSize, crs = anOutline,
                 ymn = yorigin, xmn = xorigin, ymx = yextent, xmx = xextent)
  values(rTemp) = 1L
  rDrained = mask(rTemp, mcOutline)

  ## Randomly sample drained cells
  if(sweetFlag == TRUE){
    rSweet = mask(rTemp, mcSweet)
    sweetCells = which(rSweet@data@values >= 0L)
    rDrained[sweetCells] = 2L

    nonSweetCells = which(rDrained@data@values == 1L)
    excessSweet = nPadsSweet - length(sweetCells)
    if(excessSweet <= 0){
      sweetDrainedPlaces = sample(sweetCells, size = nPadsSweet, replace = FALSE)
      nonSweetDrainedPlaces = sample(nonSweetCells, size = nPadsNonSweet, replace = FALSE)
    }else{
      bulkTotal = nPadsNonSweet + excessSweet
      bulkNonSweet = sample(nonSweetCells, size = bulkTotal, replace = FALSE)
      nonSweetDrainedPlaces = bulkNonSweet[1:nPadsNonSweet]
      sweetDrainedPlaces = c(sweetCells, bulkNonSweet[(nPadsNonSweet+1):bulkTotal])
    }
    rDrained[nonSweetDrainedPlaces] = 3L
    rDrained[sweetDrainedPlaces] = 4L
    drainedPlaces = c(nonSweetDrainedPlaces, sweetDrainedPlaces)
  }else{
    within = which(rDrained@data@values >= 1L)
    drainedPlaces = sample(within, size = nPads, replace = FALSE)
    rDrained[drainedPlaces] = 4L
  }

  ## Sample local upper-left origin of pad location on surface raster
  rangeXY = floor((drainCellSize)/deltay) - halfBuffer
  # halfBuffer is to prevent right/bottom edge placement
  localRectRange = seq(from = 1, to = rangeXY, by = 1)
  xyPadCells = sample(localRectRange, size = 2*nPads, replace = TRUE)
  xyCells = matrix(xyPadCells, nrow = nPads, ncol = 2, byrow = TRUE)

  ## Local to real cell coordinates
  xPadLocal = xyCells[,1]*deltax
  yPadLocal = xyCells[,2]*deltay

  drainLocs = xyFromCell(rDrained, cell = drainedPlaces)
  xPad = drainLocs[,1] - floor(drainCellSize/2) + xPadLocal
  # Don't need to do floor here - since it's real coordinates, but ensures within cell interior
  yPad = drainLocs[,2] - floor(drainCellSize/2) + yPadLocal

  ## Check/correct upper-left origin doesn't cause pad to be placed outside raster
  xLarge = which(xPad >= (rastBoundBox@xmax - cellBufferDist))
  if(length(xLarge) > 0){ xPad[xLarge] = rastBoundBox@xmax - cellBufferDist}
  xSmall = which(xPad <= rastBoundBox@xmin + cellBufferDist)
  if(length(xSmall) > 0){ xPad[xSmall] = rastBoundBox@xmin + cellBufferDist}
  yLarge = which(yPad >= (rastBoundBox@ymax - cellBufferDist))
  if(length(yLarge) > 0){ yPad[yLarge] = rastBoundBox@ymax - cellBufferDist}
  ySmall = which(yPad <= (rastBoundBox@ymin + cellBufferDist))
  if(length(ySmall) > 0){ yPad[ySmall] = rastBoundBox@ymin + cellBufferDist}

  originCells = cellFromXY(ecoGridInput, cbind(xPad, yPad))
  rowCells = rowFromCell(ecoGridInput, originCells)
  colCells = colFromCell(ecoGridInput, originCells)

  ## Create binary matrix of pad placement
  for(jj in 1:nPads){
    yo = rowCells[jj]
    xo = colCells[jj]
    mPads[yo:(yo + cellBuffer - 1),xo:(xo + cellBuffer - 1)] = 1L
  }

  xyPadCenter = cbind((xPad + (cellBufferDist)/2), (yPad - (cellBufferDist)/2))

  if(sweetFlag == TRUE){
    return(list(rDrained = rDrained, mPads = mPads, xyPadCenter = xyPadCenter,
                mcOutline = mcOutline, mcOutlineSweet = mcSweet))
  }else{
    return(list(rDrained = rDrained, mPads = mPads, xyPadCenter = xyPadCenter,
                mcOutline = mcOutline))
  }
  #  save(rDrained,mPads,xyPadCenter,mcOutline,mcSweet,file=paste(tracking,'spatial'))

}


########################## makeRoads K-D Tree ##############################
#' Place multiple linear roads
#' @description Default returns total 'road length' by connecting given new locations (\code{xyStarts}) and
#' known existing locations (\code{roadNodes}) using Euclidean distance. Wrapper function for a call to
#' \code{nn2} which uses a K-D tree from the \code{\link{RANN}} package.
#' @param xyStarts Matrix of easting (1st column) and northing (2nd column) of new pad center locations.
#' Corresponding eastings and northings assumed to be in the same row.
#' @param roadNodes Matrix of existing vertices of current road network. Easting (1st column) and
#' northing (2nd column) pairs. Corresponding eastings and northings assumed to be in the same row.
#' @param totalLength Default value is \code{TRUE} and only the total length of roads is returned.
#' If \code{FALSE}, a list \code{xyEnds} containing the ending xy coordinates from the
#' \code{roadNodes} are returned for each pad in addition to the total length of roads.
#' @return If \code{totalLength} is \code{TRUE}, the total length of road using a
#' Euclidean distance between starting and ending locations is returned.
#' If \code{totalLength} is \code{FALSE}, a list \code{xyEnds} containing the ending xy coordinates
#' from the \code{roadNodes} are returned for each pad in addition to the total length of roads.
#' @note Edited by CDMartinez  15 Mar 16
#' @author Created by CDMartinez 15 Mar 16
#' @details A K-D tree, or K-dimensional tree, is a binary search tree data structure used for organizing points
#' that can be used for quick and efficient sorting and nearest neighbor searching of large sets of points.
#' @import RANN
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
#' spatialPrep <- prepareSimSpatial(surfaceRaster = rBase, shape, OGasmt)
#' distributionPrep <- prepareSimDistributions(spatialPrep,wellsPerPad = 3,
#' padArea = 500, EA = OGasmt, numIterations=5)
#'
#' # Create road network
#' nVertices <- 500
#' road1 <- cbind(seq(0, 2000, length.out = nVertices),
#' seq(0, 100, length.out = nVertices)*sin(seq(-pi, 1.5*pi, length.out = nVertices)) + 600)
#' road2 <- cbind(200*cos(seq(-pi, 1.5*pi, length.out = nVertices)) +
#' seq(200, 1800, length.out = nVertices), seq(0, 2000, length.out = nVertices))
#' # Prepare road input: a two-column matrix of (Easting, Northing)
#' prepRoads <- rbind(road1, road2, cbind(road1[,1],rev(road1[,2]) + 700))
#'
#' pads <- placePads(distributionPrep, 5)
#' roadLength <- makeRoads(xyStarts = pads$xyPadCenter, roadNodes = prepRoads)
#' @export
makeRoads = function(xyStarts,roadNodes,totalLength = TRUE){
  roadKD = nn2(roadNodes,query = xyStarts,k=1,treetype='kd',searchtype = 'standard',eps=0.0)

  if(totalLength == TRUE){
    roadReturn = sum(roadKD$nn.dists)
  }else{
    roadReturn = list(xyEnds = roadNodes[roadKD$nn.idx,], roadLength = sum(roadKD$nn.dists))
  }
  roadReturn
}

########################## makeRoadsD - from Distribution ##############################
#' Place multiple linear roads
#' @description Default returns total 'road length' by drawing a segment length from a user provided
#' distribution (\code{roadNodes}) for each simulated pad (\code{xyStarts}) 
#' @param xyStarts Matrix of easting (1st column) and northing (2nd column) of new pad center locations.
#' Corresponding eastings and northings assumed to be in the same row.
#' @param roadNodes Matrix of existing vertices of current road network. Easting (1st column) and
#' northing (2nd column) pairs. Corresponding eastings and northings assumed to be in the same row.
#' @param totalLength Default value is \code{TRUE} and only the total length of roads is returned.
#' If \code{FALSE}, a list \code{xyEnds} containing the ending xy coordinates from the
#' \code{roadNodes} are returned for each pad in addition to the total length of roads.
#' @return If \code{totalLength} is \code{TRUE}, the total length of road using a
#' Euclidean distance between starting and ending locations is returned.
#' If \code{totalLength} is \code{FALSE}, a list \code{xyEnds} containing the ending xy coordinates
#' from the \code{roadNodes} are returned for each pad in addition to the total length of roads.
#' @note Edited by CDMartinez  15 Mar 16
#' @author Created by CDMartinez 15 Mar 16
#' @details A K-D tree, or K-dimensional tree, is a binary search tree data structure used for organizing points
#' that can be used for quick and efficient sorting and nearest neighbor searching of large sets of points.
#' @import RANN
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
#' spatialPrep <- prepareSimSpatial(surfaceRaster = rBase, shape, OGasmt)
#' distributionPrep <- prepareSimDistributions(spatialPrep,wellsPerPad = 3,
#' padArea = 500, EA = OGasmt, numIterations=5)
#'
#' # Create road network
#' nVertices <- 500
#' road1 <- cbind(seq(0, 2000, length.out = nVertices),
#' seq(0, 100, length.out = nVertices)*sin(seq(-pi, 1.5*pi, length.out = nVertices)) + 600)
#' road2 <- cbind(200*cos(seq(-pi, 1.5*pi, length.out = nVertices)) +
#' seq(200, 1800, length.out = nVertices), seq(0, 2000, length.out = nVertices))
#' # Prepare road input: a two-column matrix of (Easting, Northing)
#' prepRoads <- rbind(road1, road2, cbind(road1[,1],rev(road1[,2]) + 700))
#'
#' pads <- placePads(distributionPrep, 5)
#' roadLength <- makeRoads(xyStarts = pads$xyPadCenter, roadNodes = prepRoads)
#' @export
makeRoadsD = function(xyStarts, roadDist, totalLength = TRUE){
  # roadDist = input distribution sample from user (similar to prepareSimDistributions)
  
  # draw npads times
  nPads = nrow(xyStarts) # need to check and make sure this provides the correct number of pads
  roadSegs = sample(roadDist, size = nPads, replace=TRUE) 
  
  # sum total length
  if(totalLength == TRUE){
    roadReturn = sum(roadSegs)
  }else{
    roadReturn = list(roadSegLengths = roadSegs, roadLength = sum(roadSegs))
  }
  roadReturn
}

#' Calculate soil loss according to RUSLE
#' @description Calculates raster (cell) based soil loss due to either pads or roads.
#' @param rusleIn List returned from \code{\link{prepareRusle}}
#' @param padsIn List output from \code{\link{placePads}}
#' @param roadsIn Optional. List returned from \code{\link{placeRoads}}.
#' @details If argument is provided to \code{roadsIn}, output will be for roads only.
#' To get output for pads, \code{roadsIn} must be omitted. Vegetative cover \code{C} will
#' be changed from input classification to 0.2 to represent development. Soil loss will
#' be calculated for three variations of management practice \code{P}: P = 0.1 for seeded
#' and sediment dentention pond; P = 0.18 for seeded and rough surface; P = 0.26 for
#' seeded and smooth surface. Changes in \code{C} and \code{P} are based on values suggested
#' by Linard et. al. (2014) in USGS Open-File Report 2014-1158 (https://dx.doi.org/10.3133/ofr20141158).
#' @return Total soil loss change compared to baseline.
#' @author Created by CDMartinez 21 Oct 16
#' @note Edited by CDMartinez 2 Feb 17
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
#' distributionPrep <- prepareSimDistributions(prepSpatial,wellsPerPad = 3,
#' padArea = 500, EA = OGasmt, numIterations=5)
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
#' rC <- raster(matrix(sample(seq(.001, .05, length.out = 7), size = dxdy,
#' replace = TRUE), nrow = sqrt(dxdy)), xmn = -500, xmx = 2500, ymn = -500, ymx = 2500)
#'
#' # Soil erodibility factor
#' rK <- raster(matrix(sample(seq(.01, .7, length.out = 15), size = dxdy,
#' replace = TRUE), nrow = sqrt(dxdy)), xmn = -500, xmx = 2500, ymn = -500, ymx = 2500)
#'
#' # Rainfall erosivity
#' rR <- raster(matrix(rep(10, dxdy), nrow = sqrt(dxdy)), xmn = -500,
#' xmx = 2500, ymn = -500, ymx = 2500)
#'
#' # Soil Loss will be aggregated across entire area
#' prepRusle <- prepareRusle(prepSpatial, R = rR, K = rK, LS = rLS, C = rC)
#'
#' pads <- placePads(distributionPrep, 5)
#' soilLoss <- rusle(rusleIn = prepRusle, padsIn = pads)
#' @export
rusle = function(rusleIn, padsIn, roadsIn = NULL){
  if(missing(rusleIn) == TRUE){ stop('Missing input RUSLE parameters')}
  if(missing(padsIn) == TRUE){ stop('Missing input simulation parameters')}

  # Determine pads or roads...
  if(missing(roadsIn) == TRUE){
    inds2replace = which(padsIn$mPads == 1) # pads
    scaleOutput = 1
  }else{
    inds2replace = roadsIn$roadCells # roads
    scaleOutput = roadsIn$cellProportion
  }
  # mRKLSC = rusleIn$mRKLSC
  # mC = rusleIn$mC
  # rGrid = rusleIn$rGrid
  # mPads = rusleIn$mPads

  #practice = 'smooth'
  #rConstant = R*K*LS

  # vegetative cover is now developed
  C = 0.2
  # sediment detention pond, rough surface, or smooth surface
  #if(practice == 'pond'){P = 0.1}
  #if(practice == 'rough'){P = 0.18}
  #if(practice == 'smooth'){P = 0.26}
  #develop = C*P

  ppond = 0.1*C
  prough = 0.18*C
  psmooth = 0.26*C

  # mRusle = mRKSLC*(mPads-1)*(-1)
  # mPond = mRusle
  # mRough = mRusle
  # mSmooth = mRusle

  rusleBase = rusleIn$mRKLSC[inds2replace]
  repValues = (rusleBase/rusleIn$mC[inds2replace])
  # print(rusleBase)

  diffPond = ((repValues)*ppond - rusleBase)*scaleOutput
  diffRough = ((repValues)*prough - rusleBase)*scaleOutput
  diffSmooth = ((repValues)*psmooth - rusleBase)*scaleOutput

  if(is.null(rusleIn$shapeList) == TRUE){
    soilLoss = cbind(sum(diffPond, na.rm = TRUE),
                     sum(diffRough, na.rm = TRUE),
                     sum(diffSmooth, na.rm = TRUE))
  }else{
    mShapes = rusleIn$shapeList$mEco[inds2replace]
    values = rusleIn$shapeList$ecoItems
    soilLoss = matrix(NA, nrow = length(values), ncol = 3)
    for(i in 1:length(values)){
      inds = which(mShapes == values[i])
      soilLoss[i,1] = sum(diffPond[inds], na.rm = TRUE)
      soilLoss[i,2] = sum(diffRough[inds], na.rm = TRUE)
      soilLoss[i,3] = sum(diffSmooth[inds], na.rm = TRUE)
    }
  }

  colnames(soilLoss) = c('pond', 'rough', 'smooth')
  return(soilLoss)
}


########################## placeRoads ################################
#' Place roads for surface impacts
#' @description Places roads randomly for a single Monte Carlo iteration.
#' @param padsIn List returned from call to \code{\link{placePads}}
#' @param simList List returned from \code{\link{prepareSimDistributions}}.
#' @param totalRoadLength The total length of road to distribute on the landscape. Direct
#' output from \code{\link{makeRoads}} can be used as an input here.
#' @param roadWidth Optional. Number specifying road width. Default is 10 meters.
#' @param cellProportion Optional. The proportion of a single surface raster cell that can be
#' covered by road. Default is 2/3 of the cell area.
#' @details The function \code{\link{makeRoads}} simply returns the total length and
#' optionally the coordinates of the road. This function will take that total road length and
#' randomly distribute it across the landscape at a set road width for comparison to any
#' surface raster. Road segments are randomly distributed and there is no connectivity or
#' continuity implied. If optional parameters are omitted, roads may cover a maximum of 2/3
#' of the area of a single surface cell at specified resolution.
#' @return List containing road cells and area to be used in \code{\link{impactsRoads}} or \code{\link{rusle}}
#' @note Edited by CDMartinez  27 Feb 17
#' @author Created by CDMartinez 27 Feb 17
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
#' spatialPrep <- prepareSimSpatial(surfaceRaster = rBase, shape, OGasmt)
#' distributionPrep <- prepareSimDistributions(spatialPrep,wellsPerPad = 3,
#' padArea = 500, EA = OGasmt, numIterations=5)
#'
#' # Create road network
#' nVertices <- 500
#' road1 <- cbind(seq(0, 2000, length.out = nVertices),
#' seq(0, 100, length.out = nVertices)*sin(seq(-pi, 1.5*pi, length.out = nVertices)) + 600)
#' road2 <- cbind(200*cos(seq(-pi, 1.5*pi, length.out = nVertices)) +
#' seq(200, 1800, length.out = nVertices), seq(0, 2000, length.out = nVertices))
#' # Prepare road input: a two-column matrix of (Easting, Northing)
#' prepRoads <- rbind(road1, road2, cbind(road1[,1],rev(road1[,2]) + 700))
#'
#' pads <- placePads(distributionPrep, 5)
#' roadLength <- makeRoads(xyStarts = pads$xyPadCenter, roadNodes = prepRoads)
#' roads <- placeRoads(padsIn = pads, simList = distributionPrep, totalRoadLength = roadLength)
#'
#' @export
placeRoads = function(padsIn, simList, totalRoadLength, roadWidth = NULL, cellProportion = NULL){

  if(is.null(roadWidth) == TRUE){
    roadWidth = 10
  }
  if(is.null(cellProportion) == TRUE){
    cellProportion = 2/3 # max road in a cell
  }

  cellArea = simList$rastRes[1]*simList$rastRes[2]
  roadLengthPerCell = (cellArea*cellProportion)/roadWidth
  roadSegmentPerCell = round(roadLengthPerCell/max(simList$rastRes[1],simList$rastRes[2]))
  roadArea = (roadLengthPerCell/roadSegmentPerCell)*roadWidth

  mBase = simList$constantVars$mPads
  nRoads = round(totalRoadLength/roadLengthPerCell)
  mMask = mBase + padsIn$mPads
  emptyCells = which(mMask == 0)
  emptyCellsDivided = rep(emptyCells, roadSegmentPerCell) # two 30m roads in a cell

  roadCells = sample(emptyCellsDivided, size = nRoads)

  return(list(roadCells = roadCells, cellProportion = cellProportion, roadArea = roadArea))
}

