############################ Import Shapefile ############################
#'  Imports a shapefile
#' @description Wrapper function for calling \code{rgdal} function to import a shapefile.
#' @param pathLoc The shapefile directory path.
#' @param filename The file name of the shapefile to import
#' @details Omit extension in filename; expects a filename to have \code{.shp} extension.
#' \code{myShape <- importShape(pathLoc = 'path', filename = 'name')}
#' @return A shapefile
#' @exportClass SpatialPolygonsDataFrame
#' @note Edited by CDMartinez 14 Dec 15
#' @author Created by CDMartinez 03 Dec 15
#' @importFrom rgdal readOGR
#' @export
importShape = function(pathLoc, filename) {
  return(readOGR(pathLoc, filename))
}


############################ Scale Shapefile #############################
#' Scales shapefile.
#' @description Scales input shapefile to desired area assuming stationary centroid.
#' Units must be the same.
#' @param inShape A shapefile
#' @param newArea A number for the desired area
#' @details For more complex shape scaling and transformations,
#' see \code{\link{scaleShapeTransformation}}
#' @return A scaled shapefile with specified area and projection of inShape
#' @note Edited by CDMartinez 14 Dec 15
#' @author Created by CDMartinez 03 Dec 15
#' @importFrom sp SpatialPolygons Polygons Polygon
#' @export
#' @examples library(sp)
#' points <- rbind(c(250,250),c(250,1750),c(1750,1750),c(1750,250),c(250,250))
#' shape <- SpatialPolygons(list(Polygons(list(Polygon(points)), 'auOutline')))
#' shape2 <- scaleShape(shape, 2000000)
#' plot(shape, axes = TRUE)
#' plot(shape2, add = TRUE)
scaleShape = function(inShape, newArea){

  shapeArea = inShape@polygons[[1]]@area
  shapePoints = inShape@polygons[[1]]@Polygons[[1]]@coords
  if(length(inShape@polygons[[1]]@Polygons) > 1){
    stop('More than one polygon found, please parse out single polygon.')
  }

  shapeCenter = inShape@polygons[[1]]@labpt

  xLocal = shapePoints[,1] - shapeCenter[1]
  yLocal = shapePoints[,2] - shapeCenter[2]
  shapeDist = sqrt(xLocal^2 + yLocal^2)
  theta = asin(yLocal/shapeDist)

  pScale = sqrt((newArea/shapeArea))
  newDist = pScale*shapeDist
  xBelowCenter = which(shapePoints[,1] < shapeCenter[1])
  yBelowCenter = which(shapePoints[,2] < shapeCenter[2])
  xNew = newDist*cos(theta) + shapeCenter[1]
  yNew = newDist*sin(theta) + shapeCenter[2]
  xNew[xBelowCenter] = newDist[xBelowCenter]*(-1)*cos(theta[xBelowCenter]) + shapeCenter[1]

  coordsNew = cbind(xNew, yNew)
  newShape = SpatialPolygons(list(Polygons(list(Polygon(coordsNew)), 'mcOutline')),
                             proj4string = inShape@proj4string)

  return(newShape)
}

############################ transform Shapefile #############################
#' Transforms shapefile using an affine transformation.
#' @description Scales input shapefile to desired area. Units must be the same.
#' @param inShape A shapefile
#' @param newArea A number for the desired area of the transformed shape
#' @param percent The decimal percent ratio to transform the shapefile
#' @param orientation The direction to transform the shapefile
#' @param affine A 2D affine matrix to apply to shapefile
#' @return A scaled shapefile with specified area and projection of \code{inShape}
#' @details Supply either the \code{affine} matrix or the \code{percent} and
#' \code{orientation} that will be used to transform the shapefile. If \code{orientation}
#' is one of the cardinal directions then \code{percent} is the decimal ratio for
#' shapefile elongation in the specified direction relative to the perpendicular
#' direction axis. If \code{orientation} is one of the ordinal directions, then
#' \code{percent} is the decimal ratio for shear transformation of the shapefile
#' in the specified direction relative to the perpendiuclar direction axis.
#' @note Edited by CDMartinez 14 June 16
#' @author Created by CDMartinez 11 May 16
#' @importFrom sp SpatialPolygons Polygons Polygon
#' @export
#' @examples library(sp)
#' points <- rbind(c(250,250),c(250,1750),c(1750,1750),c(1750,250),c(250,250))
#' shape <- SpatialPolygons(list(Polygons(list(Polygon(points)), 'auOutline')))
#' shape2 <- scaleShapeTransformation(inShape = shape,newArea = 1500000,
#' percent = .05, orientation='SW')
#' shape3 <- scaleShapeTransformation(inShape = shape,newArea = 1500000,
#' affine = rbind(c(.2,1),c(1,.2)))
#' plot(shape, axes = TRUE)
#' plot(shape2, add = TRUE, border = 'red')
#' plot(shape3, add = TRUE, border = 'blue')
scaleShapeTransformation = function(inShape, newArea, percent = NULL,
                                    orientation = NULL, affine = NULL){

  shapeArea = inShape@polygons[[1]]@area
  shapePoints = inShape@polygons[[1]]@Polygons[[1]]@coords
  if(length(inShape@polygons[[1]]@Polygons) > 1){
    stop('More than one polygon found, please parse out single polygon.')
  }
  shapeCenter = inShape@polygons[[1]]@labpt

  xLocal = shapePoints[,1] - shapeCenter[1]
  yLocal = shapePoints[,2] - shapeCenter[2]
  xyLocal = cbind(xLocal, yLocal)

  if(missing(affine)){
    if(missing(percent)){stop('Missing percent and orientation combination OR affine matrix')}
    if(missing(orientation)){stop('Missing percent and orientation combination OR affine matrix')}
    circle = rbind('NE', 'N', 'NW', 'W' , 'SW', 'S', 'SE', 'E')
    radpi = 22.5*pi/180
    ind = which(circle == orientation)
    if(percent == 0){percent = 0.00000001}
  }else{
    tAffine = affine
    ind = 0
  }

  if(ind == 1) tAffine = rbind(c(1, percent*tan(pi/2 - radpi)), c(0, 1)) #NE
  if(ind == 2) tAffine = rbind(c(1, 0), c(0, 1 + percent)) #N
  if(ind == 3) tAffine = rbind(c(1, 0), c(percent*tan(pi/2 + radpi), 1)) #NW
  if(ind == 4) tAffine = rbind(c(1 + percent, 0), c(0, 1)) #W
  if(ind == 5) tAffine = rbind(c(1, percent*tan(pi/2 - radpi)), c(0, 1)) #SW
  if(ind == 6) tAffine = rbind(c(1, 0), c(0, 1 + percent)) #S
  if(ind == 7) tAffine = rbind(c(1, 0), c(percent*tan(pi/2 + radpi), 1)) #SE
  if(ind == 8) tAffine = rbind(c(1 + percent, 0), c(0, 1)) #E

  #First transform the shape in easting and northing
  tCoords = xyLocal%*%tAffine
  tShape = SpatialPolygons(list(Polygons(list(Polygon(tCoords)), 'Transformed Shape')),
                           proj4string = inShape@proj4string)

  #Then scale the shape for area
  shapeArea = tShape@polygons[[1]]@area
  pScale = sqrt((newArea/shapeArea))
  pAffine = rbind(c(pScale, 0), c(0, pScale))
  coordsNew = tCoords%*%pAffine
  coordsNew[,1] = coordsNew[,1] + shapeCenter[1]
  coordsNew[,2] = coordsNew[,2] + shapeCenter[2]

  newShape = SpatialPolygons(list(Polygons(list(Polygon(coordsNew)),'mcOutline')),
                             proj4string = inShape@proj4string)

  return(newShape)
}

############################ Scale square #############################
#' Scales square shapefile
#' @description Scales input square to desired area. Units must be the same.
#' Suggest using \code{\link{scaleShape}} instead.
#' @param inShape A square shapefile
#' @param newArea A number for the desired area
#' @return A scaled square shapefile with specified area and projection of \code{inShape}
#' @note Edited by CDMartinez 14 Dec 15
#' @author Created by CDMartinez 03 Dec 15
#' @importFrom sp SpatialPolygons Polygons Polygon
#' @export
#' @examples library(sp)
#' points <- rbind(c(250,250),c(250,1750),c(1750,1750),c(1750,250),c(250,250))
#' shape <- SpatialPolygons(list(Polygons(list(Polygon(points)), 'auOutline')))
#' shape2 <- scaleSquareShape(inShape = shape,newArea = 1500000)
#' plot(shape, axes = TRUE)
#' plot(shape2, add = TRUE, border = 'red')
scaleSquareShape = function(inShape, newArea){
  shapeCenter = inShape@polygons[[1]]@labpt
  halfside = sqrt(newArea)/2
  xnew = c(shapeCenter[1] - halfside, shapeCenter[1] - halfside, shapeCenter[1] + halfside,
           shapeCenter[1] + halfside, shapeCenter[1] - halfside)
  ynew = c(shapeCenter[2] - halfside, shapeCenter[2] + halfside, shapeCenter[2] + halfside,
           shapeCenter[2] - halfside, shapeCenter[2] - halfside)
  cnew = cbind(xnew, ynew)
  outOutline = SpatialPolygons(list(Polygons(list(Polygon(cnew)),'rectOutline')),
                               proj4string = inShape@proj4string)
  return(outOutline)
}

########################## Scale Shape 2 Shape ###########################
#' Scale shapefile within shapefile
#' @description Scales input shape to desired area within extent of outer shape.
#' Units must be the same.
#' @param innerShape A shapefile to be scaled
#' @param innerArea A number for the desired area of shape output
#' @param outerShape A shapefile limiting extent of shape scaling
#' @return A scaled shapefile with specified area and same projection as input.
#' @note Edited by CDMartinez 15 Mar 16
#' @author Created by CDMartinez 11 Mar 16
#' @importFrom raster crop
#' @export
#' @examples library(sp)
#' points <- rbind(c(250,250),c(250,1750),c(1750,1750),c(1750,250),c(250,250))
#' ShapeOuter <- SpatialPolygons(list(Polygons(list(Polygon(points)), 'auOutline')))
#' triPoints <- rbind(c(500,500),c(1000,1500),c(1500,500),c(500,500))
#' shape <- SpatialPolygons(list(Polygons(list(Polygon(triPoints)), 'maxOutline')))
#' shape2 <- scaleShape2Shape(innerShape = shape, outerShape = ShapeOuter, innerArea = 1500000)
#' plot(ShapeOuter, axes = TRUE)
#' plot(shape, add = TRUE, border = 'black')
#' plot(shape2, add = TRUE, border = 'red')
scaleShape2Shape = function(innerShape, outerShape, innerArea){
  if (!requireNamespace("rgeos", quietly = TRUE)) {
    stop("rgeos package is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  diffarea = 100
  newInnerShape = innerShape
  # w = 0
  while(diffarea > 0.005){
    testShape = scaleShape(newInnerShape, innerArea)
    newInnerShape = crop(testShape, outerShape)
    diffarea = innerArea - newInnerShape@polygons[[1]]@area
  }
  return(newInnerShape)
}

######################### Get shapefile outline ##########################
#' Get maximum productive area shapefile outline
#' @description Returns maximum productive area outline using basic shape
#' scaling using AU shape centroid.
#' @param auOutline A shapefile
#' @param EA Oil/gas assessment list
#' @return A shapefile with same inShape projection
#' @note Edited by CDMartinez 14 Dec 15
#' @author Created by CDMartinez 03 Dec 15
#' @export
#' @examples library(sp)
#' points <- rbind(c(250,250),c(250,1750),c(1750,1750),c(1750,250),c(250,250))
#' shape <- SpatialPolygons(list(Polygons(list(Polygon(points)), 'auOutline')))
#'
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
#' maxOutline <- getMaxOutline(OGasmt, auOutline = shape)
#' plot(maxOutline, axes=TRUE, border = 'red')
#' plot(shape, add = TRUE, border = 'black')
getMaxOutline = function(EA, auOutline){

  maxArea = max(EA$risked$ProductiveArea)
  maxOutline = scaleShape(auOutline, maxArea)
  return(maxOutline)
}


############################# Import Raster ##############################
#' Imports raster file
#' @description Wrapper function for importing a raster given the path and file name.
#' @param pathLoc The raster path directory
#' @param filename The file name of the raster to import
#' @return Raster
#' @details Include extension in file name. \code{myRaster <-
#' importRaster(pathLoc = 'path', filename = 'name.tiff')}
#' @note Edited by CDMartinez 17 Dec 15
#' @author Created by CDMartinez 17 Dec 15
#' @importFrom raster raster
#' @export
importRaster = function(pathLoc, filename) {
  def = getwd()
  setwd(pathLoc)
  inRast = raster(x = filename)
  setwd(def)

  return(inRast)
}

