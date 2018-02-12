## ----Set MC Iterations---------------------------------------------------
# Start by importing the package and setting MC iterations
library(energySim)
nMC <- 1000 # Number of spatial energy Monte Carlo iterations to perform (nMC <= ogMC)
ogMC <- 10000 # Number of USGS Petroleum Assessment iterations to perform (nMC <= ogMC)

## ------------------------------------------------------------------------
set.seed(1999) # for reproducibility
# Run USGS Petroleum Assessment; these numbers are fictitious and not based on a real AU
OGasmt <- continuousAssessment(auMC = ogMC,
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


## ------------------------------------------------------------------------
# Summary of risked outputs
round(continuousAssessmentSummary(OGasmt$risked))

## ------------------------------------------------------------------------
# This example will operate with meters as the distance measure
OGm <- convertAcre2sqMeter(OGasmt) # Convert Gas Assessment to SI units

## ------------------------------------------------------------------------
## Generate 'AU Outline' shape for example

library(sp)
points <- rbind(c(250,250),c(250,1750),c(1750,1750),c(1750,250),c(250,250))
shapeOutline <- SpatialPolygons(list(Polygons(list(Polygon(points)), 'auOutline')))

## A realistic import example

# tpsShapes <- importShape(pathLoc = '/Users/myProject/UTM/',
#                          filename = '2016MancosAUs')
# mainOutline <- tpsShapes[tpsShapes$ASSESSNAME == 'Upper Mancos Tight Gas',]
# auOutline <- SpatialPolygons(list(Polygons(
#   mainOutline@polygons[[1]]@Polygons[4], auName)),
#   proj4string = mainOutline@proj4string)


# Get maximum outline
maxOutline <- getMaxOutline(OGm, shapeOutline)
plot(maxOutline,  border = 'red', axes = TRUE,
     xlab = 'Easting (m)', ylab = 'Northing (m)')
plot(shapeOutline, border = 'black', add = TRUE)
legend('topleft', legend = c('Max Outline','AU Outline'),
       col = c('red','black'), lty = c(1,1), bty = "n")

## ------------------------------------------------------------------------
# Create/import base raster with desired resolution
library(raster)
baseRes <- 25 # raster resolution for our example
# Here, we create rHabitat which in reality could be an imported GAP or NLCD raster
rHabitat <- raster(matrix(sample(1:2, size = (3000/baseRes)^2, replace = TRUE),
                          nrow = (3000/baseRes)),
                   xmn = -500, xmx = 2500, ymn = -500, ymx = 2500)

plot(rHabitat)
plot(shapeOutline, add = TRUE)

# Call to prepareSimSpatial is a must! Output will be used later
prepSpatial <- prepareSimSpatial(surfaceRaster = rHabitat, shapeOutline, OGm)

plot(prepSpatial$rGrid)
plot(shapeOutline, add = TRUE)

## ------------------------------------------------------------------------
summary(prepSpatial)

## ------------------------------------------------------------------------
# Generate joint distribution of wells per pad and pad area
library(MASS)
sigmaIn <- rbind(c(0.06471663, 0.02990047), c(0.02990047, 0.04779028)) # covariance
barsIn  <- c(5.8083691, 1.167067) # means
jointDist  <- exp(mvrnorm(n = nMC, mu = barsIn,Sigma = sigmaIn, empirical = TRUE))
jointDist[,2]  <- ceiling(jointDist[,2]) - 1 # from continuous to discrete
plot(jointDist[,1], jointDist[,2], pch = '.',
     xlab = 'Pad Area (m^2)', ylab = 'Wells Per Pad')
jointDist[which(jointDist[,2] <= 0),2] = 1 # can't have negative numbers
detach("package:MASS", unload = TRUE)

# Call to prepareSimSpatial is a must! Output will be used later
prepDistribution <- prepareSimDistributions(prepSpatial, wellsPerPad = jointDist[,2], 
                                            padArea = jointDist[,1], OGm,
                                            numIterations = nMC)

## ------------------------------------------------------------------------
summary(prepDistribution)

## ------------------------------------------------------------------------
# Prepare surface inputs: example via a shapefile and raster

# Generate example Shapefile
criticalArea <- SpatialPolygons(list(Polygons(list(Polygon(rbind(c(0,500), c(1000,2000),
                                                        c(1500,2000), c(0,0),
                                                        c(0,500)))), 'A'),
                      Polygons(list(Polygon(rbind(c(1300,700),c(1100,950),c(900,700),
                                                  c(1100,450),c(1300,700)))), 'B')))
criticalAreaShape <- SpatialPolygonsDataFrame(criticalArea,
                                              data.frame(Years = rbind(3,5),
                                                         Designation =
                                                           rbind('Wildlife Corridor',
                                                                 'Nest Site'),
                                                         row.names = c('A','B')), 
                                              match.ID = TRUE)

plot(criticalAreaShape, axes = TRUE, main = 'Critical Areas')

## ------------------------------------------------------------------------
# Prepare shapefile for simulation - must call to use input dataset!
prepCriticalArea <- prepareChangeShape(shape = criticalAreaShape,
                                       shapeAttribute = 'Designation',
                                       spatialList = prepSpatial)

## ------------------------------------------------------------------------
# Generate example Raster
rOwnerTypes <- raster(matrix(sample(1:4, 225, replace = TRUE), nrow = 15),
                      xmn = -500, xmx = 2500, ymn = -500, ymx = 2500)

# A realistic example
# # rNLCD <- importRaster('/Users/myProject/UTM/', 'nlcd_utm.tif')

plot(rOwnerTypes, axes = TRUE, main = 'Ownership Type')

## ------------------------------------------------------------------------
# Prepare raster for simulations - must call to use input dataset!
prepOwnerTypes <- prepareChangeRaster(rRaster = rOwnerTypes, spatialList = prepSpatial)

## ------------------------------------------------------------------------
prepHabitat <- prepareChangeRaster(rHabitat,prepSpatial)

## ------------------------------------------------------------------------
# Finally, all surface imports MUST be gathered in a list - this will be used later
prepSurfaces <- list(habitat = prepHabitat, criticalConcern = prepCriticalArea,
                     ownership = prepOwnerTypes)

## ------------------------------------------------------------------------
# Create existing road network (for roads option 1)
nVertices <- 500
road1 <- cbind(seq(0, 2000, length.out = nVertices),
               seq(0, 100, length.out = nVertices)*sin(seq(-pi, 1.5*pi,
                                                           length.out = nVertices)) + 600)
road2 <- cbind(200*cos(seq(-pi, 1.5*pi, length.out = nVertices)) +
                 seq(200, 1800, length.out = nVertices), seq(0, 2000,
                                                             length.out = nVertices))

# Prepare road input: a two-column matrix of (Easting, Northing)
prepRoads <- rbind(road1, road2, cbind(road1[,1],rev(road1[,2]) + 700))
plot(prepRoads, axes = TRUE, main = 'Existing Road Network', xlab = 'Easting (m)',
     ylab = 'Northing (m)', pch = '.')

## ------------------------------------------------------------------------
# Create road segment distribution (for roads option 2)

# Creat a distribution of road segment lengths in meters
# Values must be >= 0
prepRoadDist <- rnorm(nMC, mean = 500, sd = 100)
prepRoadDist[prepRoadDist < 0] <- 0

## ------------------------------------------------------------------------
# RUSLE - create dummy factor rasters

dxdy <- 400

# Slope-length factor
tempI <- matrix(complex( real = rep(seq(0.4, .47, length.out = dxdy), each = dxdy ),
              imag = rep(seq(.3, .42, length.out = dxdy), dxdy)), ncol= dxdy, nrow = dxdy)
tempZ <- 0
for(k in 1:20){tempZ <- tempZ^2+tempI}
rLS <- raster(exp(-abs(tempZ))*20)
extent(rLS) <- extent(prepSpatial$rGrid)
plot(rLS, main = 'LS factor')

# Vegetative cover factor
rC <- raster(matrix(sample(seq(.001, .05, length.out = 7), size = dxdy, replace = TRUE), 
                    nrow = sqrt(dxdy)), xmn = -500, xmx = 2500, ymn = -500, ymx = 2500)
plot(rC, main = 'C factor')

# Soil erodibility factor
rK <- raster(matrix(sample(seq(.01, .7, length.out = 15), size = dxdy, replace = TRUE), 
                    nrow = sqrt(dxdy)), xmn = -500, xmx = 2500, ymn = -500, ymx = 2500)
plot(rK, main = 'K factor')

# Rainfall erosivity
rR <- raster(matrix(rep(10, dxdy), nrow = sqrt(dxdy)), xmn = -500, xmx = 2500,
             ymn = -500, ymx = 2500)
plot(rR, main = 'R factor')

# Soil Loss will be aggregated across entire area
prepRusle <- prepareRusle(prepSpatial, R = rR, K = rK, LS = rLS, C = rC)

# Clean up
rm(dxdy, tempI, tempZ, k)

## ------------------------------------------------------------------------
# Serial example of running simulation on one core, one iteration at a time
nTest <- nMC # number of test iterations, nMC for full simulation
timeA <- Sys.time()
simSerial <- lapply(X = 1:nTest, FUN = simIteration,
                          simDistList = prepDistribution,
                          rasterSurfaceInputs = prepSurfaces,
                          roads = prepRoads,
                          soilLoss = prepRusle)
timeB <- Sys.time()
simSerialAgg <- simAggregate(simSerial, rasterSurfaceInputs =  prepSurfaces,
                             soilLoss = prepRusle)
timeB - timeA # time taken

## ------------------------------------------------------------------------
## Parallel example of running simulation on unix/linux based machine using 2 cores
## Uncomment lines below to run
# library(parallel)
# nTest <- nMC # number of test iterations, nMC for full simulation
# timeA <- Sys.time()
# simOutputMC <- mclapply(X = 1:nTest, FUN = simIteration, mc.cores = 2,
#                           simDistList = prepDistribution,
#                           rasterSurfaceInputs = prepSurfaces,
#                           roads = prepRoads,
#                           soilLoss = prepRusle)
# timeB <- Sys.time()
# simMCagg <- simAggregate(simOutputMC, rasterSurfaceInputs =  prepSurfaces,
                             # soilLoss = prepRusle))
# timeB - timeA # time taken

## ------------------------------------------------------------------------
## Parallel example of running simulation on windows machine using 2 cores
## Uncomment lines below to run
# nTest <- nMC # number of test iterations, nMC for full simulation
# 
# library(parallel)
# nCores <- 2 # use detectCores() to see how many cores are available; suggest using ...
# # nCores <- detectCores()/2
# myCluster <- makeCluster(nCores) # MUST: make a cluster with 'nCores'
# clusterSetRNGStream(myCluster,iseed=5280) # MUST: initialize random number streams
# clusterEvalQ(myCluster,library(energySim)) # MUST: export needed libraries
# timeA <- Sys.time()
# simOutputPar <- parLapply(cl = myCluster, X = 1:nTest, fun = simIteration,
#                           simDistList = prepDistribution,
#                           rasterSurfaceInputs = prepSurfaces,
#                           roads = prepRoads,
#                           soilLoss = prepRusle)
# timeB <- Sys.time()
# stopCluster(myCluster) # MUST: stop the cluster
# # include optional surface or soil loss arguments to transfer naming conventions
# simParAgg <- simAggregate(simOutputPar, rasterSurfaceInputs =  prepSurfaces,
                             # soilLoss = prepRusle)
# timeB - timeA # time taken

## ------------------------------------------------------------------------
# the number of iterations
cat(nTest)
# an item for each iteration, nTest == nMC
length(simSerial)
# a list relating to the a single case (the first iteration) of randomly distributed wells
summary(simSerial[[1]])

# area (in square meters) for each classification converted to oil/gas pads
simSerial[[1]]$surfacePadStats
# total road length, in meters
simSerial[[1]]$roadStats
# total soil loss due to pads (for three types of management practices)
simSerial[[1]]$ruslePadStats 

## ------------------------------------------------------------------------
summary(simSerialAgg)

# Frequency distribution of total road length 
boxplot(simSerialAgg$statsRoads, main = 'Potential Total Road Needed',
        xlab = 'Road length (m)', range = 1.5, outpch = '.', outcol = 'gray',
        horizontal = TRUE)

# Frequency distribution of surface attributes
summary(simSerialAgg$statsSurfacePads)

# habitat
nHabOutputs <- dim(simSerialAgg$statsSurfacePads$habitat)[2]
for(i in 1:nHabOutputs){
  print(paste('Boxplot of Habitat',
                    colnames(simSerialAgg$statsSurfacePads$habitat)[i], 'Conversion'))
  boxplot(simSerialAgg$statsSurfacePads$habitat[,i], 
       main = paste('Habitat',
                    colnames(simSerialAgg$statsSurfacePads$habitat)[i], 'Conversion'),
       xlab = 'Converted area (m^2)', range = 1.5, outpch = '.',
       outcol = 'gray', horizontal = TRUE)
}
# criticalConcern
nCritOutputs <- dim(simSerialAgg$statsSurfacePads$criticalConcern)[2]
for(i in 1:nCritOutputs){
  print(paste('Boxplot of Critical Area',
                    colnames(simSerialAgg$statsSurfacePads$habitat)[i], 'Conversion'))
  boxplot(simSerialAgg$statsSurfacePads$criticalConcern[,i],
       main = paste('Critical Area', 
                    colnames(simSerialAgg$statsSurfacePads$criticalConcern)[i],
                    'Conversion'), range = 1.5, outpch = '.', 
       outcol = 'gray', horizontal = TRUE, xlab = 'Converted area (m^2)')
}

# ownership
nOwnOutputs <- dim(simSerialAgg$statsSurfacePads$ownership)[2]
for(i in 1:nOwnOutputs){
  print(paste('Boxplot of Ownership',
                    colnames(simSerialAgg$statsSurfacePads$habitat)[i], 'Conversion'))
  boxplot(simSerialAgg$statsSurfacePads$ownership[,i],
       main = paste('Ownership',
                    colnames(simSerialAgg$statsSurfacePads$ownership)[i], 
                    'Conversion'), range = 1.5, outpch = '.',
       outcol = 'gray', horizontal = TRUE, xlab = 'Converted area (m^2)')
}

# RUSLE
boxplot(simSerialAgg$statsRuslePads$Pond, main = 'RUSLE from pads: ponding',
        xlab = 'Potential Soil loss (tons/acre/year)', range = 1.5, outpch = '.',
        outcol = 'gray', horizontal = TRUE)
boxplot(simSerialAgg$statsRuslePads$Rough, main = 'RUSLE from pads: rough surface',
        xlab = 'Potential Soil loss (tons/acre/year)', range = 1.5, outpch = '.',
        outcol = 'gray', horizontal = TRUE)
boxplot(simSerialAgg$statsRuslePads$Smooth, main = 'RUSLE from pads: smooth surface',
        xlab = 'Potential Soil loss (tons/acre/year)', range = 1.5, outpch = '.',
        outcol = 'gray', horizontal = TRUE)

