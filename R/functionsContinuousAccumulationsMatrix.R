#' USGS Continuous Oil and Gas Assessment
#' @description Performs USGS assessment summarized in Charpentier (2010) as a vectorized Monte Carlo.
#' @param auMC number of MC iterations to perform
#' @param auType type of assessment, string specifier of 'oil' or 'gas'
#' @param auProbability Probability: as a decimal percentage
#' @param auAreaProductive Productive area of accumulation: \code{c(min,mode,max)} [acres]
#' @param auAreaDrainage Uncertainty about average drainage area of wells: \code{c(min,mode,max)} [acres]
#' @param auPercAreaUntested Percentage of total assessment-unit area that is untested: \code{c(min,mode,max)} [percentage]
#' @param auPercAreaSweet Percentage of untested assessment-unit area in sweet spots: \code{c(min,mode,max)} [percentage]
#' @param auPercFutureSS Future success ratio for sweet spots: \code{c(min,mode,max)} [percentage]
#' @param auEURss Uncertainty about sweet spot average EUR (MMBO for oil; BCFG for gas): \code{c(min,med,max)} [MMBO or BCFG]
#' @param auPercFutureNS Future success ratio for non-sweet spots: \code{c(min,mode,max)} [percentage]
#' @param auEURns Uncertainty about non-sweet spot average EUR (MMBO for oil; BCFG for gas): \code{c(min,med,max)} [MMBO or BCFG]
#' @param auGOR Gas/oil ratio (CFG/BO): \code{c(min,mode,max)} [CFG/BO]
#' @param auNGLGR NGL/gas ratio (BNGL/MMCFG): \code{c(min,mode,max)} [BNGL/MMCFG]
#' @param auLGR Liquids/gas ratio (BLIQ/MMCFG): \code{c(min,mode,max)} [BLIQ/MMCFG]
#' @param year Year [XXXX] of factsheet publication of assessment numbers.
#' @details If the AU is to be assessed as entirely sweet spot, simply omit the non-sweet spot variables.
#' Due to the nature of random sampling, the simulation is vectorized and limited only by the available storage.
#' @return A list where each row represents an MC iteration variable combination while the columns contain
#'  the sample distributions. The column names are as follows
#' \itemize{
#' \item\code{DrainageArea} \item\code{EURsweet} \item\code{EURnonSweet} \item\code{ProductiveArea} \item\code{PercAreaUntested}
#' \item\code{PercAreaSweet} \item\code{PercFutureSweet} \item\code{PercFutureNonSweet} \item\code{untestedArea} \item\code{SweetArea}
#' \item\code{NonSweetArea} \item\code{NumSweetWells} \item\code{NumNonSweetWells} \item\code{SweetAccumulation} \item\code{NonSweetAccumulation}
#' \item\code{totalAccumulation}
#' }
#' Depending on the type of assessment, also returns (for 'Gas')
#' \itemize{
#' \item\code{LGR} \item\code{SweetNGLinGas} \item\code{NonSweetNGLinGas} \item\code{totalNGLinGas}
#' }
#' or returns (for 'Oil')
#' \itemize{
#' \item\code{GOR} \item\code{NGLGR} \item\code{SweetGasinOil} \item\code{NonSweetGasinOil}
#' \item\code{totalGasInOil} \item\code{SweetNGLinOil} \item\code{NonSweetNGLinOil} \item\code{totalNGLinOil}
#' }
#' @author Created by CDMartinez 10 Nov 15
#' @note Edited by CDMartinez 24 Nov 15.
#' As of July 2016, the standard 'z-score' for the EUR distribution is implemented as 2.326 (99\%) in
#'  \code{\link{enforceRankCorrelation}}. Assigning \code{year} <= 1 will assume zscore of 3.09 (99.9\%)
#'  as used in \code{\link{conventionalAssessment}} EUR distributions since approximately 2013.
#' @references Charpentier, Ronald R, and Troy Cook, 2010, Improved USGS Methodology for Assessing Continuous Petroleum
#' Resources. U.S. Geological Survey Data Series 547, 2: 22.
#' @examples OGasmt <- continuousAssessment(auMC = 50000,
#' auType = 'Oil',
#' auProbability = 1,
#' auAreaProductive = c(2800000,3100000,3400000),
#' auAreaDrainage = c(320,400,600),
#' auPercAreaUntested = c(80,87,91),
#' auPercAreaSweet = c(24,29,70),
#' auPercFutureSS = c(98,99,100),
#' auEURss = c(0.225,0.25,0.325),
#' auPercFutureNS = c(80,90,95),
#' auEURns = c(0.075,0.15,0.25),
#' auGOR = c(500,1000,1500),
#' auNGLGR = c(35,85,115),
#' year = 2013)
#' round(continuousAssessmentSummary(OGasmt$risked))
#' @importFrom triangle rtriangle
#' @importFrom stats runif
#' @export
continuousAssessment = function(auMC,auType,auProbability,auAreaProductive,auAreaDrainage,auPercAreaUntested,
                                auPercAreaSweet,auPercFutureSS,auEURss,auPercFutureNS=NULL,auEURns=NULL,auGOR=NULL,auNGLGR=NULL,auLGR=NULL,year){

  if(missing(auEURns)){
    cat('Non-Sweet parameters missing. Assessing as ALL SWEET...\n')
    MCorder = enforceRankCorrelation(auMC=auMC,auAreaDrainage=auAreaDrainage,auEURss=auEURss,year=year)
  }else{
    MCorder = enforceRankCorrelation(auMC=auMC,auAreaDrainage=auAreaDrainage,auEURss=auEURss,auEURns=auEURns,year=year)
  }

  ########## Determine coproduct types ##########

  ogFlag = 1
  if(auType =='Oil'){
    ogFlag = 2 # flag for oil
    if(missing(auGOR)){stop("Coproduct inputs missing")}
    if(missing(auNGLGR)){stop("Coproduct inputs missing")}
  }else{
    if(missing(auLGR)){stop("Coproduct inputs missing")}
  }
  cat('Assessing continuous',auType,'...\n')

  ########## Create Storage Variables ##########

  contCalc = matrix(0,nrow=auMC,ncol=16)
  colnames(contCalc) = c("DrainageArea", "EURsweet", "EURnonSweet",
                         "ProductiveArea", "PercAreaUntested", "PercAreaSweet",
                         "PercFutureSweet", "PercFutureNonSweet", "untestedArea",
                         "SweetArea", "NonSweetArea", "NumSweetWells",
                         "NumNonSweetWells", "SweetAccumulation",
                         "NonSweetAccumulation","totalAccumulation")

  contCalc[,1:3] = MCorder
  rm(MCorder)

  contCalc[,4] = rtriangle(auMC, a=auAreaProductive[1], c=auAreaProductive[2], b=auAreaProductive[3])
  contCalc[,5] = rtriangle(auMC, a=auPercAreaUntested[1]/100, c=auPercAreaUntested[2]/100, b=auPercAreaUntested[3]/100)
  contCalc[,6] = rtriangle(auMC, a=auPercAreaSweet[1]/100, c=auPercAreaSweet[2]/100, b=auPercAreaSweet[3]/100)
  contCalc[,7] = rtriangle(auMC, a=auPercFutureSS[1]/100, c=auPercFutureSS[2]/100, b=auPercFutureSS[3]/100)

  if(missing(auPercFutureNS)){
    contCalc[,8] = matrix(0L,nrow=auMC,ncol=1)
  }else{
    contCalc[,8] = rtriangle(auMC, a=auPercFutureNS[1]/100, c=auPercFutureNS[2]/100, b=auPercFutureNS[3]/100)
  }

  cat('Calculating AU',auType,'accumulation ...\n')

  contCalc[,9] = contCalc[,4]*contCalc[,5] # untestedArea = ProductiveArea*PercAreaUntested
  contCalc[,10] = contCalc[,9]*contCalc[,6] # SweetArea = untestedArea*PercAreaSweet
  contCalc[,11] = contCalc[,9]-contCalc[,10] # NonSweetArea = untestedArea-SweetArea
  contCalc[,12] = (contCalc[,10]/contCalc[,1])*contCalc[,7] # SweetWells = (SweetArea/DrainageArea)*PercFutureSweetSuccess
  contCalc[,13] = (contCalc[,11]/contCalc[,1])*contCalc[,8] # NonSweetWells = (NonSweetArea/DrainageArea)*PercFutureNonSweetSuccess

  contCalc[,14] = contCalc[,12]*contCalc[,2] # SweetAccumulation = SweetWells*EURsweet
  contCalc[,15] = contCalc[,13]*contCalc[,3] # NonSweetAccumulation = NonSweetWells*EURnonSweet
  contCalc[,16] = contCalc[,14]+contCalc[,15] # totalAccumulation = SweetAcc + NonSweetAccu

  cat('Calculating AU coproducts ...\n')

  if(ogFlag == 1){
    coproduct = matrix(0L,nrow=auMC,ncol=4)
    colnames(coproduct) = c("LGR","SweetNGLinGas","NonSweetNGLinGas",
                            "totalNGLinGas")

    coproduct[,1] = rtriangle(auMC, a=auLGR[1], c=auLGR[2], b=auLGR[3])
    coproduct[,2] = contCalc[,14]*(coproduct[,1]/1000) # SweetNGL = SweetAccumulation*LGR
    coproduct[,3] = contCalc[,15]*(coproduct[,1]/1000) # NonSweetNGL = NonSweetAccumulation*LGR
    coproduct[,4] = coproduct[,2]+coproduct[,3] #totalNGLinGas = SweetNGL + NonSweetNGL

  }else{
    coproduct = matrix(0L,nrow=auMC,ncol=8)
    colnames(coproduct) = c("GOR","NGLGR", "SweetGasinOil","NonSweetGasinOil",
                            "totalGasInOil","SweetNGLinOil","NonSweetNGLinOil",
                            "totalNGLinOil")

    coproduct[,1] = rtriangle(auMC, a=auGOR[1], c=auGOR[2], b=auGOR[3])
    coproduct[,2] = rtriangle(auMC, a=auNGLGR[1], c=auNGLGR[2], b=auNGLGR[3])

    coproduct[,3] = contCalc[,14]*(coproduct[,1]/1000) # SweetGO = SweetAccumulation*GOR
    coproduct[,4] = contCalc[,15]*(coproduct[,1]/1000) # NonSweetGO = NonSweetAccumulation*GOR
    coproduct[,5] = coproduct[,3]+coproduct[,4] # totalGasinOil = SweetGO + NonSweetGO

    coproduct[,6] = coproduct[,3]*(coproduct[,2]/1000) # SweetNGLO = SweetAccumulation*NGLGR
    coproduct[,7] = coproduct[,4]*(coproduct[,2]/1000) # NonSweetNGLO = NonSweetAccumulation*NGLGR
    coproduct[,8] = coproduct[,6]+coproduct[,7] # totalNGLO = SweetNGLO + NonSweetNGLO
  }

  unrisk = cbind(contCalc,coproduct)
  risk = unrisk
  mcProbability = runif(auMC,min=0,max=1)
  risk[which(mcProbability > auProbability),] = 0L

  cat('Compiling risked and unrisked results ...\n')
  unrisked = as.data.frame(unrisk)
  risked = as.data.frame(risk)
  results = list(risked=risked,unrisked=unrisked)
  #   results = new.env()
  #   results$risked = risked
  #   results$unrisked = unrisked

  cat('Done!\n')
  return(results)
} #endfctn

#' Rank Correlation
#' @description Implements rank correlation method proposed by Iman & Conover (1982).
#' @inheritParams continuousAssessment
#' @details Assumes correlation matrix \code{C} required by USGS Continuous Assessment
#' methodology between mean drainage area, sweet spot EUR, and non-sweet spot EUR
#' \code{C = rbind(c(1,.5,.5),c(0.5,1,0),c(0.5,0,1))}.
#' As of July 2016, the standard 'z-score' for the EUR distribution is implemented as 2.326 (99\%).
#' Assigning \code{year} <= 1 will assume zscore of 3.09 (99.9\%) as used in
#' \code{\link{conventionalAssessment}} EUR distributions since approximately 2013.
#' @return Matrix of correlated variables in each row. Individual columns are \code{DrainageArea},
#' \code{EURsweetSpot},\code{EURnonSweetSpot}.
#' @references Iman, Ronald L, and W J Conover, 1982, A Distribution-free approach to inducing rank
#' correlation among input variables. Communications in statistics - Simulation and Computation 11(3), 311-34.
#' @author Created by CDMartinez 10 Nov 15
#' @note Edited by CDMartinez 08 Dec 15
#' @examples enforceRankCorrelation(auMC = 50000,
#' auAreaDrainage = c(10,20,40),
#' auEURss=c(0.15,0.4,0.65),year=2016)
#' @importFrom EnvStats rlnormTrunc
#' @importFrom triangle rtriangle
#' @importFrom stats cor qnorm
#' @export
enforceRankCorrelation = function(auMC,auAreaDrainage,auEURss,auEURns=NULL,year){

  if(missing(auEURns)){
    nVars = 2 # num of variables to correlate
    C = rbind(c(1,0.5),c(0.5,1)) # correlation to enforce
  }else{
    nVars = 3 # num of variables to correlate
    C = rbind(c(1,0.5,0.5),c(0.5,1,0),c(0.5,0,1)) # correlation to enfore
  }

  P=t(chol(C)) # Cholesky decomposition; transpose to lower-triangular

  r=matrix(nrow=auMC,ncol=1) # generate van der Waerden scores
  for(i in 1:auMC){
    r[i]=qnorm(i/(auMC+1))
  }

  R=matrix(nrow=auMC,ncol=nVars) # randomly sample scores for each var
  for(i in 1:nVars){
    R[,i]=sample(r,auMC,replace=FALSE)
  }

  cat('Generating correlation matrix ...\n')
  ########## Generic implementation ##########
  # RstarA=R%*%t(P)
  # MA=cor(Rstar,method="spearman") # Verify enforced correlation

  ########## More exact forcing of M=C ##########
  T=cor(R)
  Q=t(chol(T)) # cholesky decomposition of sample corr matrix
  S=P%*%solve(Q)

  Rstar=R%*%t(S) # matrix of scores with desired correlation
  # M=cor(Rstar,method="spearman") # Verify enforced correlation

  #  rm(T,Q,S,P,R,r) # Delete unused variables


  ########## Drainage Area Distribution ##########
  mcAreaDrainage = rtriangle(auMC, a=auAreaDrainage[1], c=auAreaDrainage[2], b=auAreaDrainage[3])

  ############# EUR Distributions ###############

  zscore = 2.326 # 3.09 is 99.9% rather than 99% (2.326)
  if(year <= 1){zscore=3.09}

  barEURss = log(auEURss[2]-auEURss[1])
  sigmaEURss = (log(auEURss[3]-auEURss[1])-barEURss)/zscore #3.09
  mcEURss = rlnormTrunc(auMC, meanlog=barEURss, sdlog=sigmaEURss, min=0, max=(auEURss[3]-auEURss[1]))+auEURss[1]

  if(nVars == 3){
    barEURns = log(auEURns[2]-auEURns[1])
    sigmaEURns = (log(auEURns[3]-auEURns[1])-barEURns)/zscore #3.09
    mcEURns = rlnormTrunc(auMC, meanlog=barEURns, sdlog=sigmaEURns, min=0, max=(auEURns[3]-auEURns[1]))+auEURns[1]
  }

  ########## Linking variable values to scores ##########

  cat('Correlating variables ...\n')
  MCorder=matrix(nrow=auMC,ncol=nVars) # matrix for Monte Carlo iterations

  crank=rank(Rstar[,1]) # rank of column scores
  corder=sort(mcAreaDrainage) # sort drainage random values
  MCorder[,1]=corder[crank] # assign MC order based on score ranks

  srank=rank(Rstar[,2])
  sorder=sort(mcEURss)
  MCorder[,2]=sorder[srank]

  if(nVars == 3){
    nrank=rank(Rstar[,3])
    norder=sort(mcEURns)
    MCorder[,3]=norder[nrank]
  }else{
    MCorder = cbind(MCorder,matrix(0L,nrow=auMC,ncol=1))
  }

  colnames(MCorder) = c("DrainageArea","EURsweetSpot","EURnonSweetSpot")

  return(MCorder)
}


#' Unit conversion
#' @description Convert the oil or gas assessment outputs to square meter areas instead of acres.
#' @param EA Oil or gas assessment list
#' @return Oil or gas assessment list with area variables in square meter units.
#' @note Edited by CDMartinez  15 Aug 16
#' @author Created by CDMartinez 15 Aug 16
#' @examples OGasmt <- continuousAssessment(auMC = 50000,
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
#' round(continuousAssessmentSummary(OGasmt$risked))
#' OGmeters = convertAcre2sqMeter(OGasmt)
#' round(continuousAssessmentSummary(OGmeters$risked))
#' @export
convertAcre2sqMeter = function(EA){
  ##### Conversion factors #####

  acre2metersq = 4046.86

  EA$risked$DrainageArea = acre2metersq*EA$risked$DrainageArea
  EA$risked$ProductiveArea = acre2metersq*EA$risked$ProductiveArea
  EA$risked$untestedArea = acre2metersq*EA$risked$untestedArea
  EA$risked$SweetArea = acre2metersq*EA$risked$SweetArea
  EA$risked$NonSweetArea = acre2metersq*EA$risked$NonSweetArea

  EA$unrisked$DrainageArea = acre2metersq*EA$risked$DrainageArea
  EA$unrisked$ProductiveArea = acre2metersq*EA$risked$ProductiveArea
  EA$unrisked$untestedArea = acre2metersq*EA$risked$untestedArea
  EA$unrisked$SweetArea = acre2metersq*EA$risked$SweetArea
  EA$unrisked$NonSweetArea = acre2metersq*EA$risked$NonSweetArea

  EA
}

#' Print USGS Continuous Assessment Summary
#' @description Returns a matrix table of the accumulation mean and fractiles typically published in USGS Factsheets.
#' @param dataList The list returned from \code{\link{continuousAssessment}}
#' @return Printed table in the console
#' @note Edited by CDMartinez  1 Nov 15
#' @author Created by CDMartinez  1 Nov 15
#' @examples OGasmt <- continuousAssessment(auMC = 50000,
#' auType = 'Oil',
#' auProbability = 1,
#' auAreaProductive = c(2800000,3100000,3400000),
#' auAreaDrainage = c(320,400,600),
#' auPercAreaUntested = c(80,87,91),
#' auPercAreaSweet = c(24,29,70),
#' auPercFutureSS = c(98,99,100),
#' auEURss = c(0.225,0.25,0.325),
#' auPercFutureNS = c(80,90,95),
#' auEURns = c(0.075,0.15,0.25),
#' auGOR = c(500,1000,1500),
#' auNGLGR = c(35,85,115),
#' year = 2013)
#' round(continuousAssessmentSummary(OGasmt$risked))
#' @export
continuousAssessmentSummary = function(dataList){

  qInput = c(0.05,0.5,0.95)

  if('GOR' %in% colnames(dataList)){

    accumulation = matrix(0,nrow=3,ncol=4)

    colnames(accumulation) = c('F95','F50','F5','Mean')
    rownames(accumulation) = c('Oil(MMBO)','Gas(BCFG)','NGL(MMBNGL)')

    accumulation[2,1:3] = quantile(dataList$totalGasInOil,probs=qInput)
    accumulation[2,4] = mean(dataList$totalGasInOil)
    accumulation[3,1:3] = quantile(dataList$totalNGLinOil,probs=qInput)
    accumulation[3,4] = mean(dataList$totalNGLinOil)

  }else{

    accumulation = matrix(0,nrow=2,ncol=4)

    colnames(accumulation) = c('F95','F50','F5','Mean')
    rownames(accumulation) = c('Gas(BCFG)','NGL(MMBNGL)')

    accumulation[2,1:3] = quantile(dataList$totalNGLinGas,probs=qInput)
    accumulation[2,4] = mean(dataList$totalNGLinGas)

  }

  accumulation[1,1:3] = quantile(dataList$totalAccumulation,probs=qInput)
  accumulation[1,4] = mean(dataList$totalAccumulation)

  return(accumulation)
}

#' Plot Triangle Distribution
#' @description Plots a density histogram of random sample with specified triangle distribution overlaid.
#' @param plotName Name of data being plotted (e.g. 'Drainage Area')
#' @param plotVarSample Sampled data vector
#' @param plotVarParams Triangle distribution parameters c(min,mode,max)
#' @param bins Number of bins to use for distribution
#' @author Created by CDMartinez 30 Nov 15
#' @note Edited by CDMartinez 30 Nov 15
#' @details Sample distribution bin width = 1 for unrisked or (max-min)/4 for risked
#' @return Plot in default viewer
#' @examples OGasmt = continuousAssessment(auMC = 50000,
#' auType = 'Gas',
#' auProbability = 1,
#' auAreaProductive = c(100,1000,1000000),
#' auAreaDrainage = c(10,20,40),
#' auPercAreaUntested = c(93,96,99),
#' auPercAreaSweet = c(100,100,100),
#' auPercFutureSS = c(80,92,94),
#' auEURss = c(0.15,0.4,0.65),
#' auLGR = c(.08,.5,1),
#' year = 2016)
#'
#' plotTriDistribution('Drainage Area',OGasmt$risked$DrainageArea, c(10,20,40),50)
#' @importFrom triangle dtriangle
#' @importFrom graphics legend
#' @export
plotTriDistribution = function(plotName,plotVarSample,
                            plotVarParams,bins){

  vals=seq(plotVarParams[1],plotVarParams[3],length=100)

  if(missing(bins) == TRUE){
    bins = 1000
    # if(plotType == 'Risked'){
    #   bins=bins/4
    # }
  }

  hist(plotVarSample, breaks=bins, freq=FALSE, xlab=plotName,
       main=c('Distribution of',plotName))
  lines(vals,dtriangle(vals,a=plotVarParams[1],c=plotVarParams[2],b=plotVarParams[3]),col='red',lwd=2)
  legend('topright',legend=c(paste('Sample (nbin=',bins,')'),'Theoretical'),col=c('black','red'),lty=c(1,1),lwd=c(1,2))

}

#' Plot Log-normal Distribution
#' @description Plots a density histogram of random sample with specified log-normal distribution overlaid.
#' @param plotName Name of data being plotted (e.g. 'Drainage Area')
#' @param plotVarSample Sampled data vector
#' @param plotVarParams Log-normal distribution parameters c(min,mode,max)
#' @param bins Number of bins to use for distribution
#' @author Created by CDMartinez 30 Nov 15
#' @note Edited by CDMartinez 30 Nov 15
#' @details Sample distribution bin = 1000 by default
#' @return Plot in default viewer
#' @examples OGasmt = continuousAssessment(auMC = 50000,
#' auType = 'Gas',
#' auProbability = 1,
#' auAreaProductive = c(100,1000,1000000),
#' auAreaDrainage = c(10,20,40),
#' auPercAreaUntested = c(93,96,99),
#' auPercAreaSweet = c(100,100,100),
#' auPercFutureSS = c(80,92,94),
#' auEURss = c(0.15,0.4,0.65),
#' auLGR = c(.08,.5,1),
#' year = 2016)
#'
#' plotLogNormalDistribution('Average EUR',OGasmt$risked$EURsweet, c(0.15,0.4,0.65),100)
#' @importFrom EnvStats dlnormTrunc
#' @importFrom graphics legend
#' @export
plotLogNormalDistribution = function(plotName,plotVarSample,
                               plotVarParams,bins){

  vals=seq(plotVarParams[1],plotVarParams[3],length=1000)-plotVarParams[1]

  barEUR = log(plotVarParams[2]-plotVarParams[1])
  sigmaEUR = (log(plotVarParams[3]-plotVarParams[1])-barEUR)/2.326

  ## Truncated lognormal distribution
  theory = dlnormTrunc(vals, meanlog=barEUR, sdlog=sigmaEUR, min=0, max=(plotVarParams[3]-plotVarParams[1]))#+plotVarParams[1]

  if(missing(bins) == TRUE){
    bins = 1000
    # if(plotType == 'Risked'){
    #   bins=bins/4
    # }
  }

  hist(plotVarSample, breaks=bins, freq=FALSE, xlab=plotName,
       main=c('Distribution of',plotName))
  lines(vals+plotVarParams[1],theory,col='red',lwd=2)
  legend('topright',legend=c(paste('Sample (nbin=',bins,')'),'Theoretical'),col=c('black','red'),lty=c(1,1),lwd=c(1,2))

}

#
# continuousAssessmentAllSweet = function(auMC,auType,auProbability,auAreaProductive,auAreaDrainage,auPercAreaUntested,
#                                 auPercAreaSweet,auPercFutureSS,auEURss,auGOR=NULL,auNGLGR=NULL,auLGR=NULL,year){
#
#   MCorder = enforceRankCorrelation(auMC,auAreaDrainage,auEURss,auEURss,year)
#
#   ########## Determine coproduct types ##########
#
#   ogFlag = 1
#   if(auType =='Oil'){
#     ogFlag = 2 # flag for oil
#     if(missing(auGOR)){stop("Coproduct inputs missing")}
#     if(missing(auNGLGR)){stop("Coproduct inputs missing")}
#   }else{
#     if(missing(auLGR)){stop("Coproduct inputs missing")}
#   }
#   cat('Assessing continuous',auType,'...\n')
#
#   ########## Create Storage Variables ##########
#
#   contCalc = matrix(0,nrow=auMC,ncol=16)
#   colnames(contCalc) = c("DrainageArea", "EURsweet", "EURnonSweet",
#                          "ProductiveArea", "PercAreaUntested", "PercAreaSweet",
#                          "PercFutureSweet", "PercFutureNonSweet", "untestedArea",
#                          "SweetArea", "NonSweetArea", "NumSweetWells",
#                          "NumNonSweetWells", "SweetAccumulation",
#                          "NonSweetAccumulation","totalAccumulation")
#
#
#   contCalc[,1:3] = MCorder
#   rm(MCorder)
#
#   contCalc[,4] = rtriangle(auMC, a=auAreaProductive[1], c=auAreaProductive[2], b=auAreaProductive[3])
#   contCalc[,5] = rtriangle(auMC, a=auPercAreaUntested[1]/100, c=auPercAreaUntested[2]/100, b=auPercAreaUntested[3]/100)
#   contCalc[,6] = rtriangle(auMC, a=auPercAreaSweet[1]/100, c=auPercAreaSweet[2]/100, b=auPercAreaSweet[3]/100)
#   contCalc[,7] = rtriangle(auMC, a=auPercFutureSS[1]/100, c=auPercFutureSS[2]/100, b=auPercFutureSS[3]/100)
#
#   cat('Calculating AU',auType,'accumulation ...\n')
#
#   contCalc[,9] = contCalc[,4]*contCalc[,5] # untestedArea = ProductiveArea*PercAreaUntested
#   contCalc[,10] = contCalc[,9]*contCalc[,6] # SweetArea = untestedArea*PercAreaSweet
#   contCalc[,12] = (contCalc[,10]/contCalc[,1])*contCalc[,7] # SweetWells = (SweetArea/DrainageArea)*PercFutureSweetSuccess
#
#   contCalc[,14] = contCalc[,12]*contCalc[,2] # SweetAccumulation = SweetWells*EURsweet
#   contCalc[,16] = contCalc[,14]+contCalc[,15] # totalAccumulation = SweetAcc + NonSweetAccu
#
#   cat('Calculating AU coproducts ...\n')
#
#   if(ogFlag == 1){
#     coproduct = matrix(0,nrow=auMC,ncol=4)
#     colnames(coproduct) = c("LGR","SweetNGLinGas","NonSweetNGLinGas",
#                             "totalNGLinGas")
#
#     coproduct[,1] = rtriangle(auMC, a=auLGR[1], c=auLGR[2], b=auLGR[3])
#     coproduct[,2] = contCalc[,14]*(coproduct[,1]/1000) # SweetNGL = SweetAccumulation*LGR
#     coproduct[,4] = coproduct[,2]+coproduct[,3] #totalNGLinGas = SweetNGL + NonSweetNGL
#
#   }else{
#     coproduct = matrix(0,nrow=auMC,ncol=8)
#     colnames(coproduct) = c("GOR","NGLGR", "SweetGasinOil","NonSweetGasinOil",
#                             "totalGasInOil","SweetNGLinOil","NonSweetNGLinOil",
#                             "totalNGLinOil")
#
#     coproduct[,1] = rtriangle(auMC, a=auGOR[1], c=auGOR[2], b=auGOR[3])
#     coproduct[,2] = rtriangle(auMC, a=auNGLGR[1], c=auNGLGR[2], b=auNGLGR[3])
#
#     coproduct[,3] = contCalc[,14]*(coproduct[,1]/1000) # SweetGO = SweetAccumulation*GOR
#     coproduct[,5] = coproduct[,3]+coproduct[,4] # totalGasinOil = SweetGO + NonSweetGO
#
#     coproduct[,6] = coproduct[,3]*(coproduct[,2]/1000) # SweetNGLO = SweetAccumulation*NGLGR
#     coproduct[,8] = coproduct[,6]+coproduct[,7] # totalNGLO = SweetNGLO + NonSweetNGLO
#   }
#
#   unrisk = cbind(contCalc,coproduct)
#   risk = unrisk
#   mcProbability = runif(auMC,min=0,max=1)
#   risk[which(mcProbability > auProbability),] = 0
#
#   cat('Compiling risked and unrisked results ...\n')
#   unrisked = as.data.frame(unrisk)
#   risked = as.data.frame(risk)
#   results = list(risked=risked,unrisked=unrisked)
#   #   results = new.env()
#   #   results$risked = risked
#   #   results$unrisked = unrisked
#
#   cat('Done!\n')
#   return(results)
# } #endfctn
