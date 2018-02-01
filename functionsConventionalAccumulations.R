#' USGS Conventional Oil and Gas Assessment
#' @description Performs USGS assessment summarized in Schmoker (2005) as a vectorized Monte Carlo.
#' @param auMC number of simulations
#' @param auType 'Oil' or 'Gas'
#' @param auCharge charge probability [between 0 and 1]
#' @param auRocks rock probability [between 0 and 1]
#' @param auTiming timing probability [between 0 and 1]
#' @param auFieldNum Number of possible undiscovered accumulations; \code{c(min,mode,max)} [unitless]
#' @param auFieldSize Size of undiscovered accumulations; \code{c(min,mode,max)} [unitless]
#' @param auGOR Gas to oil ratio; \code{c(min,mode,max)} [CFG/BO]
#' @param auNGLGR Natural gas liquids to gas ratio; \code{c(min,mode,max)} [BNGL/MMCFG]
#' @param auLGR Liquids to gas ratio; \code{c(min,mode,max)} [BLIQ/MMCFG]
#' @param auOGR Oil to gas ratio; \code{c(min,mode,max)} [BO/MMCFG]
#' @param year The year assessment was published. Mathematical implementation varies based on year.
#' @details Due to the nature of random sampling, the simulation is vectorized and limited only by the available storage. The assessment year determines slightly different implementations:
#' \enumerate{
#' \item \code{year < 2014} uses shifted-truncated lognormal distributions for size of fields only. Triangular for number of fields and coproducts.
#' \item \code{year >= 2014} uses shifted-truncated lognormal distributions for size of fields and number of fields. Triangular for coproducts.
#' \item \code{year < 2013} uses distribution z-score of 2.326 (99\% of the distribution)
#' \item \code{year >= 2013}  uses distribution z-score of 3.09 (99.9\% of the distribution)
#' }
#' @return A list containing three items: list of \code{risked} outputs, list of \code{unrisked}
#' outputs, and vector of Field Sizes (\code{mcFieldSize}). Each row in the
#' \code{risked} and \code{unrisked} lists represent a singel Monte Carlo iteration variable
#' combination while the columns contain the sample distributions.
#'
#' The \code{risked} and \code{unrisked} column names for 'Oil' assessment are as follows:
#' \itemize{
#' \item\code{mcProbability} \item\code{totalAccumulation} \item\code{totalGOR}
#' \item\code{totalNGLGR} \item\code{mcNumOil} \item\code{mcGOR} \item\code{mcNGLGR} \item\code{mcMEL}
#' }
#' The \code{risked} and \code{unrisked} column names for 'Gas' assessment are as follows:
#' \itemize{
#' \item\code{mcProbability} \item\code{totalAccumulation} \item\code{totalOGR}
#' \item\code{totalLGR} \item\code{mcNumGas} \item\code{mcOGR} \item\code{mcLGR} \item\code{mcMEL}
#' }
#' @note Edited by CDMartinez 29 Dec 15
#' @author Created by CDMartinez 28 Dec 15
#' @importFrom triangle rtriangle
#' @importFrom EnvStats rlnormTrunc
#' @references Schmoker, James W, and T R Klett, 2005, U.S. Geological Survey Assessment Concepts
#' for Conventional Petroleum Accumulations. In Petroleum Systems and Geologic Assessment of Oil and Gas in the
#' Uinta-Piceance Province, Utah and Colorado, U.S. Geological Survey Digital Data Series DDS-69-D, 9 pages.
#' @export
#' @examples
#' OGasmt <- conventionalAssessment(auMC = 5000, auType = "Oil", auCharge = 1, auRocks = 1,
#' auTiming = 1, auFieldNum = c(1, 2, 10), auFieldSize = c(.5, .8, 10), auGOR = c(200, 400, 600),
#' auNGLGR = c(35, 85, 115), year = 2013)
conventionalAssessment = function(auMC,auType,auCharge,auRocks,auTiming,auFieldNum,auFieldSize,
                                  auGOR=NULL,auNGLGR=NULL,auLGR=NULL,auOGR=NULL,year){

########## Check for coproduct types ##########

  if(auType =='Oil'){
    if(missing(auGOR)){stop("Coproduct inputs missing")}
    if(missing(auNGLGR)){stop("Coproduct inputs missing")}
  }else{
    if(missing(auLGR)){stop("Coproduct inputs missing")}
    if(missing(auOGR)){stop("Coproduct inputs missing")}
  }

  auProbability = auCharge*auRocks*auTiming
  mcProbability = runif(auMC,min=0,max=1)

############# Assessment ###############

  if(year >= 2013){zscore = 3.09}
  else{zscore = 2.326}

  cat('Assessing with distribution z-scores = ',zscore,'\n')

  barSize = log(auFieldSize[2]-auFieldSize[1])
  sigmaSize = (log(auFieldSize[3]-auFieldSize[1])-barSize)/zscore

  if(year >= 2014){
    barNum = log(auFieldNum[2]-auFieldNum[1])
    sigmaNum = (log(auFieldNum[3]-auFieldNum[1])-barNum)/zscore #99.9% rather than 99% (2.326)

    mcFieldNum = round(rlnormTrunc(auMC, meanlog=barNum, sdlog=sigmaNum, min=0, max=(auFieldNum[3]-auFieldNum[1]))+auFieldNum[1])
    mcMEL = matrix(nrow=auMC,ncol=1)

    cat('Assessing with field number and size distributions as truncated log-normal.\n')
  }else{
    mcFieldNum = round(rtriangle(auMC, a=auFieldNum[1], c=auFieldNum[2], b=auFieldNum[3]))
    mcMEL = matrix(nrow=auMC,ncol=1)

      cat('Assessing with field number distribution as triangular and size distribution as truncated log-normal.\n')
  }
############# Total Accumulation ###############

  totalFields = sum(mcFieldNum)
  mcFieldSize = rlnormTrunc(totalFields, meanlog=barSize, sdlog=sigmaSize, min=0, max=(auFieldSize[3]-auFieldSize[1]))+auFieldSize[1]
  totalAccumulation = matrix(nrow=auMC,ncol=1)
  totalAccumulation[1] = sum(mcFieldSize[1:mcFieldNum[1]])
  mcMEL[1] = max(mcFieldSize[1:mcFieldNum[1]])
  index = mcFieldNum[1] + 1

  for(i in 2:auMC){
    lastind = index + mcFieldNum[i] - 1
    totalAccumulation[i] = sum(mcFieldSize[index:lastind])
    mcMEL[i] = max(mcFieldSize[index:lastind])
    index = lastind + 1
  }

############# Coproduct Calculation ###############
  if(auType == 'Oil'){

    mcGOR = rtriangle(auMC, a=auGOR[1], c=auGOR[2], b=auGOR[3])
    mcNGLGR = rtriangle(auMC, a=auNGLGR[1], c=auNGLGR[2], b=auNGLGR[3])
    totalGOR = totalAccumulation*(mcGOR)/1000
    totalNGLGR = totalGOR*(mcNGLGR)/1000

    unrisk = cbind(mcProbability,totalAccumulation,totalGOR,
                   totalNGLGR,mcFieldNum,mcGOR,mcNGLGR,mcMEL)
    colnames(unrisk) = c('mcProbability','totalAccumulation','totalGOR',
                         'totalNGLGR','mcNumOil','mcGOR','mcNGLGR','mcMEL')
    risk = unrisk
    risk[which(mcProbability > auProbability),] = 0
    risked =  as.data.frame(risk)
    unrisked = as.data.frame(unrisk)

  }
  if(auType == 'Gas'){

    mcLGR = rtriangle(auMC, a=auLGR[1], c=auLGR[2], b=auLGR[3])
    mcOGR = rtriangle(auMC, a=auOGR[1], c=auOGR[2], b=auOGR[3])
    totalLGR = totalAccumulation*(mcLGR)/1000
    totalOGR = totalAccumulation*(mcOGR)/1000

    unrisk = cbind(mcProbability,totalAccumulation,totalOGR,
                   totalLGR,mcFieldNum,mcOGR,mcLGR,mcMEL)
    colnames(unrisk) = c('mcProbability','totalAccumulation','totalOGR',
                         'totalLGR','mcNumGas','mcOGR','mcLGR','mcMEL')

    risk = unrisk
    risk[which(mcProbability > auProbability),] = 0
    risked =  as.data.frame(risk)
    unrisked = as.data.frame(unrisk)
  }
  results = list(risked=risked,unrisked=unrisked,mcFieldSize=mcFieldSize)
  return(results)
}

#' Print USGS Conventional Assessment Summary
#' @description Returns a matrix table of the accumulation mean and fractiles typically published in USGS Factsheets.
#' @param dataList The \code{risked} or \code{unrisked} list returned from \code{\link{conventionalAssessment}}
#' @return Printed table in the console
#' @note Edited by CDMartinez  29 Dec 15
#' @author Created by CDMartinez  1 Nov 15
#' @export
#' @examples
#' OGasmt <- conventionalAssessment(auMC = 5000, auType = "Oil", auCharge = 1, auRocks = 1,
#' auTiming = 1, auFieldNum = c(1, 2, 10), auFieldSize = c(.5, .8, 10), auGOR = c(200, 400, 600),
#' auNGLGR = c(35, 85, 115), year = 2013)
#' round(conventionalAssessmentSummary(OGasmt$risked))
#' round(conventionalAssessmentSummary(OGasmt$unrisked))
conventionalAssessmentSummary = function(dataList){

  qInput = c(0.05,0.5,0.95)

  if('totalGOR' %in% colnames(dataList)){

    accumulation = matrix(0,nrow=3,ncol=4)

    colnames(accumulation) = c('F95','F50','F5','Mean')
    rownames(accumulation) = c('Oil(MMBO)','Gas(BCFG)','NGL(MMBNGL)')

    accumulation[2,1:3] = quantile(dataList$totalGOR,probs=qInput)
    accumulation[2,4] = mean(dataList$totalGOR)
    accumulation[3,1:3] = quantile(dataList$totalNGLGR,probs=qInput)
    accumulation[3,4] = mean(dataList$totalNGLGR)

  }else{

    accumulation = matrix(0,nrow=2,ncol=4)

    colnames(accumulation) = c('F95','F50','F5','Mean')
    rownames(accumulation) = c('Gas(BCFG)','NGL(MMBNGL)')

    accumulation[2,1:3] = quantile(dataList$totalLGR,probs=qInput)
    accumulation[2,4] = mean(dataList$totalLGR)

  }

  accumulation[1,1:3] = quantile(dataList$totalAccumulation,probs=qInput)
  accumulation[1,4] = mean(dataList$totalAccumulation)

  return(accumulation)
}
