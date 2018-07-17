
buildSamplesFlows <- function(i){
  dummy <- t(sapply(1:sizeGrid, function(x) sample(0:500,size = 1000,prob = Y[[i]][x,],replace = T) ))
  return(dummy)
}
  
# Delay Node Functions

expDelayNodeLoad <- function(population,outFlow,samplesFlow){
  dummy <- population+samplesFlow-outFlow
  dummy[dummy<0] <- 0
  dummy[dummy>population] <- population
  return(rowMeans(dummy))
}

expDelayNodeLoadMatrix <- function(population,samplesFlow){
  return(lambda* sapply(0:500, function(x) expDelayNodeLoad(population,x,samplesFlow) ) ) 
}

expDelayNodeLogLoad <- function(population,outFlow,samplesFlow){
  dummy <- population+samplesFlow-outFlow
  dummy[dummy<0] <- 0
  dummy[dummy>population] <- population
  return(rowMeans(log(0.0001+lambda * dummy)))
}

expDelayNodeLogLoadMatrix <- function(population,samplesFlow){
  return(sapply(0:500, function(x) expDelayNodeLogLoad(population,x,samplesFlow) ) ) 
}


# FCFS Node Functions

expFCFSNodeLoad <- function(processors,outFlow,samplesFlow){
  dummy <- samplesFlow-outFlow
  dummy[dummy<0] <- 0
  dummy[dummy>processors] <- processors
  return(rowMeans(dummy))
}

expFCFSNodeLoadMatrix <- function(processors,samplesFlow){
  return(alphaMu/betaMu* sapply(0:500, function(x) expFCFSNodeLoad(processors,x,samplesFlow) ) ) 
}

expFCFSNodeLogLoad <- function(processors,outFlow,samplesFlow,randomGammas){
  dummy <- samplesFlow-outFlow
  dummy[dummy<0] <- 0
  dummy[dummy>processors] <- processors
  return(rowMeans(log(0.0001+randomGammas * dummy)))
}

expFCFSNodeLogLoadMatrix <- function(processors,samplesFlow,randomGammas){
  return(sapply(0:500, function(x) expFCFSNodeLogLoad(processors,x,samplesFlow,randomGammas) ) ) 
}


