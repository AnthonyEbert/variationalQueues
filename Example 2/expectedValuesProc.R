
buildSamplesFlows <- function(i,sizeSims){
  dummy <- t(sapply(1:sizeGrid, function(x) sample(0:N,size = sizeSims,prob = Y[[i]][x,],replace = T) ))
  return(dummy)
}
  
