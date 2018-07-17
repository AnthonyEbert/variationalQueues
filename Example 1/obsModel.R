
obsFirstDirection <- function(y,observations){

  probsDelay <- (observations[1] == population - y + 0:N) * (1-epsilon) + ( (observations[1] != population - y + 0:N) & (population - y + 0:N >= 0) & (0:N <= y) ) * epsilon/N
  probsFCFS <- (observations[2] == y - 0:N) * (1-epsilon) + ( (observations[2] != y - 0:N) & (y >= 0:N) & (y - 0:N <= population) ) * epsilon/N
  prodProbs <- probsDelay * probsFCFS
  return(sum( log(0.0001 + prodProbs) * Y[[2]][indexGrid+1,] ))
  
}

obsFirstDirectionVector <- function(observations){
  
  return( exp(sapply(0:N, function(y) obsFirstDirection(y,observations))) )
  
}


###############################################################

obsSecondDirection <- function(y,observations){
  
  probsDelay <- (observations[1] == population - 0:N + y) * (1-epsilon) + ( (observations[1] != population - 0:N + y) & (population - 0:N + y >= 0) & (0:N >= y) ) * epsilon/N
  probsFCFS <- (observations[2] == 0:N - y) * (1-epsilon) + ( (observations[2] != 0:N - y) & (y <= 0:N) & (0:N - y <= population) ) * epsilon/N
  prodProbs <- probsDelay * probsFCFS
  return(sum( log(0.0001 + prodProbs) * Y[[1]][indexGrid+1,] ))
  
}

obsSecondDirectionVector <- function(observations){
  
  return( exp(sapply(0:N, function(y) obsSecondDirection(y,observations))) )
  
}