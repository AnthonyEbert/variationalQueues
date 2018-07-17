
# Specification dynamics lagrangian

# Stuff that feeds into differential equation
mat1 <- diag(0,N+1); for (i in 1:N) mat1[i,i+1] <- -1; 

# Model equations
get.R<-function(t, R, parameters) {
  
  timeIndex <- max(tail(which(gridTimes<=t),1),1)
  
  # rate of change - Use vector notation
  dummy <- k[[direction]][timeIndex,]/pmax(0.0001,Y[[direction]][timeIndex,])
  dummy <- ( 1+ dummy ) / exp( dummy )
  dR <- pmax(0,R * expected1[timeIndex,] + (mat1 %*% R) * expected2[timeIndex,] * dummy)
  
  # return list
  list(c(dR))
  
}


# Equation in log form
get.R2<-function(t, R, parameters) {
  
  timeIndex <- max(tail(which(gridTimes<=t),1),1)
  
  # rate of change - Use vector notation
  dummy <- k[[direction]][timeIndex,]/pmax(0.0001,Y[[direction]][timeIndex,])
  dummy <- ( 1+ dummy ) / exp( dummy )
  dR <- expected1[timeIndex,] - exp( - R - mat1 %*% R ) * expected2[timeIndex,] * dummy
  
  # return list
  list(c(dR))
  
}
