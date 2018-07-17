
# Master equation
get.Q.start<-function(t, Q, jumpRate) {
  
  # rate of change - Use vector notation
  mat2 <- diag(0,N+1); for (i in 1:N) mat2[i+1,i] <- 1; 
  dQ <- - c(rep(jumpRate,N),0) * Q + jumpRate * mat2 %*% Q
  
  # return list
  list(c(dQ))
  
}

# Functions to interpolate values of nu
get.nu <- function(t,direction){
  downIdx <- tail(which(t-gridTimes >= 0),1)
  return(nu[[direction]][downIdx,])
}

# Model equations
get.Q<-function(t, Q, parameters) {
  
  #print(t)
  # Matrix for vector notation in system of diff. equations
  mat2 <- diag(0,N+1); for (i in 1:N) mat2[i+1,i] <- 1; 
  
  # rate of change - Use vector notation
  dummyNu <- get.nu(t,direction)
  dQ <- mat2 %*% (Q * dummyNu) - Q * dummyNu

  #return list
  #print(t)
  list(c(dQ))
}