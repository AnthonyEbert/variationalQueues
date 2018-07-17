

# Solving ODEs
library(deSolve)

# Parameterise closed queue
mu <- 3; lambda <- 0.1

# Space to Store Parameters across iterations
toSaveAlphaMu <- c()
toSaveBetaMu <- c()
toSaveLowerBound <- c()


# Prepare, We are after mu here, this is old code
#############################################################

load("obsStates.Rdata");load("obsTimes.Rdata")
population=50; lag= 0.1; N= 500; epsilon <- 0.2 # Large N is because jumps in each direction could be a lot!
nuCap <- 100

# Define grid
gridTimes <- seq(0, 102,lag)
sizeGrid <- length(gridTimes)

# Save space for r,k and Y, most are not reusable, so do all unique just in case
directions <- 2
r <- vector(mode="list",length = directions); k <- vector(mode="list",length = directions); Y <- vector(mode="list",length = directions); nu <-vector(mode="list",length = directions);
for (i in 1:directions){
  r[[i]] <- matrix(0,sizeGrid,N+1); colnames(r[[i]])<-c(0:N); rownames(r[[i]])<- gridTimes
  k[[i]] <- matrix(0,sizeGrid,N+1); colnames(k[[i]])<-c(0:N); rownames(k[[i]])<- gridTimes
  nu[[i]] <- matrix(0,sizeGrid,N+1); colnames(nu[[i]])<-c(0:N); rownames(nu[[i]])<- gridTimes
  Y[[i]] <- matrix(0,sizeGrid,N+1); colnames(Y[[i]])<-c(0:N); rownames(Y[[i]])<- gridTimes
}

# Define priors, lambda (stands for arrivals in this code), is "fixed" 
alphaMu <- 5 
betaMu <- 2 

# Build the initial estimates across counts in Y
source("masterEq.R")
source("expectedValuesProc.R")
source("difEqLagrangian.R")
source("obsModel.R")
dummyRates <- c(mean(obsStates2) * lambda, alphaMu/betaMu) #DelayToFCFS, FCFStoDelay

for (i in 1:directions){
  
  # Parameters and grid
  Q <- c(1,rep(0,N))
  state <- c(Q=Q) 
  
  # Solve it and store it
  out <- ode(state, gridTimes, get.Q.start, parms = dummyRates[i])
  Y[[i]][,1:(N+1)] <- out[,2:(N+2)]
  Y[[i]][which(Y[[i]]<0)] <- 0
  
}

rm(dummyRates)


# Variational iterations start
###################################################

toSaveAlphaMu <- c(toSaveAlphaMu, alphaMu)
toSaveBetaMu <- c(toSaveBetaMu, betaMu)

for (iterations in 1:20){
  
  print("Current Iteration:")
  print(iterations)
  
  for (direction in 1:directions){
    
    # Assign Expected Values to Formula in r (eq.14) - across directions
    ########################################################################
    
    if (direction == 1){
      
      # Expected Delay Load
      sampleOtherTaskDirections<- buildSamplesFlows(2)
      sampleOtherTaskDirections <- sampleOtherTaskDirections - c(rowMeans(sampleOtherTaskDirections) - c(Y[[2]]%*%(0:500)))  # Smoothen sample
      expected1<- expDelayNodeLoadMatrix(population,sampleOtherTaskDirections)
      colnames(expected1) <- as.character(0:N); rownames(expected1) <- as.character(gridTimes)
      #ts.plot(sapply(1:sizeGrid, function(x)expected1[x,]%*%Y[[1]][x,]))
      
      # Expected Log Delay Load
      expected2<- expDelayNodeLogLoadMatrix(population,sampleOtherTaskDirections)
      colnames(expected2) <- as.character(0:N); rownames(expected2) <- as.character(gridTimes)
      expected2 <- exp(expected2)
      
    } else if (direction == 2){
      
      # Expected FCFS Load
      sampleOtherTaskDirections<- buildSamplesFlows(1)
      sampleOtherTaskDirections <- sampleOtherTaskDirections - c(rowMeans(sampleOtherTaskDirections) - c(Y[[1]]%*%(0:500)))  # Smoothen sample
      expected1<- expFCFSNodeLoadMatrix(1,sampleOtherTaskDirections)
      colnames(expected1) <- as.character(0:N); rownames(expected1) <- as.character(gridTimes)
      #ts.plot(sapply(1:sizeGrid, function(x)expected1[x,]%*%Y[[2]][x,]))
      
      # Expected Log FCFS Load
      randomGammas <- matrix(rgamma(1021000,shape = alphaMu, rate = betaMu),1021,1000)
      expected2<- expFCFSNodeLogLoadMatrix(1,sampleOtherTaskDirections,randomGammas)
      colnames(expected2) <- as.character(0:N); rownames(expected2) <- as.character(gridTimes)
      expected2 <- exp(expected2)
      
    } 
    
    
    ###########################
    # After all observations...
    ###########################
    
    # Start/End Parameters
    R <- rep(1,N+1) 
    state <- c(R=R) 
    
    # Solve it
    times <- sort(gridTimes[gridTimes>=tail(obsTimes,1)],decreasing = T)
    #out <- ode(state, times, get.R, parms = 0)  
    out <- ode(log(state), times, get.R2, parms = 0)  
    out[,2:502] <- exp(out[,2:502])
    
    # Store it
    r[[direction]][gridTimes>=tail(obsTimes,1),] <- out[seq(nrow(out),1,-1),2:(N+2)]
    
    
    ##############################
    # Loop between observations...
    ##############################
    
    for (i in seq(length(obsTimes)-1,0,-1)){
      
      # Manually set value of at jump...
      indexGrid <- tail(which(gridTimes<obsTimes[i+1]),1)
      
      # ObsStates is for FCFS, ObsStates2 is for Delay
      observations <- c(obsStates2[i+1],obsStates[i+1])
      
      if (direction == 1){
        
        weightsObs <- obsFirstDirectionVector(observations)
        dummy <- k[[direction]][indexGrid+1,]/pmax(0.0001,Y[[direction]][indexGrid+1,])
        dummy <- ( 1+ dummy ) / exp( dummy )
        R <- ( r[[direction]][indexGrid+1,] - pmax(0,( r[[direction]][indexGrid+1,] * expected1[indexGrid+1,] + (mat1 %*% r[[direction]][indexGrid+1,]) * expected2[indexGrid+1,] * dummy )) * lag ) * 
          weightsObs
        
        # Increase those figures (does not alter result)
        R <- R / mean(R)
        
      } else if (direction == 2){
        
        weightsObs <- obsSecondDirectionVector(observations)
        dummy <- k[[direction]][indexGrid+1,]/pmax(0.0001,Y[[direction]][indexGrid+1,])
        dummy <- ( 1+ dummy ) / exp( dummy )
        R <- ( r[[direction]][indexGrid+1,] - pmax(0,( r[[direction]][indexGrid+1,] * expected1[indexGrid+1,] + (mat1 %*% r[[direction]][indexGrid+1,]) * expected2[indexGrid+1,] * dummy )) * lag ) * 
          weightsObs
        
        # Increase those figures (does not alter result)
        R <- R / mean(R)
        
      }
      
      state <- c(R=R) 
      
      # Solve it
      times <- sort(gridTimes[gridTimes>=max(0,obsTimes[i]) & gridTimes<=obsTimes[i+1]],decreasing = T)
      # We then need to remove the first row of the output as it is the jump
      #out <- ode(state, times, get.R, parms = 0)
      out <- ode(log(state), times, get.R2, parms = 0)  
      out[,2:502] <- exp(out[,2:502])
      
      # Store it
      r[[direction]][gridTimes>=max(0,obsTimes[i]) & gridTimes<obsTimes[i+1],] <- 
        out[seq(nrow(out),1,-1),2:(N+2)][-nrow(out),]
      
      if (min(r[[1]]) <0 ) {
        print("error")
        print(i)
      }
      
      # Counter
      #print(i)
      
    }
    
    
    ########################
    # Update rates and slack
    ########################
    
    nu[[direction]][,1:N] <-  r[[direction]][,-1] / r[[direction]][,-N] * expected2[,-500]
    dummyIndex <- which(nu[[1]]>nuCap,arr.ind = T)
    if(length(dummyIndex)>0){
      
      print("RATE WENT ABOVE CAP!")
      
      k[[direction]][which(nu[[direction]]>nuCap,arr.ind = T)] <- (Y[[direction]] * log(nu[[direction]]/nuCap) )[which(nu[[direction]]>nuCap,arr.ind = T)]
      nu[[direction]][which(nu[[direction]]>nuCap,arr.ind = T)] <- nuCap
      
    }
    
    #summary(c(nu[[direction]]))
    
    
    #######################
    # Update Y in Direction 
    #######################
    
    # Parameters and grid
    Q <- c(1,rep(0,N))
    state <- c(Q=Q) 
    
    # Solve it
    times <- gridTimes
    out <- ode(state, times, get.Q, parms = 0)
    #summary(rowSums(out[,-1]))
    
    # Store it
    Y[[direction]][,] <- out[,2:(N+2)]; 
    Y[[direction]][which(Y[[direction]]<0,arr.ind = T)] <- 0
    
  }
  
  # Plot mean-average dynamics
  #plot(gridTimes, Y[[1]] %*% 0:N,type="l")
  #lines(gridTimes, Y[[2]] %*% 0:N,col="red")
  plot(gridTimes, Y[[1]] %*% 0:N - Y[[2]] %*% 0:N,type="l",ylim=c(0,50))
  points(obsTimes,obsStates,pch=20)
  points(obsTimes+0.35,population-obsStates2,pch=20,col="red")
  iepe <- matrix(0,sizeGrid,3) # Median and 95% conf interval
  for(i in 1: sizeGrid){
    ue <- sample(x = 0:500,size = 10000,replace = T,prob = Y[[1]][i,]) - sample(x = 0:500,size = 10000,replace = T,prob = Y[[2]][i,])
    iepe[i,] <- c(sort(ue)[250],sort(ue)[5000],sort(ue)[9750])
  }
  lines(gridTimes, iepe[,2],col="blue"); lines(gridTimes, iepe[,1],col="blue"); lines(gridTimes, iepe[,3],col="blue")
  
  
  
  #######################
  # Update Rates 
  #######################
  
  expIntensity <- rowSums(nu[[2]]*Y[[2]])
  expLoadFCFS <- rep(0,sizeGrid)
  for(i in 1: sizeGrid){
    ue <- sample(x = 0:500,size = 10000,replace = T,prob = Y[[1]][i,]) - sample(x = 0:500,size = 10000,replace = T,prob = Y[[2]][i,])
    ue[ue>1] <- 1
    ue[ue<0] <- 0
    expLoadFCFS[i] <- mean(ue)
  }
  
  alphaMu <- 5 + sum( expIntensity[-sizeGrid] * lag )
  betaMu <- 2 + sum( expLoadFCFS[-sizeGrid] * lag ) 
  #plot(seq(0,5,0.01),dgamma(seq(0,5,0.01),shape=alphaMu,rate=betaMu),type="l")
  
  
  #######################
  # Lower Bound to LogLik
  ####################### 
  
  # Build samples...
  samplesMu <- rgamma(1000,shape=alphaMu,rate=betaMu)
  samplesY01 <- t(sapply(1:sizeGrid, function(x) sample(0:500,size = 1000,prob = Y[[1]][x,],replace = T) ))
  samplesY10 <- t(sapply(1:sizeGrid, function(x) sample(0:500,size = 1000,prob = Y[[2]][x,],replace = T) ))
  
  # Observations
  timesObsIndex <- which(gridTimes %in% obsTimes)
  bound <- 0
  for (i in 1:50){
    observations <- c(obsStates2[i],obsStates[i])
    statesSimulated <- samplesY01[timesObsIndex[i],] - samplesY10[timesObsIndex[i],]
    probsDelay <- (observations[1] == population - statesSimulated) * (1-epsilon) + ( (observations[1] != population - statesSimulated) & (population - statesSimulated >= 0) & (0 <= statesSimulated) ) * epsilon/N
    probsFCFS <- (observations[2] == statesSimulated) * (1-epsilon) + ( (observations[2] != statesSimulated) & (statesSimulated >= 0) & (statesSimulated <= population) ) * epsilon/N
    prodProbs <- probsDelay * probsFCFS
    bound <- bound + mean(log(0.0001 + prodProbs))
  }  
    
  # Rate densities
  bound <- bound - mean(log(dgamma(samplesMu,shape=alphaMu,rate=betaMu)/dgamma(samplesMu,shape=5,rate=2)))
  
  # Path
  bound <- bound + sum( lag * sapply(1:(sizeGrid-1), function(tIndex)
  mean(
  nu[[1]][tIndex,samplesY01[tIndex,]+1] + nu[[2]][tIndex,samplesY10[tIndex,]+1] - 
    samplesMu * (samplesY01[tIndex,] - samplesY10[tIndex,] > 0) -
    lambda * (population - samplesY01[tIndex,] + samplesY10[tIndex,]) -
    nu[[1]][tIndex,samplesY01[tIndex,]+1] * log( nu[[1]][tIndex,samplesY01[tIndex,]+1] / (0.001+ lambda * pmax(0,population - samplesY01[tIndex,] + samplesY10[tIndex,])) ) -
    nu[[2]][tIndex,samplesY10[tIndex,]+1] * log( nu[[2]][tIndex,samplesY10[tIndex,]+1] / (0.001+ samplesMu * (samplesY01[tIndex,] - samplesY10[tIndex,] > 0)) )  )
  )  )
 
  
  ##################
  # Save values
  ##################
  
  toSaveAlphaMu <- c(toSaveAlphaMu, alphaMu)
  toSaveBetaMu <- c(toSaveBetaMu, betaMu)
  toSaveLowerBound <- c(toSaveLowerBound, bound)
  
}


# Store all that data in a results file
save(toSaveAlphaMu, toSaveBetaMu, toSaveLowerBound, gridTimes, Y, file = "results.Rdata")



