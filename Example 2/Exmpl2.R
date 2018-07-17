

# Solving ODEs
library(deSolve)

# Parameterise network rates for reference (Keep matrix form as it just makes it easy to interpret)
mu <- matrix(0,11,11); colnames(mu) <- c(seq(0,5),seq(1,5)); rownames(mu) <- c(seq(0,5),seq(1,5));
mu[1,2] <-0.5; mu[2,4] <- 0.25; mu[4,6] <-0.25; mu[6,1] <- 0.5; mu[1,3] <-0.5; mu[3,5] <- 1.5; mu[5,6] <- 1.5;
mu[1,7] <-1.5; mu[7,9] <- 0.5; mu[9,11] <-0.5; mu[11,1] <- 1; mu[1,8] <-1.5; mu[8,10] <- 4; mu[10,11] <- 4;

# Space to Store Parameters across iterations
toSaveAlphaMu <- c()
toSaveBetaMu<- c()

# Prepare
#############################################################

load("obsStates.Rdata");load("obsTimes.Rdata")
lag= 0.1; N= 400; epsilon <- 1e-1 # 0.01  # total amount of jobs transitioning network were 500
nuCap <- 50 
sizeSims <- 500

# Define grid
gridTimes <- seq(0, 102,lag)
sizeGrid <- length(gridTimes)

# WE HAVE A TOTAL OF 7x2 TRANSITION TYPES, that 14 values of r,k,nu,Y
directions <- 14
r <- vector(mode="list",length = directions); k <- vector(mode="list",length = directions); Y <- vector(mode="list",length = directions); nu <-vector(mode="list",length = directions);
for (i in 1:directions){
  r[[i]] <- matrix(0,sizeGrid,N+1); colnames(r[[i]])<-c(0:N); rownames(r[[i]])<- gridTimes
  k[[i]] <- matrix(0,sizeGrid,N+1); colnames(k[[i]])<-c(0:N); rownames(k[[i]])<- gridTimes
  nu[[i]] <- matrix(0,sizeGrid,N+1); colnames(nu[[i]])<-c(0:N); rownames(nu[[i]])<- gridTimes
  Y[[i]] <- matrix(0,sizeGrid,N+1); colnames(Y[[i]])<-c(0:N); rownames(Y[[i]])<- gridTimes
}

# Define exponential priors over the 10 rates that are not fixed, summary((mu[-1,])[mu[-1,]>0])
alphaMu <- rep(1,10) 
betaMu <- rep(0.3,10) 

# Build the initial estimates across counts in Y
source("masterEq.R")
source("expectedValuesProc.R")
source("difEqLagrangian.R")

arrivalRates <- c(0.5,0.5,1.5,1.5) # these are the entry rates
0.5-lm(obsStates[,1]~obsTimes-1)$coef
0.5-lm(obsStates[,2]~obsTimes-1)$coef
0.4419802-lm(obsStates[,3]~obsTimes-1)$coef
0.4956436-lm(obsStates[,4]~obsTimes-1)$coef
0.3905416 + 0.4890274 - lm(obsStates[,5]~obsTimes-1)$coef
1.5-lm(obsStates[,6]~obsTimes-1)$coef
1.5-lm(obsStates[,7]~obsTimes-1)$coef
1.38862-lm(obsStates[,8]~obsTimes-1)$coef
1.481316-lm(obsStates[,9]~obsTimes-1)$coef
1.26523 + 1.457286 - lm(obsStates[,5]~obsTimes-1)$coef

dummyRates <- c(arrivalRates, 0.4419802 , 0.4956436, 0.3905416, 0.4890274, 0.86106,
                1.38862 , 1.481316, 1.26523, 1.457286, 2.704007)

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

for (iterations in 21:23){
  
  if (iterations>20) epsilon <- 1e-10
  
  print("Current Iteration:")
  print(iterations)
  
  # First high priorities, so that it stays well calibrated for low priotiy expecteds computation
  
  for (direction in c(5,6,7,8,9,1,2,10,11,12,13,14,3,4)){

    print(direction)

    # Assign Expected Values to Formula in r (eq.14) - across directions
    ########################################################################
    
    if (direction == 5){ # FLow of high priority jobs from servers 1 to 3
      queueLengthOtherPrio <- buildSamplesFlows(3,sizeSims) - buildSamplesFlows(10,sizeSims);
      inputFlowThisPrio <- buildSamplesFlows(1,sizeSims);
    } 
    if (direction == 7){ # FLow of high priority jobs from servers 3 to 5
      queueLengthOtherPrio <- buildSamplesFlows(10,sizeSims) - buildSamplesFlows(12,sizeSims);
      inputFlowThisPrio <- buildSamplesFlows(5,sizeSims);
    }    
    if (direction == 9){ # FLow of high priority jobs from servers 5 to 0
      queueLengthOtherPrio <- buildSamplesFlows(12,sizeSims) + buildSamplesFlows(13,sizeSims) - buildSamplesFlows(14,sizeSims);
      inputFlowThisPrio <- buildSamplesFlows(7,sizeSims) + buildSamplesFlows(8,sizeSims);
    } 
    if (direction == 6){ # FLow of high priority jobs from servers 2 to 4
      queueLengthOtherPrio <- buildSamplesFlows(4,sizeSims) - buildSamplesFlows(11,sizeSims);
      inputFlowThisPrio <- buildSamplesFlows(2,sizeSims);
    }
    if (direction == 8){ # FLow of high priority jobs from servers 4 to 5
      queueLengthOtherPrio <- buildSamplesFlows(11,sizeSims) - buildSamplesFlows(13,sizeSims);
      inputFlowThisPrio <- buildSamplesFlows(6,sizeSims);
    }
    if (direction == 10){ # FLow of low priority jobs from servers 1 to 3
      queueLengthOtherPrio <- buildSamplesFlows(1,sizeSims) - buildSamplesFlows(5,sizeSims);
      inputFlowThisPrio <- buildSamplesFlows(3,sizeSims);
    } 
    if (direction == 12){ # FLow of low priority jobs from servers 3 to 5
      queueLengthOtherPrio <- buildSamplesFlows(5,sizeSims) - buildSamplesFlows(7,sizeSims);
      inputFlowThisPrio <- buildSamplesFlows(10,sizeSims);
    }    
    if (direction == 14){ # FLow of low priority jobs from servers 5 to 0
      queueLengthOtherPrio <- buildSamplesFlows(7,sizeSims) + buildSamplesFlows(8,sizeSims) - buildSamplesFlows(9,sizeSims);
      inputFlowThisPrio <- buildSamplesFlows(12,sizeSims) + buildSamplesFlows(13,sizeSims);
    } 
    if (direction == 11){ # FLow of low priority jobs from servers 2 to 4
      queueLengthOtherPrio <- buildSamplesFlows(2,sizeSims) - buildSamplesFlows(6,sizeSims);
      inputFlowThisPrio <- buildSamplesFlows(4,sizeSims);
      #epsilon <- 1e-250
    }
    if (direction == 13){ # FLow of low priority jobs from servers 4 to 5
      queueLengthOtherPrio <- buildSamplesFlows(6,sizeSims) - buildSamplesFlows(8,sizeSims);
      inputFlowThisPrio <- buildSamplesFlows(11,sizeSims);
    }

    # Functions for effective loads and log Loads in node
    if (direction %in% c(5,7,10,12)){
      queueLengthOtherPrio <- apply(queueLengthOtherPrio, 2, function(x) pmax(0,x))
      dummyFunc <- function(y) {# For each value of y...
        queueLengthThisPrio <- apply(inputFlowThisPrio -y,2,function(x) pmax(0,x))
        totals <- apply(queueLengthOtherPrio + queueLengthThisPrio,2,function(x) pmin(5/pmax(1,x),1))
        loadDummy <- queueLengthThisPrio * totals
        toReturn1 <- alphaMu[direction-4]/betaMu[direction-4] * rowMeans(loadDummy)
        toReturn2 <- rowMeans(log(apply(matrix(rgamma(1021*sizeSims,shape = alphaMu[direction-4],rate = betaMu[direction-4]),1021,sizeSims) * loadDummy,2,function(x) pmax(0.00000001,x))))
        return(cbind(toReturn1,toReturn2))
      }
    } else if (direction %in% c(9,14)){
      dummyFunc <- function(y) {# For each value of y...
        queueLengthThisPrio <- apply(inputFlowThisPrio -y,2,function(x) pmax(0,x))
        loadDummy <- queueLengthThisPrio
        toReturn1 <- alphaMu[direction-4]/betaMu[direction-4] * rowMeans(loadDummy)
        toReturn2 <- rowMeans(log(apply(matrix(rgamma(1021*sizeSims,shape = alphaMu[direction-4],rate = betaMu[direction-4]),1021,sizeSims) * loadDummy,2,function(x) pmax(0.00000001,x))))
        return(cbind(toReturn1,toReturn2))
      }
    } else if (direction %in% c(6,8)){
      dummyFunc <- function(y) {# For each value of y...
        queueLengthThisPrio <- apply(inputFlowThisPrio -y,2,function(x) pmax(0,x))
        loadDummy <- queueLengthThisPrio >0
        toReturn1 <- alphaMu[direction-4]/betaMu[direction-4] * rowMeans(loadDummy)
        toReturn2 <- rowMeans(log(apply(matrix(rgamma(1021*sizeSims,shape = alphaMu[direction-4],rate = betaMu[direction-4]),1021,sizeSims) * loadDummy,2,function(x) pmax(0.01,x))))
        return(cbind(toReturn1,toReturn2))
      }
    } else if (direction %in% c(11,13)){
      dummyFunc <- function(y) {# For each value of y...
        queueLengthThisPrio <- apply(inputFlowThisPrio -y,2,function(x) pmax(0,x))
        loadDummy <- queueLengthThisPrio >0 & queueLengthOtherPrio <=0
        toReturn1 <- alphaMu[direction-4]/betaMu[direction-4] * rowMeans(loadDummy)
        toReturn2 <- rowMeans(log(apply(matrix(rgamma(1021*sizeSims,shape = alphaMu[direction-4],rate = betaMu[direction-4]),1021,sizeSims) * loadDummy,2,function(x) pmax(0.01,x))))
        return(cbind(toReturn1,toReturn2))
      }
    }
 
    if(direction >4){
      expected <- sapply(0:N, function(x) dummyFunc(x) );
      expected1 <- expected[1:1021,]
      expected2 <- exp(expected[1022:2042,])
      # ts.plot(sapply(1:sizeGrid, function(x)expected1[x,]%*%Y[[direction]][x,])); lines(1:sizeGrid,sapply(1:sizeGrid, function(x)expected2[x,]%*%Y[[direction]][x,]),col="red")
    } else {
      expected1 <- matrix(arrivalRates[direction],1021,N+1)
      expected2 <- matrix(arrivalRates[direction],1021,N+1)
    }
       
  
    
    ###########################
    # After all observations... Equations are in log-form to ensure tractability
    ###########################
    
    # Start/End Parameters
    R <- rep(0,N+1) # Log form input
    state <- c(R=R) 
    
    # Solve it
    times <- sort(gridTimes[gridTimes>=tail(obsTimes,1)],decreasing = T)
    out <- ode(state, times, get.R2, parms = 0)  
    
    # Store it
    r[[direction]][gridTimes>=tail(obsTimes,1),] <- out[seq(nrow(out),1,-1),2:(N+2)]
    
    
    ##############################
    # Loop between observations...
    ##############################
    
    if(direction==5){
      serversToCheck <- c(1,3)
      fromOtherDirection <- 1
      toOtherDirection <- 7
    } else if(direction==10){
      serversToCheck <- c(6,8)
      fromOtherDirection <- 3
      toOtherDirection <- 12
    } else if(direction==6){
      serversToCheck <- c(2,4)
      fromOtherDirection <- 2
      toOtherDirection <- 8
    } else if(direction==11){
      serversToCheck <- c(7,9)
      fromOtherDirection <- 4
      toOtherDirection <- 13
    } else if(direction==7){
      serversToCheck <- c(3,5)
      fromOtherDirection <- 5
      toOtherDirection <- 9
      toExtraInput <- 8
    } else if(direction==12){
      serversToCheck <- c(8,10)
      fromOtherDirection <- 10
      toOtherDirection <- 14
      toExtraInput <- 13
    } else if(direction==8){
      serversToCheck <- c(4,5)
      fromOtherDirection <- 6
      toOtherDirection <- 9
      toExtraInput <- 7
    } else if(direction==13){
      serversToCheck <- c(9,10)
      fromOtherDirection <- 11
      toOtherDirection <- 14
      toExtraInput <- 12
    } else if(direction==9){
      serversToCheck <- c(5)
      fromOtherDirection1 <- 7
      fromOtherDirection2 <- 8
    } else if(direction==14){
      serversToCheck <- c(10)
      fromOtherDirection1 <- 12
      fromOtherDirection2 <- 13
    } else if(direction==1){
      serversToCheck <- c(1)
      toOtherDirection <- 5
    } else if(direction==2){
      serversToCheck <- c(2)
      toOtherDirection <- 6
    } else if(direction==3){
      serversToCheck <- c(6)
      toOtherDirection <- 10
    } else if(direction==4){
      serversToCheck <- c(7)
      toOtherDirection <- 11
    }
    

    if (direction %in% c(9,14)){
      otherDirectionSampler1 <- buildSamplesFlows(fromOtherDirection1,sizeSims)
      otherDirectionSampler2 <- buildSamplesFlows(fromOtherDirection2,sizeSims)
      dummyFunc <- function(y,observations){
        # Server from
        probs1 <- (observations[serversToCheck[1]] == otherDirectionSampler1[indexGrid+1,] + otherDirectionSampler2[indexGrid+1,] - y ) * (1-epsilon) +
          (observations[serversToCheck[1]] != otherDirectionSampler1[indexGrid+1,] + otherDirectionSampler2[indexGrid+1,] - y) * 
          (epsilon/N)^log(1+abs(  observations[serversToCheck[1]] - (otherDirectionSampler1[indexGrid+1,] + otherDirectionSampler2[indexGrid+1,] - y)  ))
        probs1 <- mean(log(probs1))
        # Combine
        return(exp(probs1))
      }
    } else if (direction %in% c(7,12,8,13)){
       extraInputSampler <- buildSamplesFlows(toExtraInput,sizeSims)
       otherDirectionSampler <- buildSamplesFlows(toOtherDirection,sizeSims)
       dummyFunc <- function(y,observations){
         # Server from
         probs1 <- (observations[serversToCheck[1]] == 0:N - y ) * (1-epsilon) +
           (observations[serversToCheck[1]] != 0:N - y)  * (epsilon/N)^log(1+abs(observations[serversToCheck[1]]-(0:N - y)))
         probs1 <- sum(log(probs1) * Y[[fromOtherDirection]][indexGrid+1,])
         # Server to   
         probs2 <- (observations[serversToCheck[2]] == y+extraInputSampler[indexGrid+1,]-otherDirectionSampler[indexGrid+1,]) * (1-epsilon) +
           (observations[serversToCheck[2]] != y+extraInputSampler[indexGrid+1,]-otherDirectionSampler[indexGrid+1,])  * 
           (epsilon/N)^log(1+abs( observations[serversToCheck[2]] - (y+extraInputSampler[indexGrid+1,]-otherDirectionSampler[indexGrid+1,])  ))
         probs2 <- mean(log(probs2))
         if(iterations == 1){
           probs2 <- 0
         }
         # Combine
         return(exp(probs2+probs1))
       }
     } else if (direction %in% c(1,2,3,4)){
       dummyFunc <- function(y,observations){
         # Server to        
         probs2 <- (observations[serversToCheck] == y - 0:N) * (1-epsilon) +
           (observations[serversToCheck] != y - 0:N)  * (epsilon/N)^log(1+abs(observations[serversToCheck]-(y-0:N)))
         probs2 <- sum(log(probs2) * Y[[toOtherDirection]][indexGrid+1,])
         # Combine
         return(exp(probs2))
       }
     } else {
       dummyFunc <- function(y,observations){
         # Server from
         probs1 <- (observations[serversToCheck[1]] == 0:N - y ) * (1-epsilon) +
           (observations[serversToCheck[1]] != 0:N - y)  * (epsilon/N)^log(1+abs(observations[serversToCheck[1]]-(0:N - y)))
         probs1 <- sum(log(probs1) * Y[[fromOtherDirection]][indexGrid+1,])
         # Server to        
         probs2 <- (observations[serversToCheck[2]] == y - 0:N) * (1-epsilon) +
           (observations[serversToCheck[2]] != y - 0:N)  * (epsilon/N)^log(1+abs(observations[serversToCheck[2]]-(y-0:N)))
         probs2 <- sum(log(probs2) * Y[[toOtherDirection]][indexGrid+1,])
         if(iterations == 1){
         probs2 <- 0
         }
         # Combine
         return(exp(probs2+probs1))
       }
     }
    
    for (i in seq(length(obsTimes)-1,0,-1)){
      
      # Manually set value of at jump...
      indexGrid <- tail(which(gridTimes<obsTimes[i+1]),1)
      
      # Update observations
      observations <- obsStates[i+1,]

      # Compute weights and limit
      weightsObs <- sapply(0:N, function(y) dummyFunc(y,observations))
      R <- r[[direction]][indexGrid+1,] + log(weightsObs) # Approximately
      # Increase those figures (does not alter result), sums with logarithms, products with bare numbers
      R <- R - mean(R)
      #ts.plot(R)
      
      state <- c(R=R) 
      # Solve it
      times <- sort(gridTimes[gridTimes>=max(0,obsTimes[i]) & gridTimes<=obsTimes[i+1]],decreasing = T)
      # We then need to remove the first row of the output as it is the jump
      out <- ode(state, times, get.R2, parms = 0)
     
      # Store it
      r[[direction]][gridTimes>=max(0,obsTimes[i]) & gridTimes<obsTimes[i+1],] <- 
        out[seq(nrow(out),1,-1),2:(N+2)][-nrow(out),]
      
      if (min(r[[direction]]) == -Inf ) {
        print("error")
        print(i)
      }
      
    }
    
    
    ########################
    # Update rates and slack
    ########################
    
    nu[[direction]][,1:N] <-  exp(r[[direction]][,-1] - r[[direction]][,-N]) * expected2[,-N]
    dummyIndex <- which(nu[[direction]]>nuCap,arr.ind = T)
    if(length(dummyIndex)>0){
      
      print("RATE WENT ABOVE CAP")
      print(direction)
      
      k[[direction]][,]<-0
      k[[direction]][which(nu[[direction]]>nuCap,arr.ind = T)] <- (Y[[direction]] * log(nu[[direction]]/nuCap) )[which(nu[[direction]]>nuCap,arr.ind = T)]
      nu[[direction]][which(nu[[direction]]>nuCap,arr.ind = T)] <- nuCap
      
    } else{
      k[[direction]][,]<-0
    }
    
    #summary(c(nu[[direction]])); ts.plot(nu[[direction]][,2])
    
    
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
  
  #######################
  # Update Rates 
  #######################
  
  # Expected intensities

  expIntensities <- t(sapply(5:14,function(x) rowSums(nu[[x]]*Y[[x]])))
  expLoadFCFS <- matrix(0,10,sizeGrid)
  # Direction 5
  thisDirectionSim <- buildSamplesFlows(1,sizeSims) - buildSamplesFlows(5,sizeSims)
  thisDirectionSim[thisDirectionSim<0]<-0
  thisDirectionSim2 <- buildSamplesFlows(3,sizeSims) - buildSamplesFlows(10,sizeSims)
  thisDirectionSim2[thisDirectionSim2<0]<-0
  totals<- apply( thisDirectionSim2 +  thisDirectionSim,2, function(x) pmin(5/pmax(1,x),1) )
  expLoadFCFS[1,] <- rowMeans(thisDirectionSim * totals)
  # Direction 10
  thisDirectionSim <- buildSamplesFlows(3,sizeSims) - buildSamplesFlows(10,sizeSims)
  thisDirectionSim[thisDirectionSim<0]<-0
  thisDirectionSim2 <- buildSamplesFlows(1,sizeSims) - buildSamplesFlows(5,sizeSims)
  thisDirectionSim2[thisDirectionSim2<0]<-0
  totals<- apply( thisDirectionSim2 +  thisDirectionSim,2, function(x) pmin(5/pmax(1,x),1) )
  expLoadFCFS[6,] <- rowMeans(thisDirectionSim * totals)
  # Direction 7
  thisDirectionSim <- buildSamplesFlows(5,sizeSims) - buildSamplesFlows(7,sizeSims)
  thisDirectionSim[thisDirectionSim<0]<-0
  thisDirectionSim2 <- buildSamplesFlows(10,sizeSims) - buildSamplesFlows(12,sizeSims)
  thisDirectionSim2[thisDirectionSim2<0]<-0
  totals<- apply(thisDirectionSim2 +  thisDirectionSim,2, function(x) pmin(5/pmax(1,x),1) )
  expLoadFCFS[3,] <- rowMeans(thisDirectionSim * totals)
  # Direction 12
  thisDirectionSim <- buildSamplesFlows(10,sizeSims) - buildSamplesFlows(12,sizeSims)
  thisDirectionSim[thisDirectionSim<0]<-0
  thisDirectionSim2 <- buildSamplesFlows(5,sizeSims) - buildSamplesFlows(7,sizeSims)
  thisDirectionSim2[thisDirectionSim2<0]<-0
  totals<- apply(thisDirectionSim2 +  thisDirectionSim,2, function(x) pmin(5/pmax(1,x),1) )
  expLoadFCFS[8,] <- rowMeans(thisDirectionSim * totals)
  # Direction 6
  thisDirectionSim <- buildSamplesFlows(2,sizeSims) - buildSamplesFlows(6,sizeSims)
  expLoadFCFS[2,] <- rowMeans(thisDirectionSim > 0)     
  # Direction 8
  thisDirectionSim <- buildSamplesFlows(6,sizeSims) - buildSamplesFlows(8,sizeSims)
  expLoadFCFS[4,] <- rowMeans(thisDirectionSim > 0)  
  # Direction 11
  thisDirectionSim <- buildSamplesFlows(4,sizeSims) - buildSamplesFlows(11,sizeSims)
  expLoadFCFS[7,] <- rowMeans(thisDirectionSim > 0 & (buildSamplesFlows(2,sizeSims) - buildSamplesFlows(6,sizeSims) <= 0))     
  # Direction 13
  thisDirectionSim <- buildSamplesFlows(11,sizeSims) - buildSamplesFlows(13,sizeSims)
  expLoadFCFS[9,] <- rowMeans(thisDirectionSim > 0 & (buildSamplesFlows(6,sizeSims) - buildSamplesFlows(8,sizeSims) <= 0))  
  # Direction 9
  thisDirectionSim <- buildSamplesFlows(7,sizeSims) + buildSamplesFlows(8,sizeSims) - buildSamplesFlows(9,sizeSims)
  expLoadFCFS[5,] <- rowMeans(thisDirectionSim)
  # Direction 14
  thisDirectionSim <- buildSamplesFlows(12,sizeSims) + buildSamplesFlows(13,sizeSims) - buildSamplesFlows(14,sizeSims)
  expLoadFCFS[10,] <- rowMeans(thisDirectionSim)
 
  # NEW RATES
  alphaMu <- 1 + rowSums( expIntensities[,-sizeGrid] * lag )
  betaMu <- 0.5 + rowSums( expLoadFCFS[,-sizeGrid] * lag ) 
  #plot(seq(0,5,0.01),dgamma(seq(0,5,0.01),shape=alphaMu,rate=betaMu),type="l")

  ##################
  # Save values
  ##################
  
  toSaveAlphaMu <- rbind(toSaveAlphaMu, alphaMu)
  toSaveBetaMu <- rbind(toSaveBetaMu, betaMu)
  print(toSaveAlphaMu/toSaveBetaMu)
}

# Store all that data in a results file
#######################################
#save(Y,nu,r,k,alphaMu,betaMu,toSaveAlphaMu,toSaveBetaMu,file="resultsFinal2.Rdata")

(toSaveAlphaMu/toSaveBetaMu)[23,]
sqrt((toSaveAlphaMu/toSaveBetaMu^2))[23,]
qgamma(0.025,shape = toSaveAlphaMu[23,],rate=toSaveBetaMu[23,])
qgamma(0.25,shape = toSaveAlphaMu[23,],rate=toSaveBetaMu[23,])
qgamma(0.5,shape = toSaveAlphaMu[23,],rate=toSaveBetaMu[23,])
qgamma(0.75,shape = toSaveAlphaMu[23,],rate=toSaveBetaMu[23,])
qgamma(0.975,shape = toSaveAlphaMu[23,],rate=toSaveBetaMu[23,])


# TRANSITION PLOTS
###########################

library(ggplot2)

toPlotData <- data.frame(t = obsTimes,x = obsStates[,6])
meanPath <- data.frame(gridTimes, mean = Y[[3]] %*% 0:400 - Y[[10]] %*% 0:400)
iepe <- data.frame(matrix(0,1021,2)); colnames(iepe) <- c('low','up')
for(i in 1: 1021){
  ue <- sample(x = 0:400,size = 20000,replace = T,prob = Y[[3]][i,]) - sample(x = 0:400,size = 20000,replace = T,prob = Y[[10]][i,])
  iepe[i,] <- c(sort(ue)[500],sort(ue)[19500])
}
meanPath <- cbind(meanPath,iepe)
meanPath$low[meanPath$low<0] <- 0

p <- ggplot(meanPath, aes(x=gridTimes, y=mean))
p1<- p +  geom_ribbon(aes(ymin = low, ymax = up, x=gridTimes),fill="gray95",colour="gray70", size=0.05) + geom_point(data=toPlotData,aes(t,x),size=1.2,shape=18) 

toPlotData <- data.frame(t = obsTimes,x = obsStates[,1])
meanPath <- data.frame(gridTimes, mean = Y[[1]] %*% 0:400 - Y[[5]] %*% 0:400)
iepe <- data.frame(matrix(0,1021,2)); colnames(iepe) <- c('low','up')
for(i in 1: 1021){
  ue <- sample(x = 0:400,size = 20000,replace = T,prob = Y[[1]][i,]) - sample(x = 0:400,size = 20000,replace = T,prob = Y[[5]][i,])
  iepe[i,] <- c(sort(ue)[500],sort(ue)[19500])
}
meanPath <- cbind(meanPath,iepe)
meanPath$low[meanPath$low<0] <- 0

p1<- p1 +  geom_ribbon(data=meanPath,aes(ymin = low, ymax = up, x=gridTimes),fill="gray80",colour="gray70", size=0.05) + 
  geom_point(data=toPlotData,aes(t,x),size=0.6) +
  scale_x_continuous(name=element_blank(), limits = c(0, 102),breaks=NULL, expand = c(0.01,0)) + 
  scale_y_continuous(name="Service station 1") + 
  theme_minimal() +
  theme(axis.title.x=element_text(size=12),legend.position="none",axis.title.y=element_text(size=13,margin=margin(l=0,r=10)))

######################

toPlotData <- data.frame(t = obsTimes,x = obsStates[,7])
meanPath <- data.frame(gridTimes, mean = Y[[4]] %*% 0:400 - Y[[11]] %*% 0:400)
iepe <- data.frame(matrix(0,1021,2)); colnames(iepe) <- c('low','up')
for(i in 1: 1021){
  ue <- sample(x = 0:400,size = 20000,replace = T,prob = Y[[4]][i,]) - sample(x = 0:400,size = 20000,replace = T,prob = Y[[11]][i,])
  iepe[i,] <- c(sort(ue)[500],sort(ue)[19500])
}
meanPath <- cbind(meanPath,iepe)
meanPath$low[meanPath$low<0] <- 0

p <- ggplot(meanPath, aes(x=gridTimes, y=mean))
p2<- p +  geom_ribbon(aes(ymin = low, ymax = up, x=gridTimes),fill="gray95",colour="gray70", size=0.05) + geom_point(data=toPlotData,aes(t,x),size=1.2,shape=18) 

toPlotData <- data.frame(t = obsTimes,x = obsStates[,2])
meanPath <- data.frame(gridTimes, mean = Y[[2]] %*% 0:400 - Y[[6]] %*% 0:400)
iepe <- data.frame(matrix(0,1021,2)); colnames(iepe) <- c('low','up')
for(i in 1: 1021){
  ue <- sample(x = 0:400,size = 20000,replace = T,prob = Y[[2]][i,]) - sample(x = 0:400,size = 20000,replace = T,prob = Y[[6]][i,])
  iepe[i,] <- c(sort(ue)[500],sort(ue)[19500])
}
meanPath <- cbind(meanPath,iepe)
meanPath$low[meanPath$low<0] <- 0

p2<- p2 +  geom_ribbon(data=meanPath,aes(ymin = low, ymax = up, x=gridTimes),fill="gray80",colour="gray70", size=0.05) + 
  geom_point(data=toPlotData,aes(t,x),size=0.6) +
  scale_x_continuous(name=element_blank(), limits = c(0, 102),breaks=NULL, expand = c(0.01,0)) + 
  scale_y_continuous(name="Service station 2") + 
  theme_minimal() +
  theme(axis.title.x=element_text(size=12),legend.position="none",axis.title.y=element_text(size=13,margin=margin(l=0,r=10)))

######################

toPlotData <- data.frame(t = obsTimes,x = obsStates[,8])
meanPath <- data.frame(gridTimes, mean = Y[[10]] %*% 0:400 - Y[[12]] %*% 0:400)
iepe <- data.frame(matrix(0,1021,2)); colnames(iepe) <- c('low','up')
for(i in 1: 1021){
  ue <- sample(x = 0:400,size = 20000,replace = T,prob = Y[[10]][i,]) - sample(x = 0:400,size = 20000,replace = T,prob = Y[[12]][i,])
  iepe[i,] <- c(sort(ue)[500],sort(ue)[19500])
}
meanPath <- cbind(meanPath,iepe)
meanPath$low[meanPath$low<0] <- 0

p <- ggplot(meanPath, aes(x=gridTimes, y=mean))
p3<- p +  geom_ribbon(aes(ymin = low, ymax = up, x=gridTimes),fill="gray95",colour="gray70", size=0.05) + geom_point(data=toPlotData,aes(t,x),size=1.2,shape=18) 

toPlotData <- data.frame(t = obsTimes,x = obsStates[,3])
meanPath <- data.frame(gridTimes, mean = Y[[5]] %*% 0:400 - Y[[7]] %*% 0:400)
iepe <- data.frame(matrix(0,1021,2)); colnames(iepe) <- c('low','up')
for(i in 1: 1021){
  ue <- sample(x = 0:400,size = 20000,replace = T,prob = Y[[5]][i,]) - sample(x = 0:400,size = 20000,replace = T,prob = Y[[7]][i,])
  iepe[i,] <- c(sort(ue)[500],sort(ue)[19500])
}
meanPath <- cbind(meanPath,iepe)
meanPath$low[meanPath$low<0] <- 0

p3<- p3 +  geom_ribbon(data=meanPath,aes(ymin = low, ymax = up, x=gridTimes),fill="gray80",colour="gray70", size=0.05) + 
  geom_point(data=toPlotData,aes(t,x),size=0.6) +
  scale_x_continuous(name=element_blank(), limits = c(0, 102),breaks=NULL, expand = c(0.01,0)) + 
  scale_y_continuous(name="Service station 3") + 
  theme_minimal() +
  theme(axis.title.x=element_text(size=12),legend.position="none",axis.title.y=element_text(size=13,margin=margin(l=0,r=10)))


######################

toPlotData <- data.frame(t = obsTimes,x = obsStates[,9])
meanPath <- data.frame(gridTimes, mean = Y[[11]] %*% 0:400 - Y[[13]] %*% 0:400)
iepe <- data.frame(matrix(0,1021,2)); colnames(iepe) <- c('low','up')
for(i in 1: 1021){
  ue <- sample(x = 0:400,size = 20000,replace = T,prob = Y[[11]][i,]) - sample(x = 0:400,size = 20000,replace = T,prob = Y[[13]][i,])
  iepe[i,] <- c(sort(ue)[500],sort(ue)[19500])
}
meanPath <- cbind(meanPath,iepe)
meanPath$low[meanPath$low<0] <- 0

p <- ggplot(meanPath, aes(x=gridTimes, y=mean))
p4<- p +  geom_ribbon(aes(ymin = low, ymax = up, x=gridTimes),fill="gray95",colour="gray70", size=0.05) + geom_point(data=toPlotData,aes(t,x),size=1.2,shape=18) 

toPlotData <- data.frame(t = obsTimes,x = obsStates[,4])
meanPath <- data.frame(gridTimes, mean = Y[[6]] %*% 0:400 - Y[[8]] %*% 0:400)
iepe <- data.frame(matrix(0,1021,2)); colnames(iepe) <- c('low','up')
for(i in 1: 1021){
  ue <- sample(x = 0:400,size = 20000,replace = T,prob = Y[[6]][i,]) - sample(x = 0:400,size = 20000,replace = T,prob = Y[[8]][i,])
  iepe[i,] <- c(sort(ue)[500],sort(ue)[19500])
}
meanPath <- cbind(meanPath,iepe)
meanPath$low[meanPath$low<0] <- 0

p4<- p4 +  geom_ribbon(data=meanPath,aes(ymin = low, ymax = up, x=gridTimes),fill="gray80",colour="gray70", size=0.05) + 
  geom_point(data=toPlotData,aes(t,x),size=0.6) +
  scale_x_continuous(name=element_blank(), limits = c(0, 102),breaks=NULL, expand = c(0.01,0)) + 
  scale_y_continuous(name="Service station 4") + 
  theme_minimal() +
  theme(axis.title.x=element_text(size=12),legend.position="none",axis.title.y=element_text(size=13,margin=margin(l=0,r=10)))


######################

toPlotData <- data.frame(t = obsTimes,x = obsStates[,10])
meanPath <- data.frame(gridTimes, mean = Y[[12]] %*% 0:400 + Y[[13]] %*% 0:400 - Y[[14]] %*% 0:400)
iepe <- data.frame(matrix(0,1021,2)); colnames(iepe) <- c('low','up')
for(i in 1: 1021){
  ue <- sample(x = 0:400,size = 20000,replace = T,prob = Y[[12]][i,]) + sample(x = 0:400,size = 20000,replace = T,prob = Y[[13]][i,]) -
    sample(x = 0:400,size = 20000,replace = T,prob = Y[[14]][i,])
  iepe[i,] <- c(sort(ue)[500],sort(ue)[19500])
}
meanPath <- cbind(meanPath,iepe)
meanPath$low[meanPath$low<0] <- 0

p <- ggplot(meanPath, aes(x=gridTimes, y=mean))
p5<- p +  geom_ribbon(aes(ymin = low, ymax = up, x=gridTimes),fill="gray95",colour="gray70", size=0.05) + geom_point(data=toPlotData,aes(t,x),size=1.2,shape=18) 

toPlotData <- data.frame(t = obsTimes,x = obsStates[,5])
meanPath <- data.frame(gridTimes, mean = Y[[7]] %*% 0:400 + Y[[8]] %*% 0:400 - Y[[9]] %*% 0:400)
iepe <- data.frame(matrix(0,1021,2)); colnames(iepe) <- c('low','up')
for(i in 1: 1021){
  ue <- sample(x = 0:400,size = 20000,replace = T,prob = Y[[7]][i,]) + sample(x = 0:400,size = 20000,replace = T,prob = Y[[8]][i,]) -
    sample(x = 0:400,size = 20000,replace = T,prob = Y[[9]][i,])
  iepe[i,] <- c(sort(ue)[500],sort(ue)[19500])
}
meanPath <- cbind(meanPath,iepe)
meanPath$low[meanPath$low<0] <- 0

p5<- p5 +  geom_ribbon(data=meanPath,aes(ymin = low, ymax = up, x=gridTimes),fill="gray80",colour="gray70", size=0.05) + 
  geom_point(data=toPlotData,aes(t,x),size=0.6) +
  scale_x_continuous(name=element_blank(), limits = c(0, 102),breaks=NULL, expand = c(0.01,0)) + 
  scale_y_continuous(name="Service station 5") + 
  theme_minimal() +
  theme(axis.title.x=element_text(size=12),legend.position="none",axis.title.y=element_text(size=13,margin=margin(l=0,r=10)))



dat1 <- data.frame(t=rep(gridTimes,2), x= c(expIntensities[1,],expLoadFCFS[1,]) , type=rep(c("Expected Intensity", "Expected Load"), each=1021))
p6 <- ggplot(dat1) + geom_line(aes(y=x, x=t),colour="gray40") +  facet_grid(type ~ .) +
  scale_x_continuous(name=element_blank(), limits = c(0, 102),breaks=NULL, expand = c(0.01,0)) + 
  scale_y_continuous(name="Jump direction (0,1,1)") + 
  theme_minimal() +
  theme(axis.title.x=element_text(size=12),legend.position="none",axis.title.y=element_text(size=13,margin=margin(l=0,r=10)))


library(grid)
library(gridExtra)
grid.arrange(p1,p2,p3,p4,p5,p6, ncol=2, nrow=3, top=textGrob("Trajectory estimate across stations and job priorities", gp=gpar(fontsize=15)))
#1000x900
















