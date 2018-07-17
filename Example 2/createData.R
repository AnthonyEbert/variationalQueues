
# Total jobs accessing the network
N=500 # 500

# Parameterise network rates (does not take queue-lengths into account)
mu <- matrix(0,11,11); colnames(mu) <- c(seq(0,5),seq(1,5)); rownames(mu) <- c(seq(0,5),seq(1,5));
mu[1,2] <-0.5; mu[2,4] <- 0.25; mu[4,6] <-0.25; mu[6,1] <- 0.5; mu[1,3] <-0.5; mu[3,5] <- 1.5; mu[5,6] <- 1.5;
mu[1,7] <-1.5; mu[7,9] <- 0.5; mu[9,11] <-0.5; mu[11,1] <- 1; mu[1,8] <-1.5; mu[8,10] <- 4; mu[10,11] <- 4;
# First 5 are High priority, next are Low

# Simulate network until jobs die
jumpTimes <- rep(0,N*4+1)
jumpStates <- matrix(0,length(jumpTimes),11); colnames(jumpStates) <- c(as.character(0:5),as.character(1:5)) 
jumpStates[1,]<-c(N,0,0,0,0,0,0,0,0,0,0)
i=2
while(N>0){
  
  # Re-scale those intensities based on current loads
  weightedMu <- mu
  weightedMu[1,] <- mu[1,] * (jumpStates[i-1,1]>0)
  
  # Nodes case by case unfortunately, many service types, could do in loop, but lots of coding involved
  
    # NODE 1 , 5-PS
    ratios <- max(1, sum(jumpStates[i-1,c(2,7)]) )
    ratios <- 5/ratios
    ratios <- min(1,ratios)
    weightedMu[2,] <- weightedMu[2,] * jumpStates[i-1,2] * ratios
    weightedMu[7,] <- weightedMu[7,] * jumpStates[i-1,7] * ratios
    # NODE 2 , 1-FCFS Priority
    ratios <- c(jumpStates[i-1,3]>0,jumpStates[i-1,3]==0 & jumpStates[i-1,8]>0)
    weightedMu[3,] <- weightedMu[3,] * ratios[1]
    weightedMu[8,] <- weightedMu[8,] * ratios[2]
    # NODE 3 , 5-PS
    ratios <- max(1, sum(jumpStates[i-1,c(4,9)]) )
    ratios <- 5/ratios
    ratios <- min(1,ratios)
    weightedMu[4,] <- weightedMu[4,] * jumpStates[i-1,4] * ratios
    weightedMu[9,] <- weightedMu[9,] * jumpStates[i-1,9] * ratios
    # NODE 4 , 1-FCFS Priority
    ratios <- c(jumpStates[i-1,5]>0,jumpStates[i-1,5]==0 & jumpStates[i-1,10]>0)
    weightedMu[5,] <- weightedMu[5,] * ratios[1]
    weightedMu[10,] <- weightedMu[10,] * ratios[2]
    # NODE 5 , Inf-PS
    weightedMu[6,] <- weightedMu[6,] * jumpStates[i-1,6] 
    weightedMu[11,] <- weightedMu[11,] * jumpStates[i-1,11]
    
  # Cumulative nodes for jump in each node
  nodeRates <- rowSums(weightedMu)
  
  # Total transition to rate
  rateJump <- sum(nodeRates)
  
  # Probabilities for jump being in each node
  jumpProbNode <- nodeRates / rateJump
  
  # Next time epoch
  delay <- rexp(1,rateJump)
  jumpTimes[i] <- jumpTimes[i-1] + delay
  
  # Add state change
  indexJump <- sample(1:11,1,prob=jumpProbNode)
  
  # Remove from node departure
  jumpStates[i,] <- jumpStates[i-1,]
  jumpStates[i,indexJump] <- jumpStates[i,indexJump] - 1
  
  # Simulate destination node
  destProbs <- weightedMu[indexJump,] / nodeRates[indexJump]
  indexJump2 <- sample(1:11,1,prob=destProbs)
  
  # Check whether this was a departure, update state vector, and update counter N
  if (!(indexJump %in% c(6,11))) {
    jumpStates[i,indexJump2] <- jumpStates[i,indexJump2] + 1
  } else {
    N <- N-1
  }
  
  # Update counter i
  i <- i+1
}

indexes <- which(jumpTimes<100)
jumpTimes <- jumpTimes[indexes]
jumpStates <- jumpStates[indexes,]

colMeans(jumpStates)
tail(jumpTimes)
plot(jumpTimes,jumpStates[,2],pch=20,cex=0.2,ylim=c(0,30)); points(jumpTimes,jumpStates[,7],pch=20,cex=0.2,col="red")
plot(jumpTimes,jumpStates[,3],pch=20,cex=0.2,ylim=c(0,30)); points(jumpTimes,jumpStates[,8],pch=20,cex=0.2,col="red")
plot(jumpTimes,jumpStates[,4],pch=20,cex=0.2,ylim=c(0,30)); points(jumpTimes,jumpStates[,9],pch=20,cex=0.2,col="red")
plot(jumpTimes,jumpStates[,5],pch=20,cex=0.2,ylim=c(0,30)); points(jumpTimes,jumpStates[,10],pch=20,cex=0.2,col="red")
plot(jumpTimes,jumpStates[,6],pch=20,cex=0.2,ylim=c(0,30)); points(jumpTimes,jumpStates[,11],pch=20,cex=0.2,col="red")


# check mean rates based on observation, to ensure good results...

# NODE 1 , 5-PS
taskHighServiced <- sum(diff(jumpStates[,2]) == -1)
taskLowServiced <- sum(diff(jumpStates[,7]) == -1)
timesInBetween <- diff(jumpTimes)
timeOccupiedTotal <- (jumpStates[,2] + jumpStates[,7])[-nrow(jumpStates)]
timeLoadHigh <- jumpStates[,2][-nrow(jumpStates)]
timeLoadLow <- jumpStates[,7][-nrow(jumpStates)]
timeLoadHigh[timeOccupiedTotal>5] <- 5/timeOccupiedTotal[timeOccupiedTotal>5] * timeLoadHigh[timeOccupiedTotal>5]
timeLoadLow[timeOccupiedTotal>5] <- 5/timeOccupiedTotal[timeOccupiedTotal>5] * timeLoadLow[timeOccupiedTotal>5]
timeLoadHigh <- sum(timeLoadHigh*timesInBetween)
timeLoadLow <- sum(timeLoadLow*timesInBetween)
taskHighServiced/timeLoadHigh #0.25
taskLowServiced/timeLoadLow #0.5
# NODE 3 , 5-PS
taskHighServiced <- sum(diff(jumpStates[,4]) == -1)
taskLowServiced <- sum(diff(jumpStates[,9]) == -1)
timesInBetween <- diff(jumpTimes)
timeOccupiedTotal <- (jumpStates[,4] + jumpStates[,9])[-nrow(jumpStates)]
timeLoadHigh <- jumpStates[,4][-nrow(jumpStates)]
timeLoadLow <- jumpStates[,9][-nrow(jumpStates)]
timeLoadHigh[timeOccupiedTotal>5] <- 5/timeOccupiedTotal[timeOccupiedTotal>5] * timeLoadHigh[timeOccupiedTotal>5]
timeLoadLow[timeOccupiedTotal>5] <- 5/timeOccupiedTotal[timeOccupiedTotal>5] * timeLoadLow[timeOccupiedTotal>5]
timeLoadHigh <- sum(timeLoadHigh*timesInBetween)
timeLoadLow <- sum(timeLoadLow*timesInBetween)
taskHighServiced/timeLoadHigh #0.25
taskLowServiced/timeLoadLow #0.5
# NODE 5 , INF-PS
taskHighServiced <- sum(diff(jumpStates[,6]) == -1)
taskLowServiced <- sum(diff(jumpStates[,11]) == -1)
timesInBetween <- diff(jumpTimes)
timeLoadHigh <- jumpStates[,6][-nrow(jumpStates)]
timeLoadLow <- jumpStates[,11][-nrow(jumpStates)]
timeLoadHigh <- sum(timeLoadHigh*timesInBetween)
timeLoadLow <- sum(timeLoadLow*timesInBetween)
taskHighServiced/timeLoadHigh #0.5
taskLowServiced/timeLoadLow #1
# NODE 2 , 1-FCFS Priority
taskHighServiced <- sum(diff(jumpStates[,3]) == -1)
taskLowServiced <- sum(diff(jumpStates[,8]) == -1)
timesInBetween <- diff(jumpTimes)
timeLoadHigh <- jumpStates[,3][-nrow(jumpStates)]>0
timeLoadLow <- jumpStates[,8][-nrow(jumpStates)] >0 & !timeLoadHigh
timeLoadHigh <- sum(timeLoadHigh*timesInBetween)
timeLoadLow <- sum(timeLoadLow*timesInBetween)
taskHighServiced/timeLoadHigh #1.5
taskLowServiced/timeLoadLow #4
# NODE 4 , 1-FCFS Priority
taskHighServiced <- sum(diff(jumpStates[,5]) == -1)
taskLowServiced <- sum(diff(jumpStates[,10]) == -1)
timesInBetween <- diff(jumpTimes)
timeLoadHigh <- jumpStates[,5][-nrow(jumpStates)]>0
timeLoadLow <- jumpStates[,10][-nrow(jumpStates)] >0 & !timeLoadHigh
timeLoadHigh <- sum(timeLoadHigh*timesInBetween)
timeLoadLow <- sum(timeLoadLow*timesInBetween)
taskHighServiced/timeLoadHigh #1.5
taskLowServiced/timeLoadLow #4


#### Pepare to save the data!
jumpTimes2 <- rep(jumpTimes,each=2)
jumpTimes2 <- jumpTimes2[-1]
jumpTimes2 <- c(jumpTimes2,100)

# Extract some partial observations
obsAmount = 50
obsTimes <- seq(100/obsAmount,100,100/obsAmount)
obsStates <- matrix(0,50,10); colnames(obsStates) <- rep(1:5,2)
for(i in 1:obsAmount) obsStates[i,] <- jumpStates[tail(which(jumpTimes-obsTimes[i]<0),1),-1]
  #No randomness in these observations, will need some low variance gaussian or indicator likelihood approximation
plot(jumpTimes2,rep(jumpStates[,3],each=2),type="l",ylim=c(0,20),xlim=c(0,100))
points(obsTimes,obsStates[,2],pch=20)
points(obsTimes,obsStates[,7],pch=20,col="red")


######################################

# Save data and make a plot with ggplot for paper
#save(jumpTimes2, jumpStates, obsTimes, file="obsTimes.Rdata")
#save(obsStates, file = "obsStates.Rdata")

