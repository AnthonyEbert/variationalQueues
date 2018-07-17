
load("results.Rdata")
load("obsStates.Rdata");load("obsTimes.Rdata")
N <- 50; epsilon <- 0.2
library(ggplot2)

# Get MCMC output to compare first:
lambda <- 0.1
mu <- 2.5
sampleStationary <-  c(obsStates[-(1:20)], N-obsStates2[-(1:20)])
sampleStationary <- sampleStationary[sample(1:length(sampleStationary),25)]
totalObs <- length(sampleStationary)
obsProb <- function(obs,state){
  return((obs==state)*(1-epsilon) + (obs!=state)*(epsilon/N))
}
obsProbMu <- function(obs,mu){
  return( sum(obsProb(obs,0:50) * (lambda/mu)^(0:N) / factorial(N-0:N) ) )
}
logLikelihood <- function(mu){
  base <- - totalObs * log( sum( (lambda/mu)^(0:N) / factorial(N-0:N) ) )
  top1 <-  sum(log( sapply(sampleStationary,function(o) obsProbMu(o,mu))  ))
  return(base+top1+log(dgamma(mu,shape = 5, rate = 2)) )
}

sampleMCMC <- rep(2.5,20000)
for (i in 2:20000){
  proposal <- rnorm(1,mean=sampleMCMC[i-1],sd = 0.25)
  ratio <- exp(logLikelihood(proposal) - logLikelihood(sampleMCMC[i-1]))
  a <- runif(1)
  if (a <= ratio){
    sampleMCMC[i] <- proposal
  } else {
    sampleMCMC[i] <- sampleMCMC[i-1]
  }
}
toDraw <- density(sampleMCMC[1000:20000])
toDraw <- data.frame(x=toDraw$x,y=toDraw$y)

toPlotData <- data.frame(mu = seq(0,6,0.001),densMu = dgamma(seq(0,6,0.001),shape = 5,rate = 2))
toPlotData2 <- data.frame(mu = seq(0,6,0.001),densMu = dgamma(seq(0,6,0.001),shape = 281.75974,rate = 101.494610065))

p <- ggplot(toPlotData, aes(x=mu, y=densMu))
ue2 <- p + geom_ribbon(aes(ymin=rep(0,6001), ymax=densMu),size=0.1,colour="black",fill=I(gray.colors(1,start = 0.9))) +
  geom_ribbon(data=toPlotData2,aes(ymin=rep(0,6001), ymax=densMu),size=0.1,colour="black",fill=I(gray.colors(1,start = 0.75))) +
  geom_point(data = data.frame(t=3,x=0),aes(t,x), size=1.5) +
  scale_x_continuous(name="Service rate", limits = c(0.5, 5), expand = c(0,0)) + scale_y_continuous(name="Density") + 
  theme_minimal() + #ggtitle("Service rate posterior density") +
  theme(axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=13,margin=margin(l=0,r=10)),plot.margin = unit(c(0.4,0.4,0.2,0.2), "cm"))

ue2 <- ue2 + geom_line(data = toDraw,aes(x,y),linetype = 2, size=0.35)

# Now do the evolution of these iterations in variational algo

toPlotData <- data.frame(t = 1:20,x = toSaveLowerBound, type="Lower bound")
toPlotData <- rbind(toPlotData, data.frame(t = 1:20,x = (toSaveAlphaMu[-1]/toSaveBetaMu[-1] + 0.278769)*300 - 1225, type="Mean service rate"))
 
((toSaveAlphaMu[-1]/toSaveBetaMu[-1])*300)[1] + 446.6599

p <- ggplot(toPlotData, aes(x = t, y = x, group = type))
p1<- p + geom_line(size=0.25, aes(linetype=type)) + geom_point(aes(shape=type),size=1.5) +
  scale_x_continuous(name="Iteration",expand = c(0.01,0.01)) + 
  scale_y_continuous(name="Lower bound", sec.axis = sec_axis(~./300 + 1225/300, name = "Mean service rate")) +
  theme_minimal() +
  theme(axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=13,margin=margin(l=0,r=10)),plot.margin = unit(c(0.4,0.2,0.2,0.4), "cm"),legend.position = c(0.8, 0.2)) 


library(grid)
library(gridExtra)
grid.arrange(ue2,p1, ncol=2, nrow=1, top=textGrob("Density estimates, lower bound and mean service rate", gp=gpar(fontsize=15)))
#1000x300
# 
# cairo_ps(filename='filename.eps', width=10, height=4)
# grid.arrange(ue2,ue3, ncol=2, nrow=1, top=textGrob("Prior and posterior densities", gp=gpar(fontsize=15)))
# dev.off()



















