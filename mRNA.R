########################################################
#   alcohol dependent mammalian circadian oscillator   #
########################################################

library(deSolve)

if(!"pracma"%in%rownames(installed.packages())){install.packages("pracma", dependencies=T)}
library(pracma)

############
#   mRNA   #
############

##################
#   Parameters   #
##################

v1b <- 9
k1b <- 1
k1i <- 0.56
c <- 0.01
p <- 8
k1d <- 0.12
k2b <- 0.3
q <- 2
k2d <- 0.05
k2t <- 0.24
k3t <- 0.02
k3d <- 0.12
v4b <- 3.6
k4b <- 2.16
r <- 3
k4d <- 0.75
k5b <- 0.24
k5d <- 0.06
k5t <- 0.45
k6t <- 0.06
k6d <- 0.12
k6a <- 0.09
k7a <- 0.003
k7d <- 0.09
#########
a <- 1
#a0 <- 1
#a1 <- 0.5
#########
#t.switch <- 25
t.switch1 <- 5000
t.switch2 <- 5000
#########


(parms <- c(v1b, k1b, k1i, c, p, k1d, k2b, q, k2d, k2t, k3t, k3d, v4b, k4b, r, k4d, k5b, k5d, k5t, k6t, k6d, k6a, k7a, k7d, a, t.switch1, t.switch2))


times <- seq(0, 5000, by=0.1)


#############
#   Model   #
#############

LVmod <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    if(t > t.switch1){
      a <- 1 ## alcohol value switch
    } 
    if(t > t.switch2){
      a <- 1
    }
    
    # concentration of Per2/Cry mRNA   
    dY1 <- ((a)*(v1b)*(Y7 + c))/((k1b)*(1 + ((Y3/k1i)^p) + (Y7 + c))) - (k1d)*Y1
    # concentration of PER2/CRY complex in the cytoplasm
    dY2 <- (k2b)*((Y1)^q) - (k2d)*(Y2) - (k2t)*(Y2) + (k3t)*(Y3)
    # concentration of PER2/CRY complex in the nucleus
    dY3 <- (k2t)*(Y2) - (k3t)*(Y3) - (k3d)*(Y3)
    # concentration of Bmal1 mRNA
    dY4 <- ((v4b)*((Y3)^r))/(((k4b)^r) + ((Y3)^r)) - (k4d)*(Y4)
    # concentration of BMAL1 protein in the cytoplasm
    dY5 <- (k5b)*(Y4) - (k5d)*(Y5) - (k5t)*(Y5) + (k6t)*(Y6)
    # concentration of BMAL1 protein in the nucleus
    dY6 <- (k5t)*(Y5) - (k6t)*(Y6) - (k6d)*(Y6) + (k7a)*(Y7) - (k6a)*(Y6)
    # concentration of transcriptionally active form BMAL1
    dY7 <- (k6a)*(Y6) - (k7a)*(Y7) - (k7d)*(Y7)
    # the output
    list(c(dY1, dY2, dY3, dY4, dY5, dY6, dY7)) 
  })
}

#########################
#   Intial Conditions   #
#########################

(xstart <- c(Y1 = 0.25, Y2 = 0.3, Y3 = 1.15, Y4 = 0.85, Y5 = 0.72, Y6 = 1.35, Y7 = 1.075))

##############
#   Solver   #
##############

out2  <- rk(xstart, times, LVmod, parms, method = "rk4")
out2 <- as.data.frame(out2)

############
#   Plot   #
############

plot(out2$time, out2$Y1, type='l', col='gray', lty=1, lwd=3, ylab="nM", xlab="Time", ylim=c(0, 3), xlim=c(4900,5000))
lines(out2$time, out2$Y4, type='l', col='black', lty=2, lwd=3)

#lines(out2$time, out2$Y1, type='l', col='red', lty=1, lwd=3)
#lines(out2$time, out2$Y4, type='l', col='red3', lty=2, lwd=3)

grid()

legend("topleft", c("Per2/Cry mRNA", "Bmal1 mRNA", "Per2/Cry mRNA Adjusted", "Bmal1 mRNA Adjusted"), col=c('gray','black','red','red3'), lty=c("solid","dashed","solid","dashed"), lwd=3)
#legend("topleft", c("Per2/Cry mRNA", "Bmal1 mRNA"), col=c('gray','black'), lty=c("solid","dashed"), lwd=3)
title(main="mRNA Oscillations with Alcohol Dependency of 0.1")

#############
#   Peaks   #
#############

findpeaks(out2$Y1, npeaks=50000, threshold=0)
findpeaks(out2$Y4, npeaks=50000, threshold=0)
#                         [,1]              [,2]
#peak 1 [1,]    peak 1 y-value    peak 1 x-value
#peak 2 [2,]    peak 2 y-value    peak 2 x-value

#Per2/Cry mRNA
maxPCm <- findpeaks(out2$Y1, npeaks=50000, threshold=0)
points(0.1*maxPCm[,2]-0.1, maxPCm[,1], pch=0, col="gray", cex=2, lwd=3)
#points(0.1*maxPCm[,2]-0.1, maxPCm[,1], pch=0, col="red", cex=2, lwd=3)
peakPCm <- maxPCm[1,1]
peakPCm

#Bmal1 mRNA
maxBm <- findpeaks(out2$Y4, npeaks=50000, threshold=0)
points(0.1*maxBm[,2]-0.1, maxBm[,1], pch=1, col="black", cex=2, lwd=3)#
points(0.1*maxBm[,2]-0.1, maxBm[,1], pch=1, col="red3", cex=2, lwd=3)
peakBm <- maxBm[1,1]
peakBm

###############
#   Periods   #
###############

#I can go the distance
periodPCm <- 0.1*maxPCm[2,2]-0.1*maxPCm[1,2]
periodPCm

periodBm <- 0.1*maxBm[2,2]-0.1*maxBm[1,2]
periodBm 

#################
#   Averages    #
#################

for(i in 1:210){
crazy <- maxPCm[,1]
crazyagain <- sum(crazy)
averagepeaksPCm <- crazyagain/210
}
averagepeaksPCm

for(i in 1:210){
  crazy <- maxBm[,1]
  crazyagain <- sum(crazy)
  averagepeaksBm <- crazyagain/210
}
averagepeaksBm

for(i in 1:209){
  crazy <- 0.1*maxPCm[i+1,2] - 0.1*maxPCm[i,2]
  crazyagain <- sum(crazy)
  averageperiodPCm <- crazyagain#/209
}
averageperiodPCm

for(i in 1:209){
  crazy <- 0.1*maxBm[i+1,2] - 0.1*maxBm[i,2]
  crazyagain <- sum(crazy)
  averageperiodBm <- crazyagain#/209
}
averageperiodBm

##############################################################
#   Attempt at Graphing Periods and Peaks Against Alcohol    #
##############################################################

#############
#   Peaks   #
#############


#Per2/Cry mRNA
#plot(a, peakPCm, type='p', pch=0, cex=2, col='gray', lwd=3, xlab="Alcohol", ylab="Per2/Cry mRNA Peaks", xlim=c(0,1), ylim=c(0,2))
#points(a, peakPCm, type='p', pch=0, cex=2, col='gray', lwd=3, xlim=c(0,1), ylim=c(0,2))
#grid()
#title(main="Per2/Cry mRNA Peaks vs. Alcohol")

#Bmal1 mRNA
#plot(a, peakBm, type='p', pch=1, cex=2, col='black', lwd=3, xlab="Alcohol", ylab="Bmal1 mRNA Peaks", xlim=c(0,1), ylim=c(0,2))
#points(a, peakBm, type='p', pch=1, cex=2, col='black', lwd=3, xlim=c(0,1), ylim=c(0,2))
#grid()
#title(main="Bmal1 mRNA Peaks vs. Alcohol")

###############
#   Period   #
###############

#Per2/Cry mRNA
#plot(a, periodPCm, type='p', pch=0, cex=2, col='gray', lwd=3, xlab="Alcohol", ylab="Per2/Cry mRNA Period", xlim=c(0,1), ylim=c(15,30))
#points(a, periodPCm, type='p', pch=0, cex=2, col='gray', lwd=3, xlim=c(0,1), ylim=c(20,30))
#grid()
#title(main="Per2/Cry mRNA Period vs. Alcohol")

#Bmal1 mRNA
#plot(a, periodBm, type='p', pch=1, cex=2, col='black', lwd=3, xlab="Alcohol", ylab="Bmal1 mRNA Period", xlim=c(0,1), ylim=c(15,30))
#points(a, periodBm, type='p', pch=1, cex=2, col='black', lwd=3, xlim=c(0,1), ylim=c(15,30))
#grid()
#title(main="Bmal1 mRNA Period vs. Alcohol")
