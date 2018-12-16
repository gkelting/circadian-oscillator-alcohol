########################################################
#   alcohol dependent mammalian circadian oscillator   #
########################################################

library(deSolve)

if(!"pracma"%in%rownames(installed.packages())){install.packages("pracma", dependencies=T)}
library(pracma)

###############
#   Protein   #
###############

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

plot(out2$time, out2$Y3, type='l', col='gray', lty=1, lwd=3, ylab="nM", xlab="Time", ylim=c(0, 5), xlim=c(0,250))
lines(out2$time, out2$Y5 + out2$Y6 + out2$Y7, type='l', col='black', lty=2, lwd=3)

#lines(out2$time, out2$Y3, type='l', col='red', lty=1, lwd=3)
#lines(out2$time, out2$Y5 + out2$Y6 + out2$Y7, type='l', col='red3', lty=2, lwd=3)

grid()

#legend("topleft", c("PER2/CRY Protein", "BMAL1 Protein", "PER2/CRY Protein Adjusted", "BMAL1 Protein Adjusted"), col=c('gray','black','red','red3'), lty=c("solid","dashed","solid","dashed"), lwd=3)
legend("topleft", c("PER2/CRY Protein", "BMAL1 Protein"), col=c('gray','black'), lty=c("solid","dashed"), lwd=3)
title(main="Protein Oscillations with No Alcohol Dependency")

#############
#   Peaks   #
#############

findpeaks(out2$Y3, npeaks=50000, threshold=0)
findpeaks(out2$Y5 + out2$Y6 + out2$Y7, npeaks=50000, threshold=0)
#                         [,1]              [,2]
#peak 1 [1,]    peak 1 y-value    peak 1 x-value
#peak 2 [2,]    peak 2 y-value    peak 2 x-value

#Per2/Cry Protein
maxPCp <- findpeaks(out2$Y3, npeaks=50000, threshold=0)
points(0.1*maxPCp[,2]-0.1, maxPCp[,1], pch=15, col="gray", cex=2, lwd=3)
#points(0.1*maxPCp[,2]-0.1, maxPCp[,1], pch=15, col="red", cex=2, lwd=3)
peakPCp <- maxPCp[1,1]

#Bmal1 Protein
maxBp <- findpeaks(out2$Y5 + out2$Y6 + out2$Y7, npeaks=50000, threshold=0)
points(0.1*maxBp[,2]-0.1, maxBp[,1], pch=16, col="black", cex=2, lwd=3)
#points(0.1*maxBp[,2]-0.1, maxBp[,1], pch=16, col="red3", cex=2, lwd=3)
peakBp <- maxBp[1,1]

###############
#   Periods   #
###############

#I can go the distance
periodPCp <- 0.1*maxPCp[2,2]-0.1*maxPCp[1,2]
periodPCp

periodBp <- 0.1*maxBp[2,2]-0.1*maxBp[1,2]
periodBp 

#################
#   Averages    #
#################

for(i in 1:210){
  crazy <- maxPCp[,1]
  crazyagain <- sum(crazy)
  averagepeaksPCp <- crazyagain/210
}
averagepeaksPCp

for(i in 1:209){
  crazy <- maxBp[,1]
  crazyagain <- sum(crazy)
  averagepeaksBp <- crazyagain/209
}
averagepeaksBp

for(i in 1:209){
  crazy <- 0.1*maxPCp[i+1,2] - 0.1*maxPCp[i,2]
  crazyagain <- sum(crazy)
  averageperiodPCp <- crazyagain#/209
}
averageperiodPCp

for(i in 1:208){
  crazy <- 0.1*maxBp[i+1,2] - 0.1*maxBp[i,2]
  crazyagain <- sum(crazy)
  averageperiodBp <- crazyagain#/209
}
averageperiodBp

##############################################################
#   Attempt at Graphing Periods and Peaks Against Alcohol    #
##############################################################
# 
# #############
# #   Peaks   #
# #############
# 
# #Per2/Cry protein
# plot(a, peakPCp, type='p', pch=15, cex=2, col='gray', lwd=3, xlab="Alcohol", ylab="PER2/CRY Protein Peaks", xlim=c(0,1), ylim=c(0,4))
# points(a, peakPCp, type='p', pch=15, cex=2, col='gray', lwd=3, xlim=c(0,1), ylim=c(0,4))
# grid()
# title(main="PER2/CRY Protein Peaks vs. Alcohol")
# 
# #Bmal1 protein
# plot(a, peakBp, type='p', pch=16, cex=2, col='black', lwd=3, xlab="Alcohol", ylab="BMAL1 Protein Peaks", xlim=c(0,1), ylim=c(0,4))
# points(a, peakBp, type='p', pch=16, cex=2, col='black', lwd=3, xlim=c(0,1), ylim=c(0,4))
# grid()
# title(main="BMAL1 Protein Peaks vs. Alcohol")
# 
# ###############
# #   Period   #
# ###############
# 
# #Per2/Cry protein
# plot(a, periodPCp, type='p', pch=15, cex=2, col='gray', lwd=3, xlab="Alcohol", ylab="PER2/CRY Protein Period", xlim=c(0,1), ylim=c(23,25))
# points(a, periodPCp, type='p', pch=15, cex=2, col='gray', lwd=3, xlim=c(0,1), ylim=c(23,25))
# grid()
# title(main="PER2/CRY Protein Period vs. Alcohol")
# 
# #Bmal1 protein
# plot(a, periodBp, type='p', pch=16, cex=2, col='black', lwd=3, xlab="Alcohol", ylab="BMAL1 Protein Period", xlim=c(0,1), ylim=c(23,25))
# points(a, periodBp, type='p', pch=16, cex=2, col='black', lwd=3, xlim=c(0,1), ylim=c(23,25))
# grid()
# title(main="BMAL1 Protein Period vs. Alcohol")
