xPC <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
yPC <- c(0,0.01857981,0.8983868,1.31689,1.478803,1.578659,1.660758,1.724122,1.783467,1.825696,1.870114)

xB <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
yB <- c(0,0.01256404,0.5073522,1.410178,1.903684,2.263564,2.547641,2.776514,2.966509,3.142601,3.282189)

plot(xPC,yPC, type='p', pch=15, cex=2, col='gray', lwd=3, xlab="Alcohol", ylab="Protein Peaks", xlim=c(0,1), ylim=c(0,4))
points(xB,yB, type='p', pch=16, cex=2, col='black', lwd=3)
grid()
legend("topleft", c("PER2/CRY Protein", "BMAL1 Protein"), col=c('gray','black'), pch=c(15,16), pt.cex=2, pt.lwd=3)
title(main="Protein Average Peak Values with Varying Concentrations of Alcohol")

xPC <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
yPC <- c(0,0,21.4,22.7,23.2,23.4,23.5,23.6,23.7,23.8,23.9)

xB <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
yB <- c(0,0,21.4,22.7,23.2,23.4,23.6,23.6,23.7,23.8,23.8)

plot(xPC,yPC, type='p', pch=15, cex=2, col='gray', lwd=3, xlab="Alcohol", ylab="Protein Periods", xlim=c(0,1), ylim=c(20,25))
points(xB,yB, type='p', pch=16, cex=2, col='black', lwd=3)
grid()
legend("topleft", c("PER2/CRY Protein", "BMAL1 Protein"), col=c('gray','black'), pch=c(15,16), pt.cex=2, pt.lwd=3)
title(main="Protein Average Period Values with Varying Concentrations of Alcohol")