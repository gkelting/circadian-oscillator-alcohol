xPC <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
yPC <- c(0,0.1,0.9013391,1.232113,1.329201,1.387558,1.423786,1.461335,1.484509,1.502413,1.523759)

xB <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
yB <- c(0,0.01,0.264212,0.7978921,1.084644,1.280267,1.435601,1.558614,1.659156,1.751552,1.832127)

plot(xPC,yPC, type='p', pch=0, cex=2, col='gray', lwd=3, xlab="Alcohol", ylab="mRNA Peaks", xlim=c(0,1), ylim=c(0,2))
points(xB,yB, type='p', pch=1, cex=2, col='black', lwd=3)
grid()
legend("topleft", c("Per2/Cry mRNA", "Bmal1 mRNA"), col=c('gray','black'), pch=c(0,1), pt.cex=2, pt.lwd=3)
title(main="mRNA Average Peak Values with Varying Concentrations of Alcohol")

xPC <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
yPC <- c(0,0,21.4,22.7,23.1,23.4,23.5,23.6,23.8,23.8,23.8)

xB <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
yB <- c(0,0,21.4,22.7,23.2,23.4,23.5,23.6,23.7,23.8,23.8)

plot(xPC,yPC, type='p', pch=0, cex=2, col='gray', lwd=3, xlab="Alcohol", ylab="mRNA Periods", xlim=c(0,1), ylim=c(20,25))
points(xB,yB, type='p', pch=1, cex=2, col='black', lwd=3)
grid()
legend("topleft", c("Per2/Cry mRNA", "Bmal1 mRNA"), col=c('gray','black'), pch=c(0,1), pt.cex=2, pt.lwd=3)
title(main="mRNA Average Period Values with Varying Concentrations of Alcohol")