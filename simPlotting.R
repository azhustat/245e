# STAT 245E, Spring 2013
# Final Project

# plots simulation results     
# took ~7.5hr to run


load(file="simOutput.RData")  # generated in sim.R 
summary(emptnvec)
summary(ordbootvec)
summary(blockbootvec)

            

gheight <- 8
#quartz(width=floor(gheight * 0.8), height= gheight, dpi=100)    

png(paste("simComparison", ".png", sep=""), res=100,
	units="in", width=floor(gheight * 0.5), height= gheight)
par(mfrow=c(3, 1))


bvec <- seq(0, 1, length.out=101)    

hemp <- hist(emptnvec, breaks=bvec, plot=FALSE)
hord <- hist(ordbootvec, breaks=bvec, plot=FALSE)
hblk <- hist(blockbootvec, breaks=bvec, plot=FALSE)  

ylim <- range(c(0, hemp$counts, hord$counts, hord$counts)) 

plot(hemp, ylim=ylim, xlim=c(0,1), main=paste("True distribution"), xlab="", cex.main=2)
plot(hord, ylim=ylim, xlim=c(0,1), main=paste("Ordinary Bootstrap"), xlab="", cex.main=2)
plot(hblk, ylim=ylim, xlim=c(0,1), main=paste("Block Bootstrap"), xlab="", cex.main=2)

par(mfrow=c(1, 1))
dev.off()    


