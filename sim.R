# STAT 245E, Spring 2013
# Final Project

# simulation



# on SCF: 
# setwd("~/Documents/13_Spring/245E/FinalProj/code")

# need to use cluster on SCF, for more details see:
# http://statistics.berkeley.edu/computing/servers/cluster
# each user can use a maximum of 12 cores for high.q

# to run this job: 
# qsub -q high.q -pe smp 12 sim.sh 

# to check status of all users: 
# qstat -u "*"
                 

# RNG with parallel computing:
# no extra step when using package "parallel"
# see http://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf


### MAIN

# cluster setup                                             
require(parallel)
nCores <- 12  # max allowed     
cl <- makeCluster(nCores) # by default this uses sockets                       

# Bickel 2010 Subsampling paper, section 5

n <- 10000
L <- 40  # block size
nsim <- 1e5

# parameters not specified in the paper
p0 <- 0.5                              
w <- 8

getProb <- function(k, w, x) {
	prob <- p0 / 2
	lower <- max(1, k - w)
	upper <- k - 1
	if (upper > 1) {
		prob <- prob + sum(x[lower:upper]) / w * (1 - p0)
	} 
	return(prob)
}
                                                        
generateSeq <- function(n, p0, w) {
	x <- vector("integer", length=n)
	for (k in 1:n) {
		prob <- getProb(k, w, x) 
		if (runif(1) <= prob) {
			x[k] <- 1
		}
	}   
	return(x)
}

getTn <- function(x, n) {
	# Feature I: the occurrence of sequence 11,100 starting at position k
	# Feature II: the occurrence of more than six 1’s in the next 10 consecutive
	# positions including the current position k.         
	ivec <- jvec <-  vector("integer", length=n)
	for (k in 1:n) { 
		if (k <= n - 4) {
			if (paste(x[k:(k+4)], collapse="") == "11100") {
				ivec[k] <- 1
			}
			
			if (k <= n - 9) {
				if (sum(x[k:(k+9)]) > 6) {
					jvec[k] <- 1
				}
			}
		}
	}                           
	
	ijprod <- ivec * jvec
	return(max(cumsum(ijprod) / (1:n)))
}

runSingleSim <- function(dummy) {
	x <-  generateSeq(n, p0, w)
	tn <- getTn(x, n)
	return(tn)
}

system.time(tn <- runSingleSim(42))
 
#  user  system elapsed 
# 0.612   0.002   0.615  


# cat(paste(paste("\"", ls(), "\"", sep=""), collapse=", "), fill=TRUE)
clusterExport(cl, c("generateSeq", "getProb", "getTn", "L", "n", 
	"nsim", "p0", "runSingleSim", "w"))   
 
# empirical distribution of estimated S from 10,000 random sequences 
# generated under model  
# emptnvec <- replicate(nsim, runSingleSim())
emptnvec <- parSapply(cl, 1:nsim, runSingleSim)         




# The Ordinary Bootstrap distribution is derived by performing a base-by-base 
# uniform sampling of the sequence x1,...,xn to construct 10,000 sequences 
# of length n.
singleOrdBootstrap <- function(dummy) {
	bootsam <- sample(xsam, size=length(xsam), replace=TRUE)
	tn <- getTn(bootsam, n)
	return(tn)
}     
# ordbootvec <- replicate(nsim, singleOrdBootstrap(xsam)) 


# The block subsampling distribution is derived by drawing independent 
# samples of blocks of length L = 40 and stringing the blocks together 
# to construct 10,000 sequences of length n.
singleBlockBootstrap <- function(dummy) {  
	index <- sample(1:(length(xsam) - L + 1), size=n/L, replace=TRUE)
	bootsam <- sapply(index, function(i) {
		return(xsam[i:(i+L-1)])
	}) 
	bootsam <- as.vector(bootsam)  # stack by columns by default
	tn <- getTn(bootsam, n)
	return(tn)
}  
# blockbootvec <- replicate(nsim, singleBlockBootstrap(xsam, L))
 

xsam <-  generateSeq(n, p0, w)
clusterExport(cl, c("generateSeq", "getProb", "getTn", "L", "n", 
	"nsim", "p0", "runSingleSim", "singleBlockBootstrap", 
	"singleOrdBootstrap", "w", "xsam"))  
ordbootvec <- parSapply(cl, 1:nsim, singleOrdBootstrap) 
blockbootvec <- parSapply(cl, 1:nsim, singleBlockBootstrap) 

### 
save(emptnvec, ordbootvec, blockbootvec, p0, L, n, w, nsim, file="simOutput.RData")
q("no")

 










