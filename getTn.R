# STAT 245E, Spring 2013
# Final Project
  

# on SCF: 
# setwd("~/Documents/13_Spring/245E/FinalProj/code")

# need to use cluster on SCF, for more details see:
# http://statistics.berkeley.edu/computing/servers/cluster
# each user can use a maximum of 12 cores for high.q

# to run this job: 
# qsub -q high.q -pe smp 12 getTn.sh 

# to check status of all users: 
# qstat -u "*"



### MAIN

# cluster setup                                             
require(parallel)
nCores <- 12  # max allowed     
cl <- makeCluster(nCores) # by default this uses sockets

# generated in getData.R
load(file="chrlen.RData")
load(file="tfbs.RData")

pairindex <- t(combn(length(tfvec), 2)) 

getOverlap <- function(i) {
	dat1 <- tfbslist[[pairindex[i, 1]]]
	dat2 <- tfbslist[[pairindex[i, 2]]]  

	ijlist <- vector("list", length=16) 
	# window size
	wsize <- 600 

	for (cnum in 1:16) {
		mat <- matrix(0, ncol=2, nrow=sum(chrlen[cnum]))
		chrname <- paste("chr", cnum, sep="")
		v1 <- which(dat1$Chromosome == chrname)
		v2 <- which(dat2$Chromosome == chrname)
		
		ijlist[[cnum]] <- rep(0, chrlen[cnum] %/% wsize + 1)   

		if (length(v1) > 0 & length(v2) > 0) {      
			# if one of the TF has no binding sites on the given chromosome,
			# then there is no need to go further 
			for (j in 1:2) {
				for (k in 1:length(get(paste("v", j, sep="")))) {
					istart <- get(paste("dat", j, sep=""))[
					get(paste("v", j, sep=""))[k], "Start"]
					iend <- get(paste("dat", j, sep=""))[
					get(paste("v", j, sep=""))[k], "End"]
					if (iend > chrlen[cnum]) {
						cat(paste("ERROR: check i =", i, "and", chrname), fill=TRUE)
					} else {
						mat[istart:iend, j] <-  1						
					}
				}
			}
			
			ivec <- jvec <-   ijlist[[cnum]]
			for (k in 0:(chrlen[cnum] %/% wsize)) {
				ivec[k+1] <- sum(mat[(k*wsize+1):min((k+1)*wsize, chrlen[cnum]), 1]) > 0
				jvec[k+1] <- sum(mat[(k*wsize+1):min((k+1)*wsize, chrlen[cnum]), 2]) > 0
			}                                 
			ijlist[[cnum]] <- ivec * jvec
		}  
		
		
	}   

	ijprod <- unlist(ijlist)
	# length(ijprod)  # 20127
	
	# needs to convert index to double to prevent integer overflow
	return(max(as.numeric(cumsum(ijprod)) / (1:length(ijprod))))
}

# cat(paste(paste("\"", ls(), "\"", sep=""), collapse=", "), fill=TRUE)
clusterExport(cl, c("chrlen", "pairindex", "tfbslist", "tfvec"))   

out <- parSapply(cl, 1:dim(pairindex)[1], getOverlap)  
save(out, file="out.RData")

###
stopCluster(cl)
q("no")
  