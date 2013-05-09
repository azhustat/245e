
getObsTn <- function(i) {
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
		
		for (j in 1:2) {
			if (length(get(paste("v", j, sep=""))) > 0) {
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
		}  
		
		ivec <- jvec <-   ijlist[[cnum]]
		for (k in 0:(chrlen[cnum] %/% wsize)) {
			ivec[k+1] <- sum(mat[(k*wsize+1):min((k+1)*wsize, chrlen[cnum]), 1]) > 0
			jvec[k+1] <- sum(mat[(k*wsize+1):min((k+1)*wsize, chrlen[cnum]), 2]) > 0
		}                                 
		ijlist[[cnum]] <- ivec - jvec
	}   

	ijdiff <- unlist(ijlist)
	# length(ijdiff)  # 20127
	nloci <- length(ijdiff)
	
	# needs to convert integer to double to prevent integer overflow
	out <- sapply(1:nloci, function(k) {
		return( max(abs(cumsum(as.numeric(ijdiff[k:nloci])))))
	})

	return(max(abs(out)) / sqrt(nloci) )
}
