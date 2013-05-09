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

source("getTnFuncDef.R")

# cat(paste(paste("\"", ls(), "\"", sep=""), collapse=", "), fill=TRUE)
clusterExport(cl, c("chrlen", "pairindex", "tfbslist", "tfvec", "getObsTn"))   

obsTn <- parSapply(cl, 1:dim(pairindex)[1], getOverlap)  
save(obsTn, file="obsTn.RData")

###
stopCluster(cl)
q("no")
  