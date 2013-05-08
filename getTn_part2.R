# STAT 245E, Spring 2013
# Final Project
  

# on SCF: 
# setwd("~/Documents/13_Spring/245E/FinalProj/code")


# generated in getData.R
load(file="chrlen.RData")
load(file="tfbs.RData")

pairindex <- t(combn(length(tfvec), 2)) 

load("out.RData")
length(out)  # 6903

summary(out)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000e+00 0.000e+00 0.000e+00 9.947e-05 0.000e+00 2.609e-02

sum(out != 0)  # 1193

# checks if the list of TF names are in alphabetical order
sum(tfvec == sort(tfvec))  # 118, check

results <- sapply(which(out != 0), function(i) {
	pair <- tfvec[pairindex[i, ]]
	return(paste(pair, collapse="-"))
})
length(results)  # 1193

# checks with validation pairs
dat <- read.table("interactingPairs.txt", colClasses = "character")
head(dat)
str(dat)                      

ipairs <- sapply(1:dim(dat)[1], function(i) {
	return(paste(as.character(dat[i, ]), collapse="-"))
})
head(ipairs)
length(ipairs)  # 265


sum(results %in% ipairs)   # 68



















