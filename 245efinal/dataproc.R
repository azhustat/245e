# data processing
raw <- read.table("interactingPairs.txt", stringsAsFactors=F)

tf <- unique(union(raw[,1], raw[,2]))

n <- length(tf)
A <- matrix(0,n,n)

for (i in 1:dim(raw)[1]) {
  one <- raw[i,1]
  one
  two <- raw[i,2]
  two
  ind1 <- which(tf==one)
  ind1
  ind2 <- which(tf==two)
  ind2
  A[ind1,ind2] <- 1
  A[ind2,ind1] <- 1
}

library(igraph)

###### clustering ######

K = 7
temp <- bnk(A, n, K)
theta <- temp$theta
rownames(theta) <- tf
save(theta, file=paste("clust",K,".Rdata",sep=""))

sink("degree.na")
degree <- rowSums(A)
for (i in 1:n) {
 cat(paste(tf[i]," = (",degree[i],") \n", sep=""))
}
sink(NULL)


for (i in 1:n) {
  theta[i,which(theta[i,] < 0.01)] <- 0
  m <- sum(theta[i,])
  theta[i,] <- theta[i,]/m
}
theta <- round(theta * 100)

sink(paste("clust",K,".na",sep=""))
cat("pie-data \n")
for (i in 1:n) {
  add <- paste(tf[i]," = (", theta[i,1],".0", sep="")
  for (j in 2:K) {
    add <- paste(add,"::",theta[i,j],".0",sep="")
  }
  add <- paste(add,") \n",sep="")
  cat(add)
}
sink(NULL)

sink(paste("cluster_labels.na"))
cat("pie-labels \n")
for (i in 1:n) {
  add <- paste(tf[i]," = (c1",sep="")
  for (k in 2:K) {
    add <- paste(add, "::c",k,sep="")
  }
  add <- paste(add,") \n",sep="")
  cat(add)
}
sink(NULL)


### this is the command for creating colored pie nodes
#nodecharts pie nodelist="all" attributelist="pie-data" labellist=".,.,.,.,.,.,." colorlist="contrasting"
#