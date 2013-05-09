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

Rprof("time")
K = 17
temp <- bnk.null(A, n, K)
Rprof(NULL)
summaryRprof("time")
theta <- temp$theta
rownames(theta) <- tf
save(theta, file=paste("theta",K,".Rdata",sep=""))
save(temp, file=paste("clust",K,".Rdata",sep=""))

sink("degree.na")
degree <- rowSums(A)
for (i in 1:n) {
 cat(paste(tf[i]," = (",degree[i],") \n", sep=""))
}
sink(NULL)

# create pie-data for nodeCharts
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

# sink(paste("cluster_labels",K,".na"))
# cat("pie-labels \n")
# for (i in 1:n) {
#   add <- paste(tf[i]," = (c1",sep="")
#   for (k in 2:K) {
#     add <- paste(add, "::c",k,sep="")
#   }
#   add <- paste(add,") \n",sep="")
#   cat(add)
# }
# sink(NULL)

# create labels for nodeCharts
for (i in 1:K) {
  cat(i,", ",sep="")
}


### this is the command for creating colored pie nodes
#nodecharts pie nodelist="all" attributelist="pie-data" labellist=".,.,.,.,.,.,." colorlist="contrasting"
#########

### compare validating complexes vs clustered complexes ###
# we are interested in comparing the TFs shared by both val and sim cplexes
interest <- intersect(tf, raw.clust[,2])

# these are the validating complexes
raw.clust <- read.table("complexesList.txt", header=F)
raw.clust <- as.matrix(raw.clust)

complex.names <- unique(raw.clust[,1])
val.cplex <- vector(length(complex.names), mode="list")
for (i in 1:length(complex.names)) {
  add <- raw.clust[which(raw.clust[,1]==complex.names[i]),2]
  add <- intersect(add, interest)
  val.cplex[[i]] <- assign(complex.names[i], add)
}
names(val.cplex) <- complex.names

sim.cplex <- vector(K, mode="list")
for (i in 1:K) {
  sim.cplex[[i]] <- intersect(tf[which(theta[,i]>=10)], interest)
}

# search through all sim.cplex to see which sim.cplex the tf is present
search <- function(tf) {
  present <- c()
  for (i in 1:K) {
    check <- (tf %in% sim.cplex[[i]])
    if (check) {
      present <- c(present, 1)
    } else {
      present <- c(present, 0)
    }
  }
  return(present)
}

# list of matrices, each matrix is tfs x K, showing each tf's presence in the clusters 
coverage.mat <- vector(length(val.cplex), mode="list")
col <- sapply(sim.cplex, length)
for (i in 1:length(complex.names)) {
  tfs <- val.cplex[[i]]
  vals <- matrix(0,length(tfs), K)
  colnames(vals) <- col
  for (j in 1:length(tfs)) {
    present <- search(tfs[j])
    vals[j,] <- present
  }
  coverage.mat[[i]] <- vals
}

fit <- lapply(coverage.mat, colSums)
name <- sapply(val.cplex, length)
names(fit) <- name
fit

save(coverage.mat, file=paste("coverage.mat",K,".Rdata",sep=""))
save(fit, file=paste("fit",K,".Rdata",sep=""))

for (i in 1:K) {
  cat("cluster ",i," contains \n")
  cat(sim.cplex[[i]])
  cat("\n")
  cat("\n")
}