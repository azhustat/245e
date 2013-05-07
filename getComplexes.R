# STAT 245E, Spring 2013
# Final Project    

# gets verified complexes curated by Taly ;)

# cd ~/Documents/Academics/Berkeley/13_Spring/245E/FinalProj/code


comp <- read.csv("../data/CLUSTERS.csv", head=FALSE, colClasses=c("character", 
	"character", rep("NULL", 4)))  
 
index <- which(comp[, 1] == "")        
comp <- comp[-index, ]

str(comp)
dim(comp)  # 135 2

lvec <- sapply(comp[, 1], function(cname) {
	out <- gsub(":", "", cname)
	out <- gsub(" ", "", out)
	return(out)
})

comp[, 1] <- lvec
comp[, 2] <- gsub("-", "", comp[, 2])


write.table(comp, file="complexes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
system("head -n 20 complexes.txt")


###

comp <- read.table("complexes.txt")
str(comp)

table(comp[, 1])

length(unique(toupper(as.character(comp[, 2]) ) ))
