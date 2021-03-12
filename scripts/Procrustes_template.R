# To be run in R, in order to conduct Procrustes analysis. Can be repeated the same for subsets of individuals or loci.

#clear workspace
rm(list=ls())

ls()

require(fields)
require(vegan)
require(sp)

setwd(".")

data <- read.csv("smartpca_lat_long.csv", na.string="NA", header=TRUE)
A <- cbind(data$eigenvec1,data5$eigenvec2)
B <- cbind(data$lat,data$long)
pro <- procrustes(B, A)

# observation of the pro object includes information about the Procrustes rotation

# generate permutation test in order to quantify significance of Procrustes similarity score
pro_test <- protest(B,A,permutations = 100000)
# print Procrustes similiarity score to be plotted on histogram
x <- pro_test$t

# generate histograms of permuted Procrustes values
hist(x,xlim=c(0,1),breaks=seq(0, 1, by=0.05),lwd=2,main="", ylim=c(0,30000))
abline(v=0.6748,lwd=3, col="black", lty="dashed")



