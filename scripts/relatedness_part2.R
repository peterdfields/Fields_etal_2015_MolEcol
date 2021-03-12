# Starting with two matrices of pair-wise distance and pair-wise relatedness, move into R


#clear workspace
rm(list=ls())

ls()

# working directory
setwd(".")

require(ade4)
require(ecodist)

# Mantel test
geo.dist <- read.csv("distance.matrix.csv", header=FALSE)
rel1.dist <- read.csv("relatedness.matrix.csv", header=FALSE)
geo.dist <- as.dist(geo.dist)
rel1.distMN <- as.dist(rel1.dist)
mantel(geo.dist~rel1.distMN,nperm = 100000)

# repeat above for subsets of data

# estimate approximate fit of pair-wise relatedness and pair-wise distance from individual target clones
datarel_MN = read.table('relatedness_target_MN.txt', header=TRUE)
plot(datarel_MN$Distance,datarel_MN$relatedness,main="Pair-wise Relatedness", xlab="Distance, Meters", ylab="Relatedness", pch=19)
abline(lm(relatedness~Distance, data=datarel_MN), col="black", lwd=2)
model_datarel_MN <- lm(relatedness~Distance, data=datarel_MN)

# repeat for different target individuals