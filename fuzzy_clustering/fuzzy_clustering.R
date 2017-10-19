# Figure 2 Fuzzy Clustering

# Lee, J., Shih, P., Schaedel, O., Rogers, A.K., Quintero-Cadena, P., 
# and Sternberg, P.W. (2017)

# Fuzzy clustering of the 8,042 differentially expressed genes into six common 
# expression profiles

# Install the required packages, then comment out these lines using '#'
install.packages("ggplot2")
install.packages("reshape2")
source("http://bioconductor.org/biocLite.R")
biocLite("Mfuzz")
source("http://bioconductor.org/biocLite.R")
biocLite(c("Biobase"))

# Load the required packages
library(ggplot2)
library("Mfuzz")
library("Biobase")

# Read in the average count data for the 8,042 DEG
deg <- read.csv("input/average_count_data_8042DEG.csv",header=T)

# Convert to matrix
mat_deg <- data.matrix(deg[,2:7]) # transform columns 2-7 into a matrix
rownames(mat_deg) <- deg[,1] # assign row names 

# Convert to ExpressionSet object
dset <- ExpressionSet(assayData=mat_deg)

# Standardize
mat_deg.s <- standardise(dset) #standardization

# Setting parameters
m1 <- mestimate(mat_deg.s)
#[1] 1.72115

# Testing empty clusters
#cplot <- cselection(mat_deg.s, m = m1, crange = seq(2,12,1), repeats = 5,
#                    visu = TRUE)

# Testing centroid distance
#Dplot <- Dmin(mat_deg.s, m = m1, crange = seq(2,12,1), repeats = 5, 
#              visu = TRUE)

# Setting c
set.seed(1992)
c1 <- mfuzz(mat_deg.s, c = 6, m = m1)

# Test plot of the clusters
#mfuzz.plot(mat_deg.s, cl = c1, mfrow = c(3,3), 
#           time.labels = c("noDA24","noDA26","noDA34","noDA60","DA26","DA34"),
#           min.mem = 0.75)

# Test cluster stability
set.seed(1992)
c1.t <- mfuzz(mat_deg.s, c = 6, m = 2)
#mfuzz.plot(mat_deg.s, cl = c1.t, mfrow = c(2,3), 
#           time.labels = c("noDA24","noDA26","noDA34","noDA60","DA26","DA34"),
#           min.mem = 0.5)

# Get list of cluster membership 
write.csv(c1$cluster, file = "output/Mfuzz 8042 DEG clusters.csv")
#(from the manual: a vector of integers containing the indices of the clusters 
#where the data points are assigned to for the closest hard clustering, 
#as obtained by assigning points to the (first) class with maximal membership.)

# Plot the clusters
pdf("output/Mfuzz 8042 DEG clusters.pdf", height = 8.5, width = 11, 
    colormodel = "srgb")
par(mar = c(5,4,4,2)+0.1, oma = c(0,0,0,0))
mfuzz.plot(mat_deg.s, cl = c1, mfrow = c(3,2), 
           time.labels = c("L2d.24","L2d.26","cD","Dauer","cL3","L4"), 
           new.window = FALSE)
dev.off()

# Plot only the core of the clusters
pdf("output/Mfuzz 8042 DEG cluster cores.pdf", height = 8.5, width = 11,
    colormodel = "srgb")
par(mar = c(5,4,4,2)+0.1, oma = c(0,0,0,0))
mfuzz.plot(mat_deg.s, cl = c1, mfrow = c(3,2), 
           time.labels = c("L2d.24","L2d.26","cD","Dauer","cL3","L4"),
           min.mem = 0.75,
           new.window = FALSE)
dev.off()

# Plot the cluster stability test
pdf("output/Mfuzz 8042 DEG cluster Stability Test.pdf", height = 8.5, 
    width = 11, colormodel = "srgb")
par(mar = c(5,4,4,2)+0.1, oma = c(0,0,0,0))
mfuzz.plot(mat_deg.s, cl = c1.t, mfrow = c(3,2), 
           time.labels = c("L2d.24","L2d.26","cD","Dauer","cL3","L4"),
           new.window = FALSE)
dev.off()

# Get cluster sizes
c1$size
#[1] 1102 1921 1025 1497 1332 1165
