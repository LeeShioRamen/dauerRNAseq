# Figure 3 Neuronal Genome Heatmap

# Lee, J., Shih, P., Schaedel, O., Rogers, A.K., Quintero-Cadena, P., 
# and Sternberg, P.W. (2017)

# Differential expression of the neuronal effector genome of C. elegans during 
# dauer and reproductive development

# Install the required packages, then comment out these lines using '#'
install.packages("GMD")
install.packages("gplots")
install.packages("gtools")
install.packages("RColorBrewer")

# Load the required packages
library(GMD)
library(gplots)
library(gtools)
library(RColorBrewer)

# Read in the average counts for the differentially expressed neuronal genes
data <- read.csv("input/neuronal_genome_counts.csv", header = T)
row.names(data) <- data$X
data$X <- NULL

# Convert to matrix
x <- as.matrix(data)

# Scale and center rows
x.rowscaled<-t(scale(t(x)))

# Set heatmapping colors
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

# Define correlation similarity metric and average linkage clustering
distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method = "average")

# Draw heatmap
pdf("output/neuronal_genome_heatmap_unclustered.pdf",
    height = 11, width = 8.5, colormodel = "srgb")
par(mar=c(0,0,0,0))
heatmap.2(x.rowscaled, Rowv = NA, Colv = NA, distfun = distCor, 
          hclustfun = hclustAvg, cexRow = 0.1, dendrogram = "none", 
          scale = "none", col = rev(cols), trace = "none", symbreaks = T, 
          keysize = 1.5, lmat = rbind(rep(0,2), 4:3, 2:1, rep(0,2)), 
          lwid = c(2,4), lhei = c(0.1,0.8,4.5,0.2))
dev.off()
