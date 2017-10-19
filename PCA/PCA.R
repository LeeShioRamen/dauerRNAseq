# Figure S1 PCA and Scree plot

# Lee, J., Shih, P., Schaedel, O., Rogers, A.K., Quintero-Cadena, P., 
# and Sternberg, P.W. (2017)

# Generation of the principal component analysis plot of the variation in gene 
# expression across the 12 sequenced samples

# This script follows parts of the DESeq protocol Anders et al. (2013) Nat. 
# Protocols

# Install the required packages, then comment out these lines using '#'
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq")

# Load the required packages
library("DESeq")
library("genefilter")

# Collect metadata of experimental design
metadata <- read.csv("input/metadata_for_DESeq.csv")

# Count reads using HTSEQ-COUNT
metadata$countf = paste(metadata$LibraryName, "count", sep=".")
samplesDESeq=with(metadata,
                  data.frame(shortname = I(shortname), 
                             countf = paste("input/", metadata$countf, sep=""),
                             condition = condition,
                             LibraryLayout = LibraryLayout))
cds = newCountDataSetFromHTSeqCount(samplesDESeq) # Create a CountDataSet

# Estimate normalization factors
cds = estimateSizeFactors(cds)
sizeFactors(cds)
head(counts(cds, normalized = TRUE))

# Inspect sample relationships
cdsB = estimateDispersions(cds,method="blind")
vsd = varianceStabilizingTransformation(cdsB)
p = plotPCA(vsd,intgroup=c("condition","LibraryLayout"))
p

# Estimate dispersion
cds = estimateDispersions(cds)
plotDispEsts(cds)

# Customize PCA plot settings
plotPCAcol <- function (x, intgroup = "condition", ntop = 500) 
{
  rv = rowVars(exprs(x))
  select = order(rv, decreasing = TRUE)[seq_len(ntop)]
  pca = prcomp(t(exprs(x)[select, ]))
  fac = factor(apply(pData(x)[, intgroup, drop = FALSE], 1, 
                     paste, collapse = " : "))
  if (length(fac) >= 3) 
    colours = c("#88D559","#46C92C","#FFB620","#EFA100","#FF5673","#DD2067")
  else colours = c("green", "blue")
  xyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$x), 
         pch = 16, cex = 2, aspect = "iso", col = colours, 
         main = draw.key(key = list(rect = list(col = colours), 
                                    text = list(levels(fac)), rep = FALSE)))
}

plotPCAcol(vsd,intgroup=c("condition"))

# Draw Scree Plot
  # Take from the plotPCAcol function:
rv = rowVars(exprs(vsd))
select = order(rv, decreasing = TRUE)[seq_len(500)]
pca = prcomp(t(exprs(vsd)[select, ]))
fac = factor(apply(pData(vsd)[, c("condition"), drop = FALSE], 1, 
                   paste, collapse = " : "))
  # To see the function for screeplot
getS3method("screeplot", "default")
  # Customize scree plot settings
myscree <- function (x, npcs = min(10, length(x$sdev)),
                     type = c("barplot", "lines"), 
                     main = deparse(substitute(x)), ...) 
{
  main
  type <- match.arg(type)
  pcs <- x$sdev^2
  xp <- seq_len(npcs)
  dev.hold()
  on.exit(dev.flush())
  if (type == "barplot") 
    barplot(pcs[xp], names.arg = names(pcs[xp]), main = "", 
            xlab = "PC", ylab = "Variances", ...)
  else {
    plot(xp, pcs[xp], type = "b", axes = FALSE, main = "", pch=16,
         xlab = "PC", ylab = "Variances", ...)
    axis(2)
    axis(1, at = xp, labels = names(pcs[xp]))
  }
  invisible()
}
  # Generate the scree plot
myscree(pca,npcs=12,"lines")

# To plot percent of variance explained by each PC
noPCs = 12
xl <- seq_len(noPCs)
var.per <- (pca$sdev^2/sum(pca$sdev^2))*100
  # Generate the scree plot with Percentage of Variance as the y-axis
par(mar=c(5.1, 5.1, 4.1, 1.1),mgp=c(3.5,1,0))
plot(xl, var.per, type = "b", axes = FALSE, main = "", pch=16,
     xlab = "PC", ylab = "Percentage of Variance",ylim=c(0,70),cex.lab=1.2)
axis(2, las=1)
axis(1, at = xl, labels = xl)
