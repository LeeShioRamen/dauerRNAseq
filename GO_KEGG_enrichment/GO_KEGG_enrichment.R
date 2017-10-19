# Figure S3 GO and KEGG enrichment analysis

# Lee, J., Shih, P., Schaedel, O., Rogers, A.K., Quintero-Cadena, P., 
# and Sternberg, P.W. (2017)

# Analysis to identify the most enriched GO (based on descending fold 
# enrichment) and KEGG (based on ascending q value) terms in clusters 1 to 6

# Install the required packages, then comment out these lines using '#'
install.packages("ggplot2")
source("http://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")

# Load the required packages
library(ggplot2)
library("clusterProfiler")

# GO Enrichment Analysis of the Clusters
# Read in GO data for Cluster 1
GO.1 <- read.csv("input/GO/Cluster1.csv", 
                 colClasses = c("character", rep("numeric", 3), 
                                rep("character", 2), "numeric"))
GO.1$Cluster <- "Cluster 1"
colnames(GO.1) <- c("GO", "In.ref", "In.input", "Expect", "Over.under", "FE",
                    "Bonferroni.pvalue", "Cluster")
GO.1 <- GO.1[GO.1$Over.under == "+",]

# Cluster 2
GO.2 <- read.csv("input/GO/Cluster2.csv", 
                 colClasses = c("character", rep("numeric", 3), 
                                rep("character", 2), "numeric"))
GO.2$Cluster <- "Cluster 2"
colnames(GO.2) <- c("GO", "In.ref", "In.input", "Expect", "Over.under", "FE",
                    "Bonferroni.pvalue", "Cluster")
GO.2 <- GO.2[GO.2$Over.under == "+",]

# Cluster 3
GO.3 <- read.csv("input/GO/Cluster3.csv", 
                 colClasses = c("character", rep("numeric", 3), 
                                rep("character", 2), "numeric"))
GO.3$Cluster <- "Cluster 3"
colnames(GO.3) <- c("GO", "In.ref", "In.input", "Expect", "Over.under", "FE",
                    "Bonferroni.pvalue", "Cluster")
GO.3 <- GO.3[GO.3$Over.under == "+",]

# Cluster 4
GO.4 <- read.csv("input/GO/Cluster4.csv", 
                 colClasses = c("character", rep("numeric", 3), 
                                rep("character", 2), "numeric"))
GO.4$Cluster <- "Cluster 4"
colnames(GO.4) <- c("GO", "In.ref", "In.input", "Expect", "Over.under", "FE",
                    "Bonferroni.pvalue", "Cluster")
GO.4 <- GO.4[GO.4$Over.under == "+",]

# Cluster 5
GO.5 <- read.csv("input/GO/Cluster5.csv", 
                 colClasses = c("character", rep("numeric", 3), 
                                rep("character", 2), "numeric"))
GO.5$Cluster <- "Cluster 5"
colnames(GO.5) <- c("GO", "In.ref", "In.input", "Expect", "Over.under", "FE",
                    "Bonferroni.pvalue", "Cluster")
GO.5 <- GO.5[GO.5$Over.under == "+",]

# Cluster 6
GO.6 <- read.csv("input/GO/Cluster6.csv", 
                 colClasses = c("character", rep("numeric", 3), 
                                rep("character", 2), "numeric"))
GO.6$Cluster <- "Cluster 6"
colnames(GO.6) <- c("GO", "In.ref", "In.input", "Expect", "Over.under", "FE",
                    "Bonferroni.pvalue", "Cluster")
GO.6 <- GO.6[GO.6$Over.under == "+",]

# Combine to one table
GO <- rbind(GO.1,GO.2,GO.3,GO.4,GO.5,GO.6)
write.csv(GO, "input/GO/Combined GO.csv", row.names = F)

# Read back in to make FE numeric
GO <- read.csv("input/GO/Combined GO.csv")

# Order GO terms by descending FE
GO.1 <- GO[GO$Cluster == "Cluster 1",]
GO.1 <- GO.1[order(-GO.1$FE, GO.1$Bonferroni.pvalue),]
GO.2 <- GO[GO$Cluster == "Cluster 2",]
GO.2 <- GO.2[order(-GO.2$FE, GO.2$Bonferroni.pvalue),]
GO.3 <- GO[GO$Cluster == "Cluster 3",]
GO.3 <- GO.3[order(-GO.3$FE, GO.3$Bonferroni.pvalue),]
GO.4 <- GO[GO$Cluster == "Cluster 4",]
GO.4 <- GO.4[order(-GO.4$FE, GO.4$Bonferroni.pvalue),]
GO.5 <- GO[GO$Cluster == "Cluster 5",]
GO.5 <- GO.5[order(-GO.5$FE, GO.5$Bonferroni.pvalue),]
GO.6 <- GO[GO$Cluster == "Cluster 6",]
GO.6 <- GO.6[order(-GO.6$FE, GO.6$Bonferroni.pvalue),]

# Combine to one table again
GO <- rbind(GO.1,GO.2,GO.3,GO.4,GO.5,GO.6)
write.csv(GO, "input/GO/Combined GO.csv", row.names = F)

# Filter by taking top 5 GO terms (based on descending FE) per Cluster
filt1 <- head(GO[GO$Cluster == "Cluster 1",],5)
filt2 <- head(GO[GO$Cluster == "Cluster 2",],5)
filt3 <- head(GO[GO$Cluster == "Cluster 3",],5)
filt4 <- head(GO[GO$Cluster == "Cluster 4",],5)
filt5 <- head(GO[GO$Cluster == "Cluster 5",],5)
filt6 <- head(GO[GO$Cluster == "Cluster 6",],5)
GO.filt <- rbind(filt1,filt2,filt3,filt4,filt5,filt6)
write.csv(GO.filt, "input/GO/Top 5 GO Terms Per Cluster.csv", row.names = F)

# Format for R in Excel

# Read in formatted data
GO.final <- read.csv("input/GO/Top 5 GO Terms Per Cluster for R.csv")
GO.final$Cluster <- factor(GO.final$Cluster, 
                           levels = c('Cluster 1','Cluster 2','Cluster 3',
                                      'Cluster 4','Cluster 5','Cluster 6'))
GO.final$GO <- factor(GO.final$GO, 
                      levels = rev(as.character(GO.final$GO[c(1:30)])))
GO.final$InBox <- paste(GO.final$In.input, GO.final$In.ref, sep="/")
GO.final$InBox[GO.final$InBox=="NA/NA"] <- NA

# Set up plot
g<- ggplot(data = GO.final, aes(x = GO, y = FE))
#must use 'stat="identity"' when y-values are used
g<- g + geom_bar(stat = "identity", aes(fill = GO))
#guide legend
g<- g + scale_fill_discrete(guide = FALSE)
#facet by 'Cluster'
g<- g + facet_wrap(~Cluster)
#flip coordinates
g<- g + coord_flip()
#add observed number to each bar
g<- g + geom_text(aes(x = GO, y = FE/2, label = In.input))
#add labels
g<- g + labs(x = "GO Term", y = "Fold Enrichment")
#move y label
g<- g + theme(axis.title.y = element_text(margin = margin(0,17,0,0)))
#plot
g


# KEGG Enrichment Analysis of the Clusters
# Read in KEGG data for Cluster 1
K1 <- read.csv("input/KEGG/Ensembl BioMart Entrez IDs Cluster 1.csv")
K1_ent <- as.character(K1$NCBI.gene.ID)
K1_e <- enrichKEGG(gene = K1_ent, organism="cel", keyType = 'ncbi-geneid',
                   pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05)

# Read in KEGG data for Cluster 2
K2 <- read.csv("input/KEGG/Ensembl BioMart Entrez IDs Cluster 2.csv")
K2_ent <- as.character(K2$NCBI.gene.ID)
K2_e <- enrichKEGG(gene = K2_ent, organism="cel", keyType = 'ncbi-geneid',
                   pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05)

# Read in KEGG data for Cluster 3
K3 <- read.csv("input/KEGG/Ensembl BioMart Entrez IDs Cluster 3.csv")
K3_ent <- as.character(K3$NCBI.gene.ID)
K3_e <- enrichKEGG(gene = K3_ent, organism="cel", keyType = 'ncbi-geneid',
                   pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05)

# Read in KEGG data for Cluster 4
K4 <- read.csv("input/KEGG/Ensembl BioMart Entrez IDs Cluster 4.csv")
K4_ent <- as.character(K4$NCBI.gene.ID)
K4_e <- enrichKEGG(gene = K4_ent, organism="cel", keyType = 'ncbi-geneid',
                   pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05)

# Read in KEGG data for Cluster 5
K5 <- read.csv("input/KEGG/Ensembl BioMart Entrez IDs Cluster 5.csv")
K5_ent <- as.character(K5$NCBI.gene.ID)
K5_e <- enrichKEGG(gene = K5_ent, organism="cel", keyType = 'ncbi-geneid',
                   pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05)

# Read in KEGG data for Cluster 6
K6 <- read.csv("input/KEGG/Ensembl BioMart Entrez IDs Cluster 6.csv")
K6_ent <- as.character(K6$NCBI.gene.ID)
K6_e <- enrichKEGG(gene = K6_ent, organism="cel", keyType = 'ncbi-geneid',
                   pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05)

# Convert to data frames
K1d <- data.frame(Cluster = "Cluster 1",
                  Description = K1_e@result$Description,
                  qvalue = K1_e@result$qvalue,
                  GeneRatio = K1_e@result$GeneRatio,
                  BgRatio = K1_e@result$BgRatio,
                  Count = K1_e@result$Count)
K1d$FE <- as.numeric(sapply(K1_e@result$GeneRatio, function(x) 
  eval(parse(text=x)))/sapply(K1_e@result$BgRatio, 
                              function(x) eval(parse(text=x))))


K2d <- data.frame(Cluster = "Cluster 2",
                  Description = K2_e@result$Description,
                  qvalue = K2_e@result$qvalue,
                  GeneRatio = K2_e@result$GeneRatio,
                  BgRatio = K2_e@result$BgRatio,
                  Count = K2_e@result$Count)
K2d$FE <- as.numeric(sapply(K2_e@result$GeneRatio, function(x) 
  eval(parse(text=x)))/sapply(K2_e@result$BgRatio, 
                              function(x) eval(parse(text=x))))

K3d <- data.frame(Cluster = "Cluster 3",
                  Description = K3_e@result$Description,
                  qvalue = K3_e@result$qvalue,
                  GeneRatio = K3_e@result$GeneRatio,
                  BgRatio = K3_e@result$BgRatio,
                  Count = K3_e@result$Count)
K3d$FE <- as.numeric(sapply(K3_e@result$GeneRatio, function(x) 
  eval(parse(text=x)))/sapply(K3_e@result$BgRatio, 
                              function(x) eval(parse(text=x))))

K4d <- data.frame(Cluster = "Cluster 4",
                  Description = K4_e@result$Description,
                  qvalue = K4_e@result$qvalue,
                  GeneRatio = K4_e@result$GeneRatio,
                  BgRatio = K4_e@result$BgRatio,
                  Count = K4_e@result$Count)
K4d$FE <- as.numeric(sapply(K4_e@result$GeneRatio, function(x) 
  eval(parse(text=x)))/sapply(K4_e@result$BgRatio, 
                              function(x) eval(parse(text=x))))

K5d <- data.frame(Cluster = "Cluster 5",
                  Description = K5_e@result$Description,
                  qvalue = K5_e@result$qvalue,
                  GeneRatio = K5_e@result$GeneRatio,
                  BgRatio = K5_e@result$BgRatio,
                  Count = K5_e@result$Count)
K5d$FE <- as.numeric(sapply(K5_e@result$GeneRatio, function(x) 
  eval(parse(text=x)))/sapply(K5_e@result$BgRatio, 
                              function(x) eval(parse(text=x))))

K6d <- data.frame(Cluster = "Cluster 6",
                  Description = K6_e@result$Description,
                  qvalue = K6_e@result$qvalue,
                  GeneRatio = K6_e@result$GeneRatio,
                  BgRatio = K6_e@result$BgRatio,
                  Count = K6_e@result$Count)
K6d$FE <- as.numeric(sapply(K6_e@result$GeneRatio, function(x) 
  eval(parse(text=x)))/sapply(K6_e@result$BgRatio, 
                              function(x) eval(parse(text=x))))

# Combine to one table
KEGG <- rbind(K1d, K2d, K3d, K4d, K5d, K6d)

# Order KEGG terms by ascending qvalue
KO.1 <- KEGG[KEGG$Cluster == "Cluster 1",]
KO.1 <- KO.1[order(KO.1$qvalue),]
KO.2 <- KEGG[KEGG$Cluster == "Cluster 2",]
KO.2 <- KO.2[order(KO.2$qvalue),]
KO.3 <- KEGG[KEGG$Cluster == "Cluster 3",]
KO.3 <- KO.3[order(KO.3$qvalue),]
KO.4 <- KEGG[KEGG$Cluster == "Cluster 4",]
KO.4 <- KO.4[order(KO.4$qvalue),]
KO.5 <- KEGG[KEGG$Cluster == "Cluster 5",]
KO.5 <- KO.5[order(KO.5$qvalue),]
KO.6 <- KEGG[KEGG$Cluster == "Cluster 6",]
KO.6 <- KO.6[order(KO.6$qvalue),]

# Combine back into one table
KEGG <- rbind(KO.1, KO.2, KO.3, KO.4, KO.5, KO.6)
write.csv(KEGG,"input/KEGG/Combined KEGG.csv", row.names = F)

# Filter by taking the top 5 KEGG pathway terms based on ascending qvalue
filtk1 <- head(KEGG[KEGG$Cluster == "Cluster 1",],5)
filtk2 <- head(KEGG[KEGG$Cluster == "Cluster 2",],5)
filtk3 <- head(KEGG[KEGG$Cluster == "Cluster 3",],5)
filtk4 <- head(KEGG[KEGG$Cluster == "Cluster 4",],5)
filtk5 <- head(KEGG[KEGG$Cluster == "Cluster 5",],5)
filtk6 <- head(KEGG[KEGG$Cluster == "Cluster 6",],5)

# Combine into one final table
KEGG.fin <- rbind(filtk1,filtk2,filtk3,filtk4,filtk5,filtk6)
write.csv(KEGG.fin, "input/KEGG/Top 5 KEGG Terms Per Cluster.csv", 
          row.names = F)

# Format for R in Excel

# Read back formatted data
data <- read.csv("input/KEGG/Top 5 KEGG Terms Per Cluster for R.csv")

# Set factor levels
data$Cluster <- factor(data$Cluster, 
                       levels = c('Cluster 1','Cluster 2','Cluster 3',
                                  'Cluster 4','Cluster 5','Cluster 6'))
BP <- c(rev(as.character(data$Description[c(1:27)])))
data$Description <- factor(data$Description, levels = BP)


# Set up the basics of the ggplot
g <- ggplot(data = data, aes(x = Description, y = FE))
#must use 'stat="identity"' when y-values are used
g <- g + geom_bar(stat = "identity", aes(fill = Description))
#guide legend
g <- g + scale_fill_discrete(guide = FALSE)
#facet by 'Cluster'
g <- g + facet_wrap(~Cluster)
#flip coordinates
g <- g + coord_flip()
#label observed number of genes
g <- g + geom_text(aes(x = Description, y = FE/2, label = Count))
#label axes
g <- g + labs(x = "Biochemical Pathway", y = "Fold Enrichment")
#space axes
g <- g + theme(axis.title.y = element_text(margin = margin(0,17,0,0)))
#plot
g
