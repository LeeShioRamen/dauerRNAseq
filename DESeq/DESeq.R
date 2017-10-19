# Figure 1 Differential Gene Expression Analysis Using DESeq

# Lee, J., Shih, P., Schaedel, O., Rogers, A.K., Quintero-Cadena, P., 
# and Sternberg, P.W. (2017)

# Using the differential gene expression analysis program DESeq, we performed 
# pairwise comparisons between 24 hph L2d and 26 hph L2d to identify gene 
# expression changes during L2d sensory integration; between L2d and 
# dauer-committed larvae for changes during dauer-commitment; and between L2d 
# and dauer larvae for changes during differentiation and maintenance of dauer 
# (Figure 1B). With our reproductive development samples, we performed pairwise 
# comparisons between L2d and L3-committing larvae for changes during commitment 
# to reproductive development, and between L2d and L4 for changes during 
# reproductive growth. In addition, our design allowed us to perform pairwise 
# comparisons between age-matched dauer- and reproductive-developing animals at 
# 26 hph (L2d versus L3-committing larvae) and 34 hph (dauer-committed larvae 
# versus L4) to identify gene expression changes specific to one developmental 
# track (Figure 1C).

# This script follows the DESeq protocol Anders et al. (2013) Nat. Protocols

# This script requires access to the raw sequencing reads, which have been
#deposited in the NCBI Sequence Read Archive (SRA) database (accession number 
#SRP116980) 

# Install the required packages, then comment out these lines using '#'
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq")

# Load the required packages
library("DESeq")

# Download the reference genome and gene model annotations in Unix
#The C. elegans reference genome and gene transfer format files were downloaded 
#from Ensembl release 75 and genome assembly WBcel235

# Here, assess sequence quality with ShortRead

# Collect metadata of experimental design
metadata <- read.csv("input/metadata_for_DESeq.csv")

# Align the reads using tophat2 to the reference genome
gf = "Caenorhabditis_elegans.WBcel235.75.gtf"
bowind = "Cel_WBcel235_75"
cmd = with(metadata,paste("nohup tophat2 -G", gf, "-p 5 -o",
                          LibraryName, bowind, fastq1, fastq2))
cmd #Run in Unix

# Organize, sort and index the BAM files and create SAM files
for(i in seq_len(nrow(metadata))) {
  lib = metadata$LibraryName[i]
  ob = file.path(lib, "accepted_hits.bam")
  #sort by name, convert to SAM for htseq-count
  cat(paste0("nohup samtools sort -n ",ob," ",lib,"_sn"),"\n")
  cat(paste0("nohup samtools view -o ",lib,"_sn.sam ",lib,"_sn.bam"),"\n")
  #sort by position and index for IGV
  cat(paste0("nohup samtools sort ",ob," ",lib,"_s"),"\n")
  cat(paste0("nohup samtools index ",lib,"_s.bam"),"\n\n")
}

# Here, insepct alignments with IGV

# Count reads using HTSEQ-COUNT
metadata$countf = paste(metadata$LibraryName, "count", sep=".")
cmd = paste0("nohup htseq-count -s no -a 10 ", metadata$LibraryName,
             "_sn.sam ", gf," > ",metadata$countf)
cmd #Run in Unix

# DESeq SIMPLE DESIGN
# Create a data frame with the metadata
samplesDESeq=with(metadata,
                  data.frame(shortname = I(shortname), 
                             countf = paste("input/", metadata$countf, sep=""),
                             condition = condition,
                             LibraryLayout = LibraryLayout))

# Create a CountDataSet
cds = newCountDataSetFromHTSeqCount(samplesDESeq)

# Estimate normalization factors
cds = estimateSizeFactors(cds)
sizeFactors(cds)
head(counts(cds, normalized = TRUE))

# Write a table of normalized counts
write.csv(counts(cds, normalized = TRUE), 
          file="output/DESeq_normalizedcounts.csv")

# Inspect sample relationships
cdsB = estimateDispersions(cds,method="blind")
vsd = varianceStabilizingTransformation(cdsB)
p = plotPCA(vsd,intgroup=c("condition","LibraryLayout"))
p

# Estimate dispersion
cds = estimateDispersions(cds)
plotDispEsts(cds)

# Perform test for differential expression
res1=nbinomTest(cds,"noDA24","noDA26")
res2=nbinomTest(cds,"noDA24","noDA34")
res3=nbinomTest(cds,"noDA24","noDA60")
res4=nbinomTest(cds,"noDA26","noDA34")
res5=nbinomTest(cds,"noDA26","noDA60")
res6=nbinomTest(cds,"noDA34","noDA60")
res7=nbinomTest(cds,"DA26","DA34")
res8=nbinomTest(cds,"noDA24","DA26")
res9=nbinomTest(cds,"noDA24","DA34")
res10=nbinomTest(cds,"DA26","noDA26")
res11=nbinomTest(cds,"DA34","noDA34")
res12=nbinomTest(cds,"DA34","noDA60")

# Plot differential expression against expression amount
plotMA(res12,
       col=ifelse(subset(res12,baseMean!=0)$padj>=0.01,"gray32","red3"),
       linecol="#ff000080",
       xlab="", ylab="",
       cex=0.45) #plotMA of res12
title(main="Genes differentially expressed between L4 and dauer", 
      cex.main=2) #title and size of title
mtext("mean of normalized counts",side=1,line=2.7,cex=1.5) #axis title
mtext(expression(log[2]~fold~change),side=2,line=2.2,cex=1.5) #axis title

# Threshold differential expression using padj < 0.01
res1Sig = res1[which(res1$padj < 0.01),]
res2Sig = res2[which(res2$padj < 0.01),]
res3Sig = res3[which(res3$padj < 0.01),]
res4Sig = res4[which(res4$padj < 0.01),]
res5Sig = res5[which(res5$padj < 0.01),]
res6Sig = res6[which(res6$padj < 0.01),]
res7Sig = res7[which(res7$padj < 0.01),]
res8Sig = res8[which(res8$padj < 0.01),]
res9Sig = res9[which(res9$padj < 0.01),]
res10Sig = res10[which(res10$padj < 0.01),]
res11Sig = res11[which(res11$padj < 0.01),]
res12Sig = res12[which(res12$padj < 0.01),]

# Tabulate the number of differentially expressed genes
table(res1$padj < 0.01)
table(res2$padj < 0.01)
table(res3$padj < 0.01)
table(res4$padj < 0.01)
table(res5$padj < 0.01)
table(res6$padj < 0.01)
table(res7$padj < 0.01)
table(res8$padj < 0.01)
table(res9$padj < 0.01)
table(res10$padj < 0.01)
table(res11$padj < 0.01)
table(res12$padj < 0.01)

# Write the results of the negative binomial testing
write.csv(res1,file="output/res1_DESeqDauer.csv")
write.csv(res2,file="output/res2_DESeqDauer.csv")
write.csv(res3,file="output/res3_DESeqDauer.csv")
write.csv(res4,file="output/res4_DESeqDauer.csv")
write.csv(res5,file="output/res5_DESeqDauer.csv")
write.csv(res6,file="output/res6_DESeqDauer.csv")
write.csv(res7,file="output/res7_DESeqDauer.csv")
write.csv(res8,file="output/res8_DESeqDauer.csv")
write.csv(res9,file="output/res9_DESeqDauer.csv")
write.csv(res10,file="output/res10_DESeqDauer.csv")
write.csv(res11,file="output/res11_DESeqDauer.csv")
write.csv(res12,file="output/res12_DESeqDauer.csv")

# Generate histograms of p-values for each comparison
hist(res1$pval,breaks=100)
hist(res2$pval,breaks=100)
hist(res3$pval,breaks=100)
hist(res4$pval,breaks=100)
hist(res5$pval,breaks=100)
hist(res6$pval,breaks=100)
hist(res7$pval,breaks=100)
hist(res8$pval,breaks=100)
hist(res9$pval,breaks=100)
hist(res10$pval,breaks=100)
hist(res11$pval,breaks=100)
hist(res12$pval,breaks=100)
