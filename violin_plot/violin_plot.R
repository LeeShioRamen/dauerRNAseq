# Figure 1D Violin plot

# Lee, J., Shih, P., Schaedel, O., Rogers, A.K., Quintero-Cadena, P., 
# and Sternberg, P.W. (2017)

# Violin plots of the significantly up- and down-regulated genes in 
# each comparison.

# Install the required packages, then comment out these lines using '#'
install.packages("ggplot2")
install.packages("reshape2")

# Load the required packages
library(ggplot2)
library(reshape2)

# Read in the 12 comparisons data
c1<-read.csv("input/comparison 1.csv",header=T)
c2<-read.csv("input/comparison 2.csv",header=T)
c3<-read.csv("input/comparison 3.csv",header=T)
c4<-read.csv("input/comparison 4.csv",header=T)
c5<-read.csv("input/comparison 5.csv",header=T)
c6<-read.csv("input/comparison 6.csv",header=T)
c7<-read.csv("input/comparison 7.csv",header=T)
c8<-read.csv("input/comparison 8.csv",header=T)
c9<-read.csv("input/comparison 9.csv",header=T)
c10<-read.csv("input/comparison 10.csv",header=T)
c11<-read.csv("input/comparison 11.csv",header=T)
c12<-read.csv("input/comparison 12.csv",header=T)

# Add a new column for the comparison label
c1$comparison<-"L2d.24 to L2d.26"
c2$comparison<-"L2d.24 to cD"
c3$comparison<-"L2d.24 to Dauer"
c4$comparison<-"L2d.26 to cD"
c5$comparison<-"L2d.26 to Dauer"
c6$comparison<-"cD to Dauer"
c7$comparison<-"L2d.24 to cL3"
c8$comparison<-"L2d.24 to L4"
c9$comparison<-"cL3 to L4"
c10$comparison<-"cL3 to L2d.26"
c11$comparison<-"L4 to cD"
c12$comparison<-"L4 to Dauer"

# Calculate the log2(fold change) using +1 pseudocounts
c1$Aplus1<-c1$baseMeanA+1
c1$Bplus1<-c1$baseMeanB+1
c1$newFC<-c1$Bplus1/c1$Aplus1
c1$newL2FC<-log2(c1$newFC)

c2$Aplus1<-c2$baseMeanA+1
c2$Bplus1<-c2$baseMeanB+1
c2$newFC<-c2$Bplus1/c2$Aplus1
c2$newL2FC<-log2(c2$newFC)

c3$Aplus1<-c3$baseMeanA+1
c3$Bplus1<-c3$baseMeanB+1
c3$newFC<-c3$Bplus1/c3$Aplus1
c3$newL2FC<-log2(c3$newFC)

c4$Aplus1<-c4$baseMeanA+1
c4$Bplus1<-c4$baseMeanB+1
c4$newFC<-c4$Bplus1/c4$Aplus1
c4$newL2FC<-log2(c4$newFC)

c5$Aplus1<-c5$baseMeanA+1
c5$Bplus1<-c5$baseMeanB+1
c5$newFC<-c5$Bplus1/c5$Aplus1
c5$newL2FC<-log2(c5$newFC)

c6$Aplus1<-c6$baseMeanA+1
c6$Bplus1<-c6$baseMeanB+1
c6$newFC<-c6$Bplus1/c6$Aplus1
c6$newL2FC<-log2(c6$newFC)

c7$Aplus1<-c7$baseMeanA+1
c7$Bplus1<-c7$baseMeanB+1
c7$newFC<-c7$Bplus1/c7$Aplus1
c7$newL2FC<-log2(c7$newFC)

c8$Aplus1<-c8$baseMeanA+1
c8$Bplus1<-c8$baseMeanB+1
c8$newFC<-c8$Bplus1/c8$Aplus1
c8$newL2FC<-log2(c8$newFC)

c9$Aplus1<-c9$baseMeanA+1
c9$Bplus1<-c9$baseMeanB+1
c9$newFC<-c9$Bplus1/c9$Aplus1
c9$newL2FC<-log2(c9$newFC)

c10$Aplus1<-c10$baseMeanA+1
c10$Bplus1<-c10$baseMeanB+1
c10$newFC<-c10$Bplus1/c10$Aplus1
c10$newL2FC<-log2(c10$newFC)

c11$Aplus1<-c11$baseMeanA+1
c11$Bplus1<-c11$baseMeanB+1
c11$newFC<-c11$Bplus1/c11$Aplus1
c11$newL2FC<-log2(c11$newFC)

c12$Aplus1<-c12$baseMeanA+1
c12$Bplus1<-c12$baseMeanB+1
c12$newFC<-c12$Bplus1/c12$Aplus1
c12$newL2FC<-log2(c12$newFC)

# Combine all tables to one
combine <- rbind(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12)

# Remove unnecessary columns
combine<-combine[,c(9,13)]
combine$comparison<-factor(combine$comparison,
                           levels=c("L2d.24 to L2d.26","L2d.24 to cD",
                                    "L2d.24 to Dauer","L2d.26 to cD",
                                    "L2d.26 to Dauer","cD to Dauer",
                                    "L2d.24 to cL3","L2d.24 to L4","cL3 to L4",
                                    "cL3 to L2d.26","L4 to cD","L4 to Dauer"))

# Read in the fold changes between the replicates of each sequenced time point
intra_fc <- read.csv("input/intra_sample_fc.csv", header=T)
melt_intra <- melt(intra_fc, variable.name="comparison", value.name="newL2FC")

# Combine all fold changes into one table
alldata <- rbind(combine,melt_intra)

# Set violin plot colors for each comparison
alldata$color<-"#FF62BA"
alldata$color[alldata$comparison=="L2d.24 to cD"]<-"#FF65AD"
alldata$color[alldata$comparison=="L2d.24 to Dauer"]<-"#FF689F"
alldata$color[alldata$comparison=="L2d.26 to cD"]<-"#FF6C90"
alldata$color[alldata$comparison=="L2d.26 to Dauer"]<-"#FC7180"
alldata$color[alldata$comparison=="cD to Dauer"]<-"#F8766F"
alldata$color[alldata$comparison=="L2d.24 to cL3"]<-"#00BECD"
alldata$color[alldata$comparison=="L2d.24 to L4"]<-"#00B4EF"
alldata$color[alldata$comparison=="cL3 to L4"]<-"#29A3FF"
alldata$color[alldata$comparison=="cL3 to L2d.26"]<-"#9091FF"
alldata$color[alldata$comparison=="L4 to cD"]<-"#B783FF"
alldata$color[alldata$comparison=="L4 to Dauer"]<-"#D376FF"
alldata$color[alldata$comparison=="noDA24"]<-"#FAA117"
alldata$color[alldata$comparison=="noDA26"]<-"#FAA117"
alldata$color[alldata$comparison=="noDA34"]<-"#FAA117"
alldata$color[alldata$comparison=="noDA60"]<-"#FAA117"
alldata$color[alldata$comparison=="DA26"]<-"#4DD831"
alldata$color[alldata$comparison=="DA34"]<-"#4DD831"

# Re-label the sequenced time points
alldata$comparison <- as.character(alldata$comparison)
alldata$comparison[alldata$comparison=="noDA24"]<-"L2d.24 (24hph)"
alldata$comparison[alldata$comparison=="noDA26"]<-"L2d.26 (26hph)"
alldata$comparison[alldata$comparison=="noDA34"]<-"cD (34hph)"
alldata$comparison[alldata$comparison=="noDA60"]<-"Dauer (60hph)"
alldata$comparison[alldata$comparison=="DA26"]<-"cL3 (26hph)"
alldata$comparison[alldata$comparison=="DA34"]<-"L4 (34hph)"

# Set factors
alldata$comparison <- factor(alldata$comparison,
                             levels=c("L2d.24 to L2d.26","L2d.24 to cD",
                                      "L2d.24 to Dauer","L2d.26 to cD",
                                      "L2d.26 to Dauer","cD to Dauer",
                                      "L2d.24 to cL3","L2d.24 to L4",
                                      "cL3 to L4",
                                      "cL3 to L2d.26","L4 to cD","L4 to Dauer",
                                      "L2d.24 (24hph)","L2d.26 (26hph)",
                                      "cD (34hph)","Dauer (60hph)",
                                      "cL3 (26hph)","L4 (34hph)"))

# Create the violin plot
g<-ggplot(data=alldata,aes(x=comparison,y=newL2FC))
g<-g+geom_point(aes(color=factor(comparison)),shape=1)
g<-g+scale_color_manual(values=c("#FF62BA","#FF65AD","#FF689F","#FF6C90",
                                 "#FC7180","#F8766F","#00BECD","#00B4EF",
                                 "#29A3FF","#9091FF","#B783FF","#D376FF",
                                 "#FAA117","#FAA117","#FAA117","#FAA117",
                                 "#4DD831","#4DD831"),guide=F)
g<-g+geom_violin(scale="width",aes(fill=factor(comparison)))
g<-g+scale_fill_manual(values=c("#FF62BA","#FF65AD","#FF689F","#FF6C90",
                                "#FC7180","#F8766F","#00BECD","#00B4EF",
                                "#29A3FF","#9091FF","#B783FF","#D376FF",
                                "#FAA117","#FAA117","#FAA117","#FAA117",
                                "#4DD831","#4DD831"),guide=F)
g<-g+scale_y_continuous(limits=c(-20,20),breaks=seq(-20,20,5))
g<-g+labs(y=expression('log'[2]*' fold change'),x="")
g<-g+theme(axis.text.x=element_text(angle=60,hjust=1),
           panel.grid.major.x=element_line(color="gray",linetype="dashed"),
           panel.grid.major.y=element_blank(),
           panel.grid.minor=element_blank(),
           axis.text=element_text(size=18,color="black"),
           axis.title=element_text(size=24),
           panel.background=element_rect(fill="white",color="black"))
g
