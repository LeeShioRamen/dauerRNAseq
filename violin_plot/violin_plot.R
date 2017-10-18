#12-05-16 Violin Plot
setwd("~/Desktop/Violin plot/")

#DATA MANIPULATION
#read in the comparisons data
c1<-read.csv("data/comparison files/comparison 1.csv",header=T)
c2<-read.csv("data/comparison files/comparison 2.csv",header=T)
c3<-read.csv("data/comparison files/comparison 3.csv",header=T)
c4<-read.csv("data/comparison files/comparison 4.csv",header=T)
c5<-read.csv("data/comparison files/comparison 5.csv",header=T)
c6<-read.csv("data/comparison files/comparison 6.csv",header=T)
c7<-read.csv("data/comparison files/comparison 7.csv",header=T)
c8<-read.csv("data/comparison files/comparison 8.csv",header=T)
c9<-read.csv("data/comparison files/comparison 9.csv",header=T)
c10<-read.csv("data/comparison files/comparison 10.csv",header=T)
c11<-read.csv("data/comparison files/comparison 11.csv",header=T)
c12<-read.csv("data/comparison files/comparison 12.csv",header=T)

#add a new column for the comparison label--useful later
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

#calculate new log2(fold change) using +1 pseudocounts
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

#combine all tables to one combined table
combine<-rbind(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12)

#remove unnecessary columns
shortcombine<-combine[,c(9,13)]
shortcombine$comparison<-factor(shortcombine$comparison,
                                   levels=c("L2d.24 to L2d.26","L2d.24 to cD","L2d.24 to Dauer","L2d.26 to cD","L2d.26 to Dauer","cD to Dauer",
                                            "L2d.24 to cL3","L2d.24 to L4","cL3 to L4",
                                            "cL3 to L2d.26","L4 to cD","L4 to Dauer"))

#read in intra-sample fold changes (version 2)
library(reshape2)
intra_fc<-read.csv("data/intra_sample_fc.csv",header=T)
melt_intra<-melt(intra_fc,variable.name="comparison",value.name="newL2FC")

#combine inter- and intra-sample fold changes into one table
alldata<-rbind(shortcombine,melt_intra) #factors are already good!

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

#write table to csv
write.csv(alldata,"data/Inter_and_intra_sample_log2fc.csv")

#PLOT VIOLIN PLOT
library(ggplot2)
g<-ggplot(data=alldata,aes(x=comparison,y=newL2FC))
g<-g+geom_point(aes(color=factor(comparison)),shape=1)
g<-g+scale_color_manual(values=c("#FF62BA","#FF65AD","#FF689F","#FF6C90","#FC7180","#F8766F",
                                 "#00BECD","#00B4EF","#29A3FF",
                                 "#9091FF","#B783FF","#D376FF",
                                 "#FAA117","#FAA117","#FAA117","#FAA117",
                                 "#4DD831","#4DD831"),guide=F)
g<-g+geom_violin(scale="width",aes(fill=factor(comparison)))
g<-g+scale_fill_manual(values=c("#FF62BA","#FF65AD","#FF689F","#FF6C90","#FC7180","#F8766F",
                                "#00BECD","#00B4EF","#29A3FF",
                                "#9091FF","#B783FF","#D376FF",
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


#BOXPLOT
alllabs=c("L2d.24 to L2d.26","L2d.24 to cD","L2d.24 to Dauer","L2d.26 to cD","L2d.26 to Dauer","cD to Dauer",
          "L2d.24 to cL3","L2d.24 to L4","cL3 to L4","cL3 to L2d.26","L4 to cD","L4 to Dauer",
          "noDA24","noDA26","noDA34","noDA60","DA26","DA34")
allcols=c("#FF62BA", "#FF65AD", "#FF689F", "#FF6C90", "#FC7180", "#F8766F","#00BECD", "#00B4EF", "#29A3FF","#9091FF", "#B783FF", "#D376FF",
          "#FAA117","#FAA117","#FAA117","#FAA117","#4DD831","#4DD831")
pdf("Boxplot of Fold Changes W Pseudocounts v2.pdf",height=8.5,width=11,colormodel="srgb")
par(mar=c(10,6,4.1,2.1))
boxplot(newL2FC~comparison,data=alldata, col="white",outcol="white",border="white",ylim=c(-20,20),
        xaxt='n',yaxt='n')
axis(1,at=seq(1,18,1),labels=F)
text(seq(1,18,1),par("usr")[3]-2,srt=60,adj=1,
     labels=alllabs,xpd=T,cex=1.5)
axis(side=2,at=seq(-20,20,5),las=2,cex.axis=1.5)
abline(v=seq(1,18,1), lty=2, col="gray")
boxplot(newL2FC~comparison,data=alldata,main="",col=allcols,outcol=allcols,cex.main=2.5,whisklty=1,
        xlab="",ylab=expression('log'[2]*' fold change'),cex.lab=2,xaxt='n',yaxt='n',ylim=c(-20,20),add=T)
dev.off()


#STUDYING THE PLOTS
#Finding the positive and negative modes
#functions from http://stackoverflow.com/questions/16255622/peak-of-the-kernel-density-estimation
#find the dominant mode (mode with the highest peak)
dmode <- function(x) {
  den <- density(x, kernel=c("gaussian"))
  ( den$x[den$y==max(den$y)] )   
}  

#find the number of modes
n.modes <- function(x) {  
  den <- density(x, kernel=c("gaussian"))
  den.s <- smooth.spline(den$x, den$y, all.knots=TRUE, spar=0.8)
  s.0 <- predict(den.s, den.s$x, deriv=0)
  s.1 <- predict(den.s, den.s$x, deriv=1)
  s.derv <- data.frame(s0=s.0$y, s1=s.1$y)
  nmodes <- length(rle(den.sign <- sign(s.derv$s1))$values)/2
  if ((nmodes > 10) == TRUE) { nmodes <- 10 }
  if (is.na(nmodes) == TRUE) { nmodes <- 0 } 
  ( nmodes )
}

#mode,mean,median,and range of Comparison 2 (L2d.24 to cD)
positive_mode<-2^dmode(alldata$newL2FC[alldata$comparison=="L2d.24 to cD" & alldata$newL2FC>0])
negative_mode<-2^dmode(alldata$newL2FC[alldata$comparison=="L2d.24 to cD" & alldata$newL2FC<0])
#(2.63,0.44)

positive_mean<-2^mean(alldata$newL2FC[alldata$comparison=="L2d.24 to cD" & alldata$newL2FC>0])
negative_mean<-2^mean(alldata$newL2FC[alldata$comparison=="L2d.24 to cD" & alldata$newL2FC<0])
#(8.01,0.22)

positive_median<-2^median(alldata$newL2FC[alldata$comparison=="L2d.24 to cD" & alldata$newL2FC>0])
negative_median<-2^median(alldata$newL2FC[alldata$comparison=="L2d.24 to cD" & alldata$newL2FC<0])
#(4.23,0.29)

positive_max<-2^max(alldata$newL2FC[alldata$comparison=="L2d.24 to cD" & alldata$newL2FC>0])
negative_min<-2^min(alldata$newL2FC[alldata$comparison=="L2d.24 to cD" & alldata$newL2FC<0])
#(199201.1,0.002) or (1.99e5,2.06e-3)


#NO PSEUDO-COUNTS (version 3)
npc<-combine[,c(9,6)]
npc$comparison<-factor(npc$comparison,
                       levels=c("L2d.24 to L2d.26","L2d.24 to cD","L2d.24 to Dauer","L2d.26 to cD","L2d.26 to Dauer","cD to Dauer",
                                "L2d.24 to cL3","L2d.24 to L4","cL3 to L4",
                                "cL3 to L2d.26","L4 to cD","L4 to Dauer"))
#read in intra-sample fold changes
npc_intra_fc<-read.csv("data/intra_sample_no_pseudocount_fc.csv",header=T)
melt_npc_intra<-melt(npc_intra_fc,variable.name="comparison",value.name="log2FoldChange")
#combine inter- and intra-sample fold changes into one table
allnpc<-rbind(npc,melt_npc_intra) #factors are already good!
#write table to csv
write.csv(allnpc,"data/Inter_and_intra_sample_no_pseudocount_log2fc.csv")

#PLOT VIOLIN PLOT
#remove Inf/-Inf values
allnpc$log2FoldChange[allnpc$log2FoldChange=="Inf"]<-NA
allnpc$log2FoldChange[allnpc$log2FoldChange=="#NAME?"]<-NA
allnpc$log2FoldChange<-as.numeric(allnpc$log2FoldChange)

library(ggplot2)
g<-ggplot(data=allnpc,aes(x=comparison,y=log2FoldChange))
g<-g+geom_point(aes(color=factor(comparison)),shape=1)
g<-g+scale_color_manual(values=c("#FF62BA","#FF65AD","#FF689F","#FF6C90","#FC7180","#F8766F",
                                 "#00BECD","#00B4EF","#29A3FF",
                                 "#9091FF","#B783FF","#D376FF",
                                 "#FAA117","#FAA117","#FAA117","#FAA117",
                                 "#4DD831","#4DD831"),guide=F)
g<-g+geom_violin(scale="width",aes(fill=factor(comparison)))
g<-g+scale_fill_manual(values=c("#FF62BA","#FF65AD","#FF689F","#FF6C90","#FC7180","#F8766F",
                                "#00BECD","#00B4EF","#29A3FF",
                                "#9091FF","#B783FF","#D376FF",
                                "#FAA117","#FAA117","#FAA117","#FAA117",
                                "#4DD831","#4DD831"),guide=F)
g<-g+scale_y_continuous(limits=c(-22,22),breaks=seq(-20,20,5))
g<-g+labs(y=expression('log'[2]*' fold change'),x="")
g<-g+theme(axis.text.x=element_text(angle=60,hjust=1),
           panel.grid.major.x=element_line(color="gray",linetype="dashed"),
           panel.grid.major.y=element_blank(),
           panel.grid.minor=element_blank(),
           axis.text=element_text(size=18,color="black"),
           axis.title=element_text(size=24),
           panel.background=element_rect(fill="white",color="black"))
g

#STUDYING THE PLOTS
#mode,mean,median,and range of Comparison 2 (L2d.24 to cD)
mode_pos_npc<-2^dmode(allnpc$log2FoldChange[allnpc$comparison=="L2d.24 to cD" & allnpc$log2FoldChange>0 & !is.na(allnpc$log2FoldChange)])
mode_neg_npc<-2^dmode(allnpc$log2FoldChange[allnpc$comparison=="L2d.24 to cD" & allnpc$log2FoldChange<0 & !is.na(allnpc$log2FoldChange)])
#(2.59,0.44)

mean_pos_npc<-2^mean(allnpc$log2FoldChange[allnpc$comparison=="L2d.24 to cD" & allnpc$log2FoldChange>0],na.rm=T)
mean_neg_npc<-2^mean(allnpc$log2FoldChange[allnpc$comparison=="L2d.24 to cD" & allnpc$log2FoldChange<0],na.rm=T)
#(7.22,0.22)

median_pos_npc<-2^median(allnpc$log2FoldChange[allnpc$comparison=="L2d.24 to cD" & allnpc$log2FoldChange>0],na.rm=T)
median_neg_npc<-2^median(allnpc$log2FoldChange[allnpc$comparison=="L2d.24 to cD" & allnpc$log2FoldChange<0],na.rm=T)
#(4.00,0.29)

max_pos_npc<-2^max(allnpc$log2FoldChange[allnpc$comparison=="L2d.24 to cD" & allnpc$log2FoldChange>0],na.rm=T)
min_neg_npc<-2^min(allnpc$log2FoldChange[allnpc$comparison=="L2d.24 to cD" & allnpc$log2FoldChange<0],na.rm=T)
#(515960.2,0.001534065) or (5.16e5,1.53e-3)


#1-18-17 Relabel sample names (v2 relabeled, resulting in v4)
setwd("~/Desktop/Violin plot/")

#DATA MANIPULATION
#read in the comparisons data
c1<-read.csv("data/comparison files/comparison 1.csv",header=T)
c2<-read.csv("data/comparison files/comparison 2.csv",header=T)
c3<-read.csv("data/comparison files/comparison 3.csv",header=T)
c4<-read.csv("data/comparison files/comparison 4.csv",header=T)
c5<-read.csv("data/comparison files/comparison 5.csv",header=T)
c6<-read.csv("data/comparison files/comparison 6.csv",header=T)
c7<-read.csv("data/comparison files/comparison 7.csv",header=T)
c8<-read.csv("data/comparison files/comparison 8.csv",header=T)
c9<-read.csv("data/comparison files/comparison 9.csv",header=T)
c10<-read.csv("data/comparison files/comparison 10.csv",header=T)
c11<-read.csv("data/comparison files/comparison 11.csv",header=T)
c12<-read.csv("data/comparison files/comparison 12.csv",header=T)

#add a new column for the comparison label--useful later
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

#calculate new log2(fold change) using +1 pseudocounts
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

#combine all tables to one combined table
combine<-rbind(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12)

#remove unnecessary columns
shortcombine<-combine[,c(9,13)]
shortcombine$comparison<-factor(shortcombine$comparison,
                                levels=c("L2d.24 to L2d.26","L2d.24 to cD","L2d.24 to Dauer","L2d.26 to cD","L2d.26 to Dauer","cD to Dauer",
                                         "L2d.24 to cL3","L2d.24 to L4","cL3 to L4",
                                         "cL3 to L2d.26","L4 to cD","L4 to Dauer"))

#read in intra-sample fold changes (version 2)
library(reshape2)
intra_fc<-read.csv("data/intra_sample_fc.csv",header=T)
melt_intra<-melt(intra_fc,variable.name="comparison",value.name="newL2FC")

#combine inter- and intra-sample fold changes into one table
alldata<-rbind(shortcombine,melt_intra) #factors are already good!

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

#for some reason, it is necessary to convert to character before replacing values
alldata$comparison <- as.character(alldata$comparison)

alldata$comparison[alldata$comparison=="noDA24"]<-"L2d.24 (24hph)"
alldata$comparison[alldata$comparison=="noDA26"]<-"L2d.26 (26hph)"
alldata$comparison[alldata$comparison=="noDA34"]<-"cD (34hph)"
alldata$comparison[alldata$comparison=="noDA60"]<-"Dauer (60hph)"

alldata$comparison[alldata$comparison=="DA26"]<-"cL3 (26hph)"
alldata$comparison[alldata$comparison=="DA34"]<-"L4 (34hph)"

alldata$comparison <- factor(alldata$comparison,
                             levels=c("L2d.24 to L2d.26","L2d.24 to cD","L2d.24 to Dauer","L2d.26 to cD","L2d.26 to Dauer","cD to Dauer",
                                      "L2d.24 to cL3","L2d.24 to L4","cL3 to L4",
                                      "cL3 to L2d.26","L4 to cD","L4 to Dauer",
                                      "L2d.24 (24hph)","L2d.26 (26hph)","cD (34hph)","Dauer (60hph)",
                                      "cL3 (26hph)","L4 (34hph)"))

#PLOT VIOLIN PLOT
library(ggplot2)
g<-ggplot(data=alldata,aes(x=comparison,y=newL2FC))
g<-g+geom_point(aes(color=factor(comparison)),shape=1)
g<-g+scale_color_manual(values=c("#FF62BA","#FF65AD","#FF689F","#FF6C90","#FC7180","#F8766F",
                                 "#00BECD","#00B4EF","#29A3FF",
                                 "#9091FF","#B783FF","#D376FF",
                                 "#FAA117","#FAA117","#FAA117","#FAA117",
                                 "#4DD831","#4DD831"),guide=F)
g<-g+geom_violin(scale="width",aes(fill=factor(comparison)))
g<-g+scale_fill_manual(values=c("#FF62BA","#FF65AD","#FF689F","#FF6C90","#FC7180","#F8766F",
                                "#00BECD","#00B4EF","#29A3FF",
                                "#9091FF","#B783FF","#D376FF",
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

#STUDYING THE PLOTS
#Finding the positive and negative modes
#functions from http://stackoverflow.com/questions/16255622/peak-of-the-kernel-density-estimation
#find the dominant mode (mode with the highest peak)
dmode <- function(x) {
  den <- density(x, kernel=c("gaussian"))
  ( den$x[den$y==max(den$y)] )   
}

#Comparison statistics
statistics <- data.frame(comparisons=c("L2d.24 to L2d.26","L2d.24 to cD","L2d.24 to Dauer","L2d.26 to cD","L2d.26 to Dauer","cD to Dauer",
                                       "L2d.24 to cL3","L2d.24 to L4","cL3 to L4",
                                       "cL3 to L2d.26","L4 to cD","L4 to Dauer"))
statistics$comparisons <- factor(statistics$comparisons,
                             levels=c("L2d.24 to L2d.26","L2d.24 to cD","L2d.24 to Dauer","L2d.26 to cD","L2d.26 to Dauer","cD to Dauer",
                                      "L2d.24 to cL3","L2d.24 to L4","cL3 to L4",
                                      "cL3 to L2d.26","L4 to cD","L4 to Dauer"))
statistics$UP_mode = c(round(2^dmode(alldata$newL2FC[alldata$comparison=="L2d.24 to L2d.26" & alldata$newL2FC>0]),digits = 2),
                       round(2^dmode(alldata$newL2FC[alldata$comparison=="L2d.24 to cD" & alldata$newL2FC>0]),digits = 2),
                       round(2^dmode(alldata$newL2FC[alldata$comparison=="L2d.24 to Dauer" & alldata$newL2FC>0]),digits = 2),
                       round(2^dmode(alldata$newL2FC[alldata$comparison=="L2d.26 to cD" & alldata$newL2FC>0]),digits = 2),
                       round(2^dmode(alldata$newL2FC[alldata$comparison=="L2d.26 to Dauer" & alldata$newL2FC>0]),digits = 2),
                       round(2^dmode(alldata$newL2FC[alldata$comparison=="cD to Dauer" & alldata$newL2FC>0]),digits = 2),
                       round(2^dmode(alldata$newL2FC[alldata$comparison=="L2d.24 to cL3" & alldata$newL2FC>0]),digits = 2),
                       round(2^dmode(alldata$newL2FC[alldata$comparison=="L2d.24 to L4" & alldata$newL2FC>0]),digits = 2),
                       round(2^dmode(alldata$newL2FC[alldata$comparison=="cL3 to L4" & alldata$newL2FC>0]),digits = 2),
                       round(2^dmode(alldata$newL2FC[alldata$comparison=="cL3 to L2d.26" & alldata$newL2FC>0]),digits = 2),
                       round(2^dmode(alldata$newL2FC[alldata$comparison=="L4 to cD" & alldata$newL2FC>0]),digits = 2),
                       round(2^dmode(alldata$newL2FC[alldata$comparison=="L4 to Dauer" & alldata$newL2FC>0]),digits = 2))

statistics$DOWN_mode = c(round(1/(2^dmode(alldata$newL2FC[alldata$comparison=="L2d.24 to L2d.26" & alldata$newL2FC<0])),digits = 2),
                         round(1/(2^dmode(alldata$newL2FC[alldata$comparison=="L2d.24 to cD" & alldata$newL2FC<0])),digits = 2),
                         round(1/(2^dmode(alldata$newL2FC[alldata$comparison=="L2d.24 to Dauer" & alldata$newL2FC<0])),digits = 2),
                         round(1/(2^dmode(alldata$newL2FC[alldata$comparison=="L2d.26 to cD" & alldata$newL2FC<0])),digits = 2),
                         round(1/(2^dmode(alldata$newL2FC[alldata$comparison=="L2d.26 to Dauer" & alldata$newL2FC<0])),digits = 2),
                         round(1/(2^dmode(alldata$newL2FC[alldata$comparison=="cD to Dauer" & alldata$newL2FC<0])),digits = 2),
                         round(1/(2^dmode(alldata$newL2FC[alldata$comparison=="L2d.24 to cL3" & alldata$newL2FC<0])),digits = 2),
                         round(1/(2^dmode(alldata$newL2FC[alldata$comparison=="L2d.24 to L4" & alldata$newL2FC<0])),digits = 2),
                         round(1/(2^dmode(alldata$newL2FC[alldata$comparison=="cL3 to L4" & alldata$newL2FC<0])),digits = 2),
                         round(1/(2^dmode(alldata$newL2FC[alldata$comparison=="cL3 to L2d.26" & alldata$newL2FC<0])),digits = 2),
                         round(1/(2^dmode(alldata$newL2FC[alldata$comparison=="L4 to cD" & alldata$newL2FC<0])),digits = 2),
                         round(1/(2^dmode(alldata$newL2FC[alldata$comparison=="L4 to Dauer" & alldata$newL2FC<0])),digits = 2))

statistics$UP_mean = c(round(2^mean(alldata$newL2FC[alldata$comparison=="L2d.24 to L2d.26" & alldata$newL2FC>0]),digits = 2),
                       round(2^mean(alldata$newL2FC[alldata$comparison=="L2d.24 to cD" & alldata$newL2FC>0]),digits = 2),
                       round(2^mean(alldata$newL2FC[alldata$comparison=="L2d.24 to Dauer" & alldata$newL2FC>0]),digits = 2),
                       round(2^mean(alldata$newL2FC[alldata$comparison=="L2d.26 to cD" & alldata$newL2FC>0]),digits = 2),
                       round(2^mean(alldata$newL2FC[alldata$comparison=="L2d.26 to Dauer" & alldata$newL2FC>0]),digits = 2),
                       round(2^mean(alldata$newL2FC[alldata$comparison=="cD to Dauer" & alldata$newL2FC>0]),digits = 2),
                       round(2^mean(alldata$newL2FC[alldata$comparison=="L2d.24 to cL3" & alldata$newL2FC>0]),digits = 2),
                       round(2^mean(alldata$newL2FC[alldata$comparison=="L2d.24 to L4" & alldata$newL2FC>0]),digits = 2),
                       round(2^mean(alldata$newL2FC[alldata$comparison=="cL3 to L4" & alldata$newL2FC>0]),digits = 2),
                       round(2^mean(alldata$newL2FC[alldata$comparison=="cL3 to L2d.26" & alldata$newL2FC>0]),digits = 2),
                       round(2^mean(alldata$newL2FC[alldata$comparison=="L4 to cD" & alldata$newL2FC>0]),digits = 2),
                       round(2^mean(alldata$newL2FC[alldata$comparison=="L4 to Dauer" & alldata$newL2FC>0]),digits = 2))

statistics$DOWN_mean = c(round(1/(2^mean(alldata$newL2FC[alldata$comparison=="L2d.24 to L2d.26" & alldata$newL2FC<0])),digits = 2),
                         round(1/(2^mean(alldata$newL2FC[alldata$comparison=="L2d.24 to cD" & alldata$newL2FC<0])),digits = 2),
                         round(1/(2^mean(alldata$newL2FC[alldata$comparison=="L2d.24 to Dauer" & alldata$newL2FC<0])),digits = 2),
                         round(1/(2^mean(alldata$newL2FC[alldata$comparison=="L2d.26 to cD" & alldata$newL2FC<0])),digits = 2),
                         round(1/(2^mean(alldata$newL2FC[alldata$comparison=="L2d.26 to Dauer" & alldata$newL2FC<0])),digits = 2),
                         round(1/(2^mean(alldata$newL2FC[alldata$comparison=="cD to Dauer" & alldata$newL2FC<0])),digits = 2),
                         round(1/(2^mean(alldata$newL2FC[alldata$comparison=="L2d.24 to cL3" & alldata$newL2FC<0])),digits = 2),
                         round(1/(2^mean(alldata$newL2FC[alldata$comparison=="L2d.24 to L4" & alldata$newL2FC<0])),digits = 2),
                         round(1/(2^mean(alldata$newL2FC[alldata$comparison=="cL3 to L4" & alldata$newL2FC<0])),digits = 2),
                         round(1/(2^mean(alldata$newL2FC[alldata$comparison=="cL3 to L2d.26" & alldata$newL2FC<0])),digits = 2),
                         round(1/(2^mean(alldata$newL2FC[alldata$comparison=="L4 to cD" & alldata$newL2FC<0])),digits = 2),
                         round(1/(2^mean(alldata$newL2FC[alldata$comparison=="L4 to Dauer" & alldata$newL2FC<0])),digits = 2))

statistics$UP_median = c(round(2^median(alldata$newL2FC[alldata$comparison=="L2d.24 to L2d.26" & alldata$newL2FC>0]),digits = 2),
                         round(2^median(alldata$newL2FC[alldata$comparison=="L2d.24 to cD" & alldata$newL2FC>0]),digits = 2),
                         round(2^median(alldata$newL2FC[alldata$comparison=="L2d.24 to Dauer" & alldata$newL2FC>0]),digits = 2),
                         round(2^median(alldata$newL2FC[alldata$comparison=="L2d.26 to cD" & alldata$newL2FC>0]),digits = 2),
                         round(2^median(alldata$newL2FC[alldata$comparison=="L2d.26 to Dauer" & alldata$newL2FC>0]),digits = 2),
                         round(2^median(alldata$newL2FC[alldata$comparison=="cD to Dauer" & alldata$newL2FC>0]),digits = 2),
                         round(2^median(alldata$newL2FC[alldata$comparison=="L2d.24 to cL3" & alldata$newL2FC>0]),digits = 2),
                         round(2^median(alldata$newL2FC[alldata$comparison=="L2d.24 to L4" & alldata$newL2FC>0]),digits = 2),
                         round(2^median(alldata$newL2FC[alldata$comparison=="cL3 to L4" & alldata$newL2FC>0]),digits = 2),
                         round(2^median(alldata$newL2FC[alldata$comparison=="cL3 to L2d.26" & alldata$newL2FC>0]),digits = 2),
                         round(2^median(alldata$newL2FC[alldata$comparison=="L4 to cD" & alldata$newL2FC>0]),digits = 2),
                         round(2^median(alldata$newL2FC[alldata$comparison=="L4 to Dauer" & alldata$newL2FC>0]),digits = 2))

statistics$DOWN_median = c(round(1/(2^median(alldata$newL2FC[alldata$comparison=="L2d.24 to L2d.26" & alldata$newL2FC<0])),digits = 2),
                           round(1/(2^median(alldata$newL2FC[alldata$comparison=="L2d.24 to cD" & alldata$newL2FC<0])),digits = 2),
                           round(1/(2^median(alldata$newL2FC[alldata$comparison=="L2d.24 to Dauer" & alldata$newL2FC<0])),digits = 2),
                           round(1/(2^median(alldata$newL2FC[alldata$comparison=="L2d.26 to cD" & alldata$newL2FC<0])),digits = 2),
                           round(1/(2^median(alldata$newL2FC[alldata$comparison=="L2d.26 to Dauer" & alldata$newL2FC<0])),digits = 2),
                           round(1/(2^median(alldata$newL2FC[alldata$comparison=="cD to Dauer" & alldata$newL2FC<0])),digits = 2),
                           round(1/(2^median(alldata$newL2FC[alldata$comparison=="L2d.24 to cL3" & alldata$newL2FC<0])),digits = 2),
                           round(1/(2^median(alldata$newL2FC[alldata$comparison=="L2d.24 to L4" & alldata$newL2FC<0])),digits = 2),
                           round(1/(2^median(alldata$newL2FC[alldata$comparison=="cL3 to L4" & alldata$newL2FC<0])),digits = 2),
                           round(1/(2^median(alldata$newL2FC[alldata$comparison=="cL3 to L2d.26" & alldata$newL2FC<0])),digits = 2),
                           round(1/(2^median(alldata$newL2FC[alldata$comparison=="L4 to cD" & alldata$newL2FC<0])),digits = 2),
                           round(1/(2^median(alldata$newL2FC[alldata$comparison=="L4 to Dauer" & alldata$newL2FC<0])),digits = 2))

statistics$UP_max = c(2^max(alldata$newL2FC[alldata$comparison=="L2d.24 to L2d.26" & alldata$newL2FC>0]),
                      2^max(alldata$newL2FC[alldata$comparison=="L2d.24 to cD" & alldata$newL2FC>0]),
                      2^max(alldata$newL2FC[alldata$comparison=="L2d.24 to Dauer" & alldata$newL2FC>0]),
                      2^max(alldata$newL2FC[alldata$comparison=="L2d.26 to cD" & alldata$newL2FC>0]),
                      2^max(alldata$newL2FC[alldata$comparison=="L2d.26 to Dauer" & alldata$newL2FC>0]),
                      2^max(alldata$newL2FC[alldata$comparison=="cD to Dauer" & alldata$newL2FC>0]),
                      2^max(alldata$newL2FC[alldata$comparison=="L2d.24 to cL3" & alldata$newL2FC>0]),
                      2^max(alldata$newL2FC[alldata$comparison=="L2d.24 to L4" & alldata$newL2FC>0]),
                      2^max(alldata$newL2FC[alldata$comparison=="cL3 to L4" & alldata$newL2FC>0]),
                      2^max(alldata$newL2FC[alldata$comparison=="cL3 to L2d.26" & alldata$newL2FC>0]),
                      2^max(alldata$newL2FC[alldata$comparison=="L4 to cD" & alldata$newL2FC>0]),
                      2^max(alldata$newL2FC[alldata$comparison=="L4 to Dauer" & alldata$newL2FC>0]))

statistics$DOWN_max = c(1/(2^min(alldata$newL2FC[alldata$comparison=="L2d.24 to L2d.26" & alldata$newL2FC<0])),
                        1/(2^min(alldata$newL2FC[alldata$comparison=="L2d.24 to cD" & alldata$newL2FC<0])),
                        1/(2^min(alldata$newL2FC[alldata$comparison=="L2d.24 to Dauer" & alldata$newL2FC<0])),
                        1/(2^min(alldata$newL2FC[alldata$comparison=="L2d.26 to cD" & alldata$newL2FC<0])),
                        1/(2^min(alldata$newL2FC[alldata$comparison=="L2d.26 to Dauer" & alldata$newL2FC<0])),
                        1/(2^min(alldata$newL2FC[alldata$comparison=="cD to Dauer" & alldata$newL2FC<0])),
                        1/(2^min(alldata$newL2FC[alldata$comparison=="L2d.24 to cL3" & alldata$newL2FC<0])),
                        1/(2^min(alldata$newL2FC[alldata$comparison=="L2d.24 to L4" & alldata$newL2FC<0])),
                        1/(2^min(alldata$newL2FC[alldata$comparison=="cL3 to L4" & alldata$newL2FC<0])),
                        1/(2^min(alldata$newL2FC[alldata$comparison=="cL3 to L2d.26" & alldata$newL2FC<0])),
                        1/(2^min(alldata$newL2FC[alldata$comparison=="L4 to cD" & alldata$newL2FC<0])),
                        1/(2^min(alldata$newL2FC[alldata$comparison=="L4 to Dauer" & alldata$newL2FC<0])))

#range of statistic values
rangecomp = data.frame(stat = c("max_UP_mode","min_UP_mode",
                                "max_DOWN_mode","min_DOWN_mode",
                                "max_UP_mean","min_UP_mean",
                                "max_DOWN_mean","min_DOWN_mean",
                                "max_UP_median","min_UP_median",
                                "max_DOWN_median","min_DOWN_median",
                                "max_UP_max","min_UP_max",
                                "max_DOWN_max","min_DOWN_max"))
rangecomp$range = c(max(statistics$UP_mode),min(statistics$UP_mode),max(statistics$DOWN_mode),min(statistics$DOWN_mode),
                    max(statistics$UP_mean),min(statistics$UP_mean),max(statistics$DOWN_mean),min(statistics$DOWN_mean),
                    max(statistics$UP_median),min(statistics$UP_median),max(statistics$DOWN_median),min(statistics$DOWN_median),
                    max(statistics$UP_max),min(statistics$UP_max),max(statistics$DOWN_max),min(statistics$DOWN_max))

write.csv(statistics,"Comparison Statistics.csv")
write.csv(rangecomp,"Range of Comparison Statistics.csv")

#Sample statistics
standarddevs <- data.frame(samples=c("L2d.24 (24hph)","L2d.26 (26hph)","cD (34hph)",
                                     "Dauer (60hph)","cL3 (26hph)","L4 (34hph)"))

standarddevs$samples <- factor(standarddevs$samples,
                               levels=c("L2d.24 (24hph)","L2d.26 (26hph)","cD (34hph)",
                                        "Dauer (60hph)","cL3 (26hph)","L4 (34hph)"))

standarddevs$std_dev = c(round(2^sd(alldata$newL2FC[alldata$comparison=="L2d.24 (24hph)"],na.rm=T),digits = 2),
                         round(2^sd(alldata$newL2FC[alldata$comparison=="L2d.26 (26hph)"],na.rm=T),digits = 2),
                         round(2^sd(alldata$newL2FC[alldata$comparison=="cD (34hph)"],na.rm=T),digits = 2),
                         round(2^sd(alldata$newL2FC[alldata$comparison=="Dauer (60hph)"],na.rm=T),digits = 2),
                         round(2^sd(alldata$newL2FC[alldata$comparison=="cL3 (26hph)"],na.rm=T),digits = 2),
                         round(2^sd(alldata$newL2FC[alldata$comparison=="L4 (34hph)"],na.rm=T),digits = 2))

standarddevs$mean = c(round(2^mean(alldata$newL2FC[alldata$comparison=="L2d.24 (24hph)"],na.rm=T),digits = 2),
                      round(2^mean(alldata$newL2FC[alldata$comparison=="L2d.26 (26hph)"],na.rm=T),digits = 2),
                      round(2^mean(alldata$newL2FC[alldata$comparison=="cD (34hph)"],na.rm=T),digits = 2),
                      round(2^mean(alldata$newL2FC[alldata$comparison=="Dauer (60hph)"],na.rm=T),digits = 2),
                      round(2^mean(alldata$newL2FC[alldata$comparison=="cL3 (26hph)"],na.rm=T),digits = 2),
                      round(2^mean(alldata$newL2FC[alldata$comparison=="L4 (34hph)"],na.rm=T),digits = 2))

standarddevs$detected = c(length(which(!is.na(alldata$newL2FC[alldata$comparison=="L2d.24 (24hph)"]))),
                          length(which(!is.na(alldata$newL2FC[alldata$comparison=="L2d.26 (26hph)"]))),
                          length(which(!is.na(alldata$newL2FC[alldata$comparison=="cD (34hph)"]))),
                          length(which(!is.na(alldata$newL2FC[alldata$comparison=="Dauer (60hph)"]))),
                          length(which(!is.na(alldata$newL2FC[alldata$comparison=="cL3 (26hph)"]))),
                          length(which(!is.na(alldata$newL2FC[alldata$comparison=="L4 (34hph)"]))))

write.csv(standarddevs,"Sample Statistics.csv")
