#Community Assembly Analyses on FINAL filtered dataset - Lindsay I. Putman


#upload necessary packages
library (phyloseq)
library('ggplot2')
library(plyr)
library(ggrepel)
library(grid)
library(gridExtra)
library(vegan)
library(ggdendro)
library(picante)
library(mcheatmaps)
library(ape)
library(phylotools)
library(wesanderson)
library(RColorBrewer)
library(edgeR)
library(tidyr)
library(DESeq2)
library(Biostrings)
library(reshape2)
library(ecodist)
library(scales)
theme_set(theme_bw())

#Upload OTU Id information 
CROMO_OTU_IDs <- read.delim("CROMO_OTU_IDs.txt", header=FALSE)
otu.ids=CROMO_OTU_IDs

#Upload Count Table Data, turn into matrix, add OTU Ids
library(readxl)
CROMO_Filtered_FINAL_counts <- read_excel("CROMO_Filtered_FINAL_counts.xlsx")
count.table=as.matrix(CROMO_Filtered_FINAL_counts)
rownames(count.table)=otu.ids[,1]

#Upload Meta data, append row names = count table sample names
CROMO_Meta.data <- read_excel("CROMO_Meta.data.xlsx")
meta.data=as.data.frame(CROMO_Meta.data)
rownames(meta.data)=colnames(count.table)

#Upload taxonomy data, combine blast match percentage with species name, turn into matrix, append row names = otu ids
CROMO_Filtered_FINAL_taxonomy <- read_excel("CROMO_Filtered_FINAL_taxonomy.xlsx")
tax.table=as.data.frame(CROMO_Filtered_FINAL_taxonomy)
rownames(tax.table)=otu.ids[,1]
tax.table=unite(tax.table, Species, Species:BLAST_Percent_Identity, sep=" (", remove=TRUE)
tax.table=as.matrix(tax.table)

#Create phyloseq object --REMOVE REFSEQ COMMAND
merged <- phyloseq(otu_table(count.table, taxa_are_rows=TRUE), tax_table(tax.table), sample_data(meta.data))

#Upload Raup Crick (RCabund) Table
RCbray_matrix = read_excel("CROMO_RCbray_matrix.xlsx")
raup.mat=as.data.frame(RCbray_matrix)
rownames(raup.mat)=raup.mat[,1]
raup.mat[,1]<-NULL
all.raup.dist=as.dist(raup.mat)


#Define color palette
cols=brewer.pal(9, "Greys")
pal=colorRampPalette(cols)

cols2=brewer.pal(11, "RdYlBu")
pal2=colorRampPalette(cols2)

#Function to extract legends from plots
#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#Define x axis for boxplots
all.x.axis.ph=c("CSW1.4", "N08C", "CSW1.2", "QV1.2", "QV1.3", "CSW1.5", "CSWold", "CSW1.3", 
              "N08A", "N08B", "QV1.1", "CSW1.1")


#Upload Raup-Crick boxplot files
RCbray_Within_Well_Comps_ONLY <- read_excel("RCbray_Within_Well_Comps_ONLY.xlsx")
rc.w.df=as.data.frame(RCbray_Within_Well_Comps_ONLY)

RCbray_Btwn_Well_Comps_ONLY <- read_excel("RCbray_Btwn_Well_Comps_ONLY.xlsx")
rc.b.df=as.data.frame(RCbray_Btwn_Well_Comps_ONLY)

#Make Within Wells Boxplot
rc.w.boxplot= ggplot(rc.w.df, aes(x=Well, y=RC_Index, fill=Well)) + geom_boxplot()
rc.w.boxplot= rc.w.boxplot + ggtitle("RCbray Index Values Within Wells") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
rc.w.boxplot= rc.w.boxplot + geom_hline(yintercept=0, linetype="dashed", colour="black") + geom_hline(yintercept=0.95, linetype="dashed", colour="red") + geom_hline(yintercept=-0.95, linetype="dashed", colour="red")
rc.w.boxplot= rc.w.boxplot + theme(plot.title = element_text(hjust=0.5)) + scale_fill_manual(values=pal2(12), limits = all.x.axis.ph) + theme(axis.text.x = element_text(angle = 90))
rc.w.boxplot = rc.w.boxplot + labs(x="", y="") + scale_x_discrete(limits=all.x.axis.ph) + theme(legend.position="none")
rc.w.boxplot #Mean RCbray = -0.2046185

#Make Between Wells Boxplot
rc.b.boxplot= ggplot(rc.b.df, aes(x=Well, y=RC_Index, fill=Well)) + geom_boxplot()
rc.b.boxplot= rc.b.boxplot + ggtitle("RCbray Index Values Between Wells") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
rc.b.boxplot= rc.b.boxplot + geom_hline(yintercept=0, linetype="dashed", colour="black") + geom_hline(yintercept=0.95, linetype="dashed", colour="red") + geom_hline(yintercept=-0.95, linetype="dashed", colour="red")
rc.b.boxplot= rc.b.boxplot + theme(plot.title = element_text(hjust=0.5)) + scale_fill_manual(values=pal2(12), limits = all.x.axis.ph) + theme(axis.text.x = element_text(angle = 90))
rc.b.boxplot = rc.b.boxplot + labs(x="", y="") + scale_x_discrete(limits=all.x.axis.ph) + theme(legend.position="none")
rc.b.boxplot #Mean RCbray = 0.2950858

#Uploading NTI data
library(readr)
CROMO_NTI_Results <- read_csv("CROMO_NTI_Results.csv")
nti.dat=as.data.frame(CROMO_NTI_Results)
#Mean NTI is 1.73

#Make boxplot of NTI dat
all.nti.boxplot= ggplot(nti.dat, aes(x=Well, y=NTI, fill=Well)) + geom_boxplot()
all.nti.boxplot= all.nti.boxplot + ggtitle("Boxplot of Nearest Taxon Index Within Individual Wells 2011-2017") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
all.nti.boxplot= all.nti.boxplot + geom_hline(yintercept=0, linetype="dashed", colour="black") + geom_hline(yintercept=2, linetype="dashed", colour="red")
all.nti.boxplot= all.nti.boxplot + theme(plot.title = element_text(hjust=0.5))  + scale_fill_manual(values=pal2(12), limits=all.x.axis.ph)
all.nti.boxplot = all.nti.boxplot + ylab("Nearest Taxon Index")+scale_x_discrete(limits=all.x.axis.ph)
all.nti.boxplot 
#alphaNTI > 0 indicates phylogenetic clustering


#Upload BNTI data
CROMO_weighted_bNTI_matrix <- read_csv("CROMO_weighted_bNTI_matrix.csv")
bnti.mat=as.data.frame(CROMO_weighted_bNTI_matrix)
row.names(bnti.mat)=bnti.mat[,1]
bnti.mat[,1]<-NULL
bnti.dist=as.dist(bnti.mat, diag=TRUE, upper=FALSE)
bnti.dist2=as.dist(bnti.mat, diag=TRUE, upper=TRUE)
bnti.mat=as.matrix(bnti.dist2)
bnti.mat=as.data.frame(bnti.mat)

#Subset BNTI Matrix
#########################################
CSW1.1.bnti=bnti.mat[1:11, 1:11]
CSW1.1.bnti=as.dist(CSW1.1.bnti)
CSW1.2.bnti=bnti.mat[12:16, 12:16]
CSW1.2.bnti=as.dist(CSW1.2.bnti)
CSW1.3.bnti=bnti.mat[17:22, 17:22]
CSW1.3.bnti=as.dist(CSW1.3.bnti)
CSW1.4.bnti=bnti.mat[23:29, 23:29]
CSW1.4.bnti=as.dist(CSW1.4.bnti)
CSW1.5.bnti=bnti.mat[30:36, 30:36]
CSW1.5.bnti=as.dist(CSW1.5.bnti)
CSWold.bnti=bnti.mat[37:49, 37:49]
CSWold.bnti=as.dist(CSWold.bnti)
N08A.bnti=bnti.mat[50:60, 50:60]
N08A.bnti=as.dist(N08A.bnti)
N08B.bnti=bnti.mat[61:70, 61:70]
N08B.bnti=as.dist(N08B.bnti)
N08C.bnti=bnti.mat[71:79, 71:79]
N08C.bnti=as.dist(N08C.bnti)
QV1.1.bnti=bnti.mat[80:91, 80:91]
QV1.1.bnti=as.dist(QV1.1.bnti)
QV1.2.bnti=bnti.mat[92:98, 92:98]
QV1.2.bnti=as.dist(QV1.2.bnti)
QV1.3.bnti=bnti.mat[99:104, 99:104]
QV1.3.bnti=as.dist(QV1.3.bnti)

#Create Vectors and Dataframes for each
all.bnti.vect=as.vector(bnti.dist)
all.bnti.df=data.frame(Well="All Wells", RC_Index=all.bnti.vect)
CSW1.1.bnti.vect=as.vector(CSW1.1.bnti)
CSW1.1.bnti.df=data.frame(Well="CSW1.1", RC_Index=CSW1.1.bnti.vect)
CSW1.2.bnti.vect=as.vector(CSW1.2.bnti)
CSW1.2.bnti.df=data.frame(Well="CSW1.2", RC_Index=CSW1.2.bnti.vect)
CSW1.3.bnti.vect=as.vector(CSW1.3.bnti)
CSW1.3.bnti.df=data.frame(Well="CSW1.3", RC_Index=CSW1.3.bnti.vect)
CSW1.4.bnti.vect=as.vector(CSW1.4.bnti)
CSW1.4.bnti.df=data.frame(Well="CSW1.4", RC_Index=CSW1.4.bnti.vect)
CSW1.5.bnti.vect=as.vector(CSW1.5.bnti)
CSW1.5.bnti.df=data.frame(Well="CSW1.5", RC_Index=CSW1.5.bnti.vect)
CSWold.bnti.vect=as.vector(CSWold.bnti)
CSWold.bnti.df=data.frame(Well="CSWold", RC_Index=CSWold.bnti.vect)
N08A.bnti.vect=as.vector(N08A.bnti)
N08A.bnti.df=data.frame(Well="N08A", RC_Index=N08A.bnti.vect)
N08B.bnti.vect=as.vector(N08B.bnti)
N08B.bnti.df=data.frame(Well="N08B", RC_Index=N08B.bnti.vect)
N08C.bnti.vect=as.vector(N08C.bnti)
N08C.bnti.df=data.frame(Well="N08C", RC_Index=N08C.bnti.vect)
QV1.1.bnti.vect=as.vector(QV1.1.bnti)
QV1.1.bnti.df=data.frame(Well="QV1.1", RC_Index=QV1.1.bnti.vect)
QV1.2.bnti.vect=as.vector(QV1.2.bnti)
QV1.2.bnti.df=data.frame(Well="QV1.2", RC_Index=QV1.2.bnti.vect)
QV1.3.bnti.vect=as.vector(QV1.3.bnti)
QV1.3.bnti.df=data.frame(Well="QV1.3", RC_Index=QV1.3.bnti.vect)
##############################################

#All wells Dataframe
all.wells.bnti.df=rbind(CSW1.4.bnti.df, N08C.bnti.df, QV1.2.bnti.df, CSW1.2.bnti.df, CSW1.1.bnti.df, 
                        QV1.1.bnti.df, CSW1.3.bnti.df, N08B.bnti.df, CSW1.5.bnti.df, QV1.3.bnti.df, 
                        N08A.bnti.df, CSWold.bnti.df)

#Upload Between Wells Boxplot Dataframe
BNTI_Btwn_Well_Comps_Only <- read_excel("BNTI_Btwn_Well_Comps_Only.xlsx")
bnti.b.df=as.data.frame(BNTI_Btwn_Well_Comps_Only)


#Make Within Wells Boxplot
all.bnti.boxplot= ggplot(all.wells.bnti.df, aes(x=Well, y=RC_Index, fill=Well)) + geom_boxplot()
all.bnti.boxplot= all.bnti.boxplot + ggtitle("BNTI Values Within Wells") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
all.bnti.boxplot= all.bnti.boxplot + geom_hline(yintercept=0, linetype="dashed", colour="black") + geom_hline(yintercept=2, linetype="dashed", colour="red") + geom_hline(yintercept=-2, linetype="dashed", colour="red")
all.bnti.boxplot= all.bnti.boxplot + theme(plot.title = element_text(hjust=0.5)) + scale_fill_manual(values=pal2(12), limits = all.x.axis.ph) + theme(axis.text.x = element_text(angle = 90))
all.bnti.boxplot = all.bnti.boxplot + labs(x="", y="") + scale_x_discrete(limits=all.x.axis.ph) + theme(legend.position = "bottom") + guides(col=guide_legend(nrow=1), byrow=TRUE)
all.bnti.boxplot #Mean BNTI = -1.55921

#BNTI Between Wells Boxplot
bnti.b.boxplot= ggplot(bnti.b.df, aes(x=Well, y=RC_Index, fill=Well)) + geom_boxplot()
bnti.b.boxplot= bnti.b.boxplot + ggtitle("BNTI Values Between Wells") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
bnti.b.boxplot= bnti.b.boxplot + geom_hline(yintercept=0, linetype="dashed", colour="black") + geom_hline(yintercept=2, linetype="dashed", colour="red") + geom_hline(yintercept=-2, linetype="dashed", colour="red")
bnti.b.boxplot= bnti.b.boxplot + theme(plot.title = element_text(hjust=0.5)) + scale_fill_manual(values=pal2(12), limits = all.x.axis.ph) + theme(axis.text.x = element_text(angle = 90))
bnti.b.boxplot = bnti.b.boxplot + labs(x="", y="") + scale_x_discrete(limits=all.x.axis.ph) + theme(legend.position = "bottom") + guides(col=guide_legend(nrow=1))
bnti.b.boxplot #Mean BNTI = -0.1857504


#Combine Boxplots together
eco.boxplots=grid.arrange(all.bnti.boxplot + theme(legend.position = "none"), rc.w.boxplot, bnti.b.boxplot + theme(legend.position = "none"), rc.b.boxplot, nrow=2, ncol=2)


#Calculate PERMANOVA and Make Bray-Curtis NMDS Plot
############################################################
merged.props=transform_sample_counts(merged, function(x) 100 * x/sum(x))

#PERMANOVA on Community Data
bray.dist=phyloseq::distance(merged.props, method="bray")
permanova=adonis(bray.dist~ Well + Days + Well_Depth_m + Temp + pH + Cond_mS 
                 + DO_mgL + ORP_mV + DIC_uM, data=meta.data)
permanova

#Calculate ordination
nmds.bray = ordinate(merged.props, method='NMDS', distance='bray') #Solution reached, stress=0.236


#Correlate with environmental data
nmds.fit.bray = envfit(nmds.bray ~ Well + Days + Well_Depth_m + Temp
                       + pH + Cond_mS + DO_mgL + ORP_mV + DIC_uM
                       , as.data.frame(meta.data),                                   
                       perm=1000, na.rm=TRUE)

#Fit vectors for the environmental variables
nmds.fit.bray.vectors = scores(nmds.fit.bray, display=c('vectors'))

#Scale significant vectors based on correlation coefficients
multiplier.bray=ordiArrowMul(nmds.fit.bray.vectors, fill=0.25)

nmds.fit.bray.vectors = data.frame(labels=rownames(nmds.fit.bray.vectors), nmds.fit.bray.vectors*multiplier.bray)

#Plot NMDS - Bray curtis - Sig. Vectors
p1.nmds=plot_ordination(merged.props, nmds.bray, type='samples', color='Well')
p1.nmds
#Adjust Scale, Colors and Shape Sizes
##Still need to figure out scale issues
p1.nmds=p1.nmds + geom_point(size=3) + scale_color_manual(values=pal2(12), limits=all.x.axis.ph)
p1.nmds=p1.nmds + ggtitle("NMDS Plot of Bray-Curtis Dissimiliarites with Significant Environmental Variables") + theme(plot.title=element_text(hjust=0.8), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1.nmds
#Add Environmental Variables
p1.nmds = p1.nmds + geom_segment(data=nmds.fit.bray.vectors, aes(x=0, y=0, xend=NMDS1, yend=NMDS2), inherit.aes=FALSE, arrow=arrow(length=unit(0.1, 'cm')), alpha=0.8)
p1.nmds
#Add Vector Labels
p1.nmds = p1.nmds + geom_text_repel(data=rbind(nmds.fit.bray.vectors), aes(NMDS1, NMDS2, label=labels), inherit.aes=FALSE, max.iter=5000, color="black", alpha=0.8)
p1.nmds

#Vectors plotted to get accurate vector direction on plot. Manually removed insignificant variables based on 
#PERMANOVA and manually re-scaled vector length based on PERMANOVA R2 values

#############################################################################################


#Code from Danczak et al. 2018 below#

### Generating percentages from bNTI and Raup-Crick data ###  DATA THAT IS SHOWN IN TABLE S5
############################################################
## Assigning a percent explained by selection/neutral processes to each sample

# Data sorting and whatnot
x = bnti.mat
y = raup.mat

# Assigning names to begin calculating percentages
model = matrix(nrow = length(x[,1]), ncol = length(x[1,]), data = as.character(NA))

model[x > 2] = "Variable Selection"
model[x < -2] = "Homogenizing Selection"
model[y < -0.95 & x > -2 & x < 2] = "Homogenizing Dispersal"
model[y > 0.95 & x > -2 & x < 2] = "Dispersal Limitation"
model[y <0.95 & y > -0.95 & x > -2 & x < 2] = "Undominated"

dimnames(model) = dimnames(x)

# Compiling model results for all comparisons
model.results = NULL
model.results = cbind(model.results, apply(model, 1, function(x) length(which(x %in% "Variable Selection"))/(length(x)-length(which(is.na(x))))))
model.results = cbind(model.results, apply(model, 1, function(x) length(which(x %in% "Homogenizing Selection"))/(length(x)-length(which(is.na(x))))))
model.results = cbind(model.results, apply(model, 1, function(x) length(which(x %in% "Homogenizing Dispersal"))/(length(x)-length(which(is.na(x))))))
model.results = cbind(model.results, apply(model, 1, function(x) length(which(x %in% "Dispersal Limitation"))/(length(x)-length(which(is.na(x))))))
model.results = cbind(model.results, apply(model, 1, function(x) length(which(x %in% "Undominated"))/(length(x)-length(which(is.na(x))))))

row.names(model.results) = row.names(model)
colnames(model.results) = c("Variable Selection", "Homogenizing Selection", "Homogenizing Dispersal", "Dispersal Limitation", "Undominated")

write.table(model.results, "eco.model.fraction.allwells.1117.results.txt", sep="\t", row.names = TRUE, col.names = TRUE, na="NA")

#Repeat for CSW Cluster
x = bnti.mat[1:49, 1:49]
y = raup.mat[1:49, 1:49]

model.csw = matrix(nrow = length(x[,1]), ncol = length(x[1,]), data = as.character(NA))

model.csw[x > 2] = "Variable Selection"
model.csw[x < -2] = "Homogenizing Selection"
model.csw[y < -0.95 & x > -2 & x < 2] = "Homogenizing Dispersal"
model.csw[y > 0.95 & x > -2 & x < 2] = "Dispersal Limitation"
model.csw[y <0.95 & y > -0.95 & x > -2 & x < 2] = "Undominated"

dimnames(model.csw) = dimnames(x)

# Compiling model.csw results for all comparisons
model.csw.results = NULL
model.csw.results = cbind(model.csw.results, apply(model.csw, 1, function(x) length(which(x %in% "Variable Selection"))/(length(x)-length(which(is.na(x))))))
model.csw.results = cbind(model.csw.results, apply(model.csw, 1, function(x) length(which(x %in% "Homogenizing Selection"))/(length(x)-length(which(is.na(x))))))
model.csw.results = cbind(model.csw.results, apply(model.csw, 1, function(x) length(which(x %in% "Homogenizing Dispersal"))/(length(x)-length(which(is.na(x))))))
model.csw.results = cbind(model.csw.results, apply(model.csw, 1, function(x) length(which(x %in% "Dispersal Limitation"))/(length(x)-length(which(is.na(x))))))
model.csw.results = cbind(model.csw.results, apply(model.csw, 1, function(x) length(which(x %in% "Undominated"))/(length(x)-length(which(is.na(x))))))

row.names(model.csw.results) = row.names(model.csw)
colnames(model.csw.results) = c("Variable Selection", "Homogenizing Selection", "Homogenizing Dispersal", "Dispersal Limitation", "Undominated")

write.table(model.csw.results, "eco.model.csw.fraction.11-17.results.txt", sep="\t", row.names = TRUE, col.names = TRUE, na="NA")

#Repeat for QV Cluster
x = bnti.mat[50:104, 50:104]
y = raup.mat[50:104, 50:104]

model.qv = matrix(nrow = length(x[,1]), ncol = length(x[1,]), data = as.character(NA))

model.qv[x > 2] = "Variable Selection"
model.qv[x < -2] = "Homogenizing Selection"
model.qv[y < -0.95 & x > -2 & x < 2] = "Homogenizing Dispersal"
model.qv[y > 0.95 & x > -2 & x < 2] = "Dispersal Limitation"
model.qv[y <0.95 & y > -0.95 & x > -2 & x < 2] = "Undominated"

dimnames(model.qv) = dimnames(x)

# Compiling model.qv results for all comparisons
model.qv.results = NULL
model.qv.results = cbind(model.qv.results, apply(model.qv, 1, function(x) length(which(x %in% "Variable Selection"))/(length(x)-length(which(is.na(x))))))
model.qv.results = cbind(model.qv.results, apply(model.qv, 1, function(x) length(which(x %in% "Homogenizing Selection"))/(length(x)-length(which(is.na(x))))))
model.qv.results = cbind(model.qv.results, apply(model.qv, 1, function(x) length(which(x %in% "Homogenizing Dispersal"))/(length(x)-length(which(is.na(x))))))
model.qv.results = cbind(model.qv.results, apply(model.qv, 1, function(x) length(which(x %in% "Dispersal Limitation"))/(length(x)-length(which(is.na(x))))))
model.qv.results = cbind(model.qv.results, apply(model.qv, 1, function(x) length(which(x %in% "Undominated"))/(length(x)-length(which(is.na(x))))))

row.names(model.qv.results) = row.names(model.qv)
colnames(model.qv.results) = c("Variable Selection", "Homogenizing Selection", "Homogenizing Dispersal", "Dispersal Limitation", "Undominated")

write.table(model.qv.results, "eco.model.qv.fraction.11-17.results.txt", sep="\t", row.names = TRUE, col.names = TRUE, na="NA")

#Neutral pH wells
x = bnti.mat[-c(1:11), -c(1:11)]
x=x[-c(6:11), -c(6:11)]
x=x[-c(13:53), -c(13:53)]
x=x[-c(22:33), -c(22:33)]
x=x[-c(29:34), -c(29:34)]

y = raup.mat[-c(1:11), -c(1:11)]
y=y[-c(6:11), -c(6:11)]
y=y[-c(13:53), -c(13:53)]
y=y[-c(22:33), -c(22:33)]
y=y[-c(29:34), -c(29:34)]

model.neut = matrix(nrow = length(x[,1]), ncol = length(x[1,]), data = as.character(NA))

model.neut[x > 2] = "Variable Selection"
model.neut[x < -2] = "Homogenizing Selection"
model.neut[y < -0.95 & x > -2 & x < 2] = "Homogenizing Dispersal"
model.neut[y > 0.95 & x > -2 & x < 2] = "Dispersal Limitation"
model.neut[y <0.95 & y > -0.95 & x > -2 & x < 2] = "Undominated"

dimnames(model.neut) = dimnames(x)

# Compiling model.sh results for all comparisons
model.neut.results = NULL
model.neut.results = cbind(model.neut.results, apply(model.neut, 1, function(x) length(which(x %in% "Variable Selection"))/(length(x)-length(which(is.na(x))))))
model.neut.results = cbind(model.neut.results, apply(model.neut, 1, function(x) length(which(x %in% "Homogenizing Selection"))/(length(x)-length(which(is.na(x))))))
model.neut.results = cbind(model.neut.results, apply(model.neut, 1, function(x) length(which(x %in% "Homogenizing Dispersal"))/(length(x)-length(which(is.na(x))))))
model.neut.results = cbind(model.neut.results, apply(model.neut, 1, function(x) length(which(x %in% "Dispersal Limitation"))/(length(x)-length(which(is.na(x))))))
model.neut.results = cbind(model.neut.results, apply(model.neut, 1, function(x) length(which(x %in% "Undominated"))/(length(x)-length(which(is.na(x))))))

row.names(model.neut.results) = row.names(model.neut)
colnames(model.neut.results) = c("Variable Selection", "Homogenizing Selection", "Homogenizing Dispersal", "Dispersal Limitation", "Undominated")

write.table(model.neut.results, "eco.model.neut.fraction.11-17.results.txt", sep="\t", row.names = TRUE, col.names = TRUE, na="NA")

#Moderate pH wells - CSW1.3, CSW1.5, CSWold, QV1.3
x = bnti.mat[-c(1:16), -c(1:16)]
x=x[-c(7:13), -c(7:13)]
x=x[-c(27:75), -c(27:75)]

y = raup.mat[-c(1:16), -c(1:16)]
y=y[-c(7:13), -c(7:13)]
y=y[-c(27:75), -c(27:75)]

model.mod = matrix(nrow = length(x[,1]), ncol = length(x[1,]), data = as.character(NA))

model.mod[x > 2] = "Variable Selection"
model.mod[x < -2] = "Homogenizing Selection"
model.mod[y < -0.95 & x > -2 & x < 2] = "Homogenizing Dispersal"
model.mod[y > 0.95 & x > -2 & x < 2] = "Dispersal Limitation"
model.mod[y <0.95 & y > -0.95 & x > -2 & x < 2] = "Undominated"

dimnames(model.mod) = dimnames(x)

# Compiling model.mod results for all comparisons
model.mod.results = NULL
model.mod.results = cbind(model.mod.results, apply(model.mod, 1, function(x) length(which(x %in% "Variable Selection"))/(length(x)-length(which(is.na(x))))))
model.mod.results = cbind(model.mod.results, apply(model.mod, 1, function(x) length(which(x %in% "Homogenizing Selection"))/(length(x)-length(which(is.na(x))))))
model.mod.results = cbind(model.mod.results, apply(model.mod, 1, function(x) length(which(x %in% "Homogenizing Dispersal"))/(length(x)-length(which(is.na(x))))))
model.mod.results = cbind(model.mod.results, apply(model.mod, 1, function(x) length(which(x %in% "Dispersal Limitation"))/(length(x)-length(which(is.na(x))))))
model.mod.results = cbind(model.mod.results, apply(model.mod, 1, function(x) length(which(x %in% "Undominated"))/(length(x)-length(which(is.na(x))))))

row.names(model.mod.results) = row.names(model.mod)
colnames(model.mod.results) = c("Variable Selection", "Homogenizing Selection", "Homogenizing Dispersal", "Dispersal Limitation", "Undominated")

write.table(model.mod.results, "eco.model.MOD.fraction.11-17.results.txt", sep="\t", row.names = TRUE, col.names = TRUE, na="NA")

#Extreme pH wells - CSW1.1, N08A, N08B, QV1.1
x = bnti.mat[-c(12:49), -c(12:49)]
x=x[-c(33:41), -c(33:41)]
x=x[-c(45:57), -c(45:57)]

y = raup.mat[-c(12:49), -c(12:49)]
y=y[-c(33:41), -c(33:41)]
y=y[-c(45:57), -c(45:57)]

model.ex = matrix(nrow = length(x[,1]), ncol = length(x[1,]), data = as.character(NA))

model.ex[x > 2] = "Variable Selection"
model.ex[x < -2] = "Homogenizing Selection"
model.ex[y < -0.95 & x > -2 & x < 2] = "Homogenizing Dispersal"
model.ex[y > 0.95 & x > -2 & x < 2] = "Dispersal Limitation"
model.ex[y <0.95 & y > -0.95 & x > -2 & x < 2] = "Undominated"

dimnames(model.ex) = dimnames(x)

# Compiling model.ex results for all comparisons
model.ex.results = NULL
model.ex.results = cbind(model.ex.results, apply(model.ex, 1, function(x) length(which(x %in% "Variable Selection"))/(length(x)-length(which(is.na(x))))))
model.ex.results = cbind(model.ex.results, apply(model.ex, 1, function(x) length(which(x %in% "Homogenizing Selection"))/(length(x)-length(which(is.na(x))))))
model.ex.results = cbind(model.ex.results, apply(model.ex, 1, function(x) length(which(x %in% "Homogenizing Dispersal"))/(length(x)-length(which(is.na(x))))))
model.ex.results = cbind(model.ex.results, apply(model.ex, 1, function(x) length(which(x %in% "Dispersal Limitation"))/(length(x)-length(which(is.na(x))))))
model.ex.results = cbind(model.ex.results, apply(model.ex, 1, function(x) length(which(x %in% "Undominated"))/(length(x)-length(which(is.na(x))))))

row.names(model.ex.results) = row.names(model.ex)
colnames(model.ex.results) = c("Variable Selection", "Homogenizing Selection", "Homogenizing Dispersal", "Dispersal Limitation", "Undominated")

write.table(model.ex.results, "eco.model.EX.fraction.11-17.results.txt", sep="\t", row.names = TRUE, col.names = TRUE, na="NA")

#CSW1.1 Well 
x = bnti.mat[1:11, 1:11]
y = raup.mat[1:11, 1:11]

model.c1 = matrix(nrow = length(x[,1]), ncol = length(x[1,]), data = as.character(NA))

model.c1[x > 2] = "Variable Selection"
model.c1[x < -2] = "Homogenizing Selection"
model.c1[y < -0.95 & x > -2 & x < 2] = "Homogenizing Dispersal"
model.c1[y > 0.95 & x > -2 & x < 2] = "Dispersal Limitation"
model.c1[y <0.95 & y > -0.95 & x > -2 & x < 2] = "Undominated"

dimnames(model.c1) = dimnames(x)

# Compiling model.c1 results for all comparisons
model.c1.results = NULL
model.c1.results = cbind(model.c1.results, apply(model.c1, 1, function(x) length(which(x %in% "Variable Selection"))/(length(x)-length(which(is.na(x))))))
model.c1.results = cbind(model.c1.results, apply(model.c1, 1, function(x) length(which(x %in% "Homogenizing Selection"))/(length(x)-length(which(is.na(x))))))
model.c1.results = cbind(model.c1.results, apply(model.c1, 1, function(x) length(which(x %in% "Homogenizing Dispersal"))/(length(x)-length(which(is.na(x))))))
model.c1.results = cbind(model.c1.results, apply(model.c1, 1, function(x) length(which(x %in% "Dispersal Limitation"))/(length(x)-length(which(is.na(x))))))
model.c1.results = cbind(model.c1.results, apply(model.c1, 1, function(x) length(which(x %in% "Undominated"))/(length(x)-length(which(is.na(x))))))

row.names(model.c1.results) = row.names(model.c1)
colnames(model.c1.results) = c("Variable Selection", "Homogenizing Selection", "Homogenizing Dispersal", "Dispersal Limitation", "Undominated")

write.table(model.c1.results, "eco.model.c1.fraction.11-17.results.txt", sep="\t", row.names = TRUE, col.names = TRUE, na="NA")

#CSW1.2
x = bnti.mat[12:16, 12:16]
y = raup.mat[12:16, 12:16]

model.c2 = matrix(nrow = length(x[,1]), ncol = length(x[1,]), data = as.character(NA))

model.c2[x > 2] = "Variable Selection"
model.c2[x < -2] = "Homogenizing Selection"
model.c2[y < -0.95 & x > -2 & x < 2] = "Homogenizing Dispersal"
model.c2[y > 0.95 & x > -2 & x < 2] = "Dispersal Limitation"
model.c2[y <0.95 & y > -0.95 & x > -2 & x < 2] = "Undominated"

dimnames(model.c2) = dimnames(x)

# Compiling model.c2 results for all comparisons
model.c2.results = NULL
model.c2.results = cbind(model.c2.results, apply(model.c2, 1, function(x) length(which(x %in% "Variable Selection"))/(length(x)-length(which(is.na(x))))))
model.c2.results = cbind(model.c2.results, apply(model.c2, 1, function(x) length(which(x %in% "Homogenizing Selection"))/(length(x)-length(which(is.na(x))))))
model.c2.results = cbind(model.c2.results, apply(model.c2, 1, function(x) length(which(x %in% "Homogenizing Dispersal"))/(length(x)-length(which(is.na(x))))))
model.c2.results = cbind(model.c2.results, apply(model.c2, 1, function(x) length(which(x %in% "Dispersal Limitation"))/(length(x)-length(which(is.na(x))))))
model.c2.results = cbind(model.c2.results, apply(model.c2, 1, function(x) length(which(x %in% "Undominated"))/(length(x)-length(which(is.na(x))))))

row.names(model.c2.results) = row.names(model.c2)
colnames(model.c2.results) = c("Variable Selection", "Homogenizing Selection", "Homogenizing Dispersal", "Dispersal Limitation", "Undominated")

write.table(model.c2.results, "eco.model.c2.fraction.11-17.results.txt", sep="\t", row.names = TRUE, col.names = TRUE, na="NA")

#CSW1.3
x = bnti.mat[17:22, 17:22]
y = raup.mat[17:22, 17:22]

model.c3 = matrix(nrow = length(x[,1]), ncol = length(x[1,]), data = as.character(NA))

model.c3[x > 2] = "Variable Selection"
model.c3[x < -2] = "Homogenizing Selection"
model.c3[y < -0.95 & x > -2 & x < 2] = "Homogenizing Dispersal"
model.c3[y > 0.95 & x > -2 & x < 2] = "Dispersal Limitation"
model.c3[y <0.95 & y > -0.95 & x > -2 & x < 2] = "Undominated"

dimnames(model.c3) = dimnames(x)

# Compiling model.c3 results for all comparisons
model.c3.results = NULL
model.c3.results = cbind(model.c3.results, apply(model.c3, 1, function(x) length(which(x %in% "Variable Selection"))/(length(x)-length(which(is.na(x))))))
model.c3.results = cbind(model.c3.results, apply(model.c3, 1, function(x) length(which(x %in% "Homogenizing Selection"))/(length(x)-length(which(is.na(x))))))
model.c3.results = cbind(model.c3.results, apply(model.c3, 1, function(x) length(which(x %in% "Homogenizing Dispersal"))/(length(x)-length(which(is.na(x))))))
model.c3.results = cbind(model.c3.results, apply(model.c3, 1, function(x) length(which(x %in% "Dispersal Limitation"))/(length(x)-length(which(is.na(x))))))
model.c3.results = cbind(model.c3.results, apply(model.c3, 1, function(x) length(which(x %in% "Undominated"))/(length(x)-length(which(is.na(x))))))

row.names(model.c3.results) = row.names(model.c3)
colnames(model.c3.results) = c("Variable Selection", "Homogenizing Selection", "Homogenizing Dispersal", "Dispersal Limitation", "Undominated")

write.table(model.c3.results, "eco.model.c3.fraction.11-17.results.txt", sep="\t", row.names = TRUE, col.names = TRUE, na="NA")

#CSW1.4
x = bnti.mat[23:29, 23:29]
y = raup.mat[23:29, 23:29]

model.c4 = matrix(nrow = length(x[,1]), ncol = length(x[1,]), data = as.character(NA))

model.c4[x > 2] = "Variable Selection"
model.c4[x < -2] = "Homogenizing Selection"
model.c4[y < -0.95 & x > -2 & x < 2] = "Homogenizing Dispersal"
model.c4[y > 0.95 & x > -2 & x < 2] = "Dispersal Limitation"
model.c4[y <0.95 & y > -0.95 & x > -2 & x < 2] = "Undominated"

dimnames(model.c4) = dimnames(x)

# Compiling model.c4 results for all comparisons
model.c4.results = NULL
model.c4.results = cbind(model.c4.results, apply(model.c4, 1, function(x) length(which(x %in% "Variable Selection"))/(length(x)-length(which(is.na(x))))))
model.c4.results = cbind(model.c4.results, apply(model.c4, 1, function(x) length(which(x %in% "Homogenizing Selection"))/(length(x)-length(which(is.na(x))))))
model.c4.results = cbind(model.c4.results, apply(model.c4, 1, function(x) length(which(x %in% "Homogenizing Dispersal"))/(length(x)-length(which(is.na(x))))))
model.c4.results = cbind(model.c4.results, apply(model.c4, 1, function(x) length(which(x %in% "Dispersal Limitation"))/(length(x)-length(which(is.na(x))))))
model.c4.results = cbind(model.c4.results, apply(model.c4, 1, function(x) length(which(x %in% "Undominated"))/(length(x)-length(which(is.na(x))))))

row.names(model.c4.results) = row.names(model.c4)
colnames(model.c4.results) = c("Variable Selection", "Homogenizing Selection", "Homogenizing Dispersal", "Dispersal Limitation", "Undominated")

write.table(model.c4.results, "eco.model.c4.fraction.11-17.results.txt", sep="\t", row.names = TRUE, col.names = TRUE, na="NA")

#CSW1.5
x = bnti.mat[30:36, 30:36]
y = raup.mat[30:36, 30:36]

model.c5 = matrix(nrow = length(x[,1]), ncol = length(x[1,]), data = as.character(NA))

model.c5[x > 2] = "Variable Selection"
model.c5[x < -2] = "Homogenizing Selection"
model.c5[y < -0.95 & x > -2 & x < 2] = "Homogenizing Dispersal"
model.c5[y > 0.95 & x > -2 & x < 2] = "Dispersal Limitation"
model.c5[y <0.95 & y > -0.95 & x > -2 & x < 2] = "Undominated"

dimnames(model.c5) = dimnames(x)

# Compiling model.c5 results for all comparisons
model.c5.results = NULL
model.c5.results = cbind(model.c5.results, apply(model.c5, 1, function(x) length(which(x %in% "Variable Selection"))/(length(x)-length(which(is.na(x))))))
model.c5.results = cbind(model.c5.results, apply(model.c5, 1, function(x) length(which(x %in% "Homogenizing Selection"))/(length(x)-length(which(is.na(x))))))
model.c5.results = cbind(model.c5.results, apply(model.c5, 1, function(x) length(which(x %in% "Homogenizing Dispersal"))/(length(x)-length(which(is.na(x))))))
model.c5.results = cbind(model.c5.results, apply(model.c5, 1, function(x) length(which(x %in% "Dispersal Limitation"))/(length(x)-length(which(is.na(x))))))
model.c5.results = cbind(model.c5.results, apply(model.c5, 1, function(x) length(which(x %in% "Undominated"))/(length(x)-length(which(is.na(x))))))

row.names(model.c5.results) = row.names(model.c5)
colnames(model.c5.results) = c("Variable Selection", "Homogenizing Selection", "Homogenizing Dispersal", "Dispersal Limitation", "Undominated")

write.table(model.c5.results, "eco.model.c5.fraction.11-17.results.txt", sep="\t", row.names = TRUE, col.names = TRUE, na="NA")

#CSWold
x = bnti.mat[37:49, 37:49]
y = raup.mat[37:49, 37:49]

model.c.old = matrix(nrow = length(x[,1]), ncol = length(x[1,]), data = as.character(NA))

model.c.old[x > 2] = "Variable Selection"
model.c.old[x < -2] = "Homogenizing Selection"
model.c.old[y < -0.95 & x > -2 & x < 2] = "Homogenizing Dispersal"
model.c.old[y > 0.95 & x > -2 & x < 2] = "Dispersal Limitation"
model.c.old[y <0.95 & y > -0.95 & x > -2 & x < 2] = "Undominated"

dimnames(model.c.old) = dimnames(x)

# Compiling model.c.old results for all comparisons
model.c.old.results = NULL
model.c.old.results = cbind(model.c.old.results, apply(model.c.old, 1, function(x) length(which(x %in% "Variable Selection"))/(length(x)-length(which(is.na(x))))))
model.c.old.results = cbind(model.c.old.results, apply(model.c.old, 1, function(x) length(which(x %in% "Homogenizing Selection"))/(length(x)-length(which(is.na(x))))))
model.c.old.results = cbind(model.c.old.results, apply(model.c.old, 1, function(x) length(which(x %in% "Homogenizing Dispersal"))/(length(x)-length(which(is.na(x))))))
model.c.old.results = cbind(model.c.old.results, apply(model.c.old, 1, function(x) length(which(x %in% "Dispersal Limitation"))/(length(x)-length(which(is.na(x))))))
model.c.old.results = cbind(model.c.old.results, apply(model.c.old, 1, function(x) length(which(x %in% "Undominated"))/(length(x)-length(which(is.na(x))))))

row.names(model.c.old.results) = row.names(model.c.old)
colnames(model.c.old.results) = c("Variable Selection", "Homogenizing Selection", "Homogenizing Dispersal", "Dispersal Limitation", "Undominated")

write.table(model.c.old.results, "eco.model.c.old.fraction.11-17.results.txt", sep="\t", row.names = TRUE, col.names = TRUE, na="NA")

#N08A
x = bnti.mat[50:60, 50:60]
y = raup.mat[50:60, 50:60]

model.nA = matrix(nrow = length(x[,1]), ncol = length(x[1,]), data = as.character(NA))

model.nA[x > 2] = "Variable Selection"
model.nA[x < -2] = "Homogenizing Selection"
model.nA[y < -0.95 & x > -2 & x < 2] = "Homogenizing Dispersal"
model.nA[y > 0.95 & x > -2 & x < 2] = "Dispersal Limitation"
model.nA[y <0.95 & y > -0.95 & x > -2 & x < 2] = "Undominated"

dimnames(model.nA) = dimnames(x)

# Compiling model.nA results for all comparisons
model.nA.results = NULL
model.nA.results = cbind(model.nA.results, apply(model.nA, 1, function(x) length(which(x %in% "Variable Selection"))/(length(x)-length(which(is.na(x))))))
model.nA.results = cbind(model.nA.results, apply(model.nA, 1, function(x) length(which(x %in% "Homogenizing Selection"))/(length(x)-length(which(is.na(x))))))
model.nA.results = cbind(model.nA.results, apply(model.nA, 1, function(x) length(which(x %in% "Homogenizing Dispersal"))/(length(x)-length(which(is.na(x))))))
model.nA.results = cbind(model.nA.results, apply(model.nA, 1, function(x) length(which(x %in% "Dispersal Limitation"))/(length(x)-length(which(is.na(x))))))
model.nA.results = cbind(model.nA.results, apply(model.nA, 1, function(x) length(which(x %in% "Undominated"))/(length(x)-length(which(is.na(x))))))

row.names(model.nA.results) = row.names(model.nA)
colnames(model.nA.results) = c("Variable Selection", "Homogenizing Selection", "Homogenizing Dispersal", "Dispersal Limitation", "Undominated")

write.table(model.nA.results, "eco.model.nA.fraction.11-17.results.txt", sep="\t", row.names = TRUE, col.names = TRUE, na="NA")

#N08B
x = bnti.mat[61:70, 61:70]
y = raup.mat[61:70, 61:70]

model.nB = matrix(nrow = length(x[,1]), ncol = length(x[1,]), data = as.character(NA))

model.nB[x > 2] = "Variable Selection"
model.nB[x < -2] = "Homogenizing Selection"
model.nB[y < -0.95 & x > -2 & x < 2] = "Homogenizing Dispersal"
model.nB[y > 0.95 & x > -2 & x < 2] = "Dispersal Limitation"
model.nB[y <0.95 & y > -0.95 & x > -2 & x < 2] = "Undominated"

dimnames(model.nB) = dimnames(x)

# Compiling model.nB results for all comparisons
model.nB.results = NULL
model.nB.results = cbind(model.nB.results, apply(model.nB, 1, function(x) length(which(x %in% "Variable Selection"))/(length(x)-length(which(is.na(x))))))
model.nB.results = cbind(model.nB.results, apply(model.nB, 1, function(x) length(which(x %in% "Homogenizing Selection"))/(length(x)-length(which(is.na(x))))))
model.nB.results = cbind(model.nB.results, apply(model.nB, 1, function(x) length(which(x %in% "Homogenizing Dispersal"))/(length(x)-length(which(is.na(x))))))
model.nB.results = cbind(model.nB.results, apply(model.nB, 1, function(x) length(which(x %in% "Dispersal Limitation"))/(length(x)-length(which(is.na(x))))))
model.nB.results = cbind(model.nB.results, apply(model.nB, 1, function(x) length(which(x %in% "Undominated"))/(length(x)-length(which(is.na(x))))))

row.names(model.nB.results) = row.names(model.nB)
colnames(model.nB.results) = c("Variable Selection", "Homogenizing Selection", "Homogenizing Dispersal", "Dispersal Limitation", "Undominated")

write.table(model.nB.results, "eco.model.nB.fraction.11-17.results.txt", sep="\t", row.names = TRUE, col.names = TRUE, na="NA")

#N08C
x = bnti.mat[71:79, 71:79]
y = raup.mat[71:79, 71:79]

model.nC = matrix(nrow = length(x[,1]), ncol = length(x[1,]), data = as.character(NA))

model.nC[x > 2] = "Variable Selection"
model.nC[x < -2] = "Homogenizing Selection"
model.nC[y < -0.95 & x > -2 & x < 2] = "Homogenizing Dispersal"
model.nC[y > 0.95 & x > -2 & x < 2] = "Dispersal Limitation"
model.nC[y <0.95 & y > -0.95 & x > -2 & x < 2] = "Undominated"

dimnames(model.nC) = dimnames(x)

# Compiling model.nC results for all comparisons
model.nC.results = NULL
model.nC.results = cbind(model.nC.results, apply(model.nC, 1, function(x) length(which(x %in% "Variable Selection"))/(length(x)-length(which(is.na(x))))))
model.nC.results = cbind(model.nC.results, apply(model.nC, 1, function(x) length(which(x %in% "Homogenizing Selection"))/(length(x)-length(which(is.na(x))))))
model.nC.results = cbind(model.nC.results, apply(model.nC, 1, function(x) length(which(x %in% "Homogenizing Dispersal"))/(length(x)-length(which(is.na(x))))))
model.nC.results = cbind(model.nC.results, apply(model.nC, 1, function(x) length(which(x %in% "Dispersal Limitation"))/(length(x)-length(which(is.na(x))))))
model.nC.results = cbind(model.nC.results, apply(model.nC, 1, function(x) length(which(x %in% "Undominated"))/(length(x)-length(which(is.na(x))))))

row.names(model.nC.results) = row.names(model.nC)
colnames(model.nC.results) = c("Variable Selection", "Homogenizing Selection", "Homogenizing Dispersal", "Dispersal Limitation", "Undominated")

write.table(model.nC.results, "eco.model.nC.fraction.11-17.results.txt", sep="\t", row.names = TRUE, col.names = TRUE, na="NA")

#QV1.1
x = bnti.mat[80:91, 80:91]
y = raup.mat[80:91, 80:91]

model.q1 = matrix(nrow = length(x[,1]), ncol = length(x[1,]), data = as.character(NA))

model.q1[x > 2] = "Variable Selection"
model.q1[x < -2] = "Homogenizing Selection"
model.q1[y < -0.95 & x > -2 & x < 2] = "Homogenizing Dispersal"
model.q1[y > 0.95 & x > -2 & x < 2] = "Dispersal Limitation"
model.q1[y <0.95 & y > -0.95 & x > -2 & x < 2] = "Undominated"

dimnames(model.q1) = dimnames(x)

# Compiling model.q1 results for all comparisons
model.q1.results = NULL
model.q1.results = cbind(model.q1.results, apply(model.q1, 1, function(x) length(which(x %in% "Variable Selection"))/(length(x)-length(which(is.na(x))))))
model.q1.results = cbind(model.q1.results, apply(model.q1, 1, function(x) length(which(x %in% "Homogenizing Selection"))/(length(x)-length(which(is.na(x))))))
model.q1.results = cbind(model.q1.results, apply(model.q1, 1, function(x) length(which(x %in% "Homogenizing Dispersal"))/(length(x)-length(which(is.na(x))))))
model.q1.results = cbind(model.q1.results, apply(model.q1, 1, function(x) length(which(x %in% "Dispersal Limitation"))/(length(x)-length(which(is.na(x))))))
model.q1.results = cbind(model.q1.results, apply(model.q1, 1, function(x) length(which(x %in% "Undominated"))/(length(x)-length(which(is.na(x))))))

row.names(model.q1.results) = row.names(model.q1)
colnames(model.q1.results) = c("Variable Selection", "Homogenizing Selection", "Homogenizing Dispersal", "Dispersal Limitation", "Undominated")

write.table(model.q1.results, "eco.model.q1.fraction.11-17.results.txt", sep="\t", row.names = TRUE, col.names = TRUE, na="NA")

#QV1.2
x = bnti.mat[92:98, 92:98]
y = raup.mat[92:98, 92:98]

model.q2 = matrix(nrow = length(x[,1]), ncol = length(x[1,]), data = as.character(NA))

model.q2[x > 2] = "Variable Selection"
model.q2[x < -2] = "Homogenizing Selection"
model.q2[y < -0.95 & x > -2 & x < 2] = "Homogenizing Dispersal"
model.q2[y > 0.95 & x > -2 & x < 2] = "Dispersal Limitation"
model.q2[y <0.95 & y > -0.95 & x > -2 & x < 2] = "Undominated"

dimnames(model.q2) = dimnames(x)

# Compiling model.q2 results for all comparisons
model.q2.results = NULL
model.q2.results = cbind(model.q2.results, apply(model.q2, 1, function(x) length(which(x %in% "Variable Selection"))/(length(x)-length(which(is.na(x))))))
model.q2.results = cbind(model.q2.results, apply(model.q2, 1, function(x) length(which(x %in% "Homogenizing Selection"))/(length(x)-length(which(is.na(x))))))
model.q2.results = cbind(model.q2.results, apply(model.q2, 1, function(x) length(which(x %in% "Homogenizing Dispersal"))/(length(x)-length(which(is.na(x))))))
model.q2.results = cbind(model.q2.results, apply(model.q2, 1, function(x) length(which(x %in% "Dispersal Limitation"))/(length(x)-length(which(is.na(x))))))
model.q2.results = cbind(model.q2.results, apply(model.q2, 1, function(x) length(which(x %in% "Undominated"))/(length(x)-length(which(is.na(x))))))

row.names(model.q2.results) = row.names(model.q2)
colnames(model.q2.results) = c("Variable Selection", "Homogenizing Selection", "Homogenizing Dispersal", "Dispersal Limitation", "Undominated")

write.table(model.q2.results, "eco.model.q2.fraction.11-17.results.txt", sep="\t", row.names = TRUE, col.names = TRUE, na="NA")

#QV1.3
x = bnti.mat[99:104, 99:104]
y = raup.mat[99:104, 99:104]

model.q3 = matrix(nrow = length(x[,1]), ncol = length(x[1,]), data = as.character(NA))

model.q3[x > 2] = "Variable Selection"
model.q3[x < -2] = "Homogenizing Selection"
model.q3[y < -0.95 & x > -2 & x < 2] = "Homogenizing Dispersal"
model.q3[y > 0.95 & x > -2 & x < 2] = "Dispersal Limitation"
model.q3[y <0.95 & y > -0.95 & x > -2 & x < 2] = "Undominated"

dimnames(model.q3) = dimnames(x)

# Compiling model.q3 results for all comparisons
model.q3.results = NULL
model.q3.results = cbind(model.q3.results, apply(model.q3, 1, function(x) length(which(x %in% "Variable Selection"))/(length(x)-length(which(is.na(x))))))
model.q3.results = cbind(model.q3.results, apply(model.q3, 1, function(x) length(which(x %in% "Homogenizing Selection"))/(length(x)-length(which(is.na(x))))))
model.q3.results = cbind(model.q3.results, apply(model.q3, 1, function(x) length(which(x %in% "Homogenizing Dispersal"))/(length(x)-length(which(is.na(x))))))
model.q3.results = cbind(model.q3.results, apply(model.q3, 1, function(x) length(which(x %in% "Dispersal Limitation"))/(length(x)-length(which(is.na(x))))))
model.q3.results = cbind(model.q3.results, apply(model.q3, 1, function(x) length(which(x %in% "Undominated"))/(length(x)-length(which(is.na(x))))))

row.names(model.q3.results) = row.names(model.q3)
colnames(model.q3.results) = c("Variable Selection", "Homogenizing Selection", "Homogenizing Dispersal", "Dispersal Limitation", "Undominated")

write.table(model.q3.results, "eco.model.q3.fraction.11-17.results.txt", sep="\t", row.names = TRUE, col.names = TRUE, na="NA")


#######################################################################################################

#End code from Danczak et al. 2018#

#Calculating distance matrices for all variables for mantel tests
########################################################

#Extract normalized count table and calculate bray curtis
counts=t(as.data.frame(as.matrix(otu_table(merged.props))))
bray.dist=distance(counts, method="bray-curtis")
bray.dist.mat=as.matrix(bray.dist) #Mean is 0.8 --> highly dissimilar compositions

#Extract other meta.data, convert to distance then vector
days=meta.data$Days
days.dist=dist(days, method="euclidean", upper=FALSE)

depth=meta.data$Well_Depth_m
depth.dist=dist(depth, method="euclidean", upper=FALSE)

temp=meta.data$Temp
temp.dist=dist(temp, method="euclidean", upper=FALSE)

ph=meta.data$pH
ph.dist=dist(ph, method="euclidean", upper=FALSE)

cond=meta.data$Cond_mS
cond.dist=dist(cond, method="euclidean", upper=FALSE)

do=meta.data$DO_mgL
do.dist=dist(do, method="euclidean", upper=FALSE)

orp=meta.data$ORP_mV
orp.dist=dist(orp, method="euclidean", upper=FALSE)

dic=meta.data$DIC_uM
dic.dist=dist(dic, method="euclidean", upper=FALSE)

################################################################


#Mantel Tests on Full Dataset - Want these values
###########################################################################
library(ade4)
#BNTI
mantel.bnti.raup=mantel.rtest(bnti.dist, raup.dist, nrepet = 9999) #R=0.342 sig.
mantel.bnti.bray=mantel.rtest(bnti.dist, bray.dist, nrepet = 9999) #R=0.368 sig.

mantel.bnti.ph=mantel.rtest(bnti.dist, ph.dist, nrepet=9999) #R=0.292, p=1e-04
mantel.bnti.orp=mantel.rtest(bnti.dist, orp.dist, nrepet=9999) #R=0.227 p=5e-04
mantel.bnti.dic=mantel.rtest(bnti.dist, dic.dist, nrepet = 9999) #R=0.135 p=0.0436

mantel.bnti.temp=mantel.rtest(bnti.dist, temp.dist, nrepet = 9999) #not sig.
mantel.bnti.depth=mantel.rtest(bnti.dist, depth.dist, nrepet = 9999) #not sig.

#RCbray
mantel.raup.bray=mantel.rtest(raup.dist, bray.dist, nrepet = 9999) #R=0.758 sig.

mantel.raup.depth=mantel.rtest(all.raup.dist, depth.dist, nrepet=9999) #R=0.315 p=1e-04
mantel.raup.temp=mantel.rtest(all.raup.dist, temp.dist, nrepet = 9999) #R=0.148 p=0.0057
mantel.raup.ph=mantel.rtest(all.raup.dist, ph.dist, nrepet=9999) #R=0.268 sig. p=1e-04
mantel.raup.cond=mantel.rtest(all.raup.dist, cond.dist, nrepet=9999) #R=0.336 p=1e-04
mantel.raup.orp=mantel.rtest(all.raup.dist, orp.dist, nrepet=9999) #R=0.299 p=1e-04
mantel.raup.dic=mantel.rtest(all.raup.dist, dic.dist, nrepet = 9999) #R=0.331 p=1e-04

####################################################################################


#Calculating richness and evenness on data and making plots 
######################################################
#Pielou's evenness plots
#Transpose data for vegan
count.table2=t(count.table)

#First calculate Shannon Index
H=diversity(count.table2, index="shannon")

#Calulate species richness
s.rich=specnumber(count.table2)

#Calculate Pielou's evenness
p.even=H/(log(s.rich))

#Now can plot evenness and richness vs. meta.data
#Create data frame and add on evenness and richness to meta.data
s.rich=as.data.frame(s.rich)
p.even=as.data.frame(p.even)

meta.data2=cbind(meta.data, s.rich, p.even)
meta.data2=as.data.frame(meta.data2)

#Now can makes plots
#Richness Plot
rich.ph.plot=ggplot(data=meta.data2, aes(x=pH, y=s.rich)) + geom_point(aes(color=Well), size=3)
rich.ph.plot=rich.ph.plot + labs(x="pH", y="Species Richness") + scale_color_manual(values=pal2(12), limits=all.x.axis.depth)
rich.ph.plot=rich.ph.plot + stat_smooth(method="lm", formula=y~x, color="red", se=FALSE)
rich.ph.plot=rich.ph.plot + ggtitle("CROMO Species Richness with pH") + theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
rich.ph.plot=rich.ph.plot + theme(legend.position = "none")
rich.ph.plot

#EVenness Plot
even.ph.plot=ggplot(data=meta.data2, aes(x=pH, y=p.even)) + geom_point(aes(color=Well), size=3)
even.ph.plot=even.ph.plot + labs(x="pH", y="Pielou's Evenness") + scale_color_manual(values=pal2(12), limits=all.x.axis.depth)
even.ph.plot=even.ph.plot + stat_smooth(method="lm", formula=y~x, color="red", se=FALSE)
even.ph.plot=even.ph.plot + ggtitle("CROMO Pielou's Evenness with pH") + theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
even.ph.plot

#Evenness and Richness Combined Plot for pH
rich.even.plots=grid.arrange(rich.ph.plot, even.ph.plot, nrow=1, widths=c(1.5, 1.75))

###########################################################################################


#Average Geochemistry Plots 
######################################################
Avg_Well_Data <- read_excel("Avg_Well_Data_for_Geochem_Plots.xlsx")
avg.data=as.data.frame(Avg_Well_Data)

#Temperature Plot
temp.plot=ggplot(data=avg.data, aes(x=depth_m, y=temp)) + geom_point(aes(color=Well), size=4)
temp.plot=temp.plot + labs(x="Depth (m)", y="Temperature (C)") + scale_color_manual(values=pal2(12), limits=all.x.axis.ph)
temp.plot=temp.plot + stat_smooth(method="lm", formula = y~x, color="black", se=FALSE)
temp.plot=temp.plot + ggtitle("Average Temperature") + theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())  
temp.plot=temp.plot + geom_errorbar(aes(ymin=temp-temp.sd, ymax=temp+temp.sd), width=0.2) + theme(legend.position="none")
temp.plot

temp.lm=lm(temp~depth_m, data=avg.data) #R2=0.6155 sig. p=0.001527

#Specific Conductance Plot
cond.plot=ggplot(data=avg.data, aes(x=depth_m, y=sp.cond_mS)) + geom_point(aes(color=Well), size=4)
cond.plot=cond.plot + labs(x="Depth (m)", y="Specific conductance (mS)") + scale_color_manual(values=pal2(12), limits=all.x.axis.ph)
cond.plot=cond.plot + stat_smooth(method="lm", formula = y~x, color="black", se=FALSE)
cond.plot=cond.plot + ggtitle("Average Specific Conductance") + theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())  
cond.plot=cond.plot + geom_errorbar(aes(ymin=sp.cond_mS-sp.cond_mS.sd, ymax=sp.cond_mS+sp.cond_mS.sd), width=0.2) + theme(legend.position="none")
cond.plot

cond.lm=lm(sp.cond_mS~depth_m, data=avg.data) #R2=0.8802 sig. p=3.952E-06

#pH Plot
ph.plot=ggplot(data=avg.data, aes(x=depth_m, y=ph)) + geom_point(aes(color=Well), size=4)
ph.plot=ph.plot + labs(x="Depth (m)", y="pH") + scale_color_manual(values=pal2(12), limits=all.x.axis.ph)
ph.plot=ph.plot + stat_smooth(method="lm", formula = y~poly(x, 2), color="black", se=FALSE)
ph.plot=ph.plot + ggtitle("Average pH") + theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())  
ph.plot=ph.plot + geom_errorbar(aes(ymin=ph-ph.sd, ymax=ph+ph.sd), width=0.2) + theme(legend.position="none")
ph.plot

ph.lm=lm(ph~poly(depth_m, 2), data=avg.data) #R2=0.1669 not sig. p=0.1782

#Dissolved Oxygen Plot
do.plot=ggplot(data=avg.data, aes(x=depth_m, y=do_mgL)) + geom_point(aes(color=Well), size=4)
do.plot=do.plot + labs(x="Depth (m)", y="Dissolved Oxygen (mg/L)") + scale_color_manual(values=pal2(12), limits=all.x.axis.ph)
do.plot=do.plot + stat_smooth(method="lm", formula = y~log(x), color="black", se=FALSE)
do.plot=do.plot + ggtitle("Average Dissolved Oxygen") + theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())  
do.plot=do.plot + geom_errorbar(aes(ymin=do_mgL-do_mgL.sd, ymax=do_mgL+do_mgL.sd), width=0.2) + theme(legend.position="none")
do.plot

do.lm=lm(do_mgL~log(depth_m), data=avg.data) #R2=0.2248 not sig. p=0.06787

#ORP Plot
orp.plot=ggplot(data=avg.data, aes(x=depth_m, y=orp_mV)) + geom_point(aes(color=Well), size=4)
orp.plot=orp.plot + labs(x="Depth (m)", y="Oxidation Reduction Potential (mV)") + scale_color_manual(values=pal2(12), limits=all.x.axis.ph)
orp.plot=orp.plot + stat_smooth(method="lm", formula = y~log(x), color="black", se=FALSE)
orp.plot=orp.plot + ggtitle("Average Oxidation Reduction Potential") + theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())  
orp.plot=orp.plot + geom_errorbar(aes(ymin=orp_mV-orp_mV.sd, ymax=orp_mV+orp_mV.sd), width=0.2) + theme(legend.position="none")
orp.plot

orp.lm=lm(orp_mV~log(depth_m), data=avg.data) #R2=0.4705 sig. p=0.008251

#DIC Plot
dic.plot=ggplot(data=avg.data, aes(x=depth_m, y=dic_uM)) + geom_point(aes(color=Well), size=4)
dic.plot=dic.plot + labs(x="Depth (m)", y="Dissolved Inorganic Carbon (uM)") + scale_color_manual(values=pal2(12), limits=all.x.axis.ph)
dic.plot=dic.plot + stat_smooth(method="lm", formula = y~log(x), color="black", se=FALSE)
dic.plot=dic.plot + ggtitle("Average Dissolved Inorganic Carbon") + theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())  
dic.plot=dic.plot + geom_errorbar(aes(ymin=dic_uM-dic_uM.sd, ymax=dic_uM+dic_uM.sd), width=0.2) + theme(legend.position="bottom") + guides(col=guide_legend(nrow=1))
dic.plot

dic.lm=lm(dic_uM~log(depth_m), data=avg.data) #R2=0.4522 sig. p=0.009899 

#Extract legend for these plots
ph.legend<-g_legend(dic.plot)

#Evenness and Richness Combined Plot for pH
chem.plots=grid.arrange(arrangeGrob(temp.plot, ph.plot, cond.plot, do.plot, orp.plot, dic.plot + theme(legend.position = "none"), nrow=2, ncol=3), ph.legend, nrow=2, heights=c(10,1))


#Ecological Modeling Results Plots
################################################################

#UPload ecological model results
Eco_Model_Data_for_R <- read_excel("Eco_Model_Percent_Data_for_R.xlsx")
eco.dat=as.data.frame(Eco_Model_Data_for_R)

#Subset second dataset by ecological process
var.sel.df=eco.dat[1:4,]
hom.sel.df=eco.dat[5:8,]
hom.dis.df=eco.dat[9:12,]
dis.lim.df=eco.dat[13:16,]
undom.df=eco.dat[17:20,]

#Create Legend Order for well clusters
clust.ord=c("All", "Neutral", "Moderate", "Extreme")

#Make updated color scale for Plots
display.brewer.all()
brewer.pal(9, "Greys")
pal=c("#DFDFDF", "#C1C1C1", "#828282", "#454545", "#000000")


#Make Variable Selection Plot
var.sel.plot= ggplot(var.sel.df, aes(x=ID, y=percent, fill=ID)) + geom_bar(stat="identity")
var.sel.plot= var.sel.plot + ggtitle("Contribution of Variable Selection") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
var.sel.plot= var.sel.plot + theme(plot.title = element_text(hjust=0.5)) + scale_fill_manual(values=pal, limits=clust.ord)
var.sel.plot = var.sel.plot + labs(x="", y="", fill="Well Cluster") + scale_x_discrete(limits=clust.ord) 
var.sel.plot = var.sel.plot + theme(legend.position = "none")
var.sel.plot

#Make Homogenizing Selection Plot
hom.sel.plot= ggplot(hom.sel.df, aes(x=ID, y=percent, fill=ID)) + geom_bar(stat="identity")
hom.sel.plot= hom.sel.plot + ggtitle("Contribution of Homogenizing Selection") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
hom.sel.plot= hom.sel.plot + theme(plot.title = element_text(hjust=0.5)) + scale_fill_manual(values=pal, limits=clust.ord)
hom.sel.plot = hom.sel.plot + labs(x="", y="", fill="Well Cluster") + scale_x_discrete(limits=clust.ord) 
hom.sel.plot = hom.sel.plot + theme(legend.position="none")
hom.sel.plot

#Make Homogenizing Dispersal Plot
hom.dis.plot= ggplot(hom.dis.df, aes(x=ID, y=percent, fill=ID)) + geom_bar(stat="identity")
hom.dis.plot= hom.dis.plot + ggtitle("Contribution of Homogenizing Dispersal") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
hom.dis.plot= hom.dis.plot + theme(plot.title = element_text(hjust=0.5)) + scale_fill_manual(values=pal, limits=clust.ord)
hom.dis.plot = hom.dis.plot + labs(x="", y="", fill="Well Cluster") + scale_x_discrete(limits=clust.ord) 
hom.dis.plot = hom.dis.plot + theme(legend.position = "none")
hom.dis.plot

#Make Dispersal Limitation Plot
dis.lim.plot= ggplot(dis.lim.df, aes(x=ID, y=percent, fill=ID)) + geom_bar(stat="identity")
dis.lim.plot= dis.lim.plot + ggtitle("Contribution of Dispersal Limitation") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dis.lim.plot= dis.lim.plot + theme(plot.title = element_text(hjust=0.5)) + scale_fill_manual(values=pal, limits=clust.ord)
dis.lim.plot = dis.lim.plot + labs(x="", y="", fill="Well Cluster") + scale_x_discrete(limits=clust.ord) 
dis.lim.plot = dis.lim.plot + theme(legend.position = "none")
dis.lim.plot

#Make Undominated Plot
undom.plot= ggplot(undom.df, aes(x=ID, y=percent, fill=ID)) + geom_bar(stat="identity")
undom.plot= undom.plot + ggtitle("Contribution of Undominated Processes") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
undom.plot= undom.plot + theme(plot.title = element_text(hjust=0.5)) + scale_fill_manual(values=pal, limits=clust.ord)
undom.plot = undom.plot + labs(x="", y="", fill="Well Cluster") + scale_x_discrete(limits=clust.ord) 
undom.plot = undom.plot + theme(legend.position = "bottom")
undom.plot 

#Extract bottom lying legend
eco.legend=g_legend(undom.plot)

eco.model.plots=grid.arrange(var.sel.plot, hom.sel.plot, dis.lim.plot, hom.dis.plot, undom.plot + theme(legend.position = "none"), nrow=2, ncol=3)

