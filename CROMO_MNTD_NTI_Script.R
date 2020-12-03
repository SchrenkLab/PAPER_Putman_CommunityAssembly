#CROMO MNTD Analysis
#Performed on Brazelton Computing Cluster --Code Pasted into sbatch job file
#Run using: 
#srun --partition highmem --mem 80GB Rscript CRO_MNTD-analysis.r &
#job ID 27051

#load librarys
library(vegan)
library(ape)
library(picante)
library(readr)

#Read in Otu.ids -- BEFORE READING IN REPLACE SPACES WITH AN "_"
otu.ids=read.delim("CROMO_TS_filtered_FINAL_otu.ids.txt", header=FALSE)
otu.ids=as.data.frame(otu.ids)

#Upload phylogenetic tree
CROMO_Tree=read_file("CROMO_TS_FINAL.tre")
phy=read.tree(text=CROMO_Tree)

#Extract reordered otu.ids
reorg.otu.ids=as.data.frame(phy$tip.label)

#Read in count table
comm=read.delim("CROMO_TS_filtered_FINAL_counts.txt", header=TRUE)
comm=as.data.frame(comm)

#Set rownames remove column 1 and transpose
row.names(comm)=otu.ids[,1]

#Reorganize count table by tree otu.ids 
comm=cbind(otu.ids, comm)
order=reorg.otu.ids[,1]
comm=comm[match(order, comm$V1),]

#Remove organizational column and transpose
comm[,1]<-NULL
comm=t(comm)

#Create combined object
combined=match.phylo.comm(phy, comm)
phy <- combined$phy
comm <- combined$comm

#Normalize Branch Lengths
#branch.sum=sum(phy$edge.length)
#phy$edge.length=phy$edge.length/branch.sum

#Calculate phylogenetic distance
#phy.dist=cophenetic.phylo(phy)
#phy.dist1=cophenetic(phy)

#Calculate MNTD
CRO_weighted.mntd=mntd(comm, cophenetic(phy), abundance.weighted = TRUE)

#Write out MNTD matrix
write.csv(CRO_weighted.mntd,'CRO_TS_FINAL_MNTD_weighted.csv', quote=F)

#Caluclate NTI
CRO_ses.mntd=ses.mntd(comm, cophenetic(phy), null.model="taxa.labels", abundance.weighted = TRUE, runs=2, iterations=1000)

#Write out NTI info
write.csv(CRO_ses.mntd, "CRO_TS_FINAL_ses.mntd.csv", quote=F)
write.csv(CRO_ses.mntd$ntaxa, "CRO_TS_FINAL_ntaxa.csv", quote=F)
write.csv(CRO_ses.mntd$mntd.obs, "CRO_TS_FINAL_mntd.obs.csv", quote=F)
write.csv(CRO_ses.mntd$mntd.rand.mean, "CRO_TS_FINAL_mntd.rand.mean.csv", quote=F)
write.csv(CRO_ses.mntd$mntd.rand.sd, "CRO_TS_FINAL_mntd.rand.sd.csv", quote=F)
write.csv(CRO_ses.mntd$mntd.obs.rank, "CRO_TS_FINAL_mntd.obs.rank.csv", quote=F)
write.csv(CRO_ses.mntd$mntd.obs.z, "CRO_TS_FINAL_mntd.obs.z.csv", quote=F)
write.csv(CRO_ses.mntd$mntd.obs.p, "CRO_TS_FINAL_mntd.obs.p.csv", quote=F)
