#Removing contaminant OTUs in the CROMO time series 16S rRNA dataset (2011-2017)
#Will be removing contaminants according to Cody Sheiks 2018 CoDL paper
#Following/modifying a script written by Kate Fullerton (UTK) to do this

#####################################Import Files and get Data Organized and Set Up#############################

#Upload neccessary packages
library(readr)
library(readxl)
library(dplyr)
library(tibble)
library(vegan) 
library (phyloseq)
library(Biostrings)
library(seqRFLP)
library('ggplot2')
theme_set(theme_bw())

#Upload Metadata, Rename meta.data, append row names
CROMO_Meta.data <- read_excel("CROMO_Meta.data.xlsx")
meta.data=as.data.frame(CROMO_Meta.data)
row.names(meta.data)=meta.data[,7]

#Upload Count Table, Rename count.table, and create as a matrix
CROMO_avg.counts <- read_excel("OTU_Count_Table - Averaged_Counts - for Contaminant Filtering.xlsx")
count.table=as.matrix(CROMO_avg.counts)
colnames(count.table)=meta.data[,7]

#Upload Taxonomy File, turn into matrix
CROMO_taxonomy <- read_excel("Original_taxonomy_file - for contaminant filtering.xlsx")
tax.table=as.matrix(CROMO_taxonomy)

#Upload Fasta file
bac_fasta <- readDNAStringSet("CROMO.raw.fasta")
seq_name = names(bac_fasta)
sequence = paste(bac_fasta)
bac_fasta_df <- data.frame(seq_name, sequence)

#Rename count.table and tax.table rownames
rownames(count.table)=paste("OTU", 1:nrow(count.table))
rownames(tax.table)=rownames(count.table)

#Reorder fasta file and export new fasta  ###IMPORTANT STEP, MOST FASTA FILES ARE NOT IN THE SAME ORDER AS YOUR COUNT TABLE
order=tax.table[,7]
bac_fasta_df=bac_fasta_df[match(order, bac_fasta_df$seq_name),]
dataframe2fas(bac_fasta_df, file = "CROMO_reorganized.fasta")

#Upload reorganized fasta file, change sequence names to OTU #s
bac_fasta2 <- readDNAStringSet("CROMO_reorganized.fasta")
names(bac_fasta2)=rownames(count.table)

############### NOW THAT DATA IS ALL ORGANIZED AND PREPARED IT CAN ALL BE MERGED INTO A PHYLOSEQ OBJECT ###################

#Merge objects into phyloseq object
merged <- phyloseq(otu_table(count.table, taxa_are_rows=TRUE), tax_table(tax.table), sample_data(meta.data), refseq(bac_fasta2))


########################################## BEGIN FIRST STAGES OF FILTERING ##############################################

#Normalize count table data
merged.norm=transform_sample_counts(merged, function(x) 100 * x/sum(x))

#First need to remove chloroplasts, mitochondria, eukarya, archaea, unknown
merged1=subset_taxa(merged, (Kingdom!="Archaea"))
#Removed 138 Archaea
merged2=subset_taxa(merged1, (Kingdom!="Eukaryota"))
#Removed 58 Eukarya
merged3=subset_taxa(merged2, (Kingdom!="unknown"))
#Removed 842 unknown
merged4=subset_taxa(merged3, (Family!="Mitochondria"))
#Removed 46 Mitochondria
merged5=subset_taxa(merged4, (Order!="Chloroplast"))
#Removed 94 Chloroplasts

#Remove from relative abundance dataset
merged1a=subset_taxa(merged.norm, (Kingdom!="Archaea"))
merged2a=subset_taxa(merged1a, (Kingdom!="Eukaryota"))
merged3a=subset_taxa(merged2a, (Kingdom!="unknown"))
merged4a=subset_taxa(merged3a, (Family!="Mitochondria"))
merged5a=subset_taxa(merged4a, (Order!="Chloroplast"))

#Extract fasta information after first round of filtering 
merged5_fasta=refseq(merged5)
seq_name = names(merged5_fasta)
sequence = paste(merged5_fasta)
merged5_fasta_df <- data.frame(seq_name, sequence)


#Could remove sequences that show up less than some threshold amount if one would like to
#merged#=filter_taxa(merged5, function(x) sum(x>3) > (0.1*length(x)), TRUE)
#sum(x>(minimum # of appearances) > (0.1*length(x)) --> percentage of samples the taxa must show up that number of times)

######################################## NOW FILTER OUT CoDL Contaminants #################################################
#***Code written by Kate Fullerton is modified here***


#List of contaminant genera + Verrucomicrobia and Chlamydiae phyla
contamination_1 <- c("Afipia", "Aquabacterium", "Asticcacaulis", "Aurantimonas", "Beijerinckia", "Bosea", "Bradyrhizobium", "Brevundimonas", "Caulobacter", "Craurococcus", "Devosia", "Hoefleae", "Mesorhizobium", "Methylobacterium", "Novosphingobium", "Ochrobactrum", "Paracoccus", "Pedomicrobium", "Phyllobacterium", "Rhizobium", "Roseomonas", "Sphingobium", "Sphingomonas", "Sphingopyxis", "Acidovorax", "Azoarcus", "Azospira", "Burkholderia", "Comamonas", "Cupriavidus", "Curvibacter", "Delftiae", "Duganella", "Herbaspirillum", "Janthinobacterium", "Kingella", "Leptothrix", "Limnobacter", "Massilia", "Methylophilus", "Methyloversatilis", "Neisseria", "Oxalobacter", "Pelomonas", "Polaromonas", "Ralstonia", "Schlegelella", "Sulfuritalea", "Undibacterium", "Variovorax", "Acinetobactera", "Enhydrobacter", "Enterobacter", "Escherichia", "Nevskia", "Pasteurella", "Pseudomonas", "Pseudoxanthomonas", "Psychrobacter", "Stenotrophomonas", "Xanthomonas", "unclassified Acidobacteria Gp2", "Aeromicrobium", "Actinomyces", "Arthrobacter", "Beutenbergia", "Brevibacterium", "Corynebacterium", "Curtobacterium", "Dietzia", "Geodermatophilus", "Janibacter", "Kocuria", "Microbacterium", "Micrococcus", "Microlunatus", "Patulibacter", "Propionibacterium", "Rhodococcus", "Tsukamurella", "Chryseobacterium", "Dyadobacter", "Flavobacterium", "Hydrotalea", "Niastella", "Olivibacter", "Parabacteroides", "Pedobacter", "Prevotella", "Wautersiella", "Deinococcus", "Abiotrophia", "Bacillus", "Brevibacillus", "Brochothrix", "Facklamia", "Lactobacillus", "Paenibacillus", "Ruminococcus", "Staphylococcus", "Streptococcus", "Veillonella", "Fusobacterium") #Kate's list of genera
contamination_2 <- c("Verrucomicrobia", "Chlamydiae") # phyla that can be problematic from CoDL contaminants


#Make dataframes of contaminant taxa IDs
merged5_bac_taxa=as(tax_table(merged5), "matrix")
merged5_bac_taxa=as.data.frame(merged5_bac_taxa)


#Pulls out contaminant genera from the whole dataset
merged5_bac_shiek_contaminant1 <- merged5_bac_taxa %>%
  rownames_to_column(var = "OTU") %>%
  #does not remove genera that are unclassified
  filter(Genus %in% contamination_1) %>%
  select(OTU, Kingdom:Species)
#485 putative contaminants after genera screening 

#Pulls out potential contaminant phyla from whole dataset
merged5_bac_shiek_contaminant2 <- merged5_bac_taxa %>%
  rownames_to_column(var = "OTU") %>%
  #does not remove genera that are unclassified
  filter(Phylum %in% contamination_2) %>%
  select(OTU, Kingdom:Species)
#746 contaminants after phylum screening


#Combine the two sets of pulled contaminants
merged5_bac_shiek_contaminant = rbind(merged5_bac_shiek_contaminant1, merged5_bac_shiek_contaminant2)

bac_contamination_otu <- merged5_bac_shiek_contaminant$OTU


#generate phyloseq object of putative contaminants from relative abundance physeq object
merged_contamination1 <- subset_taxa(merged5, Genus %in% contamination_1)
merged_contamination2 <- subset_taxa(merged5, Phylum %in% contamination_2)
merged_contamination <- merge_phyloseq(merged_contamination1, merged_contamination2)
merged_contamination

#subset from normalized data
merged_norm_contamination1 <- subset_taxa(merged5a, Genus %in% contamination_1)
merged_norm_contamination2 <- subset_taxa(merged5a, Phylum %in% contamination_2)
merged_norm_contamination <- merge_phyloseq(merged_norm_contamination1, merged_norm_contamination2)
merged_norm_contamination


#Create function to pull out otu table to analyze with vegan
veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}


#this code is to pull out the sequences of the flagged contaminants
#then you can input this into blast
bac_fasta_filter <- merged5_fasta_df %>%
  inner_join(merged5_bac_shiek_contaminant, by = c("seq_name"="OTU"))


#determine abundance of individual contaminant ASVs within each site
bac_contamination_rel <- veganotu(merged_norm_contamination) %>%
  t()%>%
  data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  #gather(-ASV, key = "Site", value = "rel") %>%
  #filter(rel > 0.001) %>%
  inner_join(bac_fasta_filter, by = c("OTU"="seq_name"))


#export fasta for blastn
bac_fasta_filter_fa <- select(bac_fasta_filter, seq_name, sequence)
dataframe2fas(bac_fasta_filter_fa, file = "CRO_TS_shiek_contaminant.fasta")


######################## IMPORT BLAST RESULTS AND APPEND THEM TO CONTAMINATION PHYLOSEQ TAXONOMY INFO #####################

#import in blast results
CROMO_contam_blast <- read.csv("CROMO_TS-Contaminant_Seqs_BLAST_Taxonomy.csv")
colnames(CROMO_contam_blast)[1] <- "OTU"

#Extract Contam_rel taxa table
contam_taxa=as(tax_table(merged_norm_contamination), "matrix")
contam_taxa=as.data.frame(contam_taxa)

#Reorder BLAST taxonomy file to match order in the phyloseq object
order2=row.names(contam_taxa)
CROMO_contam_blast=CROMO_contam_blast[match(order2, CROMO_contam_blast$OTU),]

#Replace Species name in taxa table with BLAST results
contam_taxa=contam_taxa[,-c(7)] #Deletes current species column
Species=CROMO_contam_blast[,2] 
contam_taxa=cbind(contam_taxa, Species) #Binds Kingdom:Genus with BLAST Species classifications
contam_taxa=as.matrix(contam_taxa)#Convert to matrix so this taxonomy can be used in the phyloseq object

#Replace taxa table in contamination phyloseq with the new one
tax_table(merged_norm_contamination)<-contam_taxa

#Prevent use of scientific notation (messes up sorting in the next steps)
options (scipen = 999)


##################################### INTIAL STAGES OF COMBING THROUGH CONTAMINANTS ######################################
#I had 1231 contaminants... this would be hard to manually dig through, so I did initial filtering by looking at OTU 
#abundance of each contaminant (relative to entire dataset) to get a smaller list to manually look at
#Filtered out low abundance OTUs first, and then OTUs that are present in 10% or less samples

#Extract total dataset otu sums, get relative abundance of each otu throughout dataset
otu_sums=rowSums(otu_table(merged5))
otu_sums=as.data.frame(otu_sums)
total_otu_reads=colSums(otu_sums)
otu_sums=(otu_sums/total_otu_reads)*100
all.ids=rownames(otu_sums)
otu_sums=as.data.frame(cbind.data.frame(all.ids, otu_sums[,1]))
colnames(otu_sums)=c("OTUs", "Abundance")
otu_sums=otu_sums[order(otu_sums$Abundance),]
rownames(otu_sums)=1:22816

#Get contaminant OTU ids, and subset OTU sums data
otu.ids=rownames(otu_table(merged_norm_contamination))
contam_otu_sums=otu_sums[otu_sums$OTUs %in% otu.ids,] #Now have contaminant OTU relative abundances from within entire dataset

#Reorder dataframe by abundances 
contam_otu_sums=contam_otu_sums[order(contam_otu_sums$Abundance),] #sum of abundance = 6.25% of OTU rel. abundance of all species
rownames(contam_otu_sums)=1:1231 #Greatest total abundance is 0.84%

#Remove contaminant OTUs that were singletons and averaged out to zero during tech. rep averaging
contam.otus_0=length(which(contam_otu_sums$Abundance == 0)) #314 OTUs
contam_remove1=contam_otu_sums[1:314,]
remove1_otu.ids=contam_remove1[,1]                     

#Remove contaminants with 0 total abundance (from averaging) from the full dataset
alltaxa=taxa_names(merged5)
contam1_taxa=alltaxa[!(alltaxa %in% remove1_otu.ids)]
merged6=prune_taxa(contam1_taxa, merged5)

#Remove same contaminants from contamination phyloseq and look at remaining samples
merged_contamination2=prune_taxa(contam1_taxa, merged_contamination)

#Take a look at total counts of species (relative abundance was a bit hard to work with, switched totoal counts of each OTU)
contam_counts=as(otu_table(merged_contamination2), "matrix")
contam_counts=as.data.frame(rowSums(contam_counts))
contam_ids=rownames(contam_counts)
contam_counts=contam_counts[,1]
contam_counts=cbind.data.frame(contam_ids, contam_counts)
colnames(contam_counts)=c("ID", "Counts")
rownames(contam_counts)=1:917
contam_counts=contam_counts[order(contam_counts$Counts),]

#Get OTU IDs of OTUs with 10 or less counts (indicates presence in 10 or less samples which is less than 10% of the dataset)
otus_10.or.less_counts=length(which(contam_counts$Counts <=10)) #590 with 10 or less counts
contam_remove2=contam_counts[1:590,]
remove2_out.ids=contam_remove2[,1]

#Remove the 590 potential contaminants with 10 or less total counts from total dataset
alltaxa2=taxa_names(merged6)
contam2_taxa=alltaxa2[!(alltaxa2 %in% remove2_out.ids)]
merged7=prune_taxa(contam2_taxa, merged6)

#make datset to remove high abundance ones
otus_high.abund=contam_counts[591:917,]
remove3_otu.ids=otus_high.abund[,1]
contam3.taxa=alltaxa2[!(alltaxa2 %in% remove3_otu.ids)]
low.abund_merged.contam=prune_taxa(contam3.taxa, merged_contamination2)

#Remove these from remaining contaminant OTUs
merged_contamination3=prune_taxa(contam2_taxa, merged_contamination2) #Leaves 327 OTUs

#Export taxonomy and counts of remaining 327 samples and manually curate
contam3.taxa=as(tax_table(merged_contamination3), "matrix")
contam3.count=as(otu_table(merged_contamination3), "matrix")

contam3.taxa.count=cbind.data.frame(contam3.count, contam3.taxa)

write.csv(contam3.taxa.count, "CROMO_TS_contaminants_for_manual_curation2.csv")


########################## REMOVE FINAL CONTAMINANTS AFTER MANUAL CURATION ###################################
#Investigated each individual OTU, based on published literature about an isolate and looked at their count distribution
#in the dataset to determine if higher abundance potential contaminants were likely contaminants or not

#Upload dataset that contains final decisions(contaminant or not) on last 327 potential contaminants
CROMO_TS_Contaminant_Manual_Filtering_Results <- read_excel("CROMO_TS_Contaminant_Manual_Filtering_Results.xlsx")
final.contam.coding=as.data.frame(CROMO_TS_Contaminant_Manual_Filtering_Results)

#Reorganize and get OTU IDs for contaminants to remove
final.contam.coding=final.contam.coding[order(final.contam.coding$Number),]
rownames(final.contam.coding)=1:327
final.contam.otus=final.contam.coding[1:256,]
not.contam.otus=final.contam.coding[257:327,]
final.contam.otus=final.contam.otus[,1]
not.contam.otus=not.contam.otus[,1]

#Remove final set of contaminants from dataset
alltaxa3=taxa_names(merged7)
contam_final_taxa=alltaxa3[!(alltaxa3 %in% final.contam.otus)]
not.contam_final_taxa=alltaxa3[!(alltaxa3 %in% not.contam.otus)]
merged_contamination4=prune_taxa(not.contam_final_taxa, merged_contamination3) 
merged8=prune_taxa(contam_final_taxa, merged7)


######################3#################### CONTAMINANT FILTERING DONE ###############################################

###################### NOW NEED TO FILTER OUT ~SINGLETON TAXA THAT AVERAGED OUT TO ABUNDANCE=0 ########################


#Taking a look at all the count sums for the whole dataset
all_counts=as(otu_table(merged5), "matrix")
all_counts=as.data.frame(rowSums(all_counts))
all_ids=rownames(all_counts)
all_counts=all_counts[,1]
all_counts=cbind.data.frame(all_ids, all_counts)
colnames(all_counts)=c("ID", "Counts")
all_counts=all_counts[order(all_counts$Counts),] #There are 8526 OTUs where rowSums=0 (rounded to 0 when averaging tech. reps)

#Get OTU IDs for taxa = 0
singleton.otus=all_counts[1:8526,]
singleton.otu.ids=singleton.otus[,1]

#Filter these out of the dataset for final dataset
alltaxa4=taxa_names(merged8)
singleton_taxa=alltaxa4[!(alltaxa4 %in% singleton.otu.ids)]
merged.final=prune_taxa(singleton_taxa, merged8)



#Final removal stats:
  #Started with = 23,994 OTUs --> 6,210,850 reads
  #Singletons = 8526 (314 which were identified as Shiek CoDL contaminants)
  #Archaea = 138 --> 7,271 reads
  #Eukarya = 58 --> 1,475 reads
  #Unknown = 842 --> 9,127 reads
  #Mitochondria = 46 --> 1,295 reads
  #Chloroplasts = 94 --> 7,881 reads
  #CoDL contaminants = 1231 - 314 (0 abund) - 71 (not contams) = 846
  #CoDL contaminants = 846 taxa --> 209,745 reads
    #590 contaminant OTUs with 10 or fewer reads in dataset account for 1,802 reads total
    #202 contaminant OTUs with greater abundance account for 109,774 total reads
        #combined 792 OTUs with 111,576 reads
    #Simkania nevegensis and Akkermansia muciniphila account for 54 OTUs and 98,169 reads

#Total OTUs removed = 10,550 (80% singletons, 20% contaminants)

#Final dataset = 13,444 OTUs --> 5,974,056 reads --> removed 3.8% of reads 

#CoDL contaminants were 3.4% of total reads and 3.5% of OTUs from original dataset

#Export Final Count and Taxonomy Tables for Final Data Analyses
final_counts=as(otu_table(merged.final), "matrix")
final_taxonomy=as(tax_table(merged.final), "matrix")

#Write final datasets to .csv files 
write.csv(final_counts, "CROMO_TS_filtered_FINAL_counts.csv")
write.csv(final_taxonomy, "CROMO_TS_filtered_FINAL_taxonomy.csv")

#Export Final Fasta
final_fasta=as(refseq(merged.final), "character")
final_fasta=as.data.frame(final_fasta)
dataframe2fas(final_fasta, file="CROMO_TS_filtered_FINAL.fasta")



