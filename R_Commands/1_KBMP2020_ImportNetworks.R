# KBMP2020_ImportNetworks
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Fall 2018

# Outline
# Load tools
# Import phyloseq data (sample metadata, ASV taxonomy, tree, ASV counts)
# Modify taxonomy to deal with ambiguous assignments
# Modify metadata to group samples by major growth stages
# Add statistics for interaction edges (correlation magnitude, direction)
# Adjust edge weights for plotting

# Load tools

# 16S data handling, normalization, phylogenetics
library(phyloseq)
library(ape)
library(picante)
library(DESeq2)
library(edgeR)

# File reading and writing
library(readr)
library(MASS)

# Working with dataframes
library(reshape)
library(reshape2)
library(dplyr)
library(tidyr)
library(matrixStats)

# Making plots and tables
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(directlabels)
library(ggrepel)
library(ztable)
library(gridExtra)
library(ggstance)
library(viridis)
library(UpSetR)

# Network construction, visualization, statistics
library(SpiecEasi)
library(seqtime)
library(igraph)

# Import Phyloseq

# Import OTU TABLE from .txt file and make matrix
# The table has counts for DADA2 amplicon sequence variants (ASVs) in each sample. In phyloseq it is still called the "otu" table.
otu_table <- read.table("~/Documents/KBMP2020_Networks/Dataset/ASV_table.txt", header=TRUE, sep="\t", quote="",row.names=1)
otu_table = as.matrix(otu_table)

# Import METADATA
sample_data <- read.table("~/Documents/KBMP2020_Networks/Dataset/July18_Mapping_recoded.txt", header=TRUE, sep="\t", quote="",row.names=1)

# Import TREE (rooted)
phy_tree = read_tree("~/Documents/KBMP2020_Networks/Dataset/ASV_tree.nwk")

# Import TAXONOMY from csv file and make matrix
tax_table <- read.table("~/Documents/KBMP2020_Networks/Dataset/ASV_taxonomy.tsv",header=TRUE,sep=",",quote="",row.names=1,stringsAsFactors = FALSE)

# Modify TAXONOMY

# In the taxonomy table, there are a lot of ambiguous assignments. 
# unique(tax_table$Rank6)
# To simplify grouping for plots, the ASVs which were unassigned or given a generic assignment of uncultured or ambiguous were recoded as 
# "aggregate unclassified [rank]" in the table. 
# This allowed the ambiguous ASVs to remain in the dataset for calculations of relative abundances and diversity metrics 
# but to be easily grouped in plots.

AmbiguousGenera <- c("uncultured beta proteobacterium","uncultured alpha proteobacterium","g__","uncultured bacterium","uncultured soil bacterium","uncultured Chlorobi bacterium",
                     "uncultured","uncultured Bacteroidetes bacterium","Ambiguous_taxa","possible genus 04","uncultured Acidobacteria bacterium","uncultured forest soil bacterium",
                     "uncultured delta proteobacterium","uncultured Gemmatimonadetes bacterium","uncultured actinobacterium","uncultured Verrucomicrobia subdivision 3 bacterium",
                     "uncultured Chloroflexi bacterium","uncultured Ktedonobacter sp.","uncultured Acidobacteriales bacterium","uncultured Firmicutes bacterium",
                     "uncultured proteobacterium","uncultured organism","uncultured Candidatus Saccharibacteria bacterium","uncultured Armatimonadetes bacterium",
                     "uncultured Candidatus Saccharibacteria bacterium", "uncultured Desulfuromonadales bacterium")

AmbiguousGroups <- c("p__","c__","o__","f__")

tax_table$Rank6<-if_else(tax_table$Rank6 %in% AmbiguousGenera, "unclassified genus",tax_table$Rank6)
tax_table$Rank5<-if_else(tax_table$Rank5 %in% AmbiguousGroups | grepl("uncultur",tax_table$Rank5) | grepl("Ambig",tax_table$Rank5) | grepl("Unkn",tax_table$Rank5), "unclassified family",tax_table$Rank5)
tax_table$Rank4<-if_else(tax_table$Rank4 %in% AmbiguousGroups | grepl("uncultur",tax_table$Rank4) | grepl("Ambig",tax_table$Rank4) | grepl("Unkn",tax_table$Rank4), "unclassified order",tax_table$Rank4)
tax_table$Rank3<-if_else(tax_table$Rank3 %in% AmbiguousGroups | grepl("uncultur",tax_table$Rank3) | grepl("Ambig",tax_table$Rank3) | grepl("Unkn",tax_table$Rank3), "unclassified class",tax_table$Rank3)
tax_table$Rank2<-if_else(tax_table$Rank2 %in% AmbiguousGroups | grepl("uncultur",tax_table$Rank2) | grepl("Ambig",tax_table$Rank2) | grepl("Unkn",tax_table$Rank2), "unclassified phylum",tax_table$Rank2)

tax_table = as.matrix(tax_table)

# Use phyloseq functions to turn the imported data into phyloseq components
OTU = otu_table(otu_table, taxa_are_rows = TRUE)
TAX = tax_table(tax_table)
META = sample_data(sample_data)

# Merge phyloseq components into a class 
physeq = phyloseq(TAX,OTU,META,phy_tree)

# Modify METADATA

# Multiple vegetative stages (two-leaf through eight-leaf) were sampled in the dataset. Since the earlier stages had realtively few samples, we use 
# only the latest vegetative stage sampled in each site/year (six- or eight-leaf) as the "vegetative" data. A column is added to the metadata to 
# identify the groups of samples used for network analysis at vegetative, flowering, and senescent stages.

NetworkStage <- ifelse(sample_data(physeq)[,"Stage"]=="Flowering","Flowering",
                       ifelse(sample_data(physeq)[,"Stage"]=="Senescent","Senescent",
                              ifelse(sample_data(physeq)[,"Stage"]=="SixLeaf"|sample_data(physeq)[,"Stage"]=="EightLeaf","Vegetative","NotUsed")))

sample_data(physeq)[,"NetworkStage"] <- NetworkStage

# Import networks

ASVsInNetworks <- c()
project_path = "~/Documents/KBMP2020_Networks/AssembledNetworks"
ntwk_files <- dir(project_path, pattern="_NETWORK_")
for(nf in ntwk_files){ # For each network...
  print(paste("**********",nf,"**********",sep=" "))
  ntwk_path = paste(project_path,nf,sep="/")
  # Read the text file
  dat=read.table(ntwk_path, sep='',header=TRUE,row.names=NULL,check.names=FALSE)
  print(paste(nf,"loaded",sep=" "))
  # Make the data a matrix
  m=as.matrix(dat)
  # Import igraph object
  g=graph.adjacency(m,mode="upper",weighted=TRUE)
  E(g)$direction <- sign(E(g)$weight) # assign direction of correlation to edges
  E(g)$corr <- E(g)$weight # assign magnitude of correlation to edges
  E(g)$weight <- abs(E(g)$weight) # assign absolute value of correlation for edge width in plots
  # Name graph
  assign(nf,g)
  # Print information on nodes and edges
  print(paste("ASV nodes:",length(V(g)),sep=" "))
  print(paste("Inferred interactions:",length(E(g)),sep=" "))
  # Get ASV nodes list for filtering phyloseq data
  assign(paste("ASVs",strsplit(nf,"_")[[1]][3],sep="_"),V(g)$name)
  ASVsInNetworks <- c(ASVsInNetworks,V(g)$name)
  ASVsInNetworks <- unique(ASVsInNetworks)
}
