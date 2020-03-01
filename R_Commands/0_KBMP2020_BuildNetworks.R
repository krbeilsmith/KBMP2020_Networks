# KBMP2020_BuildNetworks
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Fall 2018

# The sequencing data used here were collected by Matt Perisin as part of his 2016 dissertation "THE DYNAMICS OF BACTERIAL COMMUNITIES 
# ASSOCIATED WITH ARABIDOPSIS THALIANA" for the Committee on Microbiology at the University of Chicago. The raw data consisted of 1427 sets of
# paired 250bp Illumina MiSeq reads.

# Data pre-processing was performed with in Qiime2. Primers were removed with cutadapt and DADA2 was used to model errors and find ASVs. 
# The initial table included 1421 samples with 10,987 sequence variants. Sequences were filtered for PhiX, mitochondrial, and chloroplast 
# sequences (removed 35 sequence variants and 102 samples). Samples with notes in the metadata indicating any irregularities in the collection 
# were excluded (removed 45 samples). Sequence variants with a frequency lower than 2 counts and samples with fewer than 10 reads were
# excluded (removed 19 features and 8 samples). This left a table with 10,803 taxa and 1,272 samples.

# The microbiome dataset was handled with phyloseq, following the tutorial here: https://joey711.github.io/phyloseq/import-data.html

# ASV prevalence and abundance varied between plant tissues and developmental stages as well as across field sites and between study years (see 
# analysis in Beilsmith, Perisin, and Bergelson 2020: Natural bacterial assemblages in Arabidopsis thaliana tissues become more distinguishable and 
# diverse during host development). Root, rosette leaf, stem, and silique samples were therefore considered independently. For each tissue type, 
# samples from different developmental stages and sites were used to infer independent networks of interactions among the ASVs that were observed
# in both years of study.

# Samples in the study were collected at the following developmental stages: 2-leaf rosettes, 4-leaf rosettes, 6-leaf rosettes, 8-leaf rosettes,
# flowering plants, and senescent plants. Since the earlier vegetative stages yielded a relatively small number of samples, the latest rosette stage
# sampled (6-leaf in Year 2 and 8-leaf in Year 1) was used as the "vegetative" plant sample. Networks were assembled for microbes in these vegetative
# plants as well as flowering and senescent plants.

# Rarely observed and low-abundance ASVs were removed from the dataset before network inference because they resulted in spurious network edges. ASVs 
# present in three or more samples and with over 250 sequence counts total in the dataset were retained.

# The network of relationships among ASVs was inferred with both an inverse covariance and a correlational approach:

# Relationships among ASVs were inferred with the SPIEC-EASI pipeline (https://github.com/zdk123/SpiecEasi,
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004226). 
# *** After SPIEC-EASI, weights were assigned to the edges of the network based on Pearson correlation coefficients for the abundances of the nodes.

# Relationships among ASVs were also inferred with SparCC (https://bitbucket.org/yonatanf/sparcc/src/default/,
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002687).
# *** After SparCC, edges were filtered for correlation significance @ alpha = 0.001 and correlations with absolute value >0.3.

# Outline
# Load tools
# Import phyloseq data (sample metadata, ASV taxonomy, tree, ASV counts)
# Modify taxonomy to deal with ambiguous assignments
# Modify metadata to group samples by major growth stages (vegetative, flowering, senescent)
# Infer networks with SPIEC-EASI and with SparCC
# Filter to ASVs supported by two years of observations in a tissue/site/stage
# Filter to remove rare and low-abundance ASVs
# SPIEC-EASI + assign edge weights with correlation matrix
# Write adjacency matrix to file and print a plot of the network
# SparCC
# Write adjacency matrix to file and print a plot of the network

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
# To simplify grouping, the ASVs which were unassigned or given a generic assignment of uncultured or ambiguous were recoded as 
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

# Infer networks

# Conditions in which to infer networks 
# (not soil; plant tissues only at stages where they are observed; not flowers or cauline leaves because of the low sample count)
tissues <- c("Roots","RosLeaves","Stems","Siliques")
stages <- c("Vegetative","Flowering","Senescent")
sites <- c("ME","WW")

# Samples in the conditions of interest
physeq_select <- subset_samples(physeq, PlantPart %in% tissues & NetworkStage %in% stages & Site %in% sites)

# Get the columns of the sample data corresponding to the variables of interest and find all the combinations of values that occur in the dataset
ConditionsFrame <- data.frame(unique(sample_data(physeq_select)[,c("Site","NetworkStage","PlantPart")]),row.names=c())
ConditionsFrame <- data.frame(lapply(ConditionsFrame, as.character), stringsAsFactors=FALSE)

# Manual cleanup (removing flowering stems and flowering siliques from conditions list because of low sample counts)
ConditionsFrame <- ConditionsFrame[-c(1,3,8,11),]

# For each tissue during each developmental stage at each site...
for(cf in rownames(ConditionsFrame)){
  print(ConditionsFrame[cf,])
  # Find the ASVs observed in Year 1 samples
  subset_physeq_1 <- subset_samples(physeq, PlantPart==ConditionsFrame[cf,"PlantPart"] & NetworkStage==ConditionsFrame[cf,"NetworkStage"] & Site==ConditionsFrame[cf,"Site"] & Year=="1")
  subset_physeq_1 <- prune_taxa(taxa_sums(subset_physeq_1)>0, subset_physeq_1)
  # Filter to remove rare ASVs (observed in fewer than three samples)
  subset_physeq_1 <- transform_sample_counts(subset_physeq_1, function(x) ifelse(x>0,1,0))
  subset_physeq_1 <- prune_taxa(taxa_sums(subset_physeq_1)>2, subset_physeq_1)
  # Find the ASVs observed in Year 2 samples
  subset_physeq_2 <- subset_samples(physeq, PlantPart==ConditionsFrame[cf,"PlantPart"] & NetworkStage==ConditionsFrame[cf,"NetworkStage"] & Site==ConditionsFrame[cf,"Site"] & Year=="2")
  subset_physeq_2 <- prune_taxa(taxa_sums(subset_physeq_2)>0, subset_physeq_2)
  # Filter to remove rare ASVs (observed in fewer than three samples)
  subset_physeq_2 <- transform_sample_counts(subset_physeq_2, function(x) ifelse(x>0,1,0))
  subset_physeq_2 <- prune_taxa(taxa_sums(subset_physeq_2)>2, subset_physeq_2)
  # Retain the ASVs observed in both years of study
  supported_ASVs <- unique(c(taxa_names(subset_physeq_1),taxa_names(subset_physeq_2)))
  subset_physeq_3 <- subset_samples(physeq, PlantPart==ConditionsFrame[cf,"PlantPart"] & NetworkStage==ConditionsFrame[cf,"NetworkStage"] & Site==ConditionsFrame[cf,"Site"])
  filtered_physeq <- prune_taxa(taxa_names(subset_physeq_3) %in% supported_ASVs, subset_physeq_3)
  # Filter to remove low-count ASVs (less than 250)
  filtered_physeq <- prune_taxa(taxa_sums(filtered_physeq)>250, filtered_physeq)
  # Remove samples that drop to zero counts
  filtered_physeq <- prune_samples(sample_sums(filtered_physeq)>0, filtered_physeq)
  print(filtered_physeq)
  
  # Build Networks
  # SPIEC-EASI
  print("Building inverse covariance network...")
  se.out <- spiec.easi(filtered_physeq, method="mb", lambda.min.ratio=1e-1, icov.select.params=list(rep.num=100)) # run SPIEC-EASI
  se.graph <- adj2igraph(se.out$refit, vertex.attr=list(name=taxa_names(filtered_physeq))) # Make graph from matrix
  #assign(paste("SEnet",ConditionsFrame[cf,"PlantPart"],ConditionsFrame[cf,"NetworkStage"],ConditionsFrame[cf,"Site"],sep="_"),se.graph)
  tOTU <- t(otu_table(filtered_physeq)) # transpose OTU table
  # Get matrix of ASV correlations
  cor <- cor(tOTU, method = "pearson")
  # Add correlation information to network
  cor.graph <- graph.adjacency(cor,weighted=TRUE,diag=FALSE,mode="undirected")
  se.graph <- intersection(cor.graph,se.graph,byname=TRUE,keep.all.vertices = FALSE)
  # Export graph
  mat <- get.adjacency(se.graph, type="both", attr="weight_1", names=TRUE, sparse=FALSE)
  write.matrix(mat,file=paste("~/Documents/KBMP2020_Networks/AssembledNetworks/","SE_NETWORK_",
                              ConditionsFrame[cf,"PlantPart"],ConditionsFrame[cf,"NetworkStage"],ConditionsFrame[cf,"Site"],sep="")) 
  # Print an image of the network
  l <- layout_with_fr(se.graph)
  ppi <- 300
  png(paste("~/Documents/KBMP2020_Networks/AssembledNetworks/","SEnet",
            ConditionsFrame[cf,"PlantPart"],ConditionsFrame[cf,"NetworkStage"],ConditionsFrame[cf,"Site"],sep=""), 
      width=6*ppi, height=6*ppi, res=ppi)
  par(mar=c(0,0,0,0))
  plot(se.graph, layout=l, vertex.label=NA, edge.color="black",vertex.size=3,vertex.color = "white")
  dev.off()
  
  # SparCC
  print("Building correlation network...")
  tOTU <- t(otu_table(filtered_physeq)) # transpose OTU table
  cc <- sparccboot(tOTU,sparcc.params=list(iter=20,inner_iter=10),R=10) # run SparCC
  cctest <- pval.sparccboot(cc,sided="both") # Empirical p-values for correlations
  cc2 <- ifelse(cctest$pvals<0.001,cctest$cors,0) # P-value threshold
  cc2 <- ifelse(abs(cctest$cors)>0.3,cctest$cors,0) # Correlation threshold for absolute value
  cc3 <- diag(ncol(tOTU)) # Put correlations that are below p-val and above correlation threshold into matrix and add row/column names
  cc3[upper.tri(cc3, diag=FALSE)] <- cc2
  colnames(cc3) <- colnames(tOTU)
  rownames(cc3) <- colnames(tOTU)
  cc.graph  <- graph.adjacency(cc3,weighted=TRUE,diag=FALSE,mode="upper",add.colnames="name") # Make graph from matrix
  #assign(paste("CCnet",ConditionsFrame[cf,"PlantPart"],ConditionsFrame[cf,"NetworkStage"],ConditionsFrame[cf,"Site"],sep="_"),cc.graph)
  # Export graph
  mat <- get.adjacency(cc.graph, type="both", attr="weight", names=TRUE, sparse=FALSE)
  write.matrix(mat,file=paste("~/Documents/KBMP2020_Networks/AssembledNetworks/","CC_NETWORK_",
                              ConditionsFrame[cf,"PlantPart"],ConditionsFrame[cf,"NetworkStage"],ConditionsFrame[cf,"Site"],sep="")) 
  # Print an image of the network
  # l <- layout_with_fr(cc.graph) # Use same layout as inverse covariance graph for comparison
  ppi <- 300
  png(paste("~/Documents/KBMP2020_Networks/AssembledNetworks/","CCnet",
            ConditionsFrame[cf,"PlantPart"],ConditionsFrame[cf,"NetworkStage"],ConditionsFrame[cf,"Site"],sep=""), 
      width=6*ppi, height=6*ppi, res=ppi)
  par(mar=c(0,0,0,0))
  plot(cc.graph, layout=l, weighted=TRUE,vertex.label=NA, edge.color="black",edge.width=E(cc.graph)$weight,vertex.size=3,vertex.color = "white")
  dev.off()
}
